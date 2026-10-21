#' Calculate Absolute rDNA Copy Number
#'
#' This function calculates the **absolute rDNA copy number** for each sample 
#' in a `Meth_object` (e.g., a list of methylation data or coverage data). It 
#' compares the coverage of the rDNA region to the coverage of either the 
#' autosomes (chromosomes 1-19) or the entire genome (excluding rDNA).
#' The absolute copy number is adjusted for ploidy (diploid or haploid).
#'
#' @param Meth_object - A named list of `methylRaw` objects (from `methylKit`) 
#' or data frames, where each element represents a sample and contains columns 
#' `chr` (chromosome name) and `coverage` (coverage depth).
#' @param num_samples - An integer specifying the number of samples in 
#' `Meth_object` to process.
#' @param autosomes - A logical (`TRUE`/`FALSE`) indicating whether to use 
#' **autosomes (chromosomes 1-19)**
#'   as the reference for genome coverage. If `FALSE`, the reference is the 
#'   entire genome **excluding rDNA**.
#'   Default: `TRUE`.
#' @param rDNA_chr - A character string specifying the name of the rDNA 
#' chromosome (e.g., `"chrRDNAm"`).
#' @param diploid - A logical (`TRUE`/`FALSE`) indicating whether the genome is 
#' **diploid** (2 copies of each chromosome).
#'   If `TRUE`, the absolute rDNA copy number is multiplied by 2. Default: 
#'   `TRUE`.
#'
#' @return
#' This function **prints** a tab-separated table to the console with the 
#' following columns:
#' \describe{
#'   \item{Sample ID:}{ The sample identifier 
#'   (extracted from `Meth_object[[i]]@sample.id`).}
#'   \item{Genome Depth:}{ The mean coverage depth for the reference region 
#'   (autosomes or genome excluding rDNA).}
#'   \item{rDNA Depth:}{ The mean coverage depth for the rDNA region.}
#'   \item{absolute rDNA CN:}{ The absolute rDNA copy number, calculated as 
#'   `(cov_rDNA / cov_wgenome) * 2` (if diploid) or 
#'   `(cov_rDNA / cov_wgenome)` (if haploid).}
#' }
#' The function does **not** return a value but prints results directly.
#'
#' @examples
#' # Load the example mouse WGBS dataset
#' file_path <- system.file(
#'   "extdata",
#'   "filtered_meth_raw_chr1_chrRDNAm_small.rds",
#'   package = "myoRDNA",
#'   mustWork = TRUE
#' )
#'
#' loaded_data <- readRDS(file_path)
#'
#' # Define experimental conditions
#' # 0 = OLD, 1 = YOUNG
#' treatment <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
#' num_samples <- length(treatment)
#'
#' # Construct a methylRawList
#' loaded_data <- methylKit::methylRawList(
#'   loaded_data,
#'   treatment = treatment
#' )
#'
#' # Calculate absolute rDNA copy number
#' abs_out <- absolute_rDNA_CN(
#'   Meth_object = loaded_data,
#'   num_samples = num_samples,
#'   autosomes = TRUE,
#'   rDNA_chr = "chrRDNAm",
#'   diploid = TRUE
#' )
#'
#' # Inspect the results
#' print(abs_out)
#'
#' @note
#' - If `autosomes = TRUE`, the reference coverage (`cov_wgenome`) is 
#' calculated from chromosomes 1-19.
#'   Ensure these chromosomes exist in your data.
#' - If `autosomes = FALSE`, the reference coverage is calculated from 
#' **all chromosomes except `rDNA_chr`**.
#' - The function assumes `Meth_object` is a **named list** and that each 
#' element has 
#' a `@sample.id` slot (if using `methylRaw` objects) or can be accessed 
#' via `names(Meth_object)[i]`.
#' - If `cov_wgenome` is `0`, the function will produce `Inf` or `NA` for 
#' `abs_rDNA_CN`. 
#' Add a check for `cov_wgenome > 0` if this is a concern.
#' - The function uses `cat` to print tab-separated output. To capture the 
#' output, redirect it to a file or use `capture.output().
#'
#' @export
absolute_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
{
  # Initialise a Bioconductor DataFrame with one row per sample
  results_df <- S4Vectors::DataFrame(
    Sample_ID = character(num_samples),
    Genome_Depth = numeric(num_samples),
    rDNA_Depth = numeric(num_samples),
    Absolute_rDNA_CN = numeric(num_samples)
  )

  ploidy_factor <- if (diploid) 2 else 1

  for (i in seq_len(num_samples))
  {
    # Retrieve sample data using the methylKit accessor
    sample_data <- methylKit::getData(Meth_object[[i]])

    if (autosomes)
    {
      chr_data <- sample_data[
        sub("^chr", "", sample_data$chr) %in%
          as.character(seq_len(19)),
        ,
        drop = FALSE
      ]
    } else
    {
      chr_data <- sample_data
    }

    rDNA_data <- sample_data[
      sample_data$chr %in% rDNA_chr,
      ,
      drop = FALSE
    ]

    cov_wgenome <- mean(chr_data$coverage, na.rm = TRUE)
    cov_rDNA <- mean(rDNA_data$coverage, na.rm = TRUE)

    results_df$Sample_ID[i] <-
      methylKit::getSampleID(Meth_object[[i]])
    results_df$Genome_Depth[i] <- cov_wgenome
    results_df$rDNA_Depth[i] <- cov_rDNA
    results_df$Absolute_rDNA_CN[i] <-
      ploidy_factor * cov_rDNA / cov_wgenome
  }

  return(results_df)
}
