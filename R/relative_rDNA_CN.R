#' Calculate Relative rDNA Copy Number
#'
#' This function calculates the **relative rDNA copy number** for each sample 
#' in a `Meth_object` (e.g., a list of methylation data or coverage data).
#'  It compares the **total coverage** of the rDNA region to the 
#'  **total coverage** of either the autosomes (chromosomes 1-19) or the 
#'  entire genome.
#' The relative copy number is adjusted for ploidy (diploid or haploid).
#'
#' @param Meth_object A named list of `methylRaw` objects (from `methylKit`) 
#' or data frames, where each element represents a sample and contains 
#' columns `chr` (chromosome name)` and `coverage` (coverage depth).
#' @param num_samples An integer specifying the number of samples in 
#' `Meth_object` to process.
#' @param autosomes A logical (`TRUE`/`FALSE`) indicating whether to use 
#' **autosomes (chromosomes 1-19)** as the reference for genome coverage. 
#' If `FALSE`, the reference is the **entire genome**. Default: `TRUE`.
#' @param rDNA_chr A character string specifying the name of the rDNA 
#' chromosome (e.g., `"chrRDNAm"`).
#' @param diploid A logical (`TRUE`/`FALSE`) indicating whether the genome is 
#' **diploid** (2 copies of each chromosome). If `TRUE`, the relative rDNA 
#' copy number is multiplied by 2. Default: `TRUE`.
#'
#' @return An S4Vectors::DataFrame containing Sample_ID,
#'   Total_Number_of_Reads, Number_of_Reads_rDNA, and Relative_rDNA_CN.
#'   The read-count columns contain summed coverage across recorded sites.
#'
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
#' treatment <- c(
#'   0, 0, 0, 0, 0, 0,
#'   1, 1, 1, 1
#' )
#'
#' num_samples <- length(treatment)
#'
#' # Calculate relative rDNA copy number
#' rel_out <- relative_rDNA_CN(
#'   Meth_object = loaded_data,
#'   num_samples = num_samples,
#'   autosomes = TRUE,
#'   rDNA_chr = "chrRDNAm",
#'   diploid = TRUE
#' )
#'
#' # Display the results
#' print(rel_out)
#'
#'
#' @note
#' - If `autosomes = TRUE`, the reference coverage (`sum_cov_wgenome`) is 
#' calculated from **all chromosomes** in the sample. This is a 
#' **deviation from the function name**, as it does not restrict to autosomes 
#' (chromosomes 1-19). To fix this, update the filtering logic for `chr_data` 
#' to include only  autosomes 
#' (e.g., `chr_data <- Meth_object[[i]][Meth_object[[i]]$chr` 
#' `%in% paste0("chr", 1:19), ]`).
#' - If `autosomes = FALSE`, the reference coverage is calculated from 
#' **all chromosomes** in the sample.
#' - The function assumes `Meth_object` is a **named list** and that each 
#' element has a `@sample.id` slot (if using `methylRaw` objects) or can be 
#' accessed via `names(Meth_object)[i]`.
#' - If `sum_cov_wgenome` is `0`, the function will produce `Inf` or `NA` for 
#' `rel_rDNA_CN`. Add a check for `sum_cov_wgenome > 0` if this is a concern.
#' - The function uses `cat` to print tab-separated output. To capture the 
#' output, redirect it to a file or use `capture.output()`.
#'
#' @export
relative_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
{
  # Check sample container
  if (!(is.list(Meth_object) ||
        methods::is(Meth_object, "methylRawList")) ||
      is.data.frame(Meth_object) ||
      length(Meth_object) == 0L) {
    stop("'Meth_object' must be a non-empty list or methylRawList.",
         call. = FALSE)
  }
  
  # Check number of samples
  if (!is.numeric(num_samples) ||
      length(num_samples) != 1L ||
      !is.finite(num_samples) ||
      num_samples < 1 ||
      num_samples != floor(num_samples) ||
      num_samples > length(Meth_object)) {
    stop("'num_samples' must be a positive whole number no greater than",
         "the number of samples in 'Meth_object'.",
         call. = FALSE)
  }
  
  # Check logical arguments
  if (!is.logical(autosomes) ||
      length(autosomes) != 1L ||
      is.na(autosomes)) {
    stop("'autosomes' must be a single TRUE or FALSE.",
         call. = FALSE)
  }
  
  if (!is.logical(diploid) ||
      length(diploid) != 1L ||
      is.na(diploid)) {
    stop("'diploid' must be a single TRUE or FALSE.",
         call. = FALSE)
  }
  
  # Check chromosome name
  if (!is.character(rDNA_chr) ||
      length(rDNA_chr) != 1L ||
      is.na(rDNA_chr) ||
      !nzchar(rDNA_chr)) {
    stop("'rDNA_chr' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  results_df <- S4Vectors::DataFrame(
    Sample_ID = character(num_samples),
    Total_Number_of_Reads = numeric(num_samples),
    Number_of_Reads_rDNA = numeric(num_samples),
    Relative_rDNA_CN = numeric(num_samples)
  )
  
  for (i in seq_len(num_samples)) {
    sample_object <- Meth_object[[i]]
    
    sample_data <- methylKit::getData(sample_object)
    sample_id <- methylKit::getSampleID(sample_object)
    
    chromosomes <- as.character(sample_data$chr)
    
    if (autosomes) {
      reference_rows <- sub("^chr", "", chromosomes) %in%
        as.character(seq_len(19))
    } else {
      reference_rows <- rep(TRUE, nrow(sample_data))
    }
    
    rdna_rows <- chromosomes %in% rDNA_chr
    
    sum_cov_wgenome <- sum(
      sample_data$coverage[reference_rows],
      na.rm = TRUE
    )
    
    sum_cov_rDNA <- sum(
      sample_data$coverage[rdna_rows],
      na.rm = TRUE
    )
    
    results_df$Sample_ID[i] <- sample_id
    results_df$Total_Number_of_Reads[i] <- sum_cov_wgenome
    results_df$Number_of_Reads_rDNA[i] <- sum_cov_rDNA
    results_df$Relative_rDNA_CN[i] <-
      sum_cov_rDNA / sum_cov_wgenome
  }
  
  return(results_df)
}