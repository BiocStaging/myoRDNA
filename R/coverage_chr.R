#' Plot Average Coverage per Chromosome
#'
#' @description
#' Calculates the average coverage for each chromosome across all samples and 
#' plots it as a bar chart.
#'
#' @param meth_object A list of data frames or `methylRaw` objects containing 
#' methylation data.
#' @param num_samples Integer. Number of samples to process.
#' @param chromosome_names Character vector. Names of chromosomes to include 
#' in the analysis.
#'
#' @return A `ggplot` bar chart visualizing the average coverage per 
#' chromosome.
#' 
#' @examples
#' # Load the example methylation dataset
#' file_path <- system.file(
#'   "extdata",
#'   "filtered_meth_raw_chr1_chrRDNAm_small.rds",
#'   package = "myoRDNA",
#'   mustWork = TRUE
#' )
#'
#' loaded_data <- readRDS(file_path)
#'
#' # Define the number of samples
#' treatment <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
#' num_samples <- length(treatment)
#'
#' # Extract unique chromosome names
#' chromosome_names <- get_unique_chr_names(
#'   loaded_data,
#'   num_samples
#' )
#'
#' # Display the chromosome names
#' print(chromosome_names)
#'
#' @note
#' - The function assumes the input data frames contain a `coverage` column.
#' - The plot uses `ggplot2` for visualization.
#' 
#' @export
coverage_chr <- function(
    meth_object,
    num_samples,
    chromosome_names)
{
  # Basic input checks
  if (!(is.list(meth_object) ||
        methods::is(meth_object, "methylRawList")) ||
      length(meth_object) == 0L) {
    stop("'meth_object' must be a non-empty list or methylRawList.",
         call. = FALSE)
  }
  
  if (!is.numeric(num_samples) ||
      length(num_samples) != 1L ||
      !is.finite(num_samples) ||
      num_samples < 1 ||
      num_samples != floor(num_samples) ||
      num_samples > length(meth_object)) {
    stop(
      "'num_samples' must be a positive whole number no greater",
      "than the number of samples.",
      call. = FALSE
    )
  }
  
  if (!is.character(chromosome_names) ||
      length(chromosome_names) == 0L ||
      anyNA(chromosome_names) ||
      any(!nzchar(chromosome_names))) {
    stop("'chromosome_names' must be a non-empty character vector",
         "without missing or empty names.",
         call. = FALSE)
  }
  
  # Retrieve data using the methylKit accessor
  sample_data <- lapply(seq_len(num_samples), function(i) {
    methylKit::getData(meth_object[[i]])
  })
  
  chromosome_names <- unique(chromosome_names)
  
  # Calculate average coverage for each chromosome
  avg_coverage <- vapply(chromosome_names, function(chromo) {
    sample_means <- vapply(sample_data, function(dat) {
      keep <- dat$chr %in% chromo &
        !is.na(dat$coverage) &
        dat$coverage >= 10
      
      mean(dat$coverage[keep], na.rm = TRUE)
    }, numeric(1))
    
    mean(sample_means, na.rm = TRUE)
  }, numeric(1))
  
  return(
    S4Vectors::DataFrame(
      chromosome = chromosome_names,
      avg_coverage = unname(avg_coverage)
    )
  )
}