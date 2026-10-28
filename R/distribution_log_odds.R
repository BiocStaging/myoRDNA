#' Calculate Log-Odds of Methylation for a Chromosome
#'
#' Computes the log-odds of methylation (log(numCs/numTs)) for each CpG site
#' in a specified chromosome across samples. Useful for assessing the 
#' suitability of logistic regression models (e.g., in `methylKit` or custom 
#' analyses).
#'
#' @param meth_object A list of `methylRaw` objects (from `methylKit`),
#'   where each object contains methylation data for a sample.
#' @param treatment A vector of treatment/condition labels for each sample.
#'   Must have the same length as `meth_object`.
#' @param chromo A character string specifying the chromosome to analyze.
#'
#' @return A data frame with columns:
#'   - `log_odds`: Log-odds of methylation (log(numCs/numTs)) for each CpG 
#'   site.
#'   - `sample_id`: Identifier for the sample.
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
#' # Calculate log-odds for CpG sites on the rDNA chromosome
#' log_odds <- distribution_log_odd(
#'   meth_object = loaded_data,
#'   treatment = treatment,
#'   chromo = "chrRDNAm"
#' )
#'
#' # Display the first few results
#' head(log_odds)
#'
#' @export
distribution_log_odd <- function(
    meth_object,
    treatment,
    chromo)
{
  # Check methylation input
  if (!(is.list(meth_object) ||
        methods::is(meth_object, "methylRawList")) ||
      is.data.frame(meth_object) ||
      length(meth_object) == 0L) {
    stop("'meth_object' must be a non-empty list of methylRaw objects",
         "or a methylRawList.",
         call. = FALSE)
  }
  
  valid_objects <- vapply(
    seq_along(meth_object),
    function(i) methods::is(meth_object[[i]], "methylRaw"),
    logical(1)
  )
  
  if (!all(valid_objects)) {
    stop("Each element of 'meth_object' must be a methylRaw object.",
         call. = FALSE)
  }
  
  # Check treatment labels
  if (!(is.character(treatment) ||
        is.numeric(treatment) ||
        is.factor(treatment)) ||
      !is.null(dim(treatment)) ||
      length(treatment) == 0L ||
      anyNA(treatment) ||
      any(!nzchar(as.character(treatment)))) {
    stop("'treatment' must be a non-empty character, numeric, or factor",
         "vector without missing or empty labels.",
         call. = FALSE)
  }
  
  if (length(meth_object) != length(treatment)) {
    stop("Length of 'meth_object' must equal length of 'treatment'.",
         call. = FALSE)
  }
  
  # Check chromosome name
  if (!is.character(chromo) ||
      length(chromo) != 1L ||
      is.na(chromo) ||
      !nzchar(chromo)) {
    stop("'chromo' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  num_samples <- length(meth_object)
  log_odds_list <- vector("list", num_samples)
  
  for (i in seq_len(num_samples)) {
    # Retrieve sample data using methylKit accessors
    sample_data <- methylKit::getData(meth_object[[i]])
    sample_id <- methylKit::getSampleID(meth_object[[i]])
    
    chr_data <- sample_data[
      sample_data$chr %in% chromo,
      ,
      drop = FALSE
    ]
    
    # Calculate log-odds where both counts are positive
    lo <- rep(NA_real_, nrow(chr_data))
    
    valid_counts <- !is.na(chr_data$numCs) &
      !is.na(chr_data$numTs) &
      chr_data$numCs > 0 &
      chr_data$numTs > 0
    
    lo[valid_counts] <- log(
      chr_data$numCs[valid_counts] /
        chr_data$numTs[valid_counts]
    )
    
    log_odds_list[[i]] <- S4Vectors::DataFrame(
      log_odds = lo,
      sample_id = rep(sample_id, length(lo))
    )
  }
  
  log_odds_out <- do.call(rbind, log_odds_list)
  
  return(log_odds_out)
}