#' Extract Methylation Site Data for a Genomic Region
#'
#' Extracts methylation counts (`numCs`, `numTs`) and positions for a 
#' specified genomic region
#' across all samples. Assigns a binary treatment variable (`Group`) based on 
#' the first group in `samples`.
#' Skips samples with mismatched vector lengths and removes zero-coverage 
#' sites.
#'
#' @param region_res A nested list where each top-level element corresponds to 
#' a sample, and each sub-element corresponds to a genomic region 
#'  (e.g., `region_res[[1]]$18S`).
#' @param samples A named list of sample IDs grouped by condition 
#' (e.g., `list(OLD = c("sample1", "sample2"), YOUNG = c("sample3"))`).
#' @param region A character string specifying the genomic region to extract. 
#' Default is `"18S"`.
#'
#' @return A data frame with columns:
#'   - `sample`: Sample identifier.
#'   - `position`: Genomic position.
#'   - `numCs`: Number of methylated reads.
#'   - `numTs`: Number of unmethylated reads.
#'   - `Group`: Binary treatment variable (0 for control group, 1 otherwise).
#' 
extract_site_data <- function(
    region_res,
    samples,
    region = "18S")
{
  # Check region data
  if (!is.list(region_res) ||
      is.data.frame(region_res) ||
      length(region_res) == 0L ||
      is.null(names(region_res)) ||
      anyNA(names(region_res)) ||
      any(!nzchar(names(region_res))) ||
      anyDuplicated(names(region_res)) > 0L) {
    stop("'region_res' must be a non-empty list with unique sample names.",
         call. = FALSE)
  }
  
  # Check conditions: downstream models use a binary Group variable
  if (!is.list(samples) ||
      is.data.frame(samples) ||
      length(samples) != 2L ||
      is.null(names(samples)) ||
      anyNA(names(samples)) ||
      any(!nzchar(names(samples))) ||
      anyDuplicated(names(samples)) > 0L) {
    stop("'samples' must be a named list of two uniquely named conditions.",
         call. = FALSE)
  }
  
  valid_samples <- vapply(samples, function(x) {
    is.character(x) &&
      is.null(dim(x)) &&
      length(x) > 0L &&
      !anyNA(x) &&
      all(nzchar(x))
  }, logical(1))
  
  if (!all(valid_samples)) {
    stop("Each condition must contain a non-empty character vector 
         of sample IDs.",
         call. = FALSE)
  }
  
  sample_ids <- unlist(samples, use.names = FALSE)
  
  if (anyDuplicated(sample_ids) > 0L) {
    stop("Each sample ID must appear only once in 'samples'.",
         call. = FALSE)
  }
  
  if (!all(sample_ids %in% names(region_res))) {
    stop("All sample IDs in 'samples' must exist in names(region_res).",
         call. = FALSE)
  }
  
  # Check region
  if (!is.character(region) ||
      length(region) != 1L ||
      is.na(region) ||
      !nzchar(region)) {
    stop("'region' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  control_samples <- samples[[1]]
  
  all_data <- S4Vectors::DataFrame(
    sample = character(),
    position = numeric(),
    numCs = numeric(),
    numTs = numeric(),
    Group = numeric()
  )
  
  # Preallocate storage for each sample
  sample_results <- vector("list", length(sample_ids))
  
  for (i in seq_along(sample_ids)) {
    sam <- sample_ids[i]
    region_data <- region_res[[sam]][[region]]
    
    numCs <- region_data$numCs
    numTs <- region_data$numTs
    position <- region_data$position
    
    Group <- if (sam %in% control_samples) 0 else 1
    
    df <- S4Vectors::DataFrame(
      sample = rep(sam, length(position)),
      position = position,
      numCs = numCs,
      numTs = numTs,
      Group = rep(Group, length(position))
    )
    
    # Remove zero-coverage sites.
    # Missing counts remain for validation by diff_meth_site().
    coverage <- df$numCs + df$numTs
    keep <- is.na(coverage) | coverage > 0
    
    sample_results[[i]] <- df[keep, , drop = FALSE]
  }
  
  all_data <- do.call(rbind, sample_results)
  
  return(all_data)
}