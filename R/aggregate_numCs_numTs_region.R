#' Computes the aggregate methylation for a specified genomic region across 
#' samples, grouped by their respective conditions.
#'
#' @param region_res A nested list where each top-level element is a sample 
#' (e.g., `region_res$SRR28431752$18S$numCs`).
#' @param region A string specifying the region of interest (e.g., `"18S"`).
#' @param conditions A named list where each name is a condition 
#' (e.g., `OLD`, `YOUNG`) and each element is a vector of sample IDs belonging 
#' to that condition.
#' @param med A logical value indicating whether to compute the median 
#' (if `TRUE`) or the mean (if `FALSE`) for the aggregate methylation values. 
#' Default is `FALSE`.
#' 
#' @return A dataframe with columns:
#' \describe{
#'   \item{sample_ID:}{The sample identifier.}
#'   \item{condition:}{The condition to which the sample belongs.}
#'   \item{aggregate_meth:}{The aggregate methylation value for the 
#'   specified region.}
#' }
#'
aggregate_numCs_numTs_region <- function(
    region_res,
    region,
    conditions)
{
  if (!is.list(region_res)) {
    stop("'region_res' must be a list.", call. = FALSE)
  }
  
  if (!is.character(region) ||
      length(region) != 1L ||
      is.na(region) ||
      !nzchar(region)) {
    stop("'region' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!is.list(conditions) || length(conditions) == 0L) {
    stop("'conditions' must be a non-empty named list.",
         call. = FALSE)
  }
  
  condition_names <- names(conditions)
  
  if (is.null(condition_names) ||
      anyNA(condition_names) ||
      any(!nzchar(condition_names)) ||
      anyDuplicated(condition_names) > 0L) {
    stop("'conditions' must have unique, non-empty names.",
         call. = FALSE)
  }
  
  valid_samples <- vapply(
    conditions,
    function(x) {
      is.character(x) &&
        length(x) > 0L &&
        !anyNA(x) &&
        all(nzchar(x))
    },
    logical(1)
  )
  
  if (!all(valid_samples)) {
    stop("Each condition must contain a non-empty character vector", 
         "of sample IDs.",
         call. = FALSE)
  }
  
  n_samples <- sum(lengths(conditions))
  
  result_df <- S4Vectors::DataFrame(
    Sample_ID = character(n_samples),
    Condition = character(n_samples),
    aggregate_numCs = numeric(n_samples),
    aggregate_numTs = numeric(n_samples)
  )
  
  i <- 0L
  
  for (condition in condition_names) {
    for (sam in conditions[[condition]]) {
      region_data <- region_res[[sam]][[region]]
      
      i <- i + 1L
      
      result_df$Sample_ID[i] <- sam
      result_df$Condition[i] <- condition
      result_df$aggregate_numCs[i] <-
        sum(region_data$numCs, na.rm = TRUE)
      result_df$aggregate_numTs[i] <-
        sum(region_data$numTs, na.rm = TRUE)
    }
  }
  
  return(result_df)
}