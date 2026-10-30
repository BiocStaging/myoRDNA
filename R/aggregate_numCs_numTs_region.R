#' Aggregate methylated and unmethylated counts within a genomic region
#'
#' Sums numCs and numTs across sites within the specified region for each
#' sample, grouped by condition.
#'
#' @param region_res A nested list where each top-level element is a sample
#'   (e.g., \code{region_res[["SRR28431752"]][["18S"]][["numCs"]]}).
#' @param region A string specifying the region of interest (e.g., "18S").
#' @param conditions A named list where each name is a condition
#'   (e.g., OLD or YOUNG) and each element is a character vector of sample
#'   IDs belonging to that condition.
#'
#' @return An \code{S4Vectors::DataFrame} with columns:
#' \describe{
#'   \item{Sample_ID}{The sample identifier.}
#'   \item{Condition}{The condition to which the sample belongs.}
#'   \item{aggregate_numCs}{The sum of methylated counts within the region,
#'     excluding missing values.}
#'   \item{aggregate_numTs}{The sum of unmethylated counts within the region,
#'     excluding missing values.}
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