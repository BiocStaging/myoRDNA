#' Analyze Methylation Across Genomic Regions
#'
#' This function analyzes differential methylation across specified genomic 
#' regions for a given set of samples. It uses an external script 
#' (`diff_meth_region.R`) to compute methylation differences, p-values, and 
#' adjusted p-values for each region. Results are combined into a single data 
#' frame with multiple-testing correction applied.
#'
#' @param region_res A nested list where each top-level element corresponds to 
#' a sample, and each sub-element corresponds to a genomic region 
#'  (e.g., `region_res[[1]]$18S`). Each region must contain methylation data 
#'  (e.g., `numCs`, `coverage`).
#' @param samples A vector of sample IDs to include in the analysis.
#' @param model_type A character string specifying the statistical model to 
#' use for differential methylation analysis. Default is `"beta-binomial"`.
#' @param p_adjust_method A character string specifying the method for p-value 
#' adjustment. Default is `"BH"` (Benjamini-Hochberg).
#'
#' @return A data frame with columns:
#'   - `region`: Name of the genomic region.
#'   - `beta_0`: Intercept term from the model.
#'   - `beta_1`: Coefficient for the methylation difference.
#'   - `p_value`: Raw p-value for the methylation difference.
#'   - `p_adjust`: Adjusted p-value after multiple-testing correction.
#'   - `meth_diff`: Estimated methylation difference between conditions.
#'
analyze_methylation_regions <- function(
    region_res,
    samples,
    model_type = "beta-binomial",
    p_adjust_method = "BH")
{
  # Basic input checks
  if (!is.list(region_res) ||
      is.data.frame(region_res) ||
      length(region_res) == 0L) {
    stop("'region_res' must be a non-empty list.", call. = FALSE)
  }
  
  if (missing(samples) || length(samples) == 0L) {
    stop("'samples' must be supplied and must not be empty.",
         call. = FALSE)
  }
  
  if (!is.character(model_type) ||
      length(model_type) != 1L ||
      is.na(model_type) ||
      !nzchar(model_type)) {
    stop("'model_type' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!is.character(p_adjust_method) ||
      length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods) {
    stop("'p_adjust_method' must be a valid p-value adjustment method.",
         call. = FALSE)
  }
  
  # Define regions
  regions <- names(region_res[[1]])
  
  # Initialise the results
  out_results <- S4Vectors::DataFrame(
    region = character(),
    beta_0 = numeric(),
    beta_1 = numeric(),
    p_value = numeric(),
    p_adjust = numeric(),
    meth_diff = numeric()
  )
  
  region_results <- lapply(regions, function(reg) {
    results <- diff_meth_region(
      region_res = region_res,
      region = reg,
      samples = samples,
      model_type = model_type,
      p_adjust_method = p_adjust_method
    )
    
    results <- S4Vectors::DataFrame(results)
    results$region <- rep(reg, nrow(results))
    
    return(results)
  })
  
  if (length(region_results) > 0L) {
    out_results <- do.call(rbind, unname(region_results))
  }
  
  # Multiple-testing correction across regions
  out_results$p_adjust <- stats::p.adjust(
    out_results$p_value,
    method = p_adjust_method
  )
  
  out_results <- out_results[
    ,
    c("region", "beta_0", "beta_1",
      "p_value", "p_adjust", "meth_diff"),
    drop = FALSE
  ]
  
  return(out_results)
}