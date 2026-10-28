#' Perform Differential Methylation Analysis for a Genomic Region
#'
#' Aggregates methylation counts (`numCs` and `numTs`) for a specified genomic 
#' region, fits a logistic regression model (binomial, quasi-binomial, or 
#' beta-binomial), and computes methylation differences, p-values, and 
#' predicted probabilities.
#'
#' @param region_res A nested list where each top-level element corresponds to 
#' a sample, and each sub-element corresponds to a genomic region 
#'  (e.g., `region_res[[1]]$18S`).
#'   Each region must contain methylation data (e.g., `numCs`, `numTs`).
#' @param region A character string specifying the name of the genomic region 
#' to analyze.
#' @param samples A named vector of sample IDs, where names correspond to 
#' conditions
#' (e.g., 
#' `c(Young = c("sample1", "sample2"), Old = c("sample3", "sample4"))`).
#' @param model_type A character string specifying the statistical model to 
#' use. Options: `"binomial"`, `"quasi-binomial"`, or `"beta-binomial"` 
#' (default).
#' @param p_adjust_method A character string specifying the method for p-value 
#' adjustment.
#'   Default is `"BH"` (Benjamini-Hochberg). Note: Adjustment is applied in 
#'   the wrapper function.
#'
#' @return A data frame with columns:
#'   - `beta_0`: Intercept term from the model.
#'   - `beta_1`: Coefficient for the methylation difference.
#'   - `p_value`: Raw p-value for the methylation difference.
#'   - `meth_diff`: Estimated methylation difference between conditions 
#'   (in percentage points).
#' 
diff_meth_region <- function(
    region_res,
    region,
    samples,
    model_type = "beta-binomial",
    p_adjust_method = "BH")
{
  # Basic input checks
  if (!is.list(region_res) ||
      is.data.frame(region_res) ||
      length(region_res) == 0L) {
    stop("'region_res' must be a non-empty list.",
         call. = FALSE)
  }
  
  if (!is.character(region) ||
      length(region) != 1L ||
      is.na(region) ||
      !nzchar(region)) {
    stop("'region' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  # This model compares exactly two conditions
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
  
  valid_sample_ids <- vapply(samples, function(x) {
    is.character(x) &&
      is.null(dim(x)) &&
      length(x) > 0L &&
      !anyNA(x) &&
      all(nzchar(x))
  }, logical(1))
  
  if (!all(valid_sample_ids)) {
    stop("Each condition in 'samples' must contain a non-empty character 
      vector of sample IDs.",
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
  
  valid_models <- c(
    "binomial",
    "quasi-binomial",
    "beta-binomial"
  )
  
  if (!is.character(model_type) ||
      length(model_type) != 1L ||
      is.na(model_type) ||
      !model_type %in% valid_models) {
    stop("'model_type' must be 'binomial', 'quasi-binomial', 
         or 'beta-binomial'.",
         call. = FALSE)
  }
  
  if (!is.character(p_adjust_method) ||
      length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods) {
    stop("'p_adjust_method' must be a valid p-value adjustment method.",
         call. = FALSE)
  }
  
  # Aggregate counts and prepare data for regression
  data <- as.data.frame(
    aggregate_numCs_numTs_region(
      region_res = region_res,
      region = region,
      conditions = samples
    )
  )
  
  # First condition is the reference group
  data$Group <- ifelse(
    data$Condition == names(samples)[1],
    0,
    1
  )
  
  if (model_type == "binomial") {
    model <- stats::glm(
      cbind(aggregate_numCs, aggregate_numTs) ~ Group,
      data = data,
      family = stats::binomial(link = "logit")
    )
    
    beta_0 <- stats::coef(model)["(Intercept)"]
    beta_1 <- stats::coef(model)["Group"]
    p_value <- summary(model)$coefficients["Group", "Pr(>|z|)"]
    
  } else if (model_type == "quasi-binomial") {
    model <- stats::glm(
      cbind(aggregate_numCs, aggregate_numTs) ~ Group,
      data = data,
      family = stats::quasibinomial(link = "logit")
    )
    
    beta_0 <- stats::coef(model)["(Intercept)"]
    beta_1 <- stats::coef(model)["Group"]
    p_value <- summary(model)$coefficients["Group", "Pr(>|t|)"]
    
  } else {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package 'glmmTMB' is required for beta-binomial regression.",
           call. = FALSE)
    }
    
    model <- glmmTMB::glmmTMB(
      cbind(aggregate_numCs, aggregate_numTs) ~ Group,
      family = glmmTMB::betabinomial(link = "logit"),
      data = data
    )
    
    beta_0 <- glmmTMB::fixef(model)$cond["(Intercept)"]
    beta_1 <- glmmTMB::fixef(model)$cond["Group"]
    p_value <- summary(model)$coefficients$cond[
      "Group", "Pr(>|z|)"
    ]
  }
  
  beta_0 <- unname(beta_0)
  beta_1 <- unname(beta_1)
  
  # Predicted methylation: second condition minus first
  p_reference <- stats::plogis(beta_0)
  p_comparison <- stats::plogis(beta_0 + beta_1)
  
  meth_diff <- 100 * (p_comparison - p_reference)
  
  region_results <- S4Vectors::DataFrame(
    beta_0 = beta_0,
    beta_1 = beta_1,
    p_value = unname(p_value),
    meth_diff = meth_diff
  )
  
  return(region_results)
}