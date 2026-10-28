#' Perform Differential Methylation Analysis at Individual Sites
#'
#' This function analyzes differential methylation at individual genomic 
#' positions.
#' It fits a logistic regression model (binomial, quasi-binomial, or 
#' beta-binomial) for each position, computes odds ratios, p-values, and 
#' predicted methylation probabilities,and applies multiple-testing correction. 
#' Positions with insufficient data are skipped.
#'
#' @param all_data A data frame containing methylation data with columns:
#'   - `position`: Genomic position.
#'   - `numCs`: Number of methylated reads.
#'   - `numTs`: Number of unmethylated reads.
#'   - `Group`: Binary treatment variable (e.g., 0 for control, 1 for 
#'   treatment).
#' @param model_type A character string specifying the statistical model to 
#' use.
#'   Options: `"beta-binomial"` (default), `"quasi-binomial"`, or 
#'   `"binomial"`.
#' @param p_adjust_method A character string specifying the method for p-value
#'  adjustment. Default is `"BH"` (Benjamini-Hochberg).
#'
#' @return A data frame with columns:
#'   - `position`: Genomic position.
#'   - `beta_0`: Intercept term from the model.
#'   - `beta_1`: Coefficient for the methylation difference.
#'   - `odds_ratio`: Odds ratio for the methylation difference.
#'   - `p_value`: Raw p-value for the methylation difference.
#'   - `p_adjust`: Adjusted p-value after multiple-testing correction.
#'   - `p_cond_1`: Predicted methylation probability for the control group.
#'   - `p_cond_2`: Predicted methylation probability for the treatment group.
#'   - `meth_diff`: Methylation difference between groups (in percentage 
#'   points).
#'   - `log10_p_adjust`: Negative log10 of the adjusted p-value.
#'
diff_meth_site <- function(
    all_data,
    model_type = "beta-binomial",
    p_adjust_method = "BH")
{
  # Check input table
  if (!(is.data.frame(all_data) ||
        methods::is(all_data, "DataFrame"))) {
    stop("'all_data' must be a data.frame or S4Vectors::DataFrame.",
         call. = FALSE)
  }
  
  required_columns <- c("position", "numCs", "numTs", "Group")
  
  if (!all(required_columns %in% names(all_data))) {
    stop("'all_data' must contain position, numCs, numTs, and Group.",
         call. = FALSE)
  }
  
  # Use a base data.frame for regression
  all_data <- as.data.frame(all_data)
  
  if (!all(vapply(
    all_data[required_columns],
    is.numeric,
    logical(1)
  ))) {
    stop("position, numCs, numTs, and Group must be numeric.",
         call. = FALSE)
  }
  
  if (!all(vapply(
    all_data[required_columns],
    function(x) all(is.finite(x)),
    logical(1)
  ))) {
    stop("Required columns must not contain missing or infinite values.",
         call. = FALSE)
  }
  
  if (any(all_data$numCs < 0 |
          all_data$numCs != floor(all_data$numCs)) ||
      any(all_data$numTs < 0 |
          all_data$numTs != floor(all_data$numTs))) {
    stop("'numCs' and 'numTs' must contain non-negative whole numbers.",
         call. = FALSE)
  }
  
  if (!all(all_data$Group %in% c(0, 1))) {
    stop("'Group' must contain only 0 and 1.",
         call. = FALSE)
  }
  
  # Check model choice
  valid_models <- c(
    "beta-binomial",
    "quasi-binomial",
    "binomial"
  )
  
  if (!is.character(model_type) ||
      length(model_type) != 1L ||
      is.na(model_type) ||
      !model_type %in% valid_models) {
    stop("'model_type' must be 'beta-binomial', 'quasi-binomial',",
         "or 'binomial'.",
         call. = FALSE)
  }
  
  # Check adjustment method
  if (!is.character(p_adjust_method) ||
      length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods) {
    stop("'p_adjust_method' must be a valid p-value adjustment method.",
         call. = FALSE)
  }
  
  if (model_type == "beta-binomial" &&
      !requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required for beta-binomial regression.",
         call. = FALSE)
  }
  
  # Initialise results
  site_results <- S4Vectors::DataFrame(
    position = numeric(),
    beta_0 = numeric(),
    beta_1 = numeric(),
    odds_ratio = numeric(),
    p_value = numeric()
  )
  
  # Preallocate storage for individual results
  positions <- unique(all_data$position)
  
  site_results_list <- S4Vectors::SimpleList(
    vector("list", length(positions))
  )
  
  for (i in seq_along(positions)) {
    pos <- positions[i]
    pos_data <- all_data[
      all_data$position == pos,
      ,
      drop = FALSE
    ]
    
    # Exclude observations with no reads
    pos_data <- pos_data[
      pos_data$numCs + pos_data$numTs > 0,
      ,
      drop = FALSE
    ]
    
    if (nrow(pos_data) < 3L ||
        length(unique(pos_data$Group)) < 2L) {
      next
    }
    
    fit_result <- tryCatch({
      if (model_type == "beta-binomial") {
        model <- glmmTMB::glmmTMB(
          cbind(numCs, numTs) ~ Group,
          family = glmmTMB::betabinomial(link = "logit"),
          data = pos_data,
          control = glmmTMB::glmmTMBControl(
            optimizer = stats::optim,
            optArgs = list(method = "Nelder-Mead"),
            optCtrl = list(maxit = 5000)
          )
        )
        
        beta_0 <- glmmTMB::fixef(model)$cond["(Intercept)"]
        beta_1 <- glmmTMB::fixef(model)$cond["Group"]
        p_value <- summary(model)$coefficients$cond[
          "Group", "Pr(>|z|)"
        ]
        
      } else if (model_type == "quasi-binomial") {
        model <- stats::glm(
          cbind(numCs, numTs) ~ Group,
          family = stats::quasibinomial(link = "logit"),
          data = pos_data
        )
        
        beta_0 <- stats::coef(model)["(Intercept)"]
        beta_1 <- stats::coef(model)["Group"]
        p_value <- summary(model)$coefficients[
          "Group", "Pr(>|t|)"
        ]
        
      } else {
        model <- stats::glm(
          cbind(numCs, numTs) ~ Group,
          family = stats::binomial(link = "logit"),
          data = pos_data
        )
        
        beta_0 <- stats::coef(model)["(Intercept)"]
        beta_1 <- stats::coef(model)["Group"]
        p_value <- summary(model)$coefficients[
          "Group", "Pr(>|z|)"
        ]
      }
      
      c(
        beta_0 = unname(beta_0),
        beta_1 = unname(beta_1),
        p_value = unname(p_value)
      )
      
    }, error = function(e) {
      warning(
        "Skipping position ", pos, ": ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    })
    
    if (is.null(fit_result)) {
      next
    }
    
    site_results <- rbind(
      site_results,
      S4Vectors::DataFrame(
        position = pos,
        beta_0 = unname(fit_result["beta_0"]),
        beta_1 = unname(fit_result["beta_1"]),
        odds_ratio = unname(exp(fit_result["beta_1"])),
        p_value = unname(fit_result["p_value"])
      )
    )
  }
  
  # Remove results with missing coefficients or p-values
  invalid_rows <- is.na(site_results$beta_0) |
    is.na(site_results$beta_1) |
    is.na(site_results$p_value)
  
  if (any(invalid_rows)) {
    warning(
      sum(invalid_rows),
      " rows with missing coefficients or p-values were removed.",
      call. = FALSE
    )
    
    site_results <- site_results[!invalid_rows, , drop = FALSE]
  }
  
  # Multiple-testing correction
  site_results$p_adjust <- stats::p.adjust(
    site_results$p_value,
    method = p_adjust_method
  )
  
  # Predicted methylation probabilities
  site_results$p_cond_1 <- stats::plogis(site_results$beta_0)
  site_results$p_cond_2 <- stats::plogis(
    site_results$beta_0 + site_results$beta_1
  )
  
  site_results$meth_diff <- 100 * (
    site_results$p_cond_2 - site_results$p_cond_1
  )
  
  site_results$log10_p_adjust <- -log10(site_results$p_adjust)
  
  return(site_results)
}