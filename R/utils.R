attach_by_id <- function(
    df1,
    df2,
    id_column = "Sample_ID",
    samples,
    condition_name = "Category",
    cor_columns = NULL)
{
  # Check input tables
  if (!is.data.frame(df1) || !is.data.frame(df2)) {
    stop("'df1' and 'df2' must be data frames.",
         call. = FALSE)
  }
  
  # Check ID column
  if (!is.character(id_column) ||
      length(id_column) != 1L ||
      is.na(id_column) ||
      !nzchar(id_column)) {
    stop("'id_column' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!id_column %in% names(df1) ||
      !id_column %in% names(df2)) {
    stop("'id_column' must exist in both data frames.",
         call. = FALSE)
  }
  
  # Check condition mapping
  if (!is.list(samples) ||
      is.data.frame(samples) ||
      length(samples) == 0L ||
      is.null(names(samples)) ||
      anyNA(names(samples)) ||
      any(!nzchar(names(samples))) ||
      anyDuplicated(names(samples)) > 0L) {
    stop("'samples' must be a non-empty list with unique condition names.",
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
    stop("Each condition must contain a non-empty character",
         "vector of sample IDs.",
         call. = FALSE)
  }
  
  sample_ids <- unlist(samples, use.names = FALSE)
  
  if (anyDuplicated(sample_ids) > 0L) {
    stop("Each sample ID must appear only once in 'samples'.",
         call. = FALSE)
  }
  
  # Check output condition column name
  if (!is.character(condition_name) ||
      length(condition_name) != 1L ||
      is.na(condition_name) ||
      !nzchar(condition_name)) {
    stop("'condition_name' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  # Check optional correlation columns
  if (!is.null(cor_columns)) {
    if (!is.character(cor_columns) ||
        length(cor_columns) != 2L ||
        anyNA(cor_columns) ||
        any(!nzchar(cor_columns)) ||
        anyDuplicated(cor_columns) > 0L) {
      stop("'cor_columns' must contain two distinct column names.",
           call. = FALSE)
    }
  }
  
  # Map sample IDs to conditions
  id_conditions <- stats::setNames(
    rep(names(samples), lengths(samples)),
    sample_ids
  )
  
  merged_df <- merge(
    df1,
    df2,
    by = id_column,
    suffixes = c("_df1", "_df2")
  )
  
  if (condition_name %in% names(merged_df)) {
    stop("'condition_name' already exists in the merged data.",
         call. = FALSE)
  }
  
  merged_df[[condition_name]] <- unname(
    id_conditions[as.character(merged_df[[id_column]])]
  )
  
  # Optional Pearson correlation
  if (!is.null(cor_columns)) {
    x <- merged_df[[cor_columns[1]]]
    y <- merged_df[[cor_columns[2]]]
    
    keep <- is.finite(x) & is.finite(y)
    
    cor_test_all <- stats::cor.test(
      x[keep],
      y[keep],
      method = "pearson"
    )
  }
  
  return(merged_df)
}