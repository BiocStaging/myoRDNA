#' Computes the aggregate methylation for a specified genomic region across 
#' samples, grouped by their respective conditions.
#'
#' @param region_res A nested list where each top-level element is a sample 
#' (e.g., `region_res$SRR28431752$18S$numCs`).
#' @param region A string specifying the region of interest (e.g., `"18S"`).
#' @param conditions A named list where each name is a condition 
#' (e.g., `CONTROL`, `CONDITION 1`) and each element is a vector of sample IDs 
#' belonging to that condition.
#' @param med A logical value indicating whether to compute the median 
#' (if `TRUE`) or the mean (if `FALSE`) for the aggregate methylation values. 
#' Default is `FALSE`.
#' 
#' @return A dataframe with columns:
#' \describe{
#'   \item{sample_ID:}{The sample identifier.}
#'   \item{condition:}{The condition to which the sample belongs.}
#'   \item{aggregate_meth:}{The aggregate methylation value for the specified 
#'   region.}
#' }
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
#' # Define the experimental groups
#' samples <- list(
#'   OLD = c(
#'     "SRR28431752", "SRR28431753", "SRR28431754",
#'     "SRR28431755", "SRR28431756", "SRR28431757"
#'   ),
#'   YOUNG = c(
#'     "SRR28431748", "SRR28431749",
#'     "SRR28431750", "SRR28431751"
#'   )
#' )
#'
#' treatment <- c(
#'   0, 0, 0, 0, 0, 0,
#'   1, 1, 1, 1
#' )
#'
#' num_samples <- length(treatment)
#'
#' # Construct a methylRawList object
#' loaded_data <- methods::new(
#'   "methylRawList",
#'   loaded_data,
#'   treatment = treatment
#' )
#'
#' # Extract rDNA methylation data
#' meth_rDNA <- get_rDNA_data(
#'   loaded_data,
#'   num_samples,
#'   rDNA_chr_name = "chrRDNAm"
#' )
#'
#' # Adjust genomic coordinates
#' meth_rDNA_shifted <- shifting(
#'   meth_rdna_raw = meth_rDNA,
#'   cutoff = 35000,
#'   treatment = treatment
#' )
#'
#' # Separate methylation data into rDNA regions
#' region_res <- sort_by_rDNA_region(
#'   meth_rDNA_shifted,
#'   num_samples,
#'   species = "mm39",
#'   min_cov = 10
#' )
#'
#' # Calculate mean methylation for the 18S region
#' mean_meth <- calculate_aggregate_meth(
#'   region_res,
#'   region = "18S",
#'   conditions = samples,
#'   med = FALSE
#' )
#'
#' # Calculate median methylation for the 18S region
#' median_meth <- calculate_aggregate_meth(
#'   region_res,
#'   region = "18S",
#'   conditions = samples,
#'   med = TRUE
#' )
#'
#' # Inspect the results
#' print(mean_meth)
#' print(median_meth)
#'
#' @export
calculate_aggregate_meth <- function(
    region_res,
    region,
    conditions,
    med = FALSE)
{
  # Basic input checks
  if (!is.list(region_res) ||
      is.data.frame(region_res) ||
      length(region_res) == 0L) {
    stop("'region_res' must be a non-empty list.", call. = FALSE)
  }
  
  if (!is.character(region) ||
      length(region) != 1L ||
      is.na(region) ||
      !nzchar(region)) {
    stop("'region' must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!is.list(conditions) ||
      is.data.frame(conditions) ||
      length(conditions) == 0L ||
      is.null(names(conditions)) ||
      anyNA(names(conditions)) ||
      any(!nzchar(names(conditions)))) {
    stop("'conditions' must be a non-empty named list.",
         call. = FALSE)
  }
  
  if (!is.logical(med) ||
      length(med) != 1L ||
      is.na(med)) {
    stop("'med' must be a single TRUE or FALSE.", call. = FALSE)
  }
  
  # Preallocate the results
  n_samples <- sum(lengths(conditions))
  
  result_df <- S4Vectors::DataFrame(
    Sample_ID = character(n_samples),
    Condition = character(n_samples),
    aggregate_meth = numeric(n_samples)
  )
  
  i <- 0L
  
  for (condition in names(conditions)) {
    for (sam in conditions[[condition]]) {
      i <- i + 1L
      
      region_data <- region_res[[sam]][[region]]
      meth <- region_data$numCs / region_data$coverage
      
      if (med) {
        aggregate_meth <- stats::median(meth, na.rm = TRUE)
      } else {
        aggregate_meth <- mean(meth, na.rm = TRUE)
      }
      
      result_df$Sample_ID[i] <- sam
      result_df$Condition[i] <- condition
      result_df$aggregate_meth[i] <- aggregate_meth
    }
  }
  
  return(result_df)
}