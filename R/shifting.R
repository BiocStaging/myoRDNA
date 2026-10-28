#' Shift rDNA Methylation Data
#'
#' This function shifts the start and end positions of rDNA methylation data
#' by a calculated shift value, ensuring all coordinates remain non-negative.
#' The shift is computed as the difference between the maximum position in the 
#' data and a user-specified cutoff. Invalid coordinates (negative after 
#' shifting)are filtered out.
#'
#' @param meth_rdna_raw A list of `methylRaw` objects (from `methylKit`),
#'   where each object contains methylation data with `start` and `end` 
#'   positions.
#' @param cutoff A numeric value used to compute the shift.
#'   The shift is calculated as `max_position - cutoff`.
#' @param treatment A character vector or factor specifying the treatment 
#' group for each sample. Used to label the output `methylRawList`.
#'
#' @return A `methylRawList` object containing the shifted methylation data,
#'   with invalid coordinates (negative after shifting) removed.
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
#' num_samples <- length(treatment)
#'
#' # Define sample groups
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
#' # Shift rDNA genomic coordinates
#' meth_rDNA_shifted <- shifting(
#'   meth_rdna_raw = meth_rDNA,
#'   cutoff = 35000,
#'   treatment = treatment
#' )
#' 
#' head(meth_rDNA_shifted)
#'
#' @export
shifting <- function(
    meth_rdna_raw,
    cutoff,
    treatment)
{
  # Check methylation objects
  if (!(is.list(meth_rdna_raw) ||
        methods::is(meth_rdna_raw, "methylRawList")) ||
      length(meth_rdna_raw) == 0L) {
    stop("'meth_rdna_raw' must be a non-empty list or methylRawList.",
         call. = FALSE)
  }
  
  if (!all(vapply(meth_rdna_raw, function(x) {
    methods::is(x, "methylRaw")
  }, logical(1)))) {
    stop("Each sample must be a methylRaw object.", call. = FALSE)
  }
  
  # Check cutoff
  if (!is.numeric(cutoff) ||
      length(cutoff) != 1L ||
      !is.finite(cutoff) ||
      cutoff < 1 ||
      cutoff != floor(cutoff)) {
    stop("'cutoff' must be a positive whole number.", call. = FALSE)
  }
  
  # Check treatment codes
  if (!is.numeric(treatment) ||
      !is.null(dim(treatment)) ||
      length(treatment) != length(meth_rdna_raw) ||
      any(!is.finite(treatment)) ||
      any(treatment < 0) ||
      any(treatment != floor(treatment))) {
    stop("'treatment' must contain one non-negative whole-number 
         code per sample.",
         call. = FALSE)
  }
  
  # Extract positions
  all_positions <- unlist(lapply(
    meth_rdna_raw,
    function(x) methylKit::getData(x)$start
  ))
  
  if (length(all_positions) == 0L ||
      any(!is.finite(all_positions))) {
    stop("Start positions must be present and finite.", call. = FALSE)
  }
  
  # Calculate and apply the shift
  max_position <- max(all_positions)
  shift <- max_position - cutoff
  
  meth_rdna_shifted <- lapply(
    meth_rdna_raw,
    .shift_methylRaw,
    shift = shift
  )
  
  meth_rdna_shifted <- methylKit::methylRawList(
    meth_rdna_shifted,
    treatment = treatment
  )
  
  return(meth_rdna_shifted)
}