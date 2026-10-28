#' Sort Methylation Data by rDNA Regions
#'
#' @description
#' Extracts and organizes methylation data for specified rDNA regions across 
#' samples.
#' Supports predefined species (`mm39`, `rno7`). Validates sample bounds and 
#' required columns.
#'
#' @param rDNA_meth_object A list of `methylRaw` objects (from `methylKit`) or 
#' data frames.
#' @param num_samples Integer. Number of samples to process (must be >= 1 and 
#' <= `length(rDNA_meth_object)`).
#' @param species Character. Species reference for rDNA regions: `"mm39"` 
#' (default), `"rno7"`, or `"hg38"` (unsupported).
#' @param min_cov Numeric. Minimum coverage threshold to filter data 
#' (default: `10`).
#'
#' @return A named list of lists, where each sublist contains:
#'   - `numCs`: Number of methylated Cs in the region.
#'   - `numTs`: Number of unmethylated Ts in the region.
#'   - `coverage`: Coverage depth for the region.
#'   - `position`: Start positions of the region.
#'   Structured as `sample_id -> region -> data`.
#'   
#'   
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
#' # Construct a methylRawList object
#' loaded_data <- methods::new(
#'   "methylRawList",
#'   loaded_data,
#'   treatment = treatment
#' )
#'
#' # Extract methylation data from the rDNA chromosome
#' meth_rDNA <- get_rDNA_data(
#'   loaded_data,
#'   num_samples,
#'   rDNA_chr_name = "chrRDNAm"
#' )
#'
#' # Separate rDNA methylation data into annotated regions
#' region_res <- sort_by_rDNA_region(
#'   meth_rDNA,
#'   num_samples,
#'   species = "mm39",
#'   min_cov = 10
#' )
#'
#' # Inspect the regions available for the first sample
#' names(region_res[[1]])
#'
#' @note
#' - For `species = "hg38"`, the function stops with a warning 
#' (unsupported).
#' - Missing columns in input data will cause an error.
#'
#' @export
sort_by_rDNA_region <- function(
    rDNA_meth_object,
    num_samples,
    species = "mm39",
    min_cov = 10)
{
  # Basic input checks
  if (!(is.list(rDNA_meth_object) ||
        methods::is(rDNA_meth_object, "methylRawList")) ||
      is.data.frame(rDNA_meth_object) ||
      length(rDNA_meth_object) == 0L) {
    stop("'rDNA_meth_object' must be a non-empty list or methylRawList.",
         call. = FALSE)
  }
  
  if (!is.numeric(num_samples) ||
      length(num_samples) != 1L ||
      !is.finite(num_samples) ||
      num_samples < 1 ||
      num_samples != floor(num_samples) ||
      num_samples > length(rDNA_meth_object)) {
    stop("'num_samples' must be a positive whole number no greater", 
         "than the number of samples.",
         call. = FALSE)
  }
  
  if (!is.character(species) ||
      length(species) != 1L ||
      is.na(species) ||
      !species %in% c("mm39", "rno7", "hg38")) {
    stop("'species' must be 'mm39', 'rno7', or 'hg38'.",
         call. = FALSE)
  }
  
  if (species == "hg38") {
    stop("'hg38' is not currently supported.",
         call. = FALSE)
  }
  
  if (!is.numeric(min_cov) ||
      length(min_cov) != 1L ||
      !is.finite(min_cov) ||
      min_cov < 0) {
    stop("'min_cov' must be a single finite, non-negative number.",
         call. = FALSE)
  }
  
  # Define regions
  message("Species selected: ", species)
  
  if (species == "mm39") {
    rdna_regions <- data.frame(
      region = c(
        "18S", "5.8S", "28S", "45S", "5_ETS",
        "ITS_1", "ITS_2", "3_ETS", "r"
      ),
      start = c(
        4008, 6878, 8123, 1, 1,
        5878, 7035, 12852, 13404
      ),
      end = c(
        5877, 7034, 12852, 13403, 4007,
        6877, 8122, 13403, 30000
      )
    )
  } else {
    rdna_regions <- data.frame(
      region = c("18S", "5.8S", "28S"),
      start = c(1, 1873, 2030),
      end = c(1872, 2029, 6836)
    )
  }
  
  region_results <- list()
  
  region_results <- vector("list", num_samples)
  
  for (i in seq_len(num_samples)) {
    sample <- rDNA_meth_object[[i]]
    
    sample_data <- methylKit::getData(sample)
    sample_id <- methylKit::getSampleID(sample)
  
    # Filter coverage and exclude missing positions
    sample_data <- sample_data[
      !is.na(sample_data$coverage) &
        sample_data$coverage >= min_cov &
        !is.na(sample_data$start),
      ,
      drop = FALSE
    ]
    
    region_results[[i]] <- stats::setNames(
      vector("list", nrow(rdna_regions)),
      rdna_regions$region
    )
    
    names(region_results)[i] <- sample_id
    
    for (r in seq_len(nrow(rdna_regions))) {
      region_name <- rdna_regions$region[r]
      
      region_data <- sample_data[
        sample_data$start >= rdna_regions$start[r] &
          sample_data$start <= rdna_regions$end[r],
        ,
        drop = FALSE
      ]
      
      region_results[[i]][[region_name]] <- list(
        numCs = region_data$numCs,
        numTs = region_data$numTs,
        coverage = region_data$coverage,
        position = region_data$start
      )
    }
  }
  
  return(region_results)
}