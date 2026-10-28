#' Shift genomic coordinates in a methylRaw object
#'
#' Subtracts a fixed offset from the start and end coordinates.
#' Rows with missing or negative shifted coordinates are removed.
#' Coordinates equal to zero are retained. Sample metadata is preserved.
#'
#' @param x A methylKit \code{methylRaw} object containing
#'   \code{start} and \code{end} coordinates.
#' @param shift A single finite numeric offset to subtract from both
#'   coordinates. Use a whole number to preserve integer coordinates.
#'
#' @return A \code{methylRaw} object containing the retained rows
#'   with shifted coordinates and the original sample metadata.
#'
#' @keywords internal
.shift_methylRaw <- function(x, shift) {
  df <- methylKit::getData(x)
  
  shifted_start <- df$start - shift
  shifted_end <- df$end - shift
  
  # Preserve your existing coordinate filter
  keep <- !is.na(shifted_start) &
    !is.na(shifted_end) &
    shifted_start >= 0 &
    shifted_end >= 0
  
  # Preserve the sample metadata
  shifted_x <- x[keep, ]
  shifted_x$start <- shifted_start[keep]
  shifted_x$end <- shifted_end[keep]
  
  methods::validObject(shifted_x)
  
  return(shifted_x)
}