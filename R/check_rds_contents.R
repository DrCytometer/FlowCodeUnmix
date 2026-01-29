# check_rds_contents.r

#' @title Check RDS Contents
#'
#' @description
#' Checks the RDS file for expected elements.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom cowplot plot_grid
#'
#' @param obj The loaded RDS object.
#' @param required List of required elements.
#' @param allow.extra Logical, whether to allow additional, unexpected elements.
#' Default is `TRUE`.
#'
#' @return Returns errors if any found; otherwise, returns `TRUE`.

check.rds.contents <- function(
    obj,
    required = list(),
    allow.extra = TRUE
  ) {

  # 1. Check all required names exist
  missing <- setdiff(names(required), names(obj))
  if (length(missing) > 0) {
    stop(
      sprintf("Missing required elements: %s",
              paste(missing, collapse = ", "))
    )
  }

  # 2. Loop through each expected element and check class/dim
  for (nm in names(required)) {
    expected <- required[[nm]]
    x <- obj[[nm]]

    # expected$class: character vector of allowed classes
    if (!is.null(expected$class)) {
      if (!inherits(x, expected$class)) {
        stop(sprintf(
          "Element '%s' must be of class %s, but is %s",
          nm,
          paste(expected$class, collapse = " OR "),
          paste(class(x), collapse = "/")
        ))
      }
    }


    # expected$nrow_min: enforce minimum number of rows
    if (!is.null(expected$nrow_min)) {
      if (!is.matrix(x) && !is.data.frame(x)) {
        stop(sprintf(
          "Element '%s' must be a matrix/data.frame to check row count.",
          nm
        ))
      }
      if (nrow(x) < expected$nrow_min) {
        stop(sprintf(
          "Element '%s' must have at least %d rows, but has %d.",
          nm, expected$nrow_min, nrow(x)
        ))
      }
    }

    # expected$length_min: e.g.,  vectors
    if (!is.null(expected$length_min)) {
      if (length(x) < expected$length_min) {
        stop(sprintf(
          "Element '%s' must have length of at least %d, but has %d.",
          nm, expected$length_min, length(x)
        ))
      }
    }
  }

  # Optionally reject unexpected extra elements
  if (!allow.extra) {
    extra <- setdiff(names(obj), names(required))
    if (length(extra) > 0) {
      stop(sprintf(
        "Unexpected elements present: %s",
        paste(extra, collapse = ", ")
      ))
    }
  }

  invisible(TRUE)
}
