## Stratified 2 x 2 data ---------------------------------------------------
##
## A set of 2 x 2 tables, one per level of a third variable. This is the
## input to every stratified method in the package.

#' Build a set of stratified 2 x 2 tables
#'
#' Accepts a list of [epi2x2()] objects, a list of four-element count vectors,
#' or a 2 x 2 x K array. Stratum labels are taken from the names of the list,
#' the third dimnames of the array, or generated as `"Stratum 1"`, and so on.
#'
#' @param ... Strata: `epi2x2` objects, length-4 numeric vectors of counts, or
#'   a single list or 2 x 2 x K array containing them.
#' @param labels Optional character vector of stratum names.
#' @param exposure Length-2 character vector labelling the rows.
#' @param outcome Length-2 character vector labelling the columns.
#'
#' @return An object of class `epi_strata`.
#'
#' @examples
#' epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25), labels = c("Men", "Women"))
#' @export
epi_strata <- function(..., labels = NULL,
                       exposure = c("Exposed", "Unexposed"),
                       outcome = c("Case", "Non-case")) {

  dots <- list(...)
  if (inherits(dots[[1L]], "epi_strata")) return(dots[[1L]])

  if (length(dots) == 1L && is.array(dots[[1L]]) && length(dim(dots[[1L]])) == 3L) {
    arr <- dots[[1L]]
    if (!identical(dim(arr)[1:2], c(2L, 2L))) {
      stop("An array input must be 2 x 2 x K.", call. = FALSE)
    }
    if (is.null(labels) && !is.null(dimnames(arr)[[3L]])) {
      labels <- dimnames(arr)[[3L]]
    }
    dots <- lapply(seq_len(dim(arr)[3L]), function(k) arr[, , k])
  } else if (length(dots) == 1L && is.list(dots[[1L]]) &&
             !inherits(dots[[1L]], "epi2x2")) {
    if (is.null(labels)) labels <- names(dots[[1L]])
    dots <- unname(dots[[1L]])
  }

  if (length(dots) < 2L) {
    stop("A stratified analysis needs at least two strata.", call. = FALSE)
  }

  tables <- lapply(dots, function(c1) {
    if (inherits(c1, "epi2x2")) return(c1)
    if (is.matrix(c1) || is.table(c1)) {
      return(epi2x2(c1, exposure = exposure, outcome = outcome))
    }
    if (is.numeric(c1) && length(c1) == 4L) {
      return(epi2x2(c1[1L], c1[2L], c1[3L], c1[4L],
                    exposure = exposure, outcome = outcome))
    }
    stop("Each stratum must be an epi2x2, a 2 x 2 matrix, or four counts.",
         call. = FALSE)
  })

  if (is.null(labels)) labels <- paste("Stratum", seq_along(tables))
  if (length(labels) != length(tables)) {
    stop("`labels` must have one entry per stratum.", call. = FALSE)
  }

  structure(
    list(tables = tables, labels = as.character(labels)),
    class = "epi_strata"
  )
}

#' Collapse strata into a single 2 x 2 table
#'
#' Adds the strata cell by cell, discarding the stratifying variable. The
#' result is the *crude* table -- the one you would have had if you had never
#' stratified. Comparing a crude estimate with its adjusted counterpart is
#' how confounding is made visible.
#'
#' @param x An `epi_strata`.
#' @return An `epi2x2`.
#' @examples
#' collapse_strata(epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25)))
#' @export
collapse_strata <- function(x) {
  stopifnot(inherits(x, "epi_strata"))
  epi2x2(
    a = sum(vapply(x$tables, function(t) t$a, numeric(1L))),
    b = sum(vapply(x$tables, function(t) t$b, numeric(1L))),
    c = sum(vapply(x$tables, function(t) t$c, numeric(1L))),
    d = sum(vapply(x$tables, function(t) t$d, numeric(1L))),
    exposure = x$tables[[1L]]$exposure,
    outcome = x$tables[[1L]]$outcome
  )
}

#' Print stratified tables
#'
#' @param x An `epi_strata`.
#' @param indent Character prefix for each line.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.epi_strata <- function(x, indent = "", ...) {
  for (k in seq_along(x$tables)) {
    cat(indent, x$labels[k], "\n", sep = "")
    print(x$tables[[k]], indent = paste0(indent, "  "))
    if (k < length(x$tables)) cat("\n")
  }
  invisible(x)
}
