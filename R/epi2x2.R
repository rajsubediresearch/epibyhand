## The 2x2 table -----------------------------------------------------------
##
## Everything in this package that takes tabulated data takes an `epi2x2`.
## Standard epidemiological orientation, matching Rothman/Greenland:
##
##                Case      Non-case
##   Exposed        a          b
##   Unexposed      c          d

#' Build a 2 x 2 table
#'
#' Accepts either four cell counts or a 2 x 2 matrix/table. Orientation is
#' exposure in rows, outcome in columns, with the index level first in each
#' -- the layout used in standard epidemiology texts.
#'
#' @param a Exposed cases, or a 2 x 2 matrix or table.
#' @param b Exposed non-cases.
#' @param c Unexposed cases.
#' @param d Unexposed non-cases.
#' @param exposure Length-2 character vector labelling the rows.
#' @param outcome Length-2 character vector labelling the columns.
#'
#' @return An object of class `epi2x2`.
#'
#' @examples
#' epi2x2(36, 14, 30, 25)
#' epi2x2(matrix(c(36, 14, 30, 25), nrow = 2, byrow = TRUE))
#' @export
epi2x2 <- function(a, b = NULL, c = NULL, d = NULL,
                   exposure = c("Exposed", "Unexposed"),
                   outcome = c("Case", "Non-case")) {

  if (inherits(a, "epi2x2")) return(a)

  if (is.matrix(a) || is.table(a)) {
    m <- as.matrix(a)
    if (!identical(dim(m), c(2L, 2L))) {
      stop("A matrix input must be 2 x 2, not ",
           paste(dim(m), collapse = " x "), ".", call. = FALSE)
    }
    if (!is.null(dimnames(m))) {
      if (!is.null(rownames(m))) exposure <- rownames(m)
      if (!is.null(colnames(m))) outcome <- colnames(m)
    }
    counts <- c(m[1L, 1L], m[1L, 2L], m[2L, 1L], m[2L, 2L])
  } else {
    if (is.null(b) || is.null(c) || is.null(d)) {
      stop("Supply either four cell counts (a, b, c, d) or a 2 x 2 matrix.",
           call. = FALSE)
    }
    counts <- c(a, b, c, d)
  }

  if (!is.numeric(counts) || anyNA(counts)) {
    stop("Cell counts must be non-missing numbers.", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("Cell counts cannot be negative.", call. = FALSE)
  }

  structure(
    list(
      a = counts[1L], b = counts[2L], c = counts[3L], d = counts[4L],
      exposure = as.character(exposure),
      outcome = as.character(outcome)
    ),
    class = "epi2x2"
  )
}

#' Print a 2 x 2 table with margins
#'
#' @param x An `epi2x2`.
#' @param indent Character prefix for each line.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.epi2x2 <- function(x, indent = "", ...) {
  m <- matrix(
    c(x$a, x$b, x$a + x$b,
      x$c, x$d, x$c + x$d,
      x$a + x$c, x$b + x$d, x$a + x$b + x$c + x$d),
    nrow = 3L, byrow = TRUE
  )
  rn <- c(x$exposure, "Total")
  cn <- c(x$outcome, "Total")

  rw <- max(nchar(rn))
  cw <- pmax(nchar(cn), apply(format(m, trim = TRUE), 2L, function(z) max(nchar(z))))

  pad <- function(s, w, left = FALSE) {
    s <- as.character(s)
    sp <- strrep(" ", max(w - nchar(s), 0L))
    if (left) paste0(s, sp) else paste0(sp, s)
  }

  cat(indent, pad("", rw, left = TRUE), "  ",
      paste(mapply(pad, cn, cw), collapse = "  "), "\n", sep = "")
  for (i in 1:3) {
    cat(indent, pad(rn[i], rw, left = TRUE), "  ",
        paste(mapply(function(v, w) pad(format(v, trim = TRUE), w), m[i, ], cw),
              collapse = "  "),
        "\n", sep = "")
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
check_zero_cells <- function(x) {
  cells <- c(a = x$a, b = x$b, c = x$c, d = x$d)
  zero <- names(cells)[cells == 0]
  if (!length(zero)) return(character())
  paste0(
    "Cell ", paste(zero, collapse = " and "), " is zero, so the ratio and ",
    "its standard error are undefined. epibyhand does not silently apply a ",
    "continuity correction. If you want one, add 0.5 to every cell yourself ",
    "and note that you did so -- the corrected estimate is biased toward the ",
    "null and its confidence interval is approximate."
  )
}
