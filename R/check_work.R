#' Check a hand calculation against a derivation
#'
#' Compare a value you computed by hand against the final estimate, or against
#' a named intermediate step. When the value does not match, `check_work()`
#' searches every step of the derivation for one that does -- so instead of
#' "wrong", the student is told *where* they stopped or slipped.
#'
#' @param x An `epibyhand_derivation`.
#' @param value Numeric. The value you computed by hand.
#' @param step Optional. Which step to compare against: a step number, or the
#'   symbol of a step such as `"SE"` or `"R1"`. Defaults to the final estimate.
#' @param tol Relative tolerance for the comparison. The default is loose
#'   enough to accept a value rounded to two decimal places.
#'
#' @return Invisibly, `TRUE` if the value matched the requested target.
#'
#' @examples
#' d <- odds_ratio(epi2x2(36, 14, 30, 25))
#' check_work(d, 2.14)
#' check_work(d, 2.571)
#' @export
check_work <- function(x, value, step = NULL, tol = 0.005) {
  stopifnot(inherits(x, "epibyhand_derivation"))
  if (!is.numeric(value) || length(value) != 1L) {
    stop("`value` must be a single number.", call. = FALSE)
  }

  tbl <- steps_table(x)
  close_to <- function(target) {
    isTRUE(abs(value - target) <= tol * max(abs(target), 1e-8))
  }

  if (is.null(step)) {
    target <- x$estimate
    target_name <- paste0("the final estimate (", x$symbol, ")")
  } else if (is.numeric(step)) {
    if (step < 1 || step > nrow(tbl)) {
      stop("This derivation has ", nrow(tbl), " steps.", call. = FALSE)
    }
    target <- tbl$result[step]
    target_name <- paste0("step ", step, " (", tbl$label[step], ")")
  } else {
    i <- match(step, tbl$symbol)
    if (is.na(i)) {
      stop("No step named \"", step, "\". Available: ",
           paste(stats::na.omit(tbl$symbol), collapse = ", "), ".", call. = FALSE)
    }
    target <- tbl$result[i]
    target_name <- paste0("step ", i, " (", tbl$label[i], ")")
  }

  if (close_to(target)) {
    cat("Correct. ", fmt(value), " matches ", target_name, ".\n", sep = "")
    return(invisible(TRUE))
  }

  cat("Not a match. You gave ", fmt(value), "; ", target_name,
      " is ", fmt(target), ".\n", sep = "")

  hits <- which(vapply(tbl$result, close_to, logical(1L)))
  if (length(hits)) {
    i <- hits[1L]
    cat("\nYour value does match step ", i, ": ", tbl$label[i], ".\n", sep = "")
    cat("  ", tbl$symbol[i], " = ", tbl$formula[i], "\n", sep = "")
    if (i < nrow(tbl)) {
      cat("You may have stopped early. The next step is: ",
          tbl$label[i + 1L], ".\n", sep = "")
    }
  } else {
    cat("\nIt does not match any intermediate step either, so the slip is",
        "probably\narithmetic rather than a wrong stopping point.",
        "Print the derivation to\ncompare line by line.\n")
  }

  invisible(FALSE)
}
