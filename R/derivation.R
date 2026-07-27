## Core architecture -----------------------------------------------------
##
## Every user-facing function in epibyhand returns one object class:
## `epibyhand_derivation`. It carries the result AND the reasoning that
## produced it. All display logic lives in the print/format methods here,
## so adding a new measure means writing arithmetic and steps -- never
## writing display code.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Format a number for display
#'
#' @param x Numeric vector.
#' @param digits Number of decimal places.
#' @return Character vector.
#' @keywords internal
#' @noRd
fmt <- function(x, digits = getOption("epibyhand.digits", 4L)) {
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (!is.finite(v)) return(if (v > 0) "Inf" else "-Inf")
    trimws(format(round(v, digits), nsmall = 0L, scientific = FALSE))
  }, character(1L))
}

#' Does this console handle box-drawing characters?
#' @keywords internal
#' @noRd
utf8_ok <- function() {
  isTRUE(l10n_info()[["UTF-8"]]) && !identical(Sys.getenv("EPIBYHAND_ASCII"), "1")
}

#' Draw a section rule
#' @keywords internal
#' @noRd
rule <- function(label = NULL, width = NULL) {
  width <- width %||% min(getOption("width", 80L), 76L)
  dash <- if (utf8_ok()) "\u2500" else "-"
  if (is.null(label)) return(strrep(dash, width))
  lead <- strrep(dash, 2L)
  txt <- paste0(lead, " ", label, " ")
  pad <- max(width - nchar(txt), 0L)
  paste0(txt, strrep(dash, pad))
}

## Steps -----------------------------------------------------------------

#' Create one step of a derivation
#'
#' A step is the atomic unit of a worked solution: a human-readable label,
#' the symbolic formula, the same formula with observed values substituted
#' in, and the numeric result.
#'
#' @param label Short description of what this step computes.
#' @param formula Character. The formula in symbols, e.g. `"a / (a + b)"`.
#' @param substituted Character. The formula with numbers filled in.
#' @param result Numeric. The value this step evaluates to.
#' @param symbol Optional character. The symbol this step defines, e.g. `"R1"`.
#' @param note Optional character. A caveat or explanation attached to the step.
#' @param table Optional data frame. Per-unit working shown as a grid beneath
#'   the formula. Used when a step aggregates over strata and a single
#'   substituted line would be unreadable.
#'
#' @return An object of class `epibyhand_step`.
#' @export
derivation_step <- function(label, formula, substituted, result,
                            symbol = NULL, note = NULL, table = NULL) {
  stopifnot(is.character(label), is.character(formula))
  structure(
    list(
      label = label,
      formula = formula,
      substituted = substituted,
      result = result,
      symbol = symbol,
      note = note,
      table = table
    ),
    class = "epibyhand_step"
  )
}

## Derivation ------------------------------------------------------------

#' Create a derivation object
#'
#' @param method Character. Name of the measure, e.g. `"Odds ratio"`.
#' @param estimate Numeric. The final point estimate.
#' @param symbol Character. Symbol for the estimate, e.g. `"OR"`.
#' @param ci Numeric length-2 vector, or `NULL`.
#' @param conf_level Confidence level used for `ci`.
#' @param data Optional. The input data object, printed above the steps.
#' @param steps List of `epibyhand_step` objects.
#' @param notes Character vector of assumption notes shown after the result.
#'
#' @return An object of class `epibyhand_derivation`.
#' @export
derivation <- function(method, estimate, symbol = NULL, ci = NULL,
                       conf_level = 0.95, data = NULL, steps = list(),
                       notes = character()) {
  structure(
    list(
      method = method,
      estimate = estimate,
      symbol = symbol %||% "estimate",
      ci = ci,
      conf_level = conf_level,
      data = data,
      steps = steps,
      notes = notes
    ),
    class = "epibyhand_derivation"
  )
}

#' Print a worked derivation
#'
#' @param x An `epibyhand_derivation`.
#' @param verbose Detail level. `0` prints the result only, `1` adds the
#'   symbolic formulas, `2` (the default) shows the full worked solution with
#'   numbers substituted in. Set globally with
#'   `options(epibyhand.verbose = 1)`.
#' @param digits Decimal places. Set globally with
#'   `options(epibyhand.digits = 3)`.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.epibyhand_derivation <- function(x, verbose = NULL, digits = NULL, ...) {
  verbose <- verbose %||% getOption("epibyhand.verbose", 2L)
  digits <- digits %||% getOption("epibyhand.digits", 4L)

  cat(rule(x$method), "\n", sep = "")

  if (verbose >= 2L && !is.null(x$data)) {
    cat("\nData\n")
    print(x$data, indent = "  ")
  }

  if (verbose >= 1L && length(x$steps)) {
    cat("\n")
    for (i in seq_along(x$steps)) {
      print_step(x$steps[[i]], index = i, verbose = verbose, digits = digits)
    }
  }

  cat("\nResult\n")
  res <- paste0("  ", x$symbol, " = ", fmt(x$estimate, digits))
  if (!is.null(x$ci)) {
    res <- paste0(
      res, " (", round(x$conf_level * 100), "% CI ",
      fmt(x$ci[1L], digits), " to ", fmt(x$ci[2L], digits), ")"
    )
  }
  cat(res, "\n", sep = "")

  if (length(x$notes)) {
    cat("\nNotes\n")
    for (n in x$notes) cat(strwrap(n, prefix = "  ", width = 76L), sep = "\n")
  }

  invisible(x)
}

#' @keywords internal
#' @noRd
print_step <- function(s, index, verbose = 2L, digits = 4L) {
  head <- paste0("Step ", index, "  ", s$label)
  cat(head, "\n", sep = "")
  lhs <- s$symbol %||% ""
  pad <- strrep(" ", 8L)

  has_result <- !is.null(s$result) && !all(is.na(s$result))

  if (nzchar(lhs)) {
    cat(pad, lhs, " = ", s$formula, "\n", sep = "")
    align <- strrep(" ", nchar(lhs))
    if (verbose >= 2L && !is.null(s$substituted)) {
      cat(pad, align, " = ", s$substituted, "\n", sep = "")
    }
    if (has_result) {
      cat(pad, align, " = ", fmt(s$result, digits), "\n", sep = "")
    }
  } else {
    cat(pad, s$formula, "\n", sep = "")
    if (verbose >= 2L && !is.null(s$substituted)) {
      cat(pad, "= ", s$substituted, "\n", sep = "")
    }
    if (has_result) {
      cat(pad, "= ", fmt(s$result, digits), "\n", sep = "")
    }
  }

  if (verbose >= 2L && !is.null(s$table)) {
    cat("\n")
    print_step_table(s$table, indent = pad, digits = digits)
  }

  if (!is.null(s$note)) {
    cat(strwrap(s$note, prefix = paste0(pad, "# "), width = 76L), sep = "\n")
  }
  cat("\n")
}

#' Render a per-stratum working table beneath a step
#' @keywords internal
#' @noRd
print_step_table <- function(df, indent = "        ", digits = 4L) {
  body <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) fmt(col, digits) else as.character(col)
  }), stringsAsFactors = FALSE)
  names(body) <- names(df)

  widths <- vapply(seq_along(body), function(j) {
    max(nchar(c(names(body)[j], body[[j]])))
  }, integer(1L))

  padcell <- function(s, w) paste0(strrep(" ", max(w - nchar(s), 0L)), s)

  cat(indent, paste(mapply(padcell, names(body), widths), collapse = "  "),
      "\n", sep = "")
  cat(indent, paste(mapply(function(w) strrep("-", w), widths), collapse = "  "),
      "\n", sep = "")
  for (i in seq_len(nrow(body))) {
    cat(indent,
        paste(mapply(padcell, unlist(body[i, ], use.names = FALSE), widths),
              collapse = "  "),
        "\n", sep = "")
  }
  cat("\n")
}

## Extractors ------------------------------------------------------------

#' Extract the point estimate from a derivation
#' @param x An `epibyhand_derivation`.
#' @return Numeric scalar.
#' @export
estimate <- function(x) {
  stopifnot(inherits(x, "epibyhand_derivation"))
  x$estimate
}

#' Extract the confidence interval from a derivation
#' @param object An `epibyhand_derivation`.
#' @param parm Unused, for S3 consistency.
#' @param level Unused; the level is fixed when the derivation is built.
#' @param ... Unused.
#' @return Numeric vector of length 2.
#' @export
confint.epibyhand_derivation <- function(object, parm, level, ...) {
  object$ci
}

#' Return the steps of a derivation as a data frame
#'
#' Useful for building answer keys, or for rendering the derivation in a
#' format this package does not provide.
#'
#' @param x An `epibyhand_derivation`.
#' @return A data frame with one row per step.
#' @export
steps_table <- function(x) {
  stopifnot(inherits(x, "epibyhand_derivation"))
  if (!length(x$steps)) {
    return(data.frame(
      step = integer(), symbol = character(), label = character(),
      formula = character(), substituted = character(), result = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    step = seq_along(x$steps),
    symbol = vapply(x$steps, function(s) s$symbol %||% NA_character_, character(1L)),
    label = vapply(x$steps, function(s) s$label, character(1L)),
    formula = vapply(x$steps, function(s) s$formula, character(1L)),
    substituted = vapply(x$steps, function(s) s$substituted %||% NA_character_, character(1L)),
    result = vapply(x$steps, function(s) as.numeric(s$result), numeric(1L)),
    stringsAsFactors = FALSE
  )
}
