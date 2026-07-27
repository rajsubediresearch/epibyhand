#' epibyhand: Worked Derivations for Classical Epidemiological Measures
#'
#' Every function returns the answer together with the reasoning that
#' produced it: each intermediate quantity, the formula in symbols, and the
#' formula with the observed numbers substituted in. The result is a
#' derivation object that can be printed at three levels of detail, converted
#' to a data frame with [steps_table()], or compared against a student's own
#' arithmetic with [check_work()].
#'
#' @section Scope:
#' The package covers methods a student can compute by hand on paper. That
#' boundary is deliberate and load-bearing: it is what keeps the package small
#' enough to stay correct without continuous maintenance, and it is why there
#' is no regression modelling here. Anything requiring an iterative fit is out
#' of scope by definition rather than by omission.
#'
#' @section Options:
#' \describe{
#'   \item{`epibyhand.verbose`}{Detail level for printing: `0` result only,
#'     `1` adds symbolic formulas, `2` (default) full worked solution.}
#'   \item{`epibyhand.digits`}{Decimal places, default `4`.}
#' }
#'
#' @keywords internal
"_PACKAGE"
