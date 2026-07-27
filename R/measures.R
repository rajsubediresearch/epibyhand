## Measures of association from a 2 x 2 table ------------------------------

#' @keywords internal
#' @noRd
zcrit <- function(conf_level) stats::qnorm(1 - (1 - conf_level) / 2)

#' @keywords internal
#' @noRd
risk_steps <- function(x, digits) {
  n1 <- x$a + x$b
  n0 <- x$c + x$d
  r1 <- x$a / n1
  r0 <- x$c / n0
  list(
    n1 = n1, n0 = n0, r1 = r1, r0 = r0,
    steps = list(
      derivation_step(
        label = paste0("Risk of ", x$outcome[1L], " among the ", tolower(x$exposure[1L])),
        symbol = "R1",
        formula = "a / (a + b)",
        substituted = paste0(fmt(x$a, digits), " / (", fmt(x$a, digits), " + ",
                             fmt(x$b, digits), ")"),
        result = r1
      ),
      derivation_step(
        label = paste0("Risk of ", x$outcome[1L], " among the ", tolower(x$exposure[2L])),
        symbol = "R0",
        formula = "c / (c + d)",
        substituted = paste0(fmt(x$c, digits), " / (", fmt(x$c, digits), " + ",
                             fmt(x$d, digits), ")"),
        result = r0
      )
    )
  )
}

#' Risk ratio with worked derivation
#'
#' @param x An `epi2x2`, or the first of four cell counts.
#' @param ... Further cell counts if `x` is a bare count.
#' @param conf_level Confidence level for the interval.
#'
#' @return An `epibyhand_derivation`.
#'
#' @examples
#' risk_ratio(epi2x2(36, 14, 30, 25))
#' @export
risk_ratio <- function(x, ..., conf_level = 0.95) {
  x <- epi2x2(x, ...)
  digits <- getOption("epibyhand.digits", 4L)
  z <- zcrit(conf_level)

  rs <- risk_steps(x, digits)
  rr <- rs$r1 / rs$r0
  se <- sqrt(1 / x$a - 1 / rs$n1 + 1 / x$c - 1 / rs$n0)
  lo <- exp(log(rr) - z * se)
  hi <- exp(log(rr) + z * se)

  steps <- c(rs$steps, list(
    derivation_step(
      label = "Risk ratio",
      symbol = "RR",
      formula = "R1 / R0",
      substituted = paste0(fmt(rs$r1, digits), " / ", fmt(rs$r0, digits)),
      result = rr
    ),
    derivation_step(
      label = "Standard error of log(RR)",
      symbol = "SE",
      formula = "sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))",
      substituted = paste0("sqrt(1/", fmt(x$a, digits), " - 1/", fmt(rs$n1, digits),
                           " + 1/", fmt(x$c, digits), " - 1/", fmt(rs$n0, digits), ")"),
      result = se,
      note = "The interval is built on the log scale because RR is a ratio: its sampling distribution is skewed, but log(RR) is roughly normal."
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, lower limit"),
      symbol = "lower",
      formula = "exp(log(RR) - z * SE)",
      substituted = paste0("exp(", fmt(log(rr), digits), " - ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = lo
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, upper limit"),
      symbol = "upper",
      formula = "exp(log(RR) + z * SE)",
      substituted = paste0("exp(", fmt(log(rr), digits), " + ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = hi
    )
  ))

  derivation(
    method = "Risk ratio",
    estimate = rr, symbol = "RR",
    ci = c(lo, hi), conf_level = conf_level,
    data = x, steps = steps,
    notes = c(
      check_zero_cells(x),
      "A risk ratio requires that everyone was followed for the same period. If follow-up time varied, use a rate ratio on person-time instead."
    )
  )
}

#' Odds ratio with worked derivation
#'
#' @inheritParams risk_ratio
#' @return An `epibyhand_derivation`.
#'
#' @examples
#' odds_ratio(epi2x2(36, 14, 30, 25))
#' @export
odds_ratio <- function(x, ..., conf_level = 0.95) {
  x <- epi2x2(x, ...)
  digits <- getOption("epibyhand.digits", 4L)
  z <- zcrit(conf_level)

  odds1 <- x$a / x$b
  odds0 <- x$c / x$d
  or <- odds1 / odds0
  se <- sqrt(1 / x$a + 1 / x$b + 1 / x$c + 1 / x$d)
  lo <- exp(log(or) - z * se)
  hi <- exp(log(or) + z * se)

  r1 <- x$a / (x$a + x$b)
  rr <- r1 / (x$c / (x$c + x$d))

  steps <- list(
    derivation_step(
      label = paste0("Odds of ", x$outcome[1L], " among the ", tolower(x$exposure[1L])),
      symbol = "odds1", formula = "a / b",
      substituted = paste0(fmt(x$a, digits), " / ", fmt(x$b, digits)),
      result = odds1
    ),
    derivation_step(
      label = paste0("Odds of ", x$outcome[1L], " among the ", tolower(x$exposure[2L])),
      symbol = "odds0", formula = "c / d",
      substituted = paste0(fmt(x$c, digits), " / ", fmt(x$d, digits)),
      result = odds0
    ),
    derivation_step(
      label = "Odds ratio",
      symbol = "OR", formula = "odds1 / odds0",
      substituted = paste0(fmt(odds1, digits), " / ", fmt(odds0, digits)),
      result = or,
      note = "Equivalently OR = ad / bc, which is why the odds ratio is the same whether you condition on exposure or on outcome. That symmetry is what makes it usable in case-control studies."
    ),
    derivation_step(
      label = "Standard error of log(OR), Woolf's method",
      symbol = "SE", formula = "sqrt(1/a + 1/b + 1/c + 1/d)",
      substituted = paste0("sqrt(1/", fmt(x$a, digits), " + 1/", fmt(x$b, digits),
                           " + 1/", fmt(x$c, digits), " + 1/", fmt(x$d, digits), ")"),
      result = se
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, lower limit"),
      symbol = "lower", formula = "exp(log(OR) - z * SE)",
      substituted = paste0("exp(", fmt(log(or), digits), " - ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = lo
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, upper limit"),
      symbol = "upper", formula = "exp(log(OR) + z * SE)",
      substituted = paste0("exp(", fmt(log(or), digits), " + ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = hi
    )
  )

  rare <- paste0(
    "These data are tabulated as a cohort, where risk is estimable. The risk ratio here is ",
    fmt(rr, digits), " against an odds ratio of ", fmt(or, digits),
    " -- the odds ratio is the more extreme of the two, and always will be. ",
    "The two converge only when the outcome is rare; baseline risk here is ",
    fmt(x$c / (x$c + x$d), digits), "."
  )

  derivation(
    method = "Odds ratio",
    estimate = or, symbol = "OR",
    ci = c(lo, hi), conf_level = conf_level,
    data = x, steps = steps,
    notes = c(check_zero_cells(x), rare)
  )
}

#' Risk difference with worked derivation
#'
#' @inheritParams risk_ratio
#' @return An `epibyhand_derivation`.
#'
#' @examples
#' risk_difference(epi2x2(36, 14, 30, 25))
#' @export
risk_difference <- function(x, ..., conf_level = 0.95) {
  x <- epi2x2(x, ...)
  digits <- getOption("epibyhand.digits", 4L)
  z <- zcrit(conf_level)

  rs <- risk_steps(x, digits)
  rd <- rs$r1 - rs$r0
  se <- sqrt(rs$r1 * (1 - rs$r1) / rs$n1 + rs$r0 * (1 - rs$r0) / rs$n0)
  lo <- rd - z * se
  hi <- rd + z * se

  steps <- c(rs$steps, list(
    derivation_step(
      label = "Risk difference",
      symbol = "RD", formula = "R1 - R0",
      substituted = paste0(fmt(rs$r1, digits), " - ", fmt(rs$r0, digits)),
      result = rd
    ),
    derivation_step(
      label = "Standard error of RD",
      symbol = "SE",
      formula = "sqrt(R1(1-R1)/(a+b) + R0(1-R0)/(c+d))",
      substituted = paste0(
        "sqrt(", fmt(rs$r1, digits), "*", fmt(1 - rs$r1, digits), "/", fmt(rs$n1, digits),
        " + ", fmt(rs$r0, digits), "*", fmt(1 - rs$r0, digits), "/", fmt(rs$n0, digits), ")"
      ),
      result = se,
      note = "No log transform here. A difference can be negative, so it is already on a scale where the normal approximation applies directly."
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, lower limit"),
      symbol = "lower", formula = "RD - z * SE",
      substituted = paste0(fmt(rd, digits), " - ", fmt(z, digits), " * ", fmt(se, digits)),
      result = lo
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, upper limit"),
      symbol = "upper", formula = "RD + z * SE",
      substituted = paste0(fmt(rd, digits), " + ", fmt(z, digits), " * ", fmt(se, digits)),
      result = hi
    )
  ))

  derivation(
    method = "Risk difference",
    estimate = rd, symbol = "RD",
    ci = c(lo, hi), conf_level = conf_level,
    data = x, steps = steps,
    notes = paste0(
      "The risk difference is on the absolute scale: ", fmt(rd * 100, 2L),
      " excess cases per 100 exposed. Its reciprocal, ", fmt(1 / abs(rd), 1L),
      ", is the number needed to expose for one additional case."
    )
  )
}
