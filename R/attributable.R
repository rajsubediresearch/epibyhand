## Attributable fractions --------------------------------------------------
##
## Three formulas circulate for the population attributable fraction and
## students are rarely shown that they are the same quantity. This file
## computes all three and prints them together.

#' Attributable fraction, with worked derivation
#'
#' The proportion of disease that would be prevented by removing the exposure,
#' either among the exposed (`among = "exposed"`) or in the whole population
#' (`among = "population"`).
#'
#' For the population fraction the derivation computes the same quantity three
#' ways -- directly from risks, by Levin's formula from exposure prevalence,
#' and by Miettinen's formula from the proportion of cases exposed -- and shows
#' that they agree. They are not competing estimators; they are one quantity
#' written three ways, which is easier to believe once seen numerically.
#'
#' @param x An `epi2x2`, or the first of four cell counts.
#' @param ... Further cell counts if `x` is a bare count.
#' @param among `"exposed"` (the default) or `"population"`.
#' @param conf_level Confidence level for the interval.
#'
#' @return An `epibyhand_derivation`.
#'
#' @examples
#' attributable_fraction(epi2x2(36, 14, 30, 25))
#' attributable_fraction(epi2x2(36, 14, 30, 25), among = "population")
#' @export
attributable_fraction <- function(x, ..., among = c("exposed", "population"),
                                  conf_level = 0.95) {
  among <- match.arg(among)
  x <- epi2x2(x, ...)
  digits <- getOption("epibyhand.digits", 4L)
  z <- zcrit(conf_level)

  rs <- risk_steps(x, digits)
  r1 <- rs$r1
  r0 <- rs$r0
  rr <- r1 / r0

  se_rr <- sqrt(1 / x$a - 1 / rs$n1 + 1 / x$c - 1 / rs$n0)
  rr_lo <- exp(log(rr) - z * se_rr)
  rr_hi <- exp(log(rr) + z * se_rr)

  rr_step <- derivation_step(
    label = "Risk ratio",
    symbol = "RR", formula = "R1 / R0",
    substituted = paste0(fmt(r1, digits), " / ", fmt(r0, digits)),
    result = rr
  )

  if (among == "exposed") {
    afe <- (r1 - r0) / r1
    lo <- 1 - 1 / rr_lo
    hi <- 1 - 1 / rr_hi

    steps <- c(rs$steps, list(
      rr_step,
      derivation_step(
        label = "Attributable fraction among the exposed",
        symbol = "AFe", formula = "(R1 - R0) / R1",
        substituted = paste0("(", fmt(r1, digits), " - ", fmt(r0, digits), ") / ",
                             fmt(r1, digits)),
        result = afe
      ),
      derivation_step(
        label = "The same quantity from the risk ratio alone",
        symbol = "AFe", formula = "(RR - 1) / RR",
        substituted = paste0("(", fmt(rr, digits), " - 1) / ", fmt(rr, digits)),
        result = (rr - 1) / rr,
        note = paste(
          "Identical, because dividing through by R0 cancels it. This is why",
          "the attributable fraction among the exposed can be computed from a",
          "case-control study, where absolute risks are not available but the",
          "ratio is estimable."
        )
      ),
      derivation_step(
        label = paste0(round(conf_level * 100), "% confidence interval"),
        symbol = "CI", formula = "1 - 1/RR  applied to each limit of the RR interval",
        substituted = paste0("1 - 1/", fmt(rr_lo, digits), "  to  1 - 1/",
                             fmt(rr_hi, digits)),
        result = NA_real_,
        table = data.frame(
          limit = c("lower", "upper"),
          RR = c(rr_lo, rr_hi),
          AFe = c(lo, hi),
          stringsAsFactors = FALSE
        ),
        note = paste(
          "AFe increases monotonically with RR, so transforming the two RR",
          "limits gives the AFe limits directly. No new standard error is",
          "needed."
        )
      )
    ))

    return(derivation(
      method = "Attributable fraction among the exposed",
      estimate = afe, symbol = "AFe",
      ci = c(lo, hi), conf_level = conf_level,
      data = x, steps = steps,
      notes = c(
        check_zero_cells(x),
        paste0(
          "Read as: ", fmt(100 * afe, 1L), "% of cases among the exposed are ",
          "attributable to the exposure, if the association is causal and ",
          "unconfounded. The arithmetic cannot tell you whether it is."
        )
      )
    ))
  }

  ## Population attributable fraction
  n <- x$a + x$b + x$c + x$d
  rt <- (x$a + x$c) / n
  p <- (x$a + x$b) / n
  pc <- x$a / (x$a + x$c)

  paf <- (rt - r0) / rt
  levin <- p * (rr - 1) / (1 + p * (rr - 1))
  miettinen <- pc * (rr - 1) / rr

  paf_from_rr <- function(r) p * (r - 1) / (1 + p * (r - 1))
  lo <- paf_from_rr(rr_lo)
  hi <- paf_from_rr(rr_hi)

  steps <- c(rs$steps, list(
    rr_step,
    derivation_step(
      label = "Risk in the whole population",
      symbol = "Rt", formula = "(a + c) / n",
      substituted = paste0("(", fmt(x$a, digits), " + ", fmt(x$c, digits), ") / ",
                           fmt(n, digits)),
      result = rt
    ),
    derivation_step(
      label = "Population attributable fraction, directly from risks",
      symbol = "PAF", formula = "(Rt - R0) / Rt",
      substituted = paste0("(", fmt(rt, digits), " - ", fmt(r0, digits), ") / ",
                           fmt(rt, digits)),
      result = paf,
      note = paste(
        "The share of the population's risk that would disappear if everyone",
        "had the risk of the unexposed."
      )
    ),
    derivation_step(
      label = "Levin's formula, from exposure prevalence",
      symbol = "PAF", formula = "p(RR - 1) / (1 + p(RR - 1))",
      substituted = paste0(fmt(p, digits), "(", fmt(rr, digits), " - 1) / (1 + ",
                           fmt(p, digits), "(", fmt(rr, digits), " - 1))"),
      result = levin,
      note = paste0(
        "p = (a + b)/n = ", fmt(p, digits), " is the proportion exposed. Use ",
        "this form when you have the risk ratio from one study and exposure ",
        "prevalence from another."
      )
    ),
    derivation_step(
      label = "Miettinen's formula, from the proportion of cases exposed",
      symbol = "PAF", formula = "pc * (RR - 1) / RR",
      substituted = paste0(fmt(pc, digits), " * (", fmt(rr, digits), " - 1) / ",
                           fmt(rr, digits)),
      result = miettinen,
      note = paste0(
        "pc = a/(a + c) = ", fmt(pc, digits), " is the proportion of cases who ",
        "were exposed. All three lines above are the same number. They are one ",
        "quantity written three ways, not three estimators to choose between."
      )
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval"),
      symbol = "CI", formula = "Levin's formula applied to each limit of the RR interval",
      substituted = NULL, result = NA_real_,
      table = data.frame(
        limit = c("lower", "upper"),
        RR = c(rr_lo, rr_hi),
        PAF = c(lo, hi),
        stringsAsFactors = FALSE
      ),
      note = paste(
        "This holds exposure prevalence fixed at its observed value, so the",
        "interval is slightly too narrow. It is the version you can compute by",
        "hand; a delta-method interval that propagates uncertainty in p as well",
        "is wider."
      )
    )
  ))

  derivation(
    method = "Population attributable fraction",
    estimate = paf, symbol = "PAF",
    ci = c(lo, hi), conf_level = conf_level,
    data = x, steps = steps,
    notes = c(
      check_zero_cells(x),
      paste0(
        "Read as: ", fmt(100 * paf, 1L), "% of all cases in this population ",
        "are attributable to the exposure, if the association is causal."
      ),
      paste(
        "PAF depends on how common the exposure is, so it does not transfer",
        "between populations the way a risk ratio does. A strong risk factor",
        "that is rare has a small PAF; a weak one that is universal can have a",
        "large one. This is also why PAF from a case-control study needs the",
        "exposure prevalence of the source population, not of the controls."
      )
    )
  )
}
