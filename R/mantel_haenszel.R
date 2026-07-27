## Stratified analysis: Mantel-Haenszel -----------------------------------
##
## The weights are the point of this file. Every package returns the pooled
## estimate; almost none show the student that it is a weighted average of
## the stratum-specific estimates, or what the weights are.

#' @keywords internal
#' @noRd
stratum_frame <- function(x) {
  data.frame(
    Stratum = x$labels,
    a = vapply(x$tables, function(t) t$a, numeric(1L)),
    b = vapply(x$tables, function(t) t$b, numeric(1L)),
    c = vapply(x$tables, function(t) t$c, numeric(1L)),
    d = vapply(x$tables, function(t) t$d, numeric(1L)),
    stringsAsFactors = FALSE
  )
}

#' Mantel-Haenszel odds ratio with worked derivation
#'
#' Pools stratum-specific odds ratios into a single adjusted estimate,
#' showing the stratum weights explicitly and comparing the result against
#' the crude estimate so that confounding is visible as arithmetic.
#'
#' The confidence interval uses the Robins-Breslow-Greenland variance, which
#' is valid both for a few large strata and for many small ones.
#'
#' @param x An `epi_strata`, or arguments passed to [epi_strata()].
#' @param ... Further strata if `x` is not already an `epi_strata`.
#' @param conf_level Confidence level for the interval.
#'
#' @return An `epibyhand_derivation`.
#'
#' @examples
#' s <- epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25),
#'                 labels = c("Men", "Women"))
#' mh_odds_ratio(s)
#' @export
mh_odds_ratio <- function(x, ..., conf_level = 0.95) {
  x <- epi_strata(x, ...)
  digits <- getOption("epibyhand.digits", 4L)
  z <- zcrit(conf_level)

  df <- stratum_frame(x)
  n <- df$a + df$b + df$c + df$d
  or_i <- (df$a * df$d) / (df$b * df$c)
  Ri <- df$a * df$d / n
  Si <- df$b * df$c / n
  R <- sum(Ri)
  S <- sum(Si)
  or_mh <- R / S

  Pi <- (df$a + df$d) / n
  Qi <- (df$b + df$c) / n
  varlog <- sum(Pi * Ri) / (2 * R^2) +
    sum(Pi * Si + Qi * Ri) / (2 * R * S) +
    sum(Qi * Si) / (2 * S^2)
  se <- sqrt(varlog)
  lo <- exp(log(or_mh) - z * se)
  hi <- exp(log(or_mh) + z * se)

  crude_tab <- collapse_strata(x)
  or_crude <- (crude_tab$a * crude_tab$d) / (crude_tab$b * crude_tab$c)
  pct <- 100 * (or_crude - or_mh) / or_mh

  steps <- list(
    derivation_step(
      label = "Odds ratio within each stratum",
      symbol = "OR_i", formula = "(a_i * d_i) / (b_i * c_i)",
      substituted = NULL,
      result = NA_real_,
      table = data.frame(
        Stratum = df$Stratum, a = df$a, b = df$b, c = df$c, d = df$d,
        n_i = n, OR_i = or_i, stringsAsFactors = FALSE
      ),
      note = paste(
        "Look at these before pooling. If they disagree substantially the",
        "stratifying variable is an effect modifier, and a single pooled",
        "number hides the finding rather than reporting it. Test with",
        "homogeneity()."
      )
    ),
    derivation_step(
      label = "Stratum contributions",
      symbol = NULL, formula = "R_i = a_i*d_i/n_i     S_i = b_i*c_i/n_i",
      substituted = NULL,
      result = NA_real_,
      table = data.frame(
        Stratum = df$Stratum, R_i = Ri, S_i = Si, stringsAsFactors = FALSE
      ),
      note = paste0(
        "S_i is the weight this stratum carries in the pooled estimate. ",
        "It is largest where the stratum has the most information, so small ",
        "or unbalanced strata contribute little."
      )
    ),
    derivation_step(
      label = "Pooled odds ratio",
      symbol = "OR_MH", formula = "sum(R_i) / sum(S_i)",
      substituted = paste0(fmt(R, digits), " / ", fmt(S, digits)),
      result = or_mh,
      note = paste0(
        "Equivalently sum(S_i * OR_i) / sum(S_i): a weighted average of the ",
        "stratum odds ratios with weights S_i. The pooled value must fall ",
        "between the smallest and largest stratum estimate; if yours does ",
        "not, the arithmetic is wrong."
      )
    ),
    derivation_step(
      label = "Standard error of log(OR_MH), Robins-Breslow-Greenland",
      symbol = "SE",
      formula = "sqrt( sum(P_i R_i)/(2R^2) + sum(P_i S_i + Q_i R_i)/(2RS) + sum(Q_i S_i)/(2S^2) )",
      substituted = paste0(
        "sqrt(", fmt(sum(Pi * Ri) / (2 * R^2), digits), " + ",
        fmt(sum(Pi * Si + Qi * Ri) / (2 * R * S), digits), " + ",
        fmt(sum(Qi * Si) / (2 * S^2), digits), ")"
      ),
      result = se,
      table = data.frame(
        Stratum = df$Stratum, P_i = Pi, Q_i = Qi, stringsAsFactors = FALSE
      ),
      note = "P_i = (a_i+d_i)/n_i and Q_i = (b_i+c_i)/n_i are the concordant and discordant proportions in each stratum."
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, lower limit"),
      symbol = "lower", formula = "exp(log(OR_MH) - z * SE)",
      substituted = paste0("exp(", fmt(log(or_mh), digits), " - ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = lo
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, upper limit"),
      symbol = "upper", formula = "exp(log(OR_MH) + z * SE)",
      substituted = paste0("exp(", fmt(log(or_mh), digits), " + ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = hi
    ),
    derivation_step(
      label = "Crude odds ratio, ignoring the strata",
      symbol = "OR_crude", formula = "(A * D) / (B * C)  on the collapsed table",
      substituted = paste0("(", fmt(crude_tab$a, digits), " * ", fmt(crude_tab$d, digits),
                           ") / (", fmt(crude_tab$b, digits), " * ",
                           fmt(crude_tab$c, digits), ")"),
      result = or_crude,
      note = paste0(
        "The crude estimate differs from the adjusted one by ",
        fmt(pct, 1L), "%. A change beyond about 10% is the usual working ",
        "signal that the stratifying variable confounds the association. ",
        "This is a judgement about the data, not a hypothesis test -- do not ",
        "decide it with a p-value."
      )
    )
  )

  out <- derivation(
    method = "Mantel-Haenszel odds ratio",
    estimate = or_mh, symbol = "OR_MH",
    ci = c(lo, hi), conf_level = conf_level,
    data = x, steps = steps,
    notes = c(
      unlist(lapply(x$tables, check_zero_cells)),
      paste0(
        "Pooling assumes one common odds ratio underlies every stratum. ",
        "Check that with homogeneity() before reporting this number."
      )
    )
  )
  out$stratum_estimates <- stats::setNames(or_i, df$Stratum)
  out$crude <- or_crude
  out
}

#' Mantel-Haenszel risk ratio with worked derivation
#'
#' @inheritParams mh_odds_ratio
#' @return An `epibyhand_derivation`.
#' @examples
#' mh_risk_ratio(epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25)))
#' @export
mh_risk_ratio <- function(x, ..., conf_level = 0.95) {
  x <- epi_strata(x, ...)
  digits <- getOption("epibyhand.digits", 4L)
  z <- zcrit(conf_level)

  df <- stratum_frame(x)
  n1 <- df$a + df$b
  n0 <- df$c + df$d
  n <- n1 + n0
  rr_i <- (df$a / n1) / (df$c / n0)

  num_i <- df$a * n0 / n
  den_i <- df$c * n1 / n
  num <- sum(num_i)
  den <- sum(den_i)
  rr_mh <- num / den

  varlog <- sum((n1 * n0 * (df$a + df$c) - df$a * df$c * n) / n^2) / (num * den)
  se <- sqrt(varlog)
  lo <- exp(log(rr_mh) - z * se)
  hi <- exp(log(rr_mh) + z * se)

  crude_tab <- collapse_strata(x)
  rr_crude <- (crude_tab$a / (crude_tab$a + crude_tab$b)) /
    (crude_tab$c / (crude_tab$c + crude_tab$d))
  pct <- 100 * (rr_crude - rr_mh) / rr_mh

  steps <- list(
    derivation_step(
      label = "Risk ratio within each stratum",
      symbol = "RR_i", formula = "(a_i/(a_i+b_i)) / (c_i/(c_i+d_i))",
      substituted = NULL, result = NA_real_,
      table = data.frame(
        Stratum = df$Stratum, a = df$a, b = df$b, c = df$c, d = df$d,
        RR_i = rr_i, stringsAsFactors = FALSE
      )
    ),
    derivation_step(
      label = "Stratum contributions",
      symbol = NULL,
      formula = "num_i = a_i*n0_i/n_i     den_i = c_i*n1_i/n_i",
      substituted = NULL, result = NA_real_,
      table = data.frame(
        Stratum = df$Stratum, num_i = num_i, den_i = den_i,
        stringsAsFactors = FALSE
      ),
      note = "n1_i is the exposed total and n0_i the unexposed total in stratum i."
    ),
    derivation_step(
      label = "Pooled risk ratio",
      symbol = "RR_MH", formula = "sum(num_i) / sum(den_i)",
      substituted = paste0(fmt(num, digits), " / ", fmt(den, digits)),
      result = rr_mh
    ),
    derivation_step(
      label = "Standard error of log(RR_MH), Greenland-Robins",
      symbol = "SE",
      formula = "sqrt( sum((n1_i*n0_i*(a_i+c_i) - a_i*c_i*n_i)/n_i^2) / (sum(num_i)*sum(den_i)) )",
      substituted = paste0(
        "sqrt(", fmt(sum((n1 * n0 * (df$a + df$c) - df$a * df$c * n) / n^2), digits),
        " / (", fmt(num, digits), " * ", fmt(den, digits), "))"
      ),
      result = se
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, lower limit"),
      symbol = "lower", formula = "exp(log(RR_MH) - z * SE)",
      substituted = paste0("exp(", fmt(log(rr_mh), digits), " - ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = lo
    ),
    derivation_step(
      label = paste0(round(conf_level * 100), "% confidence interval, upper limit"),
      symbol = "upper", formula = "exp(log(RR_MH) + z * SE)",
      substituted = paste0("exp(", fmt(log(rr_mh), digits), " + ", fmt(z, digits),
                           " * ", fmt(se, digits), ")"),
      result = hi
    ),
    derivation_step(
      label = "Crude risk ratio, ignoring the strata",
      symbol = "RR_crude", formula = "(A/(A+B)) / (C/(C+D))  on the collapsed table",
      substituted = paste0("(", fmt(crude_tab$a, digits), "/",
                           fmt(crude_tab$a + crude_tab$b, digits), ") / (",
                           fmt(crude_tab$c, digits), "/",
                           fmt(crude_tab$c + crude_tab$d, digits), ")"),
      result = rr_crude,
      note = paste0("The crude estimate differs from the adjusted one by ",
                    fmt(pct, 1L), "%.")
    )
  )

  out <- derivation(
    method = "Mantel-Haenszel risk ratio",
    estimate = rr_mh, symbol = "RR_MH",
    ci = c(lo, hi), conf_level = conf_level,
    data = x, steps = steps,
    notes = c(
      unlist(lapply(x$tables, check_zero_cells)),
      "A risk ratio needs closed cohorts with equal follow-up in every stratum."
    )
  )
  out$stratum_estimates <- stats::setNames(rr_i, df$Stratum)
  out$crude <- rr_crude
  out
}

#' Breslow-Day test of homogeneity, with worked derivation
#'
#' Tests whether one common odds ratio can be assumed across strata -- the
#' assumption that makes [mh_odds_ratio()] meaningful. The statistic is a sum
#' of per-stratum contributions, and those contributions are shown, so a
#' single badly behaved stratum can be seen rather than inferred.
#'
#' Tarone's correction is applied by default; without it the statistic is
#' slightly too large.
#'
#' @param x An `epi_strata`, or arguments passed to [epi_strata()].
#' @param ... Further strata if `x` is not already an `epi_strata`.
#' @param tarone Apply Tarone's correction. Defaults to `TRUE`.
#'
#' @return An `epibyhand_derivation`.
#'
#' @examples
#' homogeneity(epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25)))
#' @export
homogeneity <- function(x, ..., tarone = TRUE) {
  x <- epi_strata(x, ...)
  digits <- getOption("epibyhand.digits", 4L)

  df <- stratum_frame(x)
  psi <- estimate(mh_odds_ratio(x))

  n1 <- df$a + df$b
  m1 <- df$a + df$c
  N <- df$a + df$b + df$c + df$d

  expected <- vapply(seq_len(nrow(df)), function(i) {
    solve_expected_a(psi, n1[i], m1[i], N[i])
  }, numeric(1L))

  varA <- 1 / (1 / expected + 1 / (n1 - expected) + 1 / (m1 - expected) +
                 1 / (N - n1 - m1 + expected))
  contrib <- (df$a - expected)^2 / varA
  stat <- sum(contrib)

  correction <- 0
  if (tarone) {
    correction <- (sum(df$a) - sum(expected))^2 / sum(varA)
    stat <- stat - correction
  }

  k <- nrow(df)
  p <- stats::pchisq(stat, df = k - 1L, lower.tail = FALSE)

  steps <- list(
    derivation_step(
      label = "Expected exposed cases in each stratum under a common OR",
      symbol = "A_i",
      formula = "root of  (1-psi)A^2 + (N - n1 - m1 + psi(n1+m1))A - psi*n1*m1 = 0",
      substituted = paste0("psi = OR_MH = ", fmt(psi, digits)),
      result = NA_real_,
      table = data.frame(
        Stratum = df$Stratum, a_obs = df$a, A_exp = expected,
        Var_A = varA, stringsAsFactors = FALSE
      ),
      note = paste(
        "A_i is what cell a would be if this stratum had exactly the pooled",
        "odds ratio, holding its margins fixed. Solving a quadratic is the",
        "one step here you would not do by hand."
      )
    ),
    derivation_step(
      label = "Contribution of each stratum to the statistic",
      symbol = "X2_i", formula = "(a_i - A_i)^2 / Var(A_i)",
      substituted = NULL, result = NA_real_,
      table = data.frame(
        Stratum = df$Stratum, contribution = contrib, stringsAsFactors = FALSE
      ),
      note = "A single large contribution means one stratum is driving the result."
    ),
    derivation_step(
      label = if (tarone) "Statistic, with Tarone's correction" else "Statistic",
      symbol = "X2",
      formula = if (tarone) "sum(X2_i) - (sum(a_i) - sum(A_i))^2 / sum(Var(A_i))" else "sum(X2_i)",
      substituted = if (tarone) {
        paste0(fmt(sum(contrib), digits), " - ", fmt(correction, digits))
      } else {
        paste(fmt(contrib, digits), collapse = " + ")
      },
      result = stat
    ),
    derivation_step(
      label = "Reference distribution",
      symbol = "p", formula = "P(chi-squared with K - 1 df > X2)",
      substituted = paste0("K = ", k, ", so df = ", k - 1L),
      result = p
    )
  )

  derivation(
    method = "Breslow-Day test of homogeneity",
    estimate = stat, symbol = "X2",
    ci = NULL, data = x, steps = steps,
    notes = c(
      paste0(
        "p = ", fmt(p, digits), " on ", k - 1L, " degrees of freedom."
      ),
      paste(
        "A large p-value is not evidence that the odds ratios are equal.",
        "This test has poor power with small strata, so it will usually fail",
        "to reject whether or not effect modification is present. Inspect the",
        "stratum-specific estimates as well; they are the more informative",
        "thing."
      )
    )
  )
}

#' Solve for the expected count in cell a under a fixed odds ratio
#' @keywords internal
#' @noRd
solve_expected_a <- function(psi, n1, m1, N) {
  lower <- max(0, n1 + m1 - N)
  upper <- min(n1, m1)

  if (abs(psi - 1) < 1e-10) return(n1 * m1 / N)

  alpha <- 1 - psi
  beta <- N - n1 - m1 + psi * (n1 + m1)
  gamma <- -psi * n1 * m1

  disc <- beta^2 - 4 * alpha * gamma
  if (disc < 0) return(n1 * m1 / N)
  roots <- c((-beta + sqrt(disc)) / (2 * alpha),
             (-beta - sqrt(disc)) / (2 * alpha))

  ok <- roots[roots >= lower - 1e-8 & roots <= upper + 1e-8]
  if (!length(ok)) return(n1 * m1 / N)
  min(max(ok[1L], lower + 1e-10), upper - 1e-10)
}
