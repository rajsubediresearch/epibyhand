# epibyhand 0.1.0

First release.

* `epi2x2()` builds a 2 x 2 table in standard epidemiological orientation.
* `risk_ratio()`, `odds_ratio()` and `risk_difference()` return the measure
  together with the full worked derivation.
* `attributable_fraction()` computes the fraction among the exposed or in the
  population, showing that the direct, Levin and Miettinen formulas agree.
* `epi_strata()` and `collapse_strata()` handle stratified data.
* `mh_odds_ratio()` and `mh_risk_ratio()` pool across strata, printing the
  stratum weights and the crude estimate alongside the adjusted one.
* `homogeneity()` implements the Breslow-Day test with Tarone's correction.
* `check_work()` compares a hand calculation against every step of a
  derivation and reports where it diverged.
* `steps_table()` returns a derivation as a data frame.
* Output detail is controlled with `options(epibyhand.verbose = )` and
  `options(epibyhand.digits = )`.
