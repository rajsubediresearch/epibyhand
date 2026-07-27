# epibyhand

<!-- badges: start -->
[![R-CMD-check](https://github.com/rajsubediresearch/epibyhand/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rajsubediresearch/epibyhand/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Classical epidemiological measures that show their work.

Every function returns the answer **and** the reasoning that produced it: each
intermediate quantity, the formula in symbols, and the formula with the
observed numbers substituted in. It is written for teaching, for checking hand
calculations, and for generating worked solutions in course materials.

## Scope

**epibyhand covers methods a student can compute by hand on paper.**

That boundary is deliberate and load-bearing. It is what keeps the package
small enough to remain correct without continuous maintenance, and it is why
there is no regression modelling here — anything requiring an iterative fit is
out of scope by definition rather than by omission. Feature requests that cross
that line will be declined, not deferred.

The package imports `stats` and nothing else.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rajsubediresearch/epibyhand")
```

## Usage

```r
library(epibyhand)

odds_ratio(epi2x2(36, 14, 30, 25))
```

```
-- Odds ratio --------------------------------------------------------------

Data
             Case  Non-case  Total
  Exposed      36        14     50
  Unexposed    30        25     55
  Total        66        39    105

Step 1  Odds of Case among the exposed
        odds1 = a / b
              = 36 / 14
              = 2.5714

Step 2  Odds of Case among the unexposed
        odds0 = c / d
              = 30 / 25
              = 1.2

Step 3  Odds ratio
        OR = odds1 / odds0
           = 2.5714 / 1.2
           = 2.1429
        # Equivalently OR = ad / bc, which is why the odds ratio is the
        # same whether you condition on exposure or on outcome. That
        # symmetry is what makes it usable in case-control studies.

...

Result
  OR = 2.1429 (95% CI 0.9493 to 4.8369)
```

### Stratified analysis

The Mantel-Haenszel functions show the stratum weights explicitly, and print
the crude estimate alongside the adjusted one so that confounding is visible
as arithmetic rather than asserted.

```r
s <- epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25),
                labels = c("Men", "Women"))
mh_odds_ratio(s)
homogeneity(s)
```

### Checking a hand calculation

`check_work()` compares your value against the final answer, and when it does
not match, searches every intermediate step for one that does — so instead of
"wrong", you are told *where* you stopped or slipped.

```r
d <- odds_ratio(epi2x2(36, 14, 30, 25))

check_work(d, 2.5714)
#> Not a match. You gave 2.5714; the final estimate (OR) is 2.1429.
#>
#> Your value does match step 1: Odds of Case among the exposed.
#>   odds1 = a / b
#> You may have stopped early. The next step is: Odds of Case among the unexposed.
```

### Controlling detail

```r
options(epibyhand.verbose = 0)   # result only
options(epibyhand.verbose = 1)   # add symbolic formulas
options(epibyhand.verbose = 2)   # full worked solution (default)
options(epibyhand.digits = 3)
```

`steps_table()` returns the derivation as a data frame, for building answer
keys or rendering the working in a format this package does not provide.

## Design

One S3 class, `epibyhand_derivation`, carries the result and the steps that
produced it. All display logic lives in its print method, so adding a measure
means writing arithmetic and steps — never writing display code.

Silent corrections are avoided on principle. A zero cell produces an
explanation of what a continuity correction would do to the estimate, and
leaves the choice to you.

## Currently implemented

| Function | Measure |
|---|---|
| `risk_ratio()` | Risk ratio, with log-scale interval |
| `odds_ratio()` | Odds ratio, Woolf interval |
| `risk_difference()` | Risk difference, number needed to expose |
| `attributable_fraction()` | AFe and PAF, with all three PAF formulas shown to agree |
| `mh_odds_ratio()` | Mantel-Haenszel OR, Robins-Breslow-Greenland interval |
| `mh_risk_ratio()` | Mantel-Haenszel RR, Greenland-Robins interval |
| `homogeneity()` | Breslow-Day test, Tarone corrected |
| `check_work()` | Locate where a hand calculation diverged |

## License

MIT
