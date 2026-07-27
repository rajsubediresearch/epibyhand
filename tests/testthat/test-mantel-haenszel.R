strata2 <- epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25),
                      labels = c("Men", "Women"))

test_that("epi_strata accepts counts, matrices, lists and arrays alike", {
  from_counts <- epi_strata(c(10, 20, 5, 40), c(30, 15, 20, 25))
  from_list <- epi_strata(list(c(10, 20, 5, 40), c(30, 15, 20, 25)))
  arr <- array(c(10, 5, 20, 40, 30, 20, 15, 25), dim = c(2, 2, 2))
  from_array <- epi_strata(arr)

  expect_equal(estimate(mh_odds_ratio(from_counts)),
               estimate(mh_odds_ratio(from_list)))
  expect_equal(estimate(mh_odds_ratio(from_counts)),
               estimate(mh_odds_ratio(from_array)))
})

test_that("a stratified analysis needs at least two strata", {
  expect_error(epi_strata(c(1, 2, 3, 4)), "at least two")
})

test_that("MH odds ratio matches stats::mantelhaen.test", {
  arr <- array(c(10, 5, 20, 40, 30, 20, 15, 25), dim = c(2, 2, 2))
  ref <- stats::mantelhaen.test(arr, correct = FALSE)
  d <- mh_odds_ratio(strata2)
  expect_equal(estimate(d), unname(ref$estimate), tolerance = 1e-10)
  expect_equal(confint(d), as.numeric(ref$conf.int), tolerance = 1e-8)
})

test_that("the pooled estimate lies between the stratum estimates", {
  d <- mh_odds_ratio(strata2)
  expect_gte(estimate(d), min(d$stratum_estimates))
  expect_lte(estimate(d), max(d$stratum_estimates))
})

test_that("MH is a weighted average with weights S_i", {
  d <- mh_odds_ratio(strata2)
  tt <- steps_table(d)
  s_step <- d$steps[[2L]]$table
  weighted <- sum(s_step$S_i * d$stratum_estimates) / sum(s_step$S_i)
  expect_equal(weighted, estimate(d), tolerance = 1e-10)
})

test_that("MH risk ratio matches the hand calculation", {
  d <- mh_risk_ratio(strata2)
  num <- 10 * 45 / 75 + 30 * 45 / 90
  den <- 5 * 30 / 75 + 20 * 45 / 90
  expect_equal(estimate(d), num / den)
})

test_that("collapsing strata reproduces the crude table", {
  crude <- collapse_strata(strata2)
  expect_equal(crude$a, 40)
  expect_equal(crude$b, 35)
  expect_equal(crude$c, 25)
  expect_equal(crude$d, 65)
  expect_equal(mh_odds_ratio(strata2)$crude, estimate(odds_ratio(crude)))
})

test_that("no confounding leaves crude and adjusted nearly equal", {
  same <- epi_strata(c(20, 20, 10, 40), c(40, 40, 20, 80))
  d <- mh_odds_ratio(same)
  expect_equal(d$crude, estimate(d), tolerance = 1e-8)
})

test_that("identical stratum odds ratios give a homogeneity statistic near zero", {
  same <- epi_strata(c(20, 20, 10, 40), c(40, 40, 20, 80))
  expect_lt(estimate(homogeneity(same)), 1e-6)
})

test_that("divergent stratum odds ratios raise the statistic", {
  split <- epi_strata(c(40, 10, 10, 40), c(10, 40, 40, 10))
  expect_gt(estimate(homogeneity(split)), estimate(homogeneity(strata2)))
})

test_that("Tarone's correction reduces the statistic", {
  raw <- estimate(homogeneity(strata2, tarone = FALSE))
  corrected <- estimate(homogeneity(strata2, tarone = TRUE))
  expect_lte(corrected, raw)
})

test_that("expected counts respect the table margins", {
  a <- epibyhand:::solve_expected_a(psi = 3, n1 = 30, m1 = 15, N = 75)
  expect_gte(a, max(0, 30 + 15 - 75))
  expect_lte(a, min(30, 15))
})

test_that("a common odds ratio of one gives independence expected counts", {
  expect_equal(epibyhand:::solve_expected_a(psi = 1, n1 = 30, m1 = 15, N = 75),
               30 * 15 / 75)
})

test_that("steps carrying only a table print without an NA result", {
  out <- capture.output(print(mh_odds_ratio(strata2)))
  expect_false(any(grepl("= NA", out, fixed = TRUE)))
})

test_that("check_work reaches into a stratified derivation", {
  d <- mh_odds_ratio(strata2)
  expect_true(check_work(d, 2.9286))
  expect_true(check_work(d, 0.3541, step = "SE"))
})
