tab <- epi2x2(36, 14, 30, 25)

test_that("the three PAF formulas give one number", {
  d <- attributable_fraction(tab, among = "population")
  tt <- steps_table(d)
  paf_rows <- tt$result[tt$symbol == "PAF"]
  expect_length(paf_rows, 3L)
  expect_equal(diff(range(paf_rows)), 0, tolerance = 1e-12)
})

test_that("PAF matches the direct definition", {
  d <- attributable_fraction(tab, among = "population")
  r0 <- 30 / 55
  rt <- 66 / 105
  expect_equal(estimate(d), (rt - r0) / rt)
})

test_that("AFe matches both of its formulas", {
  d <- attributable_fraction(tab, among = "exposed")
  r1 <- 36 / 50
  r0 <- 30 / 55
  rr <- r1 / r0
  expect_equal(estimate(d), (r1 - r0) / r1)
  expect_equal(estimate(d), (rr - 1) / rr)
})

test_that("AFe is at least as large as PAF", {
  afe <- estimate(attributable_fraction(tab, among = "exposed"))
  paf <- estimate(attributable_fraction(tab, among = "population"))
  expect_gte(afe, paf)
})

test_that("PAF approaches AFe as exposure becomes universal", {
  nearly_all <- epi2x2(3600, 1400, 3, 5)
  afe <- estimate(attributable_fraction(nearly_all, among = "exposed"))
  paf <- estimate(attributable_fraction(nearly_all, among = "population"))
  expect_lt(abs(afe - paf), 0.01)
})

test_that("PAF is near zero when exposure is rare", {
  rare_exposure <- epi2x2(2, 1, 300, 400)
  paf <- estimate(attributable_fraction(rare_exposure, among = "population"))
  expect_lt(paf, 0.02)
})

test_that("no association gives an attributable fraction of zero", {
  null_tab <- epi2x2(20, 30, 20, 30)
  expect_equal(estimate(attributable_fraction(null_tab)), 0)
  expect_equal(estimate(attributable_fraction(null_tab, among = "population")), 0)
})

test_that("a protective exposure gives a negative fraction", {
  protective <- epi2x2(10, 40, 30, 20)
  expect_lt(estimate(attributable_fraction(protective)), 0)
})

test_that("confidence limits bracket the estimate and respect the level", {
  d <- attributable_fraction(tab, among = "population")
  expect_lt(confint(d)[1], estimate(d))
  expect_gt(confint(d)[2], estimate(d))
  wide <- attributable_fraction(tab, among = "population", conf_level = 0.99)
  expect_lt(confint(wide)[1], confint(d)[1])
})

test_that("among is validated", {
  expect_error(attributable_fraction(tab, among = "everyone"))
})
