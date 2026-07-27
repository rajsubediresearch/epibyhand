tab <- epi2x2(36, 14, 30, 25)

test_that("epi2x2 accepts counts and matrices identically", {
  m <- matrix(c(36, 14, 30, 25), nrow = 2, byrow = TRUE)
  expect_equal(unclass(epi2x2(m))[c("a", "b", "c", "d")],
               unclass(tab)[c("a", "b", "c", "d")])
})

test_that("epi2x2 rejects malformed input", {
  expect_error(epi2x2(matrix(1:6, nrow = 2)), "2 x 2")
  expect_error(epi2x2(1, 2, 3), "four cell counts")
  expect_error(epi2x2(-1, 2, 3, 4), "negative")
})

test_that("odds ratio matches the closed form ad/bc", {
  d <- odds_ratio(tab)
  expect_equal(estimate(d), (36 * 25) / (14 * 30))
})

test_that("odds ratio interval matches Woolf by hand", {
  d <- odds_ratio(tab)
  se <- sqrt(1 / 36 + 1 / 14 + 1 / 30 + 1 / 25)
  z <- qnorm(0.975)
  expect_equal(confint(d),
               exp(log((36 * 25) / (14 * 30)) + c(-1, 1) * z * se))
})

test_that("risk ratio and risk difference agree with hand values", {
  expect_equal(estimate(risk_ratio(tab)), (36 / 50) / (30 / 55))
  expect_equal(estimate(risk_difference(tab)), 36 / 50 - 30 / 55)
})

test_that("odds ratio is further from the null than the risk ratio", {
  or <- estimate(odds_ratio(tab))
  rr <- estimate(risk_ratio(tab))
  expect_gt(or, rr)
  expect_gt(rr, 1)
})

test_that("a rare outcome pulls OR toward RR", {
  rare <- epi2x2(3, 997, 1, 999)
  or <- estimate(odds_ratio(rare))
  rr <- estimate(risk_ratio(rare))
  expect_lt(abs(or - rr), 0.01)
})

test_that("confidence level is respected", {
  wide <- odds_ratio(tab, conf_level = 0.99)
  narrow <- odds_ratio(tab, conf_level = 0.90)
  expect_lt(confint(wide)[1], confint(narrow)[1])
  expect_gt(confint(wide)[2], confint(narrow)[2])
})

test_that("zero cells produce a note rather than a silent correction", {
  d <- odds_ratio(epi2x2(0, 10, 5, 5))
  expect_true(any(grepl("continuity correction", d$notes)))
  expect_true(is.nan(estimate(d)) || estimate(d) == 0)
})

test_that("steps_table exposes every step numerically", {
  tt <- steps_table(odds_ratio(tab))
  expect_s3_class(tt, "data.frame")
  expect_equal(nrow(tt), 6L)
  expect_true(all(c("OR", "SE", "lower", "upper") %in% tt$symbol))
  expect_equal(tt$result[tt$symbol == "OR"], estimate(odds_ratio(tab)))
})

test_that("check_work locates a value that matches an earlier step", {
  d <- odds_ratio(tab)
  expect_true(check_work(d, 2.1429))
  expect_false(check_work(d, 2.5714))
  expect_output(check_work(d, 2.5714), "stopped early")
})

test_that("check_work targets a named step", {
  d <- odds_ratio(tab)
  expect_true(check_work(d, 0.4154, step = "SE"))
  expect_error(check_work(d, 1, step = "nope"), "No step named")
})

test_that("verbosity controls output length", {
  d <- odds_ratio(tab)
  short <- capture.output(print(d, verbose = 0))
  full <- capture.output(print(d, verbose = 2))
  expect_lt(length(short), length(full))
  expect_true(any(grepl("OR = ", short)))
})
