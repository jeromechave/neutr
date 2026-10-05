test_that("zeroremove drops zeros and missing values", {
  expect_equal(zeroremove(c(2.1, 0.0, 3.0)), c(2.1, 3.0))
  expect_equal(zeroremove(c(NA, 0, 1, NA, 2)), c(1, 2))
  expect_equal(zeroremove(c(NaN, 5)), 5)
  expect_equal(zeroremove(c(0, 0)), numeric(0))
  expect_equal(zeroremove(c(-1, 0, 1)), c(-1, 1))
  expect_identical(zeroremove(integer(0)), integer(0))
})

test_that("rmultinom_big preserves the total beyond .Machine$integer.max", {
  set.seed(1)
  for (size in c(10, 2^31 - 1, 2^31, 1e10, 1e12 + 7)) {
    x <- neutr:::rmultinom_big(size, c(0.5, 0.3, 0.2))
    expect_equal(sum(x), size)
    expect_true(all(x >= 0))
  }
  x <- neutr:::rmultinom_big(1e12, c(0.5, 0.3, 0.2))
  expect_equal(x / 1e12, c(0.5, 0.3, 0.2), tolerance = 1e-4)
})

test_that("the expected Pitman richness matches the exact recurrence", {
  ## E[K_{n+1}] = E[K_n] (1 + sigma / (theta + n)) + theta / (theta + n)
  exact <- function(theta, sigma, J) {
    E <- 0
    for (n in 0:(J - 1)) E <- E * (1 + sigma / (theta + n)) + theta / (theta + n)
    E
  }
  for (p in list(c(10, 0.3, 500), c(50, 0.6, 1000), c(1, 0.5, 10), c(1e4, 0.6, 100),
                 c(1e7, 0.6, 1000), c(0.5, 0.9, 200))) {
    expect_equal(neutr:::kest_pitman(p[1], p[2], p[3]), exact(p[1], p[2], p[3]),
                 tolerance = 1e-9, info = paste(p, collapse = " "))
  }
  ## theta = 0 is allowed (limit of the formula)
  expect_equal(neutr:::kest_pitman(0, 0.5, 100), exact(1e-12, 0.5, 100), tolerance = 1e-6)
})

test_that("input checks give informative messages naming the function", {
  expect_error(optim.ewens(c(1, -2, 3)), "optim.ewens\\(\\).*non-negative")
  expect_error(optim.ewens(c(1, NA, 3)), "optim.ewens\\(\\).*NA")
  expect_error(optim.ewens(c(1.5, 2)), "optim.ewens\\(\\).*whole numbers")
  expect_error(optim.ewens(numeric(0)), "optim.ewens\\(\\).*non-empty")
  expect_error(optim.ewens(c(0, 0)), "at least 1 species")
  expect_error(optim.ewens("a"), "optim.ewens\\(\\)")
})
