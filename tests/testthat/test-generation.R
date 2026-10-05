test_that("generators never return NULL and always preserve the sample size", {
  set.seed(1)
  ## small samples: the abundance vector must always be returned
  draws <- replicate(500, generate.hoppe.urn(1, 3)$abundance, simplify = FALSE)
  expect_false(any(vapply(draws, is.null, logical(1))))
  expect_true(all(vapply(draws, sum, numeric(1)) == 3))

  draws <- replicate(300, generate.pitman.urn(0.5, 0.5, 5)$abundance, simplify = FALSE)
  expect_false(any(vapply(draws, is.null, logical(1))))
  expect_true(all(vapply(draws, sum, numeric(1)) == 5))

  ## J = 1 is the smallest sample
  expect_equal(generate.hoppe.urn(1, 1)$abundance, 1)
  expect_equal(generate.pitman.urn(1, 0.5, 1)$abundance, 1)
  expect_equal(generate.hoppe.urn0(1, 1), 1L)
})

test_that("the total number of individuals is exact up to 1e12", {
  set.seed(2)
  for (J in c(2^31, 1e10, 1e11, 1e12)) {
    expect_equal(sum(generate.hoppe.urn(100, J)$abundance), J)
    expect_equal(sum(generate.pitman.urn(10, 0.3, J)$abundance), J)
  }
})

test_that("outputs are sorted in decreasing order and have the documented structure", {
  set.seed(3)
  h <- generate.hoppe.urn(11.3, 1e5)
  expect_named(h, c("abundance", "k"))
  expect_false(is.unsorted(rev(h$abundance)))
  expect_true(all(h$abundance > 0))
  expect_equal(h$k, kest(11.3, 1e5))
  p <- generate.pitman.urn(11.3, 0.1, 1e5)
  expect_named(p, c("abundance", "k"))
  expect_false(is.unsorted(rev(p$abundance)))
  u <- generate.hoppe.urn0(11.3, 1e5, seed = 1)
  expect_type(u, "integer")
  expect_false(is.unsorted(rev(u)))
  expect_equal(sum(u), 1e5)
})

test_that("generate.hoppe.urn0 is reproducible through `seed`", {
  a <- generate.hoppe.urn0(5, 1000, seed = 42)
  b <- generate.hoppe.urn0(5, 1000, seed = 42)
  d <- generate.hoppe.urn0(5, 1000, seed = 43)
  expect_identical(a, b)
  expect_false(identical(a, d))
})

test_that("GEM generators follow set.seed()", {
  set.seed(9); a <- generate.hoppe.urn(10, 1e4)$abundance
  set.seed(9); b <- generate.hoppe.urn(10, 1e4)$abundance
  expect_identical(a, b)
  set.seed(9); a <- generate.pitman.urn(10, 0.4, 1e4)$abundance
  set.seed(9); b <- generate.pitman.urn(10, 0.4, 1e4)$abundance
  expect_identical(a, b)
})

test_that("argument errors name the right function", {
  expect_error(generate.hoppe.urn(-1, 10), "generate.hoppe.urn\\(\\).*theta")
  expect_error(generate.hoppe.urn(0, 10), "generate.hoppe.urn\\(\\).*theta")
  expect_error(generate.hoppe.urn(1, 0), "generate.hoppe.urn\\(\\).*J")
  expect_error(generate.hoppe.urn(1, 2.5), "generate.hoppe.urn\\(\\).*whole")
  expect_error(generate.hoppe.urn0(-1, 10), "generate.hoppe.urn0\\(\\).*theta")
  expect_error(generate.hoppe.urn0(NA_real_, 10), "generate.hoppe.urn0\\(\\)")
  expect_error(generate.hoppe.urn0(1, 3e9), "generate.hoppe.urn0\\(\\).*too large")
  expect_error(generate.pitman.urn(-1, 0.5, 10), "generate.pitman.urn\\(\\).*theta")
  expect_error(generate.pitman.urn(1, 0, 10), "generate.pitman.urn\\(\\).*sigma")
  expect_error(generate.pitman.urn(1, 1, 10), "generate.pitman.urn\\(\\).*sigma")
  expect_error(generate.pitman.urn(1, 0.5, -3), "generate.pitman.urn\\(\\).*J")
})

## --- statistical correctness -----------------------------------------------------
## These compare the fast generators with the exact urn processes. They use a fixed
## seed, so they are deterministic; tolerances are several standard errors wide.

mean_se <- function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x)))
## two-sample KS p-value; the data are discrete (ties), which makes the p-value approximate
ks_p <- function(a, b) suppressWarnings(stats::ks.test(a, b)$p.value)

test_that("GEM Ewens sampler has the right law", {
  set.seed(100)
  th <- 11.3; J <- 3000; R <- 400
  gem <- replicate(R, { a <- generate.hoppe.urn(th, J)$abundance; c(length(a), max(a)) })
  urn <- replicate(R, { a <- generate.hoppe.urn0(th, J); c(length(a), max(a)) })
  ## number of species: theory
  k <- mean_se(gem[1, ])
  expect_lt(abs(k["mean"] - kest(th, J)), 4 * k["se"])
  ## number of species and largest abundance: same law as the exact urn
  expect_gt(ks_p(gem[1, ], urn[1, ]), 0.001)
  expect_gt(ks_p(gem[2, ], urn[2, ]), 0.001)
  expect_lt(abs(mean(gem[2, ]) - mean(urn[2, ])),
            4 * sqrt(stats::var(gem[2, ]) / R + stats::var(urn[2, ]) / R))
})

test_that("GEM Pitman sampler has the right law, also for heavy tails", {
  set.seed(101)
  for (p in list(c(10, 0.3), c(50, 0.6), c(5, 0.85))) {
    R <- 250; J <- 1500
    k <- replicate(R, length(generate.pitman.urn(p[1], p[2], J)$abundance))
    ks <- mean_se(k)
    expect_lt(abs(ks["mean"] - neutr:::kest_pitman(p[1], p[2], J)), 4 * ks["se"],
              label = paste("theta, sigma =", p[1], p[2]))
  }
  ## largest abundance and singletons against the independent sequential urn
  py_exact <- function(theta, sigma, J) {
    n <- numeric(J); k <- 0L
    for (j in 1:J) {
      if (j == 1 || stats::runif(1) < (theta + k * sigma) / (theta + j - 1)) { k <- k + 1L; n[k] <- 1 }
      else { s <- sample.int(k, 1, prob = n[1:k] - sigma); n[s] <- n[s] + 1 }
    }
    n[1:k]
  }
  stat <- function(a) c(max(a), sum(a == 1))
  G <- replicate(300, stat(generate.pitman.urn(10, 0.3, 800)$abundance))
  E <- replicate(120, stat(py_exact(10, 0.3, 800)))
  expect_gt(ks_p(G[1, ], E[1, ]), 0.001)
  expect_gt(ks_p(G[2, ], E[2, ]), 0.001)
})

test_that("the sequential Pitman urn used for the tail has the right law", {
  set.seed(102)
  R <- 400
  k <- replicate(R, length(neutr:::pitman_urn(20, 0.5, 600)))
  expect_lt(abs(mean(k) - neutr:::kest_pitman(20, 0.5, 600)), 4 * stats::sd(k) / sqrt(R))
  expect_true(all(replicate(50, sum(neutr:::pitman_urn(3, 0.9, 200)) == 200)))
})

test_that("heavy-tailed Pitman draws terminate and are exact", {
  set.seed(103)
  out <- generate.pitman.urn(50, 0.6, 5e4)
  expect_equal(sum(out$abundance), 5e4)
  out <- generate.pitman.urn(10, 0.9, 1e5)
  expect_equal(sum(out$abundance), 1e5)
})
