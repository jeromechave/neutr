test_that("kest is the exact expected number of species and is vectorised", {
  expect_equal(kest(11.3, 234373), 112.8284, tolerance = 1e-6)
  expect_equal(kest(c(1, 5), 100), c(kest(1, 100), kest(5, 100)))
  expect_equal(kest(1, 10), sum(1 / (1 + 0:9)))        # E[K] = sum_i theta / (theta + i)
})

test_that("kest.gt is a deterministic approximation close to the simulation", {
  set.seed(1)
  det <- kest.gt(11.3, 234373, 50)
  expect_s3_class(det, "integrate")
  sim <- kest.gt0(11.3, 234373, 50, nrep = 60)
  expect_equal(det$value, sim$k, tolerance = 0.05)
  ## monotone in N and bounded by the total richness
  vals <- vapply(c(0, 1, 10, 100), function(N) kest.gt(11.3, 234373, N)$value, numeric(1))
  expect_true(all(diff(vals) < 0))
  expect_lt(vals[1], kest(11.3, 234373) + 1e-6)
})

test_that("kest.gt0 has an nrep argument and a documented return value", {
  set.seed(2)
  out <- kest.gt0(5, 5000, 10, nrep = 10)
  expect_named(out, c("k", "sigmak"))
  expect_true(out$k > 0 && out$sigmak >= 0)
  expect_error(kest.gt0(5, 5000, 10, nrep = 1), "nrep")
})

test_that("richness functions validate their arguments", {
  expect_error(kest.gt(-1, 100, 5), "kest.gt\\(\\).*theta")
  expect_error(kest.gt(1, 100, 100), "kest.gt\\(\\).*N")
  expect_error(kest.gt(1, 100, -1), "kest.gt\\(\\).*N")
  expect_error(kest.gt0(0, 100, 5), "kest.gt0\\(\\).*theta")
  expect_error(kest.gt0(1, 100, NA_real_), "kest.gt0\\(\\).*N")
})
