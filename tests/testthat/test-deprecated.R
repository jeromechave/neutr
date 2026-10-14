# Deprecated dotted names: each one signals a deprecation and returns exactly
# what its snake_case replacement returns.

ab <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
zoo <- c(8, 5, 3, 2, 1, 1) # example of Etienne (2005), as in test-etienne.R

test_that("deprecated estimation functions match their replacements", {
  lifecycle::expect_deprecated(old <- optim.ewens(ab))
  expect_identical(old, optim_ewens(ab))

  lifecycle::expect_deprecated(old <- optim.pitman(ab))
  expect_identical(old, optim_pitman(ab))

  m2 <- rbind(c(30, 10, 4, 0), c(2, 12, 9, 5))
  lifecycle::expect_deprecated(old <- optim.multideme(m2))
  expect_identical(old, optim_multideme(m2))
})

test_that("deprecated Etienne functions match their replacements", {
  lifecycle::expect_deprecated(old <- logl.etienne(zoo, 5, 0.5))
  expect_identical(old, logl_etienne(zoo, 5, 0.5))

  lifecycle::expect_deprecated(old <- optim.etienne(zoo))
  expect_identical(old, optim_etienne(zoo))
})

test_that("deprecated generation functions match their replacements", {
  set.seed(1)
  lifecycle::expect_deprecated(old <- generate.hoppe.urn(5, 1000))
  set.seed(1)
  expect_identical(old, generate_hoppe_urn(5, 1000))

  set.seed(1)
  lifecycle::expect_deprecated(old <- generate.pitman.urn(5, 0.3, 1000))
  set.seed(1)
  expect_identical(old, generate_pitman_urn(5, 0.3, 1000))

  lifecycle::expect_deprecated(old <- generate.hoppe.urn0(5, 1000, seed = 1))
  expect_identical(old, generate_hoppe_urn0(5, 1000, seed = 1))
})

test_that("deprecated richness functions match their replacements", {
  lifecycle::expect_deprecated(old <- kest.gt(5, 1000, 10))
  expect_identical(old, kest_gt(5, 1000, 10))

  set.seed(1)
  lifecycle::expect_deprecated(old <- kest.gt0(5, 1000, 10, nrep = 5))
  set.seed(1)
  expect_identical(old, kest_gt0(5, 1000, 10, nrep = 5))
})
