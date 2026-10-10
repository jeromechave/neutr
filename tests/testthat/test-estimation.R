ab <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)

## ---- optim.ewens -------------------------------------------------------------------

test_that("optim.ewens returns the log-likelihood of the partition", {
  fit <- optim.ewens(ab)
  J <- sum(ab); k <- length(ab)
  direct <- k * log(fit$theta) + lgamma(fit$theta) - lgamma(fit$theta + J) + sum(lgamma(ab))
  expect_named(fit, c("theta", "logl"))
  expect_lt(fit$logl, 0)
  expect_equal(fit$logl, direct)
  expect_equal(fit$logl, -433.618716843, tolerance = 1e-7)
})

test_that("optim.ewens solves the score equation to full precision", {
  set.seed(1)
  for (th in c(0.5, 3, 20, 200)) for (J in c(500, 20000)) {
    a <- generate.hoppe.urn(th, J)$abundance
    fit <- optim.ewens(a)
    k <- length(a)
    expect_equal(fit$theta * (digamma(fit$theta + J) - digamma(fit$theta)), k, tolerance = 1e-10)
  }
  expect_equal(optim.ewens(ab)$theta, 2.476993036, tolerance = 1e-8)
})

test_that("optim.ewens finds the maximum of the likelihood", {
  fit <- optim.ewens(ab)
  ll <- function(th) fit$logl + (length(ab) * log(th / fit$theta) -
         (lgamma(th + sum(ab)) - lgamma(fit$theta + sum(ab))) + (lgamma(th) - lgamma(fit$theta)))
  for (f in c(0.5, 0.9, 0.99, 1.01, 1.1, 2)) expect_lt(ll(fit$theta * f), fit$logl)
})

test_that("the optim.ewens estimate is not restricted to theta <= J", {
  a <- c(rep(1, 98), 2)                   # J = 100, k = 99: the MLE is far above J
  fit <- optim.ewens(a)
  expect_gt(fit$theta, 1000)
  expect_equal(fit$theta * (digamma(fit$theta + 100) - digamma(fit$theta)), 99, tolerance = 1e-8)
})

test_that("optim.ewens ignores empty species and recovers theta", {
  expect_equal(optim.ewens(c(5, 0, 3, 0, 2)), optim.ewens(c(5, 3, 2)))
  set.seed(5)
  a <- generate.hoppe.urn(20, 1e5)$abundance
  expect_lt(abs(optim.ewens(a)$theta / 20 - 1), 0.3)
})

test_that("optim.ewens handles the boundary cases with a warning", {
  expect_warning(f1 <- optim.ewens(10), "single species")
  expect_equal(f1$theta, 0)
  expect_warning(f2 <- optim.ewens(rep(1, 20)), "single individual")
  expect_equal(f2$theta, Inf)
  expect_equal(f1$logl, 0)
  expect_equal(f2$logl, 0)
})

test_that("the Ewens log-likelihood is the sigma = 0 case of the Pitman log-likelihood", {
  fit <- optim.ewens(ab)
  th <- fit$theta; J <- sum(ab); k <- length(ab)
  u <- sort(unique(ab)); cnt <- tabulate(match(ab, u))
  expect_equal(neutr:::pitman_logl(th, 1e-9, u, cnt, J, k), fit$logl, tolerance = 1e-5)
})

## ---- optim.pitman ------------------------------------------------------------------

test_that("the Pitman likelihood matches its closed-form expression", {
  naive <- function(th, sg, n) {
    J <- sum(n); k <- length(n)
    k * log(sg) - k * lgamma(1 - sg) + lgamma(th) - lgamma(th + J) +
      lgamma(th / sg + k) - lgamma(th / sg) + sum(lgamma(n - sg))
  }
  u <- sort(unique(ab)); cnt <- tabulate(match(ab, u)); J <- sum(ab); k <- length(ab)
  for (p in list(c(10, 0.1), c(3, 0.5), c(50, 0.8), c(0.7, 0.28)))
    expect_equal(neutr:::pitman_logl(p[1], p[2], u, cnt, J, k), naive(p[1], p[2], ab))
})

test_that("the analytic gradient matches finite differences", {
  set.seed(2)
  n <- generate.pitman.urn(10, 0.3, 5000)$abundance
  u <- sort(unique(n)); cnt <- tabulate(match(n, u)); J <- sum(n); k <- length(n)
  for (p in list(c(10, 0.1), c(3, 0.5), c(50, 0.8), c(0.7, 0.28))) {
    h <- 1e-6
    num <- c((neutr:::pitman_logl(p[1] + h, p[2], u, cnt, J, k) - neutr:::pitman_logl(p[1] - h, p[2], u, cnt, J, k)) / (2 * h),
             (neutr:::pitman_logl(p[1], p[2] + h, u, cnt, J, k) - neutr:::pitman_logl(p[1], p[2] - h, u, cnt, J, k)) / (2 * h))
    expect_equal(neutr:::pitman_grad(p[1], p[2], u, cnt, J, k), num, tolerance = 1e-5)
  }
})

test_that("optim.pitman finds an interior maximum", {
  fit <- optim.pitman(ab, c(10, 0.1))
  expect_named(fit, c("theta", "sigma", "logl", "converged"))
  expect_true(fit$converged)
  expect_equal(fit$logl, -432.4740558, tolerance = 1e-7)
  expect_equal(fit$theta, 0.7110667, tolerance = 1e-4)
  expect_equal(fit$sigma, 0.2781073, tolerance = 1e-4)
  u <- sort(unique(ab)); cnt <- tabulate(match(ab, u)); J <- sum(ab); k <- length(ab)
  for (d in list(c(1.05, 1), c(0.95, 1), c(1, 1.05), c(1, 0.95)))
    expect_lt(neutr:::pitman_logl(fit$theta * d[1], fit$sigma * d[2], u, cnt, J, k), fit$logl)
})

test_that("optim.pitman does not depend on the starting values", {
  a <- optim.pitman(ab, c(10, 0.1))
  b <- optim.pitman(ab, c(100, 0.8))
  c0 <- optim.pitman(ab, c(0.5, 0.05))
  expect_equal(b$logl, a$logl, tolerance = 1e-8)
  expect_equal(c0$logl, a$logl, tolerance = 1e-8)
  expect_equal(c0$theta, a$theta, tolerance = 1e-3)
})

test_that("the Ewens model is nested in the Pitman model: their log-likelihoods are comparable", {
  set.seed(7)
  for (par in list(c(10, 0.1), c(5, 0.4), c(30, 0.0001))) {
    a <- if (par[2] > 0.001) generate.pitman.urn(par[1], par[2], 20000)$abundance
         else generate.hoppe.urn(par[1], 20000)$abundance
    e <- optim.ewens(a)
    p <- suppressWarnings(optim.pitman(a))
    expect_gte(p$logl, e$logl - 1e-6)                # nested: LRT statistic >= 0
  }
})

test_that("optim.pitman recovers known parameters", {
  set.seed(11)
  a <- generate.pitman.urn(10, 0.3, 1e5)$abundance
  fit <- optim.pitman(a)
  expect_lt(abs(fit$sigma - 0.3), 0.05)
  expect_lt(abs(fit$theta / 10 - 1), 0.5)
  a <- generate.pitman.urn(50, 0.6, 5e4)$abundance
  fit <- optim.pitman(a)
  expect_lt(abs(fit$sigma - 0.6), 0.05)
})

test_that("optim.pitman accepts theta_init < 1 and validates its arguments", {
  expect_no_error(optim.pitman(ab, c(0.5, 0.1)))      # theta_init may be below 1
  expect_error(optim.pitman(ab, c(-1, 0.1)), "init_vals")
  expect_error(optim.pitman(ab, c(10, 1.2)), "init_vals")
  expect_error(optim.pitman(ab, c(10, 0)), "init_vals")
  expect_error(optim.pitman(ab, 10), "init_vals")
  expect_error(optim.pitman(ab, c(NA, 0.1)), "init_vals")
  expect_error(optim.pitman(ab, lower = c(5, 0.1), upper = c(1, 0.5)), "lower")
  expect_error(optim.pitman(10), "at least 2 species")
})

test_that("optim.pitman warns when an estimate sits on a bound", {
  ## all singletons: both parameters run to the edge of their interval
  expect_warning(expect_warning(optim.pitman(rep(1, 50)), "`theta`.*upper bound"),
                 "`sigma`.*upper bound")
  expect_warning(optim.pitman(c(5, 3, 2, 2, 1, 1, 1)), "lower bound of the search interval")
  expect_no_warning(optim.pitman(ab))
})

## ---- optim.multideme ----------------------------------------------------------------

m2 <- rbind(c(44, 37, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1),
            c(240, 20, 48, 2, 21, 1, 3, 2, 5, 2, 0, 1, 1))

test_that("optim.multideme estimates the immigration rates", {
  fit <- optim.multideme(m2)
  expect_named(fit, c("I", "m", "J", "k", "logl", "se_I", "se_m"))
  expect_equal(fit$I, c(69.444824, 197.921875), tolerance = 1e-4)
  expect_equal(fit$m, c(0.33638443, 0.36454946), tolerance = 1e-4)
  expect_equal(fit$J, c(138, 346))
  expect_equal(fit$k, c(13, 12))
  expect_equal(fit$m, fit$I / (fit$J - 1 + fit$I))
})

test_that("optim.multideme finds the global maximum of each deme likelihood", {
  fit <- optim.multideme(m2)
  x <- colSums(m2) / sum(m2)
  for (j in 1:2) {
    n <- m2[j, ]; pr <- n > 0
    ll <- function(I) sum(lgamma(I * x[pr] + n[pr]) - lgamma(I * x[pr])) - lgamma(I + sum(n)) + lgamma(I)
    grid <- exp(seq(log(1e-3), log(sum(n)), length.out = 4000))
    expect_equal(fit$I[j], grid[which.max(vapply(grid, ll, numeric(1)))], tolerance = 0.01)
    expect_equal(fit$logl[j], ll(fit$I[j]), tolerance = 1e-8)
  }
})

test_that("optim.multideme is silent unless asked, and reports progress with `verbose`", {
  expect_silent(optim.multideme(m2))
  expect_message(optim.multideme(m2, verbose = TRUE), "Deme 1 of 2")
  expect_message(optim.multideme(m2, verbose = TRUE), "Deme 2 of 2")
})

test_that("optim.multideme accepts data frames and validates its input", {
  expect_equal(optim.multideme(as.data.frame(m2)), optim.multideme(m2))
  expect_error(optim.multideme(m2 * -1), "non-negative whole numbers")
  m <- m2; m[1, 1] <- NA
  expect_error(optim.multideme(m), "non-negative whole numbers")
  expect_error(optim.multideme(m2 + 0.5), "non-negative whole numbers")
  expect_error(optim.multideme(matrix(0, 2, 3)), "no individuals")
  expect_error(optim.multideme(matrix(numeric(0), 0, 0)), "numeric matrix")
  expect_error(optim.multideme("a"), "numeric matrix")
  expect_error(optim.multideme(m2, I_max = -1), "I_max")
  expect_error(optim.multideme(m2, I_max = c(1, 2, 3)), "I_max")
})

test_that("degenerate demes are handled explicitly", {
  m <- rbind(c(5, 0, 0), c(1, 0, 0), c(30, 10, 4), c(2, 12, 9), c(0, 3, 20))
  expect_warning(fit <- optim.multideme(m), "deme 2 has a single individual")
  expect_true(is.na(fit$I[2]) && is.na(fit$m[2]))
  expect_equal(fit$I[1], 0)                                         # single species: MLE at I = 0
  expect_equal(fit$m[1], 0)
  expect_false(anyNA(fit$I[c(1, 3, 4, 5)]))
  expect_equal(fit$J, c(5, 1, 44, 23, 23))
  expect_warning(f0 <- optim.multideme(rbind(c(0, 0, 0), c(30, 10, 4), c(2, 12, 9), c(0, 3, 20))),
                 "deme 1 has no individual")
  expect_true(is.na(f0$I[1]))
})

test_that("optim.multideme warns at the search bound and I_max moves it", {
  one <- matrix(c(5, 3, 2), 1)          # a single deme is a pure sample of the 'region'
  expect_warning(f <- optim.multideme(one), "upper bound")
  expect_equal(f$I, 10, tolerance = 1e-3)                           # default bound: J
  f2 <- suppressWarnings(optim.multideme(one, I_max = 1e3))
  expect_gt(f2$I, 500)
})

test_that("optim.multideme standard errors use the exact curvature", {
  fit <- optim.multideme(m2)
  x <- colSums(m2) / sum(m2)
  for (j in 1:2) {
    n <- m2[j, ]; pr <- n > 0
    ll <- function(I) sum(lgamma(I * x[pr] + n[pr]) - lgamma(I * x[pr])) - lgamma(I + sum(n)) + lgamma(I)
    h <- 1e-4 * fit$I[j]
    d2 <- (ll(fit$I[j] + h) - 2 * ll(fit$I[j]) + ll(fit$I[j] - h)) / h^2
    expect_equal(fit$se_I[j], 1 / sqrt(-d2), tolerance = 1e-4)
  }
  expect_equal(fit$se_m, fit$se_I * (fit$J - 1) / (fit$J - 1 + fit$I)^2)
})

test_that("optim.multideme likelihood-ratio intervals", {
  fit <- optim.multideme(m2, ci = TRUE, I_max = 1e5)     # wide enough for interior limits
  expect_named(fit, c("I", "m", "J", "k", "logl", "se_I", "se_m", "ci_I", "ci_m"))
  x <- colSums(m2) / sum(m2)
  for (j in 1:2) {
    n <- m2[j, ]; pr <- n > 0
    ll <- function(I) sum(lgamma(I * x[pr] + n[pr]) - lgamma(I * x[pr])) - lgamma(I + sum(n)) + lgamma(I)
    expect_equal(fit$logl[j] - vapply(fit$ci_I[j, ], ll, numeric(1)),
                 rep(stats::qchisq(0.95, 1) / 2, 2), tolerance = 1e-6, ignore_attr = TRUE)
    expect_true(fit$ci_I[j, 1] < fit$I[j] && fit$I[j] < fit$ci_I[j, 2])
  }
  expect_equal(fit$ci_m, fit$ci_I / (fit$J - 1 + fit$ci_I))
  wide <- optim.multideme(m2, ci = TRUE, I_max = 1e5, level = 0.99)
  expect_true(all(wide$ci_I[, 1] < fit$ci_I[, 1] & wide$ci_I[, 2] > fit$ci_I[, 2]))
  ## with the default search interval [1e-6, J], the upper limit stops at the bound
  def <- optim.multideme(m2, ci = TRUE)
  expect_equal(def$ci_I[, "upper"], def$J)
  expect_error(optim.multideme(m2, ci = NA), "`ci`")
  expect_error(optim.multideme(m2, ci = TRUE, level = 1.5), "`level`")
})

test_that("optim.multideme handles degenerate demes in the uncertainty outputs", {
  m <- rbind(c(5, 0, 0), c(30, 10, 4), c(2, 12, 9), c(0, 3, 20))
  fit <- optim.multideme(m, ci = TRUE)
  expect_true(is.na(fit$se_I[1]))                       # single species: I = 0
  expect_equal(fit$ci_I[1, "lower"], 0, ignore_attr = TRUE)
  expect_gt(fit$ci_I[1, "upper"], 0)
})

test_that("optim.multideme reproduces the Jabot et al. (2008) estimates of TeTame 2.1", {
  ## trial dataset of the TeTame 2.1 manual (option 'j'); TeTame's Std_I is not
  ## reproduced: it lacks a factor (0.01 I) in its curvature
  tt <- rbind(c(1, 8, 80, 2, 1, 20, 0), c(15, 8, 8, 20, 10, 2, 0), c(15, 81, 80, 2, 1, 2, 0))
  fit <- optim.multideme(tt)
  expect_equal(fit$I, c(7.4197, 5.738, 14.6712), tolerance = 1e-4)
  expect_equal(fit$m, c(0.0626559, 0.0847088, 0.0753638), tolerance = 1e-5)
  const <- vapply(1:3, function(j) { n <- tt[j, tt[j, ] > 0]; lfactorial(sum(n)) - sum(lfactorial(n)) },
                  numeric(1))
  expect_equal(-(fit$logl + const), c(15.9501, 20.1273, 15.9371), tolerance = 1e-5)
})

