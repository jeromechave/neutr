zoo <- c(8, 5, 3, 2, 1, 1)                 # example of Etienne (2005), supplementary material

## exact sequential sampler of Etienne's model: individual j of the local sample is an
## immigrant with probability I / (I + j - 1); immigrants follow Hoppe's urn with theta
sim_etienne <- function(theta, m, J) {
  I <- m / (1 - m) * (J - 1)
  sp <- integer(J)
  imm <- integer(0)
  k <- 0L
  for (j in seq_len(J)) {
    if (stats::runif(1) < I / (I + j - 1)) {
      r <- length(imm)
      if (stats::runif(1) < theta / (theta + r)) {
        k <- k + 1L
        s <- k
      } else {
        s <- imm[sample.int(r, 1)]
      }
      imm <- c(imm, s)
      sp[j] <- s
    } else {
      sp[j] <- sp[sample.int(j - 1, 1)]
    }
  }
  tabulate(sp)
}

## all partitions of n (as abundance vectors)
partitions_of <- function(n, max_part = n) {
  if (n == 0) return(list(numeric(0)))
  do.call(c, lapply(seq_len(min(n, max_part)), function(i)
    lapply(partitions_of(n - i, i), function(p) c(i, p))))
}

## ---- logkda -------------------------------------------------------------------------

test_that("logkda reproduces the example of Etienne (2005)", {
  expect_equal(logkda(zoo), c(0, 1.970769, 3.322156, 4.324023, 5.073586, 5.613912, 5.963897,
                              6.128735, 6.103853, 5.876264, 5.423493, 4.709155, 3.672072,
                              2.197225, 0), tolerance = 1e-6)
})

test_that("logkda matches PARI/GP to numerical precision, beyond double-precision overflow", {
  ## computed with logkda.pari() of package untb 1.7; J = 2000 overflows doubles
  ref <- readRDS(test_path("fixtures", "pari_reference.rds"))
  for (r in ref) expect_equal(logkda(r$abundances), r$logkda, tolerance = 1e-12)
})

test_that("logkda of a single species gives the Stirling numbers of the first kind", {
  s1 <- c(120, 274, 225, 85, 15, 1)        # |s(6, a)|, a = 1..6
  expect_equal(exp(logkda(6)) * factorial(5) / factorial(0:5), s1)
})

test_that("logkda handles edge cases and does not depend on species order", {
  expect_equal(logkda(rep(1, 4)), 0)                       # all singletons
  expect_equal(logkda(1), 0)
  expect_equal(logkda(c(3, 0, 1)), logkda(c(3, 1)))
  expect_length(logkda(c(7, 3, 1)), 11 - 3 + 1)
  expect_equal(logkda(c(9, 4, 4, 2, 1)), logkda(c(1, 2, 4, 9, 4)))
  expect_error(logkda(c(0, 0)), "logkda\\(\\).*at least 1 species")
  expect_error(logkda(c(2.5, 1)), "logkda\\(\\).*whole numbers")
})

## ---- logl.etienne -----------------------------------------------------------------

test_that("logl.etienne matches Etienne's equation 6 evaluated directly", {
  naive <- function(n, theta, m) {          # log-probability of the partition
    J <- sum(n); k <- length(n); A <- k:J
    I <- m / (1 - m) * (J - 1)
    k * log(theta) + lgamma(theta) - lgamma(theta + J) + sum(lgamma(n)) +
      log(sum(exp(logkda(n) + lgamma(theta + J) - lgamma(theta + A) +
                    lgamma(I) - lgamma(I + J) + A * log(I))))
  }
  for (p in list(c(7, 0.2), c(30, 0.01), c(2, 0.9), c(0.5, 1e-4)))
    expect_equal(logl.etienne(zoo, p[1], p[2]), naive(zoo, p[1], p[2]))
})

test_that("with m = 1, logl.etienne is the log-likelihood of optim.ewens", {
  ew <- optim.ewens(zoo)
  expect_equal(logl.etienne(zoo, ew$theta, 1), ew$logl)
  ## the limit m -> 1 is continuous, without loss of precision
  d <- logl.etienne(zoo, 3, 1 - 10^-(4:12)) - logl.etienne(zoo, 3, 1)
  expect_true(all(d > 0))
  expect_equal(d / 10^-(4:12), rep(d[1] / 1e-4, 9), tolerance = 1e-3)
})

test_that("logl.etienne(full = TRUE) is a probability distribution over abundance vectors", {
  for (p in list(c(1.7, 0.4), c(10, 0.05))) {
    pr <- vapply(partitions_of(7), function(n) logl.etienne(n, p[1], p[2], full = TRUE), numeric(1))
    expect_equal(sum(exp(pr)), 1)
  }
})

test_that("logl.etienne is vectorised and accepts a precomputed log_kda", {
  lk <- logkda(zoo)
  expect_equal(logl.etienne(zoo, c(2, 7), c(0.1, 0.5)),
               c(logl.etienne(zoo, 2, 0.1), logl.etienne(zoo, 7, 0.5)))
  expect_equal(logl.etienne(zoo, 7, c(0.1, 0.5, 1)),
               logl.etienne(zoo, 7, c(0.1, 0.5, 1), log_kda = lk))
  expect_length(logl.etienne(zoo, 7, seq(0.1, 0.9, by = 0.1)), 9)
})

test_that("logl.etienne validates its arguments", {
  expect_error(logl.etienne(zoo, -1, 0.5), "logl.etienne\\(\\).*theta")
  expect_error(logl.etienne(zoo, 5, 0), "logl.etienne\\(\\).*`m`")
  expect_error(logl.etienne(zoo, 5, 1.2), "`m`")
  expect_error(logl.etienne(zoo, 5, NA), "`m`")
  expect_error(logl.etienne(zoo, 5, 0.5, log_kda = 1:3), "log_kda")
  expect_error(logl.etienne(c(1, -2), 5, 0.5), "non-negative")
})

test_that("the analytic gradient matches finite differences", {
  n <- readRDS(test_path("fixtures", "etienne_hard_cases.rds"))$lnsrch
  lk <- logkda(n)
  for (p in list(c(7, 3), c(30, 0.2), c(2, 300), c(0.5, 5000), c(40, 1e8))) {
    h <- 1e-5 * p
    num <- c((neutr:::etienne_logl(p[1] + h[1], p[2], n, lk) -
                neutr:::etienne_logl(p[1] - h[1], p[2], n, lk)) / (2 * h[1]),
             (neutr:::etienne_logl(p[1], p[2] + h[2], n, lk) -
                neutr:::etienne_logl(p[1], p[2] - h[2], n, lk)) / (2 * h[2]))
    expect_equal(neutr:::etienne_grad(p[1], p[2], n, lk), num, tolerance = 1e-5)
  }
})

## ---- optim.etienne ----------------------------------------------------------------

test_that("optim.etienne reproduces the estimates of Etienne (2005)", {
  fit <- optim.etienne(zoo)
  expect_named(fit, c("theta", "m", "I", "logl", "converged", "se_theta", "se_m", "se_I",
                      "local_maxima"))
  expect_true(fit$converged)
  expect_equal(fit$theta, 7.047958, tolerance = 1e-5)
  expect_equal(fit$m, 0.22635923, tolerance = 1e-5)
  expect_equal(fit$I, fit$m / (1 - fit$m) * (sum(zoo) - 1))
  expect_equal(fit$logl, logl.etienne(zoo, fit$theta, fit$m))
})

test_that("optim.etienne finds the global maximum when the likelihood is multimodal", {
  hard <- readRDS(test_path("fixtures", "etienne_hard_cases.rds"))
  ## bimodal surface: local maximum near (51, 0.19), global one near (214, 0.0137)
  fit <- optim.etienne(hard$bimodal)
  expect_equal(fit$theta, 214.144, tolerance = 1e-3)
  expect_equal(fit$m, 0.0136808, tolerance = 1e-3)
  lk <- logkda(hard$bimodal)
  grid <- expand.grid(theta = exp(seq(log(1), log(2000), length.out = 80)),
                      m = stats::plogis(seq(-10, 10, length.out = 80)))
  expect_gte(fit$logl, max(logl.etienne(hard$bimodal, grid$theta, grid$m, log_kda = lk)))
  ## narrow peak next to a ridge running to the theta bound: missed by a profile over
  ## log(I) alone, found by TeTame 2.1 at theta = 281.85, m = 0.0014431
  fit3 <- optim.etienne(hard$ridge)
  expect_equal(fit3$theta, 281.85, tolerance = 1e-3)
  expect_equal(fit3$m, 0.0014431, tolerance = 1e-3)
  expect_equal(-logl.etienne(hard$ridge, fit3$theta, fit3$m, full = TRUE), 164.6275, tolerance = 1e-6)
  ## the optimiser stops on a line-search failure at the optimum: still converged
  fit2 <- optim.etienne(hard$lnsrch)
  expect_true(fit2$converged)
  expect_equal(fit2$logl, -948.79894455, tolerance = 1e-9)
})

test_that("optim.etienne does not depend on the starting values", {
  n <- readRDS(test_path("fixtures", "etienne_hard_cases.rds"))$lnsrch
  a <- optim.etienne(n)
  for (s in list(c(1, 0.01), c(200, 0.9), c(0.1, 0.5))) {
    b <- optim.etienne(n, init_vals = s)
    expect_equal(b$logl, a$logl, tolerance = 1e-10)
    expect_equal(b$theta, a$theta, tolerance = 1e-4)
  }
})

test_that("the Ewens model is nested in the Etienne model: their log-likelihoods are comparable", {
  set.seed(3)
  for (p in list(c(10, 0.05), c(30, 0.5), c(5, 0.9))) {
    a <- sim_etienne(p[1], p[2], 1000)
    et <- suppressWarnings(optim.etienne(a))
    ew <- optim.ewens(a)
    expect_gte(et$logl, ew$logl - 1e-8)            # nested: LRT statistic >= 0
  }
})

test_that("optim.etienne estimates are consistent with the true parameters", {
  ## the true parameters lie within the 99.9% likelihood-ratio confidence region
  set.seed(11)
  for (p in list(c(5, 0.01), c(50, 0.1), c(20, 0.5))) {
    a <- sim_etienne(p[1], p[2], 3000)
    fit <- suppressWarnings(optim.etienne(a))
    expect_true(fit$converged)
    expect_lt(2 * (fit$logl - logl.etienne(a, p[1], p[2])), stats::qchisq(0.999, df = 2))
  }
})

test_that("optim.etienne warns when m is at its upper bound (no dispersal limitation)", {
  set.seed(5)
  a <- generate.hoppe.urn(20, 2000)$abundance          # Ewens: no dispersal limitation
  expect_warning(fit <- optim.etienne(a), "`m`.*upper bound")
  expect_equal(fit$m, 1, tolerance = 1e-6)
  expect_equal(fit$theta, optim.ewens(a)$theta, tolerance = 1e-3)
  expect_no_warning(optim.etienne(zoo))
  ## the example vector of neutr shows no evidence of dispersal limitation
  ab <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
  expect_warning(fa <- optim.etienne(ab), "`m`.*upper bound")
  expect_equal(fa$logl, optim.ewens(ab)$logl, tolerance = 1e-8)
})

test_that("optim.etienne validates its arguments", {
  expect_error(optim.etienne(zoo, c(-1, 0.5)), "init_vals")
  expect_error(optim.etienne(zoo, c(10, 1)), "init_vals")
  expect_error(optim.etienne(zoo, c(10, 0)), "init_vals")
  expect_error(optim.etienne(zoo, 10), "init_vals")
  expect_error(optim.etienne(zoo, c(NA, 0.5)), "init_vals")
  expect_error(optim.etienne(zoo, lower = c(5, 0.1), upper = c(1, 0.5)), "lower")
  expect_error(optim.etienne(zoo, upper = c(100, 1)), "upper")
  expect_error(optim.etienne(10), "optim.etienne\\(\\).*at least 2 species")
})

## ---- TeTame 2.1: trial dataset, uncertainty and local maxima -------------------------

tt <- list(c(1, 8, 80, 2, 1, 20), c(15, 8, 8, 20, 10, 2), c(15, 81, 80, 2, 1, 2))

test_that("optim.etienne reproduces TeTame 2.1 on its trial dataset", {
  ## TeTame reports -log P(abundance vector), i.e. logl.etienne(full = TRUE)
  f2 <- optim.etienne(tt[[2]])
  expect_equal(c(f2$theta, f2$I), c(2.94777, 7.04503), tolerance = 1e-5)
  expect_equal(-logl.etienne(tt[[2]], f2$theta, f2$m, full = TRUE), 12.2275, tolerance = 1e-5)
  f3 <- optim.etienne(tt[[3]])
  expect_equal(c(f3$theta, f3$I), c(21.7474, 1.24466), tolerance = 1e-5)
  expect_equal(-logl.etienne(tt[[3]], f3$theta, f3$m, full = TRUE), 13.5365, tolerance = 1e-5)
  ## Ewens fits
  expect_equal(vapply(tt, function(y) optim.ewens(y)$theta, numeric(1)),
               c(1.19529, 1.43357, 1.05382), tolerance = 1e-5)
})

test_that("optim.etienne reports the local maxima reported by TeTame", {
  lm3 <- optim.etienne(tt[[3]])$local_maxima
  expect_equal(nrow(lm3), 2)
  expect_equal(c(lm3$theta[2], lm3$m[2]), c(1.0563, 0.980993), tolerance = 1e-4)  # TeTame Theta2, m2
  expect_false(any(lm3$boundary))
  expect_true(lm3$logl[1] > lm3$logl[2])
  expect_equal(nrow(optim.etienne(tt[[2]])$local_maxima), 1)
  ## a flat ridge: both ends reported, with the same likelihood
  lm1 <- suppressWarnings(optim.etienne(tt[[1]]))$local_maxima
  expect_equal(nrow(lm1), 2)
  expect_true(all(lm1$boundary))
  expect_equal(lm1$logl[1], lm1$logl[2], tolerance = 1e-6)
})

test_that("optim.etienne standard errors come from the inverse Hessian", {
  y <- tt[[2]]
  fit <- optim.etienne(y)
  J <- sum(y)
  f <- function(th, I) logl.etienne(y, th, I / (I + J - 1))
  h <- 1e-4; th <- fit$theta; I <- fit$I
  H <- matrix(0, 2, 2)
  H[1, 1] <- (f(th * (1 + h), I) - 2 * f(th, I) + f(th * (1 - h), I)) / (h * th)^2
  H[2, 2] <- (f(th, I * (1 + h)) - 2 * f(th, I) + f(th, I * (1 - h))) / (h * I)^2
  H[1, 2] <- H[2, 1] <- (f(th * (1 + h), I * (1 + h)) - f(th * (1 + h), I * (1 - h)) -
                           f(th * (1 - h), I * (1 + h)) + f(th * (1 - h), I * (1 - h))) / (4 * h^2 * th * I)
  V <- solve(-H)
  expect_equal(c(fit$se_theta, fit$se_I), sqrt(diag(V)), tolerance = 1e-4)
  expect_equal(fit$se_m, fit$se_I * (J - 1) / (J - 1 + fit$I)^2, tolerance = 1e-8)
  ## undefined on a bound
  ab <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
  fb <- suppressWarnings(optim.etienne(ab))
  expect_true(is.na(fb$se_theta) && is.na(fb$se_m) && is.na(fb$se_I))
})

test_that("optim.etienne likelihood-ratio intervals", {
  fit <- optim.etienne(zoo, ci = TRUE)
  expect_true(all(c("ci_theta", "ci_m", "ci_I") %in% names(fit)))
  lk <- logkda(zoo)
  prof_theta <- function(th) {
    max(vapply(stats::plogis(seq(-12, 12, length.out = 400)),
               function(m) logl.etienne(zoo, th, m, log_kda = lk), numeric(1)))
  }
  half <- stats::qchisq(0.95, 1) / 2
  expect_equal(fit$logl - prof_theta(fit$ci_theta[1]), half, tolerance = 0.02)
  expect_true(fit$ci_theta[1] < fit$theta && fit$theta < fit$ci_theta[2])
  expect_true(fit$ci_m[1] < fit$m && fit$m <= fit$ci_m[2])
  expect_equal(fit$ci_m, fit$ci_I / (sum(zoo) - 1 + fit$ci_I))
  wide <- optim.etienne(zoo, ci = TRUE, level = 0.99)
  expect_lt(wide$ci_theta[1], fit$ci_theta[1])
  expect_error(optim.etienne(zoo, ci = "yes"), "`ci`")
  expect_error(optim.etienne(zoo, level = 0), "`level`")
})

