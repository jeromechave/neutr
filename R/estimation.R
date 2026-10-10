# Maximum likelihood estimation ---------------------------------------------------

#' Maximum likelihood estimation of \eqn{\theta} based on Ewens sampling formula
#'
#' @description `optim.ewens()` computes the maximum likelihood estimate of Ewens
#' \eqn{\theta} parameter for a given species abundance distribution
#' \eqn{(n_1, ..., n_k)} where \eqn{n_i} is the number of organisms of species i in
#' the sample and k is the total number of species in the sample (for all
#' \eqn{i, n_i > 0}). The input should be a vector; species with zero abundance
#' are ignored.
#'
#' The maximum likelihood estimate is the root of
#' \eqn{k = \theta[\psi(\theta+J)-\psi(\theta)]} (the expected number of species
#' equals the observed one), where \eqn{\psi} is the digamma function and
#' \eqn{J=\sum n_i}. It is found by a one-dimensional root search to full numerical
#' precision, so no general-purpose optimiser is involved.
#'
#' @section Boundary cases:
#' The likelihood has no interior maximum when all individuals belong to a single
#' species (`k = 1`, the estimate is `theta = 0`) or when every species has a single
#' individual (`k = J`, the estimate is `theta = Inf`). A warning is issued and the
#' boundary value is returned.
#'
#' @section Log-likelihood:
#' `logl` is the log-probability of the partition (the exchangeable partition
#' probability function of Ewens sampling formula). Terms that depend on the data
#' only, and not on the model, are omitted. It is on the same scale as the `logl`
#' returned by [optim.pitman()]. The Ewens model is the Pitman model with
#' \eqn{\sigma = 0}, so the two models are nested and can be compared directly,
#' for instance with the likelihood-ratio statistic
#' \eqn{2(\mathrm{logl}_{Pitman} - \mathrm{logl}_{Ewens})} or with AIC.
#'
#' @param input_abundances a vector of integers (species abundances)
#'
#' @return A list with
#' \describe{
#'   \item{theta}{the \eqn{\theta} parameter of Ewens sampling formula, a real value}
#'   \item{logl}{the maximal value of the log-likelihood (see Details), a real value;
#'     it is comparable with the `logl` value of [optim.pitman()]}
#' }
#' @seealso [optim.pitman()], [optim.multideme()]
#' @export
#'
#' @examples
#' input_abundances <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
#' optim.ewens(input_abundances)
optim.ewens <- function(input_abundances) {
  fn <- "optim.ewens"
  n <- check_abundances(input_abundances, fn)
  J <- sum(n)
  k <- length(n)

  logl_at <- function(theta) {
    k * log(theta) + lgamma(theta) - lgamma(theta + J) + sum(lgamma(n))
  }

  if (k == 1L) {
    warning(fn, "(): all individuals belong to a single species; the likelihood is ",
            "maximal at the boundary theta = 0.", call. = FALSE)
    theta <- 0
    logl <- 0                           # the partition has probability 1 as theta -> 0
  } else if (k == J) {
    warning(fn, "(): every species has a single individual; the likelihood increases ",
            "with theta without bound (estimate: Inf).", call. = FALSE)
    theta <- Inf
    logl <- 0                           # probability -> 1 as theta -> Inf
  } else {
    ## score equation, on the log(theta) scale; increasing from 1 - k < 0 to J - k > 0
    score <- function(lt) {
      th <- exp(lt)
      th * (digamma(th + J) - digamma(th)) - k
    }
    lower <- log(1e-10)
    upper <- log(max(J, 10))
    while (score(upper) <= 0 && upper < 600) upper <- upper + log(10)
    if (score(upper) <= 0) {
      warning(fn, "(): the estimate of `theta` is numerically unbounded; returning Inf.",
              call. = FALSE)
      theta <- Inf
      logl <- 0
    } else {
      theta <- exp(stats::uniroot(score, c(lower, upper), tol = 1e-13)$root)
      logl <- logl_at(theta)
    }
  }
  list(theta = theta, logl = logl)
}

#' Maximum likelihood estimation of local immigration rates \eqn{(m_1,...,m_K)} based on the K-deme sampling formula
#'
#' @description `optim.multideme()` computes the maximum likelihood estimate for the
#' multi-deme parameters for a given species abundance matrix
#' \eqn{(n_{1j}, ..., n_{kj})} for deme j (\eqn{j\in\{1,...,K\}})
#' where \eqn{n_{ij}} is the number of organisms of species i in deme j
#' and k is the total number of species (for all \eqn{i, n_i \geq 0}).
#' The formulation in this description is that of species abundance distributions
#' but is valid for any partition in a subdivided setting (e.g., word frequencies in multiple books).
#'
#' **The matrix must have one row per deme and one column per species.** The regional
#' species abundance distribution is obtained by summing the demes (rows).
#'
#' @section Boundary and degenerate cases:
#' * A deme with a single individual carries no information on the immigration
#'   rate: `I` and `m` are `NA` and a warning is issued.
#' * A deme containing a single species has its likelihood maximal at `I = 0`
#'   (`m = 0`); this value is returned as is.
#' * The search interval for `I` is `[1e-6, I_max]` with `I_max` equal to the
#'   deme size `J` by default. If an estimate lies on the
#'   edge of this interval (for instance a deme that is essentially a random sample
#'   of the regional pool, for which the likelihood keeps increasing with `I`), a
#'   warning is issued; increase `I_max` to explore larger values.
#'
#' @section Uncertainty:
#' The standard errors `se_I` and `se_m` are obtained from the exact second
#' derivative of the log-likelihood in `I` (and the delta method for `m`); they are
#' `NA` when the estimate lies on a bound of the search interval or for a deme with
#' a single species. With `ci = TRUE`, likelihood-ratio confidence intervals are
#' also returned: the range of `I` over which the log-likelihood stays within
#' `qchisq(level, 1) / 2` of its maximum. A limit equal to a bound of the search
#' interval means that the interval extends at least to that bound. These replace
#' the `Std_I` and `Std_m` columns of the TeTame software.
#'
#' @param input_abundance_matrix a matrix \eqn{n_{ij}} of abundances (one row per
#'   deme (local site) j, one column per species i). A data frame of counts is accepted.
#' @param verbose if `TRUE`, report progress over the demes with [message()]
#'   (default `FALSE`)
#' @param I_max upper bound of the search interval for `I`: `NULL` (default, the
#'   deme size `J`), or a positive number, or one number per deme
#' @param ci if `TRUE`, also return likelihood-ratio confidence intervals (default
#'   `FALSE`)
#' @param level confidence level of the intervals (default 0.95)
#'
#' @return A list with
#' \describe{
#'   \item{I}{rescaled immigration rate, \eqn{I=\frac{m}{1-m}(J-1)}, a vector of real numbers}
#'   \item{m}{immigration rate, a vector of real numbers}
#'   \item{J}{local community size, a vector of integers}
#'   \item{k}{number of classes/species, a vector of integers}
#'   \item{logl}{maximal log-likelihood of each deme (up to a constant that does not
#'     depend on `I`), a vector of real numbers}
#'   \item{se_I, se_m}{standard errors of `I` and `m`, vectors of real numbers}
#'   \item{ci_I, ci_m}{only if `ci = TRUE`: confidence intervals of `I` and `m`,
#'     matrices with one row per deme and columns `lower` and `upper`}
#' }
#' @seealso [optim.ewens()], [optim.pitman()], [optim.etienne()]
#' @export
#'
#' @examples
#' input_abundances1 <- c(44, 37, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
#' input_abundances2 <- c(240, 20, 48, 2, 21, 1, 3, 2, 5, 2, 0, 1, 1)
#' input_abundance_matrix <- rbind(input_abundances1, input_abundances2)
#' optim.multideme(input_abundance_matrix)
optim.multideme <- function(input_abundance_matrix, verbose = FALSE, I_max = NULL,
                            ci = FALSE, level = 0.95) {
  fn <- "optim.multideme"
  mat <- tryCatch(as.matrix(input_abundance_matrix), error = function(e) NULL)
  if (is.null(mat) || !is.numeric(mat) || length(mat) == 0L || length(dim(mat)) != 2L)
    stop(fn, "(): `input_abundance_matrix` must be a numeric matrix (rows = demes, ",
         "columns = species).", call. = FALSE)
  if (anyNA(mat) || any(!is.finite(mat)) || any(mat < 0) ||
      any(abs(mat - round(mat)) > sqrt(.Machine$double.eps)))
    stop(fn, "(): `input_abundance_matrix` must contain non-negative whole numbers.",
         call. = FALSE)
  storage.mode(mat) <- "double"

  demes <- nrow(mat)
  if (!is.null(I_max)) {
    if (!is.numeric(I_max) || !(length(I_max) %in% c(1L, demes)) ||
        any(!is.finite(I_max)) || any(I_max <= 0))
      stop(fn, "(): `I_max` must be NULL, or positive numbers (one, or one per deme).",
           call. = FALSE)
    I_max <- rep_len(I_max, demes)
  }

  check_ci_args(ci, level, fn)

  regional <- colSums(mat)                     # regional abundance of each species
  if (sum(regional) == 0) stop(fn, "(): the matrix contains no individuals.", call. = FALSE)
  x <- regional / sum(regional)                # regional species abundance distribution

  I <- rep(NA_real_, demes)
  J <- rep(0, demes)
  k <- rep(0, demes)
  logl <- rep(NA_real_, demes)
  se_I <- rep(NA_real_, demes)
  ci_I <- matrix(NA_real_, demes, 2L, dimnames = list(NULL, c("lower", "upper")))
  I_lower <- 1e-6

  for (j in seq_len(demes)) {
    if (verbose) message(sprintf("Deme %d of %d ...", j, demes))
    local_abundance <- mat[j, ]
    present <- which(local_abundance != 0)
    n1 <- local_abundance[present]
    x1 <- x[present]
    J0 <- sum(n1)
    J[j] <- J0
    k[j] <- length(present)

    ## log-likelihood of the Dirichlet-multinomial (K-deme sampling formula)
    logL <- function(I) {
      sum(lgamma(I * x1 + n1) - lgamma(I * x1)) - lgamma(I + J0) + lgamma(I)
    }

    if (J0 <= 1) {
      warning(fn, "(): deme ", j, " has ", if (J0 == 0) "no individual" else "a single individual",
              ": the immigration rate cannot be estimated (NA).", call. = FALSE)
      next
    }
    if (k[j] == 1L) {
      I[j] <- 0                                # likelihood maximal at I = 0
      logl[j] <- logL(I_lower)
      if (ci) {
        upper <- if (is.null(I_max)) J0 else I_max[j]
        target <- -logl[j] + stats::qchisq(level, df = 1) / 2
        ci_I[j, ] <- c(0, exp(lr_interval(function(lI) -logL(exp(lI)), log(I_lower),
                                          log(I_lower), log(max(upper, I_lower)), target)[2]))
      }
      next
    }

    upper <- if (is.null(I_max)) J0 else I_max[j]
    if (upper <= I_lower) {
      I[j] <- upper
      logl[j] <- logL(upper)
      next
    }
    ## coarse scan on the log scale (guards against local maxima), then refinement
    grid <- seq(log(I_lower), log(upper), length.out = 41L)
    vals <- vapply(exp(grid), logL, numeric(1))
    best <- which.max(vals)
    lo <- grid[max(best - 1L, 1L)]
    hi <- grid[min(best + 1L, length(grid))]
    opt <- stats::optimize(function(lI) logL(exp(lI)), c(lo, hi), maximum = TRUE,
                           tol = 1e-10)
    I[j] <- exp(opt$maximum)
    logl[j] <- opt$objective
    on_bound <- warn_boundary(I[j], I_lower, upper, sprintf("I[%d]", j), fn, log_scale = TRUE)

    ## standard error from the exact second derivative of the log-likelihood in I
    if (!on_bound) {
      d2 <- sum(x1^2 * (trigamma(I[j] * x1 + n1) - trigamma(I[j] * x1))) -
        trigamma(I[j] + J0) + trigamma(I[j])
      if (is.finite(d2) && d2 < 0) se_I[j] <- 1 / sqrt(-d2)
    }
    if (ci) {
      target <- -logl[j] + stats::qchisq(level, df = 1) / 2
      ci_I[j, ] <- exp(lr_interval(function(lI) -logL(exp(lI)), log(I[j]),
                                   log(I_lower), log(upper), target))
    }
  }
  m <- I / (J - 1 + I)
  se_m <- se_I * (J - 1) / (J - 1 + I)^2

  out <- list(I = I, m = m, J = J, k = k, logl = logl, se_I = se_I, se_m = se_m)
  if (ci) {
    out$ci_I <- ci_I
    out$ci_m <- ci_I / (J - 1 + ci_I)
  }
  out
}

#' Maximum likelihood estimation of \eqn{\theta} and \eqn{\sigma} based on Pitman sampling formula
#'
#' @description `optim.pitman()` computes the maximum likelihood estimate of Pitman
#' \eqn{(\theta,\sigma)} parameters for a given species abundance distribution
#' \eqn{(n_1, ..., n_k)} where \eqn{n_i} is the number of organisms of species i in
#' the sample and k is the total number of species (for all \eqn{i, n_i > 0}).
#' The input should be a vector and initial values of \eqn{(\theta,\sigma)}
#' can be supplied with the argument `init_vals = c(theta_init, sigma_init)`
#' with the condition `theta_init > 0` and `0 < sigma_init < 1`.
#'
#' The log-likelihood is maximised with the box-constrained quasi-Newton method
#' L-BFGS-B ([stats::optim()]) using the exact gradient, on the scale
#' \eqn{(\log\theta,\sigma)}. The optimisation is run from `init_vals` and from a
#' second starting point derived from the Ewens estimate, and the better result is
#' kept.
#'
#' @section Search interval and boundary warnings:
#' By default \eqn{\theta} is searched in `[1e-6, J]` and \eqn{\sigma} in
#' `[1e-8, 1 - 1e-8]`, where `J` is the number of individuals. A warning is issued
#' if an estimate lies on the edge of its interval (the likelihood then has no
#' interior maximum, or the interval is too narrow); use `lower` and `upper`
#' to change it.
#'
#' @param input_abundances a vector of integers (species abundances); species with
#'   zero abundance are ignored
#' @param init_vals initial values of (theta, sigma), a vector of two real numbers
#' @param lower,upper bounds of the search interval for `c(theta, sigma)`
#'
#' @return A list with
#' \describe{
#'   \item{theta}{parameter \eqn{\theta}, a real value}
#'   \item{sigma}{parameter \eqn{\sigma}, a real value}
#'   \item{logl}{maximal log-likelihood (log-probability of the partition, see
#'     [optim.ewens()]), a real value; it is comparable with the `logl` value of
#'     [optim.ewens()]}
#'   \item{converged}{`TRUE` if the optimiser reported convergence}
#' }
#' @seealso [optim.ewens()], [optim.multideme()]
#' @export
#'
#' @examples
#' optim.pitman(c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1), c(10.0, 0.1))
optim.pitman <- function(input_abundances, init_vals = c(10.0, 0.1),
                         lower = c(1e-6, 1e-8), upper = NULL) {
  fn <- "optim.pitman"
  n <- check_abundances(input_abundances, fn, min_species = 2L)
  J <- sum(n)
  k <- length(n)
  if (is.null(upper)) upper <- c(max(J, 2e-6), 1 - 1e-8)

  if (!is.numeric(init_vals) || length(init_vals) != 2L || any(!is.finite(init_vals)))
    stop(fn, "(): `init_vals` must be two finite numbers, c(theta, sigma).", call. = FALSE)
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 2L || length(upper) != 2L ||
      any(!is.finite(lower)) || any(!is.finite(upper)) || any(lower >= upper) ||
      lower[1] <= 0 || lower[2] <= 0 || upper[2] >= 1)
    stop(fn, "(): `lower` and `upper` must be c(theta, sigma) with 0 < lower < upper ",
         "and upper sigma < 1.", call. = FALSE)
  if (init_vals[1] <= 0 || init_vals[2] <= 0 || init_vals[2] >= 1)
    stop(fn, "(): `init_vals` must satisfy theta > 0 and 0 < sigma < 1.", call. = FALSE)
  init_vals <- pmin(pmax(init_vals, lower), upper)

  ## Data summarised by distinct abundance values (many species share the same value)
  u <- sort(unique(n))
  cnt <- tabulate(match(n, u))

  ## minimise -logL over x = (log theta, sigma)
  f <- function(x) {
    v <- -pitman_logl(exp(x[1]), x[2], u, cnt, J, k)
    if (is.finite(v)) v else 1e300
  }
  g <- function(x) {
    gr <- -pitman_grad(exp(x[1]), x[2], u, cnt, J, k)
    gr[1] <- gr[1] * exp(x[1])
    if (all(is.finite(gr))) gr else c(0, 0)
  }
  lo <- c(log(lower[1]), lower[2])
  hi <- c(log(upper[1]), upper[2])

  run <- function(start) {
    stats::optim(c(log(start[1]), start[2]), f, g, method = "L-BFGS-B", lower = lo,
                 upper = hi, control = list(maxit = 1000, factr = 10, pgtol = 0))
  }
  fits <- list(run(init_vals))
  ## second start: Ewens estimate for theta, small sigma
  ew <- suppressWarnings(optim.ewens(n))
  if (is.finite(ew$theta) && ew$theta > 0) {
    fits[[2]] <- run(pmin(pmax(c(ew$theta, 0.02), lower), upper))
  }
  ## keep the best run; among runs that tie (to numerical precision) prefer one that
  ## the optimiser reports as converged
  vals <- vapply(fits, function(o) o$value, numeric(1))
  ok <- vapply(fits, function(o) o$convergence == 0L, logical(1))
  tied <- which(ok & vals <= min(vals) + 1e-8 * (1 + abs(min(vals))))
  best <- if (length(tied)) fits[[tied[which.min(vals[tied])]]] else fits[[which.min(vals)]]

  theta <- exp(best$par[1])
  sigma <- best$par[2]
  warn_boundary(theta, lower[1], upper[1], "theta", fn, log_scale = TRUE)
  warn_boundary(sigma, lower[2], upper[2], "sigma", fn)
  logl <- -best$value
  list(theta = theta, sigma = sigma, logl = logl,
       converged = best$convergence == 0L)
}
