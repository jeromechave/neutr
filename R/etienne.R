# Etienne sampling formula ---------------------------------------------------------
#
# Adapted, with improvements, from the etienne() and logkda() functions of package
# untb by Robin K. S. Hankin (Hankin 2007, Journal of Statistical Software 22(12)).

#' Etienne's coefficients \eqn{K(D,A)} of the sampling formula with dispersal limitation
#'
#' @description `logkda()` computes the logarithm of the coefficients \eqn{K(D,A)} of
#' Etienne's sampling formula (Etienne 2005, equation A11), for a species abundance
#' distribution \eqn{D = (n_1, ..., n_k)} and all \eqn{A} from k to \eqn{J=\sum n_i}.
#' These coefficients depend on the data only, not on the parameters, and are the
#' expensive part of [logl.etienne()] and [optim.etienne()].
#'
#' \eqn{K(D,A)} is the coefficient of \eqn{x^A} in the polynomial
#' \eqn{\prod_{i=1}^k T_{n_i}(x)}, where \eqn{T_n(x)=\sum_a T_n[a] x^a} with
#' \eqn{T_1[1]=1} and \eqn{T_n[a] = T_{n-1}[a] + T_{n-1}[a-1](a-1)/(n-1)}, i.e.
#' \eqn{T_n[a] = |s(n,a)|(a-1)!/(n-1)!} with \eqn{s} the Stirling numbers of the first
#' kind. The recursion and the polynomial products are evaluated entirely on the log
#' scale, so that the calculation cannot overflow whatever the sample size.
#'
#' @section Sample size:
#' The computing time of \eqn{K(D,A)} grows roughly as \eqn{J^2}: about 2 s for
#' \eqn{10^4} individuals and 30 s for \eqn{4 \times 10^4} on a typical computer.
#' Unlike the other estimators of the package, the Etienne model is therefore
#' practical up to about \eqn{10^5} individuals.
#'
#' @param input_abundances a vector of integers (species abundances); species with
#'   zero abundance are ignored
#'
#' @return \eqn{\log K(D,A)} for \eqn{A = k, k+1, ..., J}, a vector of `J - k + 1`
#'   real values (the first and the last are always 0)
#' @references Etienne, R. S. (2005) A new sampling formula for neutral biodiversity.
#'   *Ecology Letters* 8: 253-260.
#' @seealso [logl.etienne()], [optim.etienne()]
#' @export
#'
#' @examples
#' logkda(c(8, 5, 3, 2, 1, 1))
#'
#' ## a single species of abundance n gives the Stirling numbers |s(n, a)|
#' exp(logkda(5)) * factorial(4) / factorial(0:4)
logkda <- function(input_abundances) {
  logkda_core(check_abundances(input_abundances, "logkda"))
}

#' Log-likelihood of Etienne's sampling formula
#'
#' @description `logl.etienne()` computes the log-likelihood of Etienne's sampling
#' formula for a species abundance distribution \eqn{(n_1, ..., n_k)}, given the
#' fundamental biodiversity number \eqn{\theta} of the metacommunity and the
#' immigration rate \eqn{m} into the local community. With the rescaled immigration
#' rate \eqn{I=\frac{m}{1-m}(J-1)} and \eqn{J=\sum n_i}, the log-likelihood is
#' \deqn{k\log\theta + \log\Gamma(\theta) + \sum_i\log\Gamma(n_i) + \log\Gamma(I) -
#'   \log\Gamma(I+J) + \log\sum_{A=k}^{J} K(D,A)\frac{I^A}{\Gamma(\theta+A)}}
#' with \eqn{K(D,A)} given by [logkda()]. The function is vectorised over `theta`
#' and `m`.
#'
#' @section Log-likelihood:
#' By default `logl.etienne()` returns the log-probability of the partition, on the
#' same scale as the `logl` of [optim.ewens()] and [optim.pitman()]. For \eqn{m = 1}
#' (no dispersal limitation) Etienne's formula reduces to Ewens sampling formula and
#' the two log-likelihoods are identical, so the Ewens model is nested in the
#' Etienne model and the two can be compared by likelihood ratio or AIC.
#'
#' With `full = TRUE` the log-probability of the species abundance distribution
#' itself (Etienne 2005, equation 6) is returned instead: it adds the term
#' \eqn{\log J! - \sum_i \log n_i! - \sum_j \log \Phi_j!}, where \eqn{\Phi_j} is the
#' number of species with \eqn{j} individuals, which depends on the data only.
#'
#' @section Acknowledgement:
#' This function is adapted from the `etienne()` function of package `untb` by
#' Robin K. S. Hankin (Hankin 2007). The present version computes the coefficients
#' \eqn{K(D,A)} on the log scale in pure R (no PARI/GP), evaluates the likelihood
#' without loss of precision when \eqn{m} is close to 1, and puts it on the same
#' scale as the other estimators of `neutr`.
#'
#' @param input_abundances a vector of integers (species abundances); species with
#'   zero abundance are ignored
#' @param theta fundamental biodiversity number \eqn{\theta}, positive real values
#' @param m immigration rate, real values in \eqn{(0, 1]}
#' @param log_kda optional, the output of `logkda(input_abundances)`. It depends on
#'   the data only and is the expensive part of the calculation: supply it when the
#'   likelihood is evaluated many times on the same data
#' @param full if `TRUE`, return the log-probability of the species abundance
#'   distribution rather than of the partition (see Details); default `FALSE`
#'
#' @return the log-likelihood, a vector of real values (one per value of `theta`
#'   and `m`, recycled to a common length)
#' @references Etienne, R. S. (2005) A new sampling formula for neutral biodiversity.
#'   *Ecology Letters* 8: 253-260.
#'
#' Hankin, R. K. S. (2007) Introducing untb, an R package for simulating ecological
#'   drift under the unified neutral theory of biodiversity. *Journal of Statistical
#'   Software* 22(12).
#' @seealso [optim.etienne()], [logkda()], [optim.ewens()]
#' @export
#'
#' @examples
#' zoo <- c(8, 5, 3, 2, 1, 1)
#' logl.etienne(zoo, theta = 7, m = 0.2)
#'
#' ## profile over m, computing K(D,A) once
#' lk <- logkda(zoo)
#' logl.etienne(zoo, theta = 7, m = c(0.01, 0.1, 0.2, 0.5, 1), log_kda = lk)
#'
#' ## m = 1 is Ewens sampling formula
#' all.equal(logl.etienne(zoo, theta = 3, m = 1),
#'           length(zoo) * log(3) + lgamma(3) - lgamma(3 + sum(zoo)) + sum(lgamma(zoo)))
logl.etienne <- function(input_abundances, theta, m, log_kda = NULL, full = FALSE) {
  fn <- "logl.etienne"
  n <- check_abundances(input_abundances, fn)
  J <- sum(n)
  k <- length(n)
  if (!is.numeric(theta) || length(theta) == 0L || any(!is.finite(theta)) || any(theta <= 0))
    stop(fn, "(): `theta` must be positive finite numbers.", call. = FALSE)
  if (!is.numeric(m) || length(m) == 0L || anyNA(m) || any(m <= 0) || any(m > 1))
    stop(fn, "(): `m` must be numbers in (0, 1].", call. = FALSE)
  log_kda <- check_log_kda(log_kda, n, fn)

  len <- max(length(theta), length(m))
  theta <- rep_len(theta, len)
  m <- rep_len(m, len)
  I <- m / (1 - m) * (J - 1)
  logl <- vapply(seq_len(len), function(i) etienne_logl(theta[i], I[i], n, log_kda),
                 numeric(1))
  if (full) logl <- logl + lfactorial(J) - sum(lfactorial(n)) - sum(lfactorial(tabulate(n)))
  logl
}

#' Maximum likelihood estimation of \eqn{\theta} and \eqn{m} based on Etienne sampling formula
#'
#' @description `optim.etienne()` computes the maximum likelihood estimates of the
#' fundamental biodiversity number \eqn{\theta} and the immigration rate \eqn{m} of
#' Etienne's sampling formula (Etienne 2005) for a species abundance distribution
#' \eqn{(n_1, ..., n_k)} sampled from a single local community that receives
#' immigrants from a neutral metacommunity. The input should be a vector and
#' initial values of \eqn{(\theta, m)} can be supplied with the argument
#' `init_vals = c(theta_init, m_init)` with the condition `theta_init > 0` and
#' `0 < m_init < 1`.
#'
#' The coefficients \eqn{K(D,A)} are computed once with [logkda()]. The
#' log-likelihood is then maximised with the box-constrained quasi-Newton method
#' L-BFGS-B ([stats::optim()]) using the exact gradient, on the scale
#' \eqn{(\log\theta, \log I)} with \eqn{I=\frac{m}{1-m}(J-1)}. The optimisation is
#' The likelihood surface can be multimodal, with narrow peaks next to ridges that
#' run to the edge of the search domain, so it is first profiled in both directions:
#' on a grid of 61 values of \eqn{\log I} the likelihood is maximised over
#' \eqn{\theta}, and on a grid of 61 values of \eqn{\log\theta} it is maximised
#' over \eqn{I}. The optimisation is then run from every local maximum of the two
#' profiles, from `init_vals`, and from the Ewens estimate with \eqn{m} close to 1,
#' and the best result is kept.
#'
#' @section Search interval and boundary warnings:
#' By default \eqn{\theta} is searched in `[1e-6, 100 J]` and \eqn{m} in
#' `[1e-8, 1 - 1e-8]`, where `J` is the number of individuals. A warning is issued
#' if an estimate lies on the edge of its interval; use `lower` and `upper` to
#' change it. An estimate of \eqn{m} at its upper bound means that the data show no
#' evidence of dispersal limitation: the Ewens model ([optim.ewens()]) is then
#' adequate.
#'
#' @section Uncertainty:
#' The standard errors `se_theta`, `se_m` and `se_I` are obtained from the inverse
#' of the Hessian of the log-likelihood (on the scale \eqn{(\log\theta, \log I)},
#' then the delta method). They account for the correlation between the two
#' parameters, unlike the `Std_` columns of the TeTame software, and are `NA` when an
#' estimate lies on a bound of the search interval. As the likelihood is often
#' strongly asymmetric, the likelihood-ratio confidence intervals returned with
#' `ci = TRUE` are preferable: for each parameter, the range over which the profile
#' log-likelihood (maximised over the other parameter) stays within
#' `qchisq(level, 1) / 2` of its maximum. A limit equal to a bound of the search
#' interval means that the interval extends at least to that bound. When the
#' likelihood runs along a ridge towards \eqn{m = 1}, \eqn{\theta} and \eqn{m} are
#' weakly identified, the chi-squared approximation is conservative, and the
#' intervals are wider than their nominal level (in 200 simulations with
#' \eqn{\theta = 20}, \eqn{m = 0.1}, \eqn{J = 1000}, 95% intervals covered the true
#' values in 99.5% of the replicates).
#'
#' @section Local maxima:
#' `local_maxima` lists the distinct local maxima of the likelihood found by the
#' optimisation (TeTame reports a second maximum in the same way). The first row is
#' the global maximum. Points on a bound of the search interval are flagged
#' (`boundary`); they are kept only when they are local maxima of the profile
#' likelihood there. Two flagged points with the same likelihood typically mark the
#' two ends of a flat ridge, along which the data cannot discriminate between
#' \eqn{\theta} and \eqn{m}.
#'
#' @section Log-likelihood:
#' `logl` is on the same scale as the `logl` of [optim.ewens()] and
#' [optim.pitman()] (see [logl.etienne()]). The Ewens model is the Etienne model
#' with \eqn{m = 1}, so the two models are nested and can be compared by the
#' likelihood-ratio statistic \eqn{2(\mathrm{logl}_{Etienne} - \mathrm{logl}_{Ewens})}
#' or by AIC.
#'
#' @section Sample size:
#' The computing time of \eqn{K(D,A)} grows roughly as \eqn{J^2}: about 2 s for
#' \eqn{10^4} individuals and 30 s for \eqn{4 \times 10^4} on a typical computer.
#' Unlike the other estimators of the package, the Etienne model is therefore
#' practical up to about \eqn{10^5} individuals.
#'
#' @param input_abundances a vector of integers (species abundances); species with
#'   zero abundance are ignored
#' @param init_vals initial values of (theta, m), a vector of two real numbers
#' @param lower,upper bounds of the search interval for `c(theta, m)`
#' @param ci if `TRUE`, also return likelihood-ratio confidence intervals (default
#'   `FALSE`)
#' @param level confidence level of the intervals (default 0.95)
#'
#' @return A list with
#' \describe{
#'   \item{theta}{fundamental biodiversity number \eqn{\theta}, a real value}
#'   \item{m}{immigration rate, a real value}
#'   \item{I}{rescaled immigration rate, \eqn{I=\frac{m}{1-m}(J-1)}, a real value}
#'   \item{logl}{maximal log-likelihood, a real value; it is comparable with the
#'     `logl` value of [optim.ewens()] and [optim.pitman()]}
#'   \item{converged}{`TRUE` if the optimiser reported convergence, or stopped on a
#'     line-search failure at a point where the (projected) gradient vanishes}
#'   \item{se_theta, se_m, se_I}{standard errors of `theta`, `m` and `I`, real values}
#'   \item{local_maxima}{a data frame of the distinct local maxima, best first, with
#'     columns `theta`, `m`, `I`, `logl` and `boundary`}
#'   \item{ci_theta, ci_m, ci_I}{only if `ci = TRUE`: confidence intervals,
#'     vectors `c(lower, upper)`}
#' }
#' @references Etienne, R. S. (2005) A new sampling formula for neutral biodiversity.
#'   *Ecology Letters* 8: 253-260.
#' @seealso [logl.etienne()], [logkda()], [optim.ewens()], [optim.multideme()]
#' @export
#'
#' @examples
#' ## Example of Etienne (2005), supplementary material: theta = 7.05, m = 0.226
#' zoo <- c(8, 5, 3, 2, 1, 1)
#' et <- optim.etienne(zoo)
#' et
#'
#' ## Is there evidence of dispersal limitation? Compare with the nested Ewens model
#' ew <- optim.ewens(zoo)
#' 2 * (et$logl - ew$logl)          # likelihood-ratio statistic
#'
#' ## Confidence intervals; a sample with two local maxima (TeTame test data)
#' optim.etienne(zoo, ci = TRUE)[c("ci_theta", "ci_m")]
#' optim.etienne(c(15, 81, 80, 2, 1, 2))$local_maxima
optim.etienne <- function(input_abundances, init_vals = c(10, 0.5),
                          lower = c(1e-6, 1e-8), upper = NULL, ci = FALSE, level = 0.95) {
  fn <- "optim.etienne"
  n <- check_abundances(input_abundances, fn, min_species = 2L)
  J <- sum(n)
  k <- length(n)
  if (is.null(upper)) upper <- c(100 * J, 1 - 1e-8)

  if (!is.numeric(init_vals) || length(init_vals) != 2L || any(!is.finite(init_vals)))
    stop(fn, "(): `init_vals` must be two finite numbers, c(theta, m).", call. = FALSE)
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 2L || length(upper) != 2L ||
      any(!is.finite(lower)) || any(!is.finite(upper)) || any(lower >= upper) ||
      lower[1] <= 0 || lower[2] <= 0 || upper[2] >= 1)
    stop(fn, "(): `lower` and `upper` must be c(theta, m) with 0 < lower < upper ",
         "and upper m < 1.", call. = FALSE)
  if (init_vals[1] <= 0 || init_vals[2] <= 0 || init_vals[2] >= 1)
    stop(fn, "(): `init_vals` must satisfy theta > 0 and 0 < m < 1.", call. = FALSE)
  init_vals <- pmin(pmax(init_vals, lower), upper)
  check_ci_args(ci, level, fn)

  log_kda <- logkda_core(n)
  m_to_I <- function(m) m / (1 - m) * (J - 1)

  ## minimise -logL over x = (log theta, log I)
  f <- function(x) {
    v <- -etienne_logl(exp(x[1]), exp(x[2]), n, log_kda)
    if (is.finite(v)) v else 1e300
  }
  g <- function(x) {
    gr <- -etienne_grad(exp(x[1]), exp(x[2]), n, log_kda) * exp(x)
    if (all(is.finite(gr))) gr else c(0, 0)
  }
  lo <- c(log(lower[1]), log(m_to_I(lower[2])))
  hi <- c(log(upper[1]), log(m_to_I(upper[2])))

  refine <- function(x) {
    stats::optim(x, f, g, method = "L-BFGS-B", lower = lo, upper = hi,
                 control = list(maxit = 1000, factr = 10, pgtol = 0))
  }
  run <- function(start) refine(c(log(start[1]), log(m_to_I(start[2]))))

  ## The likelihood surface can be multimodal, with narrow peaks next to ridges that
  ## run to the edge of the search domain. Profile it in both directions: on a grid
  ## of log(I), maximise over theta, and on a grid of log(theta), maximise over I
  ## (cheap one-dimensional searches; a peak that is narrow in one direction is wide
  ## in the other). Then refine with L-BFGS-B from every local maximum of the two
  ## profiles, from `init_vals`, and from the Ewens estimate with m close to 1.
  profile_peaks <- function(grid, inner, other) {
    prof <- lapply(grid, function(z) stats::optimize(function(w) f(inner(z, w)), other, tol = 1e-4))
    pv <- vapply(prof, function(o) o$objective, numeric(1))
    pk <- which(pv <= c(Inf, pv[-length(pv)]) & pv <= c(pv[-1], Inf))
    lapply(pk, function(i) inner(grid[i], prof[[i]]$minimum))
  }
  starts <- c(profile_peaks(seq(lo[2], hi[2], length.out = 61L), function(lI, lt) c(lt, lI),
                            c(lo[1], hi[1])),
              profile_peaks(seq(lo[1], hi[1], length.out = 61L), function(lt, lI) c(lt, lI),
                            c(lo[2], hi[2])))
  fits <- lapply(starts, refine)
  fits[[length(fits) + 1L]] <- run(init_vals)
  ew <- suppressWarnings(optim.ewens(n))
  if (is.finite(ew$theta) && ew$theta > 0) {
    fits[[length(fits) + 1L]] <- run(pmin(pmax(c(ew$theta, 0.999), lower), upper))
  }
  ## keep the best run; among runs that tie (to numerical precision) prefer one that
  ## the optimiser reports as converged
  vals <- vapply(fits, function(o) o$value, numeric(1))
  ok <- vapply(fits, function(o) o$convergence == 0L, logical(1))
  tied <- which(ok & vals <= min(vals) + 1e-8 * (1 + abs(min(vals))))
  best <- if (length(tied)) fits[[tied[which.min(vals[tied])]]] else fits[[which.min(vals)]]
  ## L-BFGS-B with these tight tolerances often stops on a line-search failure at the
  ## optimum itself, when floating-point precision prevents any further decrease: such
  ## a run is counted as converged if the projected gradient vanishes there
  is_converged <- function(o) {
    gr <- g(o$par)
    gr[o$par <= lo & gr > 0] <- 0
    gr[o$par >= hi & gr < 0] <- 0
    o$convergence == 0L ||
      (isTRUE(grepl("LNSRCH", o$message)) && max(abs(gr)) < 1e-4 * (1 + abs(o$value)))
  }
  converged <- is_converged(best)

  theta <- exp(best$par[1])
  I <- exp(best$par[2])
  m <- I / (J - 1 + I)
  on_bound <- warn_boundary(theta, lower[1], upper[1], "theta", fn, log_scale = TRUE) |
    warn_boundary(m, lower[2], upper[2], "m", fn)

  ## standard errors from the inverse Hessian on the (log theta, log I) scale, then
  ## the delta method; undefined (NA) when an estimate is on a bound
  se <- se_from_gradient(g, best$par, on_bound)    # standard errors of log theta, log I
  se_theta <- theta * se[1]
  se_I <- I * se[2]
  se_m <- m * (1 - m) * se[2]

  ## distinct local maxima reached by the runs (the first is the global one). For a
  ## run that ends on, or near, a bound, two conditions are required: the gradient
  ## must point out of the domain there, and moving three log-units inwards (while
  ## re-optimising the other parameter) must not increase the likelihood. The first
  ## condition alone is not enough: as m -> 1 the likelihood flattens out, and a run
  ## can stop at a saddle point of the constrained problem, close to the bound.
  edge_of <- function(o) {                     # -1: near lower bound, 1: near upper, 0: interior
    th <- exp(o$par[1])
    mm <- exp(o$par[2]) / (J - 1 + exp(o$par[2]))
    near <- function(v, l, u, log_scale) {
      if (!at_bound(v, l, u, log_scale = log_scale)) return(0L)
      if (abs(if (log_scale) log(v / l) else v - l) < abs(if (log_scale) log(u / v) else u - v)) -1L else 1L
    }
    c(near(th, lower[1], upper[1], TRUE), near(mm, lower[2], upper[2], FALSE))
  }
  is_local_max <- function(o) {
    e <- edge_of(o)
    gr <- g(o$par)
    if (any(abs(gr[e == 0L]) >= 1e-4 * (1 + abs(o$value)))) return(FALSE)
    if (any(gr[e == -1L] < 0) || any(gr[e == 1L] > 0)) return(FALSE)
    for (i in which(e != 0L)) {
      probe <- o$par
      probe[i] <- probe[i] - 3 * e[i]
      j <- 3L - i                              # the other parameter
      val <- if (e[j] != 0L) {
        f(probe)
      } else {
        stats::optimize(function(w) f(replace(probe, j, w)),
                        c(max(lo[j], probe[j] - 2), min(hi[j], probe[j] + 2)), tol = 1e-8)$objective
      }
      if (val < o$value - 1e-13 * (1 + abs(o$value))) return(FALSE)
    }
    TRUE
  }
  same <- function(o, q) {                     # same point, or same end of a flat ridge
    max(abs(o$par - q$par)) < 1e-2 ||
      (any(edge_of(o) != 0L & edge_of(o) == edge_of(q)) && abs(o$par[1] - q$par[1]) < 1e-2 &&
         abs(o$value - q$value) < 1e-6 * (1 + abs(o$value)))
  }
  keep <- Filter(function(o) identical(o, best) || is_local_max(o),
                 fits[order(vapply(fits, function(o) o$value, numeric(1)))])
  pars <- list()
  for (o in keep) {
    if (!any(vapply(pars, function(q) same(o, q), logical(1)))) pars[[length(pars) + 1L]] <- o
  }
  lm_theta <- exp(vapply(pars, function(o) o$par[1], numeric(1)))
  lm_I <- exp(vapply(pars, function(o) o$par[2], numeric(1)))
  local_maxima <- data.frame(
    theta = lm_theta, m = lm_I / (J - 1 + lm_I), I = lm_I,
    logl = -vapply(pars, function(o) o$value, numeric(1)),
    boundary = at_bound(lm_theta, lower[1], upper[1], log_scale = TRUE) |
      at_bound(lm_I / (J - 1 + lm_I), lower[2], upper[2])
  )

  out <- list(theta = theta, m = m, I = I, logl = -best$value, converged = converged,
              se_theta = se_theta, se_m = se_m, se_I = se_I, local_maxima = local_maxima)

  if (ci) {
    ## profile likelihoods: for a fixed value of one parameter, maximise over the other
    ## (coarse grid, then refinement around the best grid point)
    profile <- function(fixed, which) {
      other <- if (which == 1L) c(lo[2], hi[2]) else c(lo[1], hi[1])
      at <- function(w) if (which == 1L) f(c(fixed, w)) else f(c(w, fixed))
      grid <- seq(other[1], other[2], length.out = 41L)
      vals <- vapply(grid, at, numeric(1))
      i <- which.min(vals)
      stats::optimize(at, c(grid[max(i - 1L, 1L)], grid[min(i + 1L, 41L)]), tol = 1e-6)$objective
    }
    target <- best$value + stats::qchisq(level, df = 1) / 2
    ci_lt <- lr_interval(function(z) profile(z, 1L), best$par[1], lo[1], hi[1], target)
    ci_lI <- lr_interval(function(z) profile(z, 2L), best$par[2], lo[2], hi[2], target)
    out$ci_theta <- exp(ci_lt)
    out$ci_I <- exp(ci_lI)
    out$ci_m <- out$ci_I / (J - 1 + out$ci_I)
  }
  out
}

# Internal helpers ----------------------------------------------------------------

# log(exp(x) + exp(y)), elementwise, without overflow; exact when x or y is -Inf.
log_add <- function(x, y) {
  mx <- pmax(x, y)
  out <- mx + log1p(exp(-abs(x - y)))
  out[mx == -Inf] <- -Inf
  out
}

# log(sum(exp(x))) without overflow.
log_sum_exp <- function(x) {
  mx <- max(x)
  mx + log(sum(exp(x - mx)))
}

# log K(D,A), A = k..J, for a validated vector `n` of positive whole numbers.
# The log T_n coefficients are kept only for the abundances present, so memory is
# proportional to the sum of the distinct abundances, not to max(n)^2. In the
# polynomial products, the R-level loop runs over the shorter of the two factors.
logkda_core <- function(n) {
  needed <- sort(unique(n))
  logT <- vector("list", length(needed))
  names(logT) <- needed
  cur <- 0                                     # log T_1
  if (needed[1] == 1) logT[["1"]] <- cur
  if (max(n) > 1) {
    for (s in 2:max(n)) {
      a <- seq_len(s)
      cur <- log_add(c(cur, -Inf), c(-Inf, cur) + log(a - 1) - log(s - 1))
      if (s %in% needed) logT[[as.character(s)]] <- cur
    }
  }
  L <- 0
  for (s in n) {
    t <- logT[[as.character(s)]]
    if (length(t) > length(L)) {
      tmp <- L
      L <- t
      t <- tmp
    }
    out <- rep(-Inf, length(L) + length(t) - 1L)
    idx <- seq_along(L)
    for (j in seq_along(t)) out[idx + j - 1L] <- log_add(out[idx + j - 1L], L + t[j])
    L <- out
  }
  L
}

# Validates a user-supplied `log_kda` (or computes it when NULL).
check_log_kda <- function(log_kda, n, fn) {
  if (is.null(log_kda)) return(logkda_core(n))
  if (!is.numeric(log_kda) || length(log_kda) != sum(n) - length(n) + 1L || anyNA(log_kda))
    stop(fn, "(): `log_kda` must be the output of logkda() for the same abundances.",
         call. = FALSE)
  log_kda
}

# Log-likelihood of Etienne's sampling formula (log-probability of the partition,
# on the scale of optim.ewens()), for theta > 0 and rescaled immigration rate I > 0:
#   k log(theta) + lgamma(theta) + sum(lgamma(n))
#     + log sum_A K(D,A) I^A Gamma(I) / (Gamma(I + J) Gamma(theta + A))
# (the lgamma(theta + J) terms of Etienne's equation 6 cancel). For large I (m close
# to 1), lgamma(I) - lgamma(I + J) + A log(I) is a difference of nearly equal large
# numbers; it is evaluated exactly as (A - J) log(I) - sum_{i<J} log1p(i / I), which
# has no cancellation. I = Inf is the Ewens limit m = 1.
etienne_logl <- function(theta, I, n, log_kda) {
  J <- sum(n)
  k <- length(n)
  ewens <- k * log(theta) + lgamma(theta) - lgamma(theta + J) + sum(lgamma(n))
  if (is.infinite(I)) return(ewens)
  A <- k:J
  t <- log_kda - lgamma(theta + A) + (A - J) * log(I)
  k * log(theta) + lgamma(theta) + sum(lgamma(n)) - sum(log1p(seq_len(J - 1) / I)) +
    log_sum_exp(t)
}

# Gradient of etienne_logl() with respect to (theta, I).
etienne_grad <- function(theta, I, n, log_kda) {
  J <- sum(n)
  k <- length(n)
  A <- k:J
  t <- log_kda - lgamma(theta + A) + (A - J) * log(I)
  w <- exp(t - log_sum_exp(t))                 # weights of the terms of the sum
  i <- seq_len(J - 1)
  c(k / theta + digamma(theta) - sum(w * digamma(theta + A)),
    (sum(w * A) - J) / I + sum(i / (I * (I + i))))
}

# Validates the `ci` and `level` arguments of the estimators.
check_ci_args <- function(ci, level, fn) {
  if (!is.logical(ci) || length(ci) != 1L || is.na(ci))
    stop(fn, "(): `ci` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop(fn, "(): `level` must be a number in (0, 1).", call. = FALSE)
  invisible(TRUE)
}

# TRUE where `value` lies on the edge of [lower, upper], with the tolerance of
# warn_boundary() (vectorised, silent).
at_bound <- function(value, lower, upper, log_scale = FALSE, tol = 1e-3) {
  v <- if (log_scale) log(value) else value
  lo <- if (log_scale) log(lower) else lower
  hi <- if (log_scale) log(upper) else upper
  is.finite(v) & ((v - lo < tol * (hi - lo)) | (hi - v < tol * (hi - lo)))
}

# Likelihood-ratio interval: the range around `est` where the profile `prof` (a
# negative log-likelihood) stays below `target`. A limit equal to `lo` or `hi` means
# that the interval extends at least to that bound of the search domain.
lr_interval <- function(prof, est, lo, hi, target) {
  side <- function(b) {
    if (abs(b - est) <= 1e-12 * (1 + abs(est)) || prof(b) <= target) return(b)
    stats::uniroot(function(z) prof(z) - target, sort(c(est, b)), tol = 1e-9)$root
  }
  c(lower = side(lo), upper = side(hi))
}

# Standard errors of the parameters x of a negative log-likelihood whose gradient is
# `g`, from the inverse of the Hessian (central differences of the exact gradient).
# NA when `on_bound` is TRUE or the Hessian is not positive definite.
se_from_gradient <- function(g, x, on_bound, h = 1e-5) {
  if (on_bound) return(rep(NA_real_, length(x)))
  H <- vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (g(x + e) - g(x - e)) / (2 * h)
  }, numeric(length(x)))
  H <- (H + t(H)) / 2
  ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA)
  if (anyNA(ev) || any(ev <= 0)) return(rep(NA_real_, length(x)))
  sqrt(diag(solve(H)))
}
