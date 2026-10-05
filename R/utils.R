# General use functions ------------------------------------------------------------

#' General use function
#'
#' @description `zeroremove()` is used to remove zeros and NAs in a vector of reals
#' or of integers. This is a simple utility function.
#'
#' @param x vector of real
#'
#' @return the vector `x` without its zero and missing (`NA`/`NaN`) elements
#' @export
#'
#' @examples
#' zeroremove(c(2.1, 0.0, 3.0, NA))
zeroremove <- function(x) x[!is.na(x) & x != 0]

# Internal helpers --------------------------------------------------------------

# Checks that `x` is a single finite number, optionally within bounds.
check_scalar <- function(x, name, fn, lower = -Inf, upper = Inf,
                         lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    (if (lower_open) x > lower else x >= lower) &&
    (if (upper_open) x < upper else x <= upper)
  if (!ok) {
    range_txt <- if (is.finite(lower) || is.finite(upper)) {
      sprintf(" in %s%s, %s%s", if (lower_open) "(" else "[", format(lower),
              format(upper), if (upper_open) ")" else "]")
    } else ""
    stop(sprintf("%s(): `%s` must be a single finite number%s.", fn, name, range_txt),
         call. = FALSE)
  }
  invisible(x)
}

# Checks that `J` is a single whole number >= 1 (stored as double: up to 2^53).
check_size <- function(J, fn) {
  check_scalar(J, "J", fn, lower = 1)
  if (J != round(J))
    stop(sprintf("%s(): `J` must be a whole number of individuals.", fn), call. = FALSE)
  if (J > 2^53)
    stop(sprintf("%s(): `J` must be at most 2^53.", fn), call. = FALSE)
  invisible(J)
}

# Validates a vector of species abundances and drops the empty species.
# Returns a numeric vector of strictly positive whole numbers.
check_abundances <- function(x, fn, min_species = 1L) {
  if (!is.numeric(x) || length(x) == 0L)
    stop(sprintf("%s(): `input_abundances` must be a non-empty numeric vector.", fn),
         call. = FALSE)
  if (anyNA(x) || any(!is.finite(x)))
    stop(sprintf("%s(): `input_abundances` must not contain NA, NaN or infinite values.", fn),
         call. = FALSE)
  if (any(x < 0))
    stop(sprintf("%s(): `input_abundances` must be non-negative.", fn), call. = FALSE)
  if (any(abs(x - round(x)) > sqrt(.Machine$double.eps)))
    stop(sprintf("%s(): `input_abundances` must be counts (whole numbers).", fn), call. = FALSE)
  x <- as.numeric(round(x[x > 0]))
  if (length(x) < min_species)
    stop(sprintf("%s(): at least %d species with positive abundance are required.",
                 fn, min_species), call. = FALSE)
  x
}

# Multinomial draw that accepts a total size above .Machine$integer.max
# (stats::rmultinom cannot). The total is split into near-equal chunks, each one
# drawn independently with the same probabilities; their sum is exactly a
# multinomial of the full size. Returns a numeric vector of length(prob).
rmultinom_big <- function(size, prob) {
  limit <- .Machine$integer.max
  if (size <= limit) return(as.numeric(stats::rmultinom(1L, size, prob)))
  n_chunks <- ceiling(size / limit)
  base <- floor(size / n_chunks)
  sizes <- rep(base, n_chunks)
  extra <- size - base * n_chunks              # 0 <= extra < n_chunks
  if (extra > 0) sizes[seq_len(extra)] <- base + 1
  rowSums(vapply(sizes, function(s) stats::rmultinom(1L, s, prob)[, 1L],
                 numeric(length(prob))))
}

# Exact abundance vector (unsorted, zeros removed) of J individuals drawn from the
# Pitman-Yor (sigma > 0) / Ewens (sigma = 0) model, using its GEM stick-breaking
# representation for the bulk of the community:
#   W_i ~ Beta(1 - sigma, theta + i * sigma),  P_i = W_i * prod_{l < i} (1 - W_l).
#
# The multinomial draw of the individuals over a block of B sticks has one extra
# "leftover" category (the mass of all later sticks). The leftover individuals are
# *exactly* a Pitman(sigma, theta + B * sigma) partition of the leftover count: the
# later sticks are independent, with Beta parameters shifted by the B sticks already
# used. The draw is therefore exact. What happens to the leftover depends on sigma:
#  * sigma = 0 (Ewens): the stick masses decay geometrically, the leftover is
#    essentially always empty; if not, the next block continues the same process.
#  * sigma > 0: the masses decay only polynomially, so the last few individuals sit
#    in sticks of astronomically large index and no block size reaches them. The
#    (small) leftover is instead drawn with the exact sequential Pitman urn.
gem_abundances <- function(theta, sigma, J, max_block = 2e6, urn_max = 2e6,
                           max_rounds = 500) {
  remaining <- J
  counts <- list()
  for (round in seq_len(max_rounds)) {
    if (round > 1 && sigma > 0 && remaining <= urn_max) {
      counts[[round]] <- pitman_urn(theta, sigma, remaining)
      remaining <- 0
      break
    }
    ## block size: a multiple of the expected number of species of the current
    ## (sub-)problem, which is a safe overestimate of the sticks that will be used
    kh <- if (sigma > 0) kest_pitman(theta, sigma, remaining) else kest_digamma(theta, remaining)
    if (!is.finite(kh) || kh < 1) kh <- 1
    block <- as.integer(min(max_block, max(64, ceiling(3 * min(kh, remaining)))))

    idx <- seq_len(block)
    W <- stats::rbeta(block, 1 - sigma, theta + idx * sigma)
    ###W <- stats::rbeta(block, 1 - sigma, theta + (idx - 1L) * sigma)
    one_minus <- 1 - W
    P <- W * c(1, cumprod(one_minus[-block]))
    draw <- rmultinom_big(remaining, c(P, prod(one_minus)))
    counts[[round]] <- draw[idx]
    remaining <- draw[block + 1L]
    if (remaining <= 0) break
    theta <- theta + block * sigma     # the tail is Pitman(sigma, theta + B * sigma)
  }
  if (remaining > 0)
    stop("the parameters give too many species to be generated exactly ",
         "(try a smaller sigma or J).", call. = FALSE)
  x <- unlist(counts, use.names = FALSE)
  x[x > 0]
}

# Exact sequential Pitman(sigma, theta) urn for L individuals; returns the species
# abundances. Individual n + 1 founds a new species with probability
# (theta + k * sigma) / (theta + n), otherwise joins species i with probability
# (n_i - sigma) / (theta + n). The latter is drawn by rejection: pick a uniformly
# chosen earlier individual (species i with probability n_i / n) and accept with
# probability (n_i - sigma) / n_i >= 1 - sigma.
pitman_urn <- function(theta, sigma, L) {
  L <- as.integer(L)
  lab <- integer(L)                 # species of each individual
  size <- integer(L)                # abundance of each species
  k <- 0L
  u_new <- stats::runif(L)
  u_pick <- stats::runif(L)
  u_acc <- stats::runif(L)
  for (n in seq_len(L) - 1L) {      # n individuals are already in the urn
    if (n == 0L || u_new[n + 1L] * (theta + n) < theta + k * sigma) {
      k <- k + 1L
      lab[n + 1L] <- k
      size[k] <- 1L
    } else {
      s <- lab[floor(u_pick[n + 1L] * n) + 1L]
      while (u_acc[n + 1L] * size[s] >= size[s] - sigma) {   # rejected: draw again
        v <- stats::runif(2L)
        s <- lab[floor(v[1L] * n) + 1L]
        u_acc[n + 1L] <- v[2L]
      }
      lab[n + 1L] <- s
      size[s] <- size[s] + 1L
    }
  }
  size[seq_len(k)]
}

# Expected number of species under the Ewens model.
kest_digamma <- function(theta, J) theta * (digamma(theta + J) - digamma(theta))

# Expected number of species under the Pitman-Yor model (sigma > 0):
#   E[K] = (theta / sigma) * ((theta + sigma)_J / (theta)_J - 1).
# Written with lbeta() and expm1() so that it stays accurate when theta is large
# compared with J (the naive form subtracts two nearly equal numbers).
kest_pitman <- function(theta, sigma, J) {
  if (theta == 0) return(exp(lgamma(J + sigma) - lgamma(sigma + 1) - lgamma(J)))
  theta / sigma * expm1(lbeta(theta, sigma) - lbeta(theta + J, sigma))
}

# Warn when an estimate sits on the edge of its search interval.
warn_boundary <- function(value, lower, upper, name, fn, log_scale = FALSE, tol = 1e-3) {
  if (!is.finite(value)) return(invisible(FALSE))
  v <- if (log_scale) log(value) else value
  lo <- if (log_scale) log(lower) else lower
  hi <- if (log_scale) log(upper) else upper
  width <- hi - lo
  side <- if (v - lo < tol * width) "lower" else if (hi - v < tol * width) "upper" else NA
  if (!is.na(side)) {
    warning(sprintf(paste0("%s(): the estimate of `%s` (%s) is at the %s bound of the search interval ",
                           "[%s, %s]; the likelihood may have no interior maximum for these data."),
                    fn, name, format(signif(value, 4)), side,
                    format(signif(lower, 4)), format(signif(upper, 4))), call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

# Log-probability of a partition under the Pitman sampling formula (EPPF), for
# species abundances summarised as distinct values `u` with multiplicities `cnt`
# (J = sum(u * cnt) individuals, k = sum(cnt) >= 2 species):
#   (k-1) log(sigma) + lgamma(theta/sigma + k) - lgamma(theta/sigma + 1)
#     - lgamma(theta + J) + lgamma(theta + 1) + sum_i [lgamma(n_i - sigma) - lgamma(1 - sigma)]
# lgamma(a + k) - lgamma(a + 1) is evaluated as lgamma(k - 1) - lbeta(a + 1, k - 1),
# which stays accurate when a = theta / sigma is very large (small sigma).
pitman_logl <- function(theta, sigma, u, cnt, J, k) {
  a <- theta / sigma
  (k - 1) * log(sigma) + (lgamma(k - 1) - lbeta(a + 1, k - 1)) -
    (lgamma(theta + J) - lgamma(theta + 1)) +
    sum(cnt * (lgamma(u - sigma) - lgamma(1 - sigma)))
}

# Gradient of pitman_logl() with respect to (theta, sigma).
pitman_grad <- function(theta, sigma, u, cnt, J, k) {
  a <- theta / sigma
  d <- digamma(a + k) - digamma(a + 1)
  c(d / sigma + digamma(theta + 1) - digamma(theta + J),
    (k - 1) / sigma - theta / sigma^2 * d -
      sum(cnt * (digamma(u - sigma) - digamma(1 - sigma))))
}
