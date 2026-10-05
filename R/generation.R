# Generation of neutral species abundance distributions ---------------------------

#' Generation of a neutral partition given parameter \eqn{\theta} using Hoppe's urn
#'
#' @description `generate.hoppe.urn0()` creates a rank-abundance distribution based on
#' Hoppe's urn scheme for a given parameter \eqn{\theta} and sample size `J`.
#' The process starts at zero species, corresponding to the single black ball with
#' weight \eqn{\theta}. If the black ball is picked, a new class is created, else the
#' picked colored ball is duplicated (its abundance increases by one unit).
#'
#' This is an *exact* simulation of the urn, but it does not loop over the
#' individuals: the genealogy of the `J` individuals is resolved in a handful of
#' vectorised passes, and the random numbers come from the fast `dqrng` generator.
#' Memory use is proportional to `J`, so it is meant for `J` up to a few hundred
#' million; for larger samples use [generate.hoppe.urn()].
#'
#' @section Reproducibility:
#' `generate.hoppe.urn0()` draws its random numbers from `dqrng`, which has its own
#' random number stream: `set.seed()` has **no effect** on it. Use the `seed`
#' argument instead. [generate.hoppe.urn()] and [generate.pitman.urn()] use R's
#' standard generator and follow `set.seed()`.
#'
#' @param theta parameter \eqn{\theta}, a positive real value
#' @param J number of individuals, a whole number
#' @param seed optional integer; if given, the `dqrng` generator is seeded with it
#'   (this changes the state of the `dqrng` generator for the whole session)
#'
#' @return `generate.hoppe.urn0()` returns the species abundances, an integer vector
#'   sorted in decreasing order and summing to `J`.
#' @seealso [generate.hoppe.urn()] for a faster approach on very large samples
#' @export
#'
#' @examples
#' generate.hoppe.urn0(11.3, 234373, seed = 1)
generate.hoppe.urn0 <- function(theta, J, seed = NULL) {
  fn <- "generate.hoppe.urn0"
  check_scalar(theta, "theta", fn, lower = 0, lower_open = TRUE)
  check_size(J, fn)
  if (J > .Machine$integer.max)
    stop(fn, "(): J > ", .Machine$integer.max, " is too large for an individual-based ",
         "simulation; use generate.hoppe.urn() instead.", call. = FALSE)
  if (!is.null(seed)) dqrng::dqset.seed(seed)

  J <- as.integer(J)
  j <- seq_len(J)

  ## Individual j founds a new species with probability theta / (theta + j - 1)
  ## (always the case for j = 1); otherwise it copies a uniformly chosen earlier one.
  is_new <- dqrng::dqrunif(J) <= theta / (theta + j - 1)
  parent <- as.integer(floor(dqrng::dqrunif(J) * (j - 1)) + 1)
  parent[is_new] <- j[is_new]

  ## Each individual belongs to the species of its oldest ancestor. Follow the parent
  ## pointers by repeated squaring: the depth of the genealogy is O(log J), so this
  ## converges in a few passes.
  root <- parent
  repeat {
    nxt <- root[root]
    if (identical(nxt, root)) break
    root <- nxt
  }

  x <- tabulate(root, nbins = J)
  sort(x[x > 0L], decreasing = TRUE)
}

#' Generation of a neutral partition given parameter \eqn{\theta} using the Griffiths-Engen-McCloskey representation
#'
#' @description `generate.hoppe.urn()` creates a rank-abundance distribution for a given
#' \eqn{\theta} and `J`. Beta-distributed random variables are drawn from
#' \eqn{Beta(1,\theta)} to build the stick-breaking weights of the Griffiths-Engen-McCloskey
#' (GEM) distribution, and the partition is a multinomial draw of `J` individuals
#' according to these weights. Its cost does not depend on `J` (up to the
#' multinomial draw), which makes it suitable for very large samples
#' (10^12 individuals and beyond).
#'
#' The draw is exact: sticks are generated in blocks until every individual has been
#' assigned. The total number of individuals is always exactly `J`.
#'
#' @param theta parameter \eqn{\theta}, a positive real value
#' @param J number of individuals, a whole number (up to 2^53)
#'
#' @return A list with
#' \describe{
#'   \item{abundance}{the species abundances, sorted in decreasing order and summing to `J`}
#'   \item{k}{the *expected* number of species for these parameters (see [kest()]);
#'     the number of species in this particular draw is `length(abundance)`}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' out <- generate.hoppe.urn(11.3, 234373)
#' length(out$abundance)  # close to out$k
generate.hoppe.urn <- function(theta, J) {
  fn <- "generate.hoppe.urn"
  check_scalar(theta, "theta", fn, lower = 0, lower_open = TRUE)
  check_size(J, fn)

  ## expected number of species
  kest <- theta * (digamma(theta + J) - digamma(theta))

  x <- gem_abundances(theta, 0, J)
  list(abundance = sort(x, decreasing = TRUE), k = kest)
}

#' Generation of a Pitman partition given parameters \eqn{\theta} and \eqn{\sigma} from the Griffiths-Engen-McCloskey representation
#'
#' @description The Pitman urn scheme is as follows: starting at zero species,
#' corresponding to the single black ball with weight \eqn{\theta}, a species i, of k
#' species total, is selected with probability \eqn{(n_i-\sigma)/(n+\theta)}, where
#' \eqn{\sigma} is strictly between 0.0 and 1.0 and n is the current number of balls;
#' the black ball is picked with probability \eqn{(k\sigma+\theta)/(n+\theta)}.
#' `generate.pitman.urn()` creates a rank-abundance distribution for given
#' \eqn{\theta}, \eqn{\sigma} and `J`. Beta-distributed random variables are drawn from
#' \eqn{Beta(1-\sigma,\theta+i\sigma)} (for the i-th stick) and the partition is the
#' multinomial draw of `J` individuals according to the resulting GEM weights.
#'
#' The draw is exact: sticks are generated in blocks until every individual has been
#' assigned, and the few individuals left in the far tail are drawn with the exact
#' sequential Pitman urn. The total number of individuals is always exactly `J`.
#'
#' @param theta parameter \eqn{\theta}, a real value \eqn{\ge 0}
#' @param sigma parameter \eqn{\sigma}, a real value strictly between 0 and 1
#'   (use [generate.hoppe.urn()] for \eqn{\sigma = 0})
#' @param J number of individuals, a whole number (up to 2^53)
#'
#' @return A list with
#' \describe{
#'   \item{abundance}{the species abundances, sorted in decreasing order and summing to `J`}
#'   \item{k}{the *expected* number of species for these parameters;
#'     the number of species in this particular draw is `length(abundance)`}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' generate.pitman.urn(11.3, 0.1, 234373)$k
generate.pitman.urn <- function(theta, sigma, J) {
  fn <- "generate.pitman.urn"
  check_scalar(theta, "theta", fn, lower = 0)
  check_scalar(sigma, "sigma", fn, lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE)
  check_size(J, fn)

  ## expected number of species under the Pitman model
  kest <- kest_pitman(theta, sigma, J)

  x <- gem_abundances(theta, sigma, J)
  list(abundance = sort(x, decreasing = TRUE), k = kest)
}
