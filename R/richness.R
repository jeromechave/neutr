# Expected number of species ------------------------------------------------------

#' Expected number of classes given parameter \eqn{\theta} and sample size J for the Hoppe urn model
#'
#' @description The expected number of species is computed using an exact expression
#' based on the digamma function, \eqn{\theta[\psi(\theta+J)-\psi(\theta)]}.
#' The function is vectorised over `theta` and `J`.
#'
#' @param theta parameter \eqn{\theta}, a positive real value
#' @param J number of individuals
#'
#' @return the expected number of classes (species), a real number
#' @export
#'
#' @examples
#' kest(11.3, 234373)
kest <- function(theta, J) {
  theta * digamma(theta + J) - theta * digamma(theta)
}

#' Expected number of classes with more than N individuals given parameters \eqn{\theta} and sample size J using the GEM representation
#'
#' @description The expected number of species with more than `N` individuals is
#' estimated by averaging `nrep` independent draws of the neutral model with the
#' Griffiths-Engen-McCloskey representation ([generate.hoppe.urn()]).
#'
#' @param theta parameter \eqn{\theta}, a positive real value
#' @param J number of individuals, a whole number
#' @param N minimal abundance (a species is counted if its abundance is `> N`)
#' @param nrep number of simulated communities, at least 2 (default 100)
#'
#' @return A list with
#' \describe{
#'   \item{k}{the expected number of classes with more than `N` individuals, a real number}
#'   \item{sigmak}{the standard deviation of that number across the `nrep` draws,
#'     a real number (a measure of the spread of single communities, not the
#'     standard error of `k`)}
#' }
#' @seealso [kest.gt()] for a deterministic approximation
#' @export
#'
#' @examples
#' set.seed(1)
#' kest.gt0(11.3, 234373, 50)
kest.gt0 <- function(theta, J, N, nrep = 100) {
  fn <- "kest.gt0"
  check_scalar(theta, "theta", fn, lower = 0, lower_open = TRUE)
  check_size(J, fn)
  check_scalar(N, "N", fn, lower = 0)
  check_scalar(nrep, "nrep", fn, lower = 2)
  nbsp <- vapply(seq_len(round(nrep)), function(i) {
    sum(generate.hoppe.urn(theta, J)$abundance > N)
  }, numeric(1))
  list(k = mean(nbsp), sigmak = stats::sd(nbsp))
}

#' Expected number of classes with more than N individuals given parameters \eqn{\theta} and sample size J using the species frequency spectrum
#'
#' @description The expected number of species \eqn{>N} is computed using the approximation
#'  based on the species frequency spectrum \eqn{\psi(x)=\theta x^{-1}(1-x)^{\theta-1}}
#'  such that the number of species \eqn{>N} is the integral from N/J to 1 of \eqn{(1-(1-x)^J)\psi(x)}.
#'  This is an approximation of the exact formula, but it turns out to be accurate.
#'
#' @param theta parameter \eqn{\theta}, a positive real value
#' @param J number of individuals, a whole number
#' @param N minimal abundance, with `0 <= N < J`
#'
#' @return An object of class `"integrate"` as returned by [stats::integrate()]:
#'   the expected number of classes \eqn{>N} is its `$value` component, and
#'   `$abs.error` is the estimated absolute error of the numerical integration.
#' @export
#'
#' @examples
#' kest.gt(11.3, 234373, 50)$value
kest.gt <- function(theta, J, N) {
  fn <- "kest.gt"
  check_scalar(theta, "theta", fn, lower = 0, lower_open = TRUE)
  check_size(J, fn)
  check_scalar(N, "N", fn, lower = 0, upper = J, upper_open = TRUE)
  f <- function(x) theta * (1 - (1 - x)^J) / x * (1 - x)^(theta - 1)
  stats::integrate(f, N / J, 1.0)
}
