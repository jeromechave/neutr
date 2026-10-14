# Deprecated function names ----------------------------------------------------
#
# In neutr 0.2.1 the exported functions were renamed from dotted names
# (e.g. `optim_ewens()`) to snake_case (e.g. `optim_ewens()`), following the
# tidyverse style guide: dots in function names are reserved for S3 methods,
# and names such as `kest_gt()` read as methods of a `kest()` generic.
# The old names remain available but are soft-deprecated.

#' Deprecated functions in neutr
#'
#' @description
#' These functions were renamed in neutr 0.2.1 to follow snake_case naming.
#' The old names still work but are deprecated and will be removed in a
#' future release. Please use the new names instead:
#'
#' | Deprecated                | Replacement               |
#' |---------------------------|---------------------------|
#' | `generate.hoppe.urn()`    | [generate_hoppe_urn()]    |
#' | `generate.hoppe.urn0()`   | [generate_hoppe_urn0()]   |
#' | `generate.pitman.urn()`   | [generate_pitman_urn()]   |
#' | `kest.gt()`               | [kest_gt()]               |
#' | `kest.gt0()`              | [kest_gt0()]              |
#' | `logl.etienne()`          | [logl_etienne()]          |
#' | `optim.etienne()`         | [optim_etienne()]         |
#' | `optim.ewens()`           | [optim_ewens()]           |
#' | `optim.multideme()`       | [optim_multideme()]       |
#' | `optim.pitman()`          | [optim_pitman()]          |
#'
#' The deprecated functions take the same arguments and return the same
#' values as their replacements. They signal a soft deprecation warning
#' with [lifecycle::deprecate_soft()]: it is shown when a deprecated function
#' is called from the global environment or from tests, but not when it is
#' called from another package.
#'
#' @param theta,sigma,J,seed,N,nrep,input_abundances,input_abundance_matrix
#'   See the documentation of the replacement function.
#' @param m,log_kda,full,init_vals,lower,upper,ci,level,verbose,I_max
#'   See the documentation of the replacement function.
#' @return See the documentation of the replacement function.
#' @name neutr-deprecated
#' @keywords internal
NULL

#' @rdname neutr-deprecated
#' @export
generate.hoppe.urn <- function(theta, J) {
  lifecycle::deprecate_soft(
    "0.2.1", "generate.hoppe.urn()", "generate_hoppe_urn()"
  )
  generate_hoppe_urn(theta = theta, J = J)
}

#' @rdname neutr-deprecated
#' @export
generate.hoppe.urn0 <- function(theta, J, seed = NULL) {
  lifecycle::deprecate_soft(
    "0.2.1", "generate.hoppe.urn0()", "generate_hoppe_urn0()"
  )
  generate_hoppe_urn0(theta = theta, J = J, seed = seed)
}

#' @rdname neutr-deprecated
#' @export
generate.pitman.urn <- function(theta, sigma, J) {
  lifecycle::deprecate_soft(
    "0.2.1", "generate.pitman.urn()", "generate_pitman_urn()"
  )
  generate_pitman_urn(theta = theta, sigma = sigma, J = J)
}

#' @rdname neutr-deprecated
#' @export
kest.gt <- function(theta, J, N) {
  lifecycle::deprecate_soft("0.2.1", "kest.gt()", "kest_gt()")
  kest_gt(theta = theta, J = J, N = N)
}

#' @rdname neutr-deprecated
#' @export
kest.gt0 <- function(theta, J, N, nrep = 100) {
  lifecycle::deprecate_soft("0.2.1", "kest.gt0()", "kest_gt0()")
  kest_gt0(theta = theta, J = J, N = N, nrep = nrep)
}

#' @rdname neutr-deprecated
#' @export
logl.etienne <- function(input_abundances, theta, m, log_kda = NULL,
                         full = FALSE) {
  lifecycle::deprecate_soft("0.2.1", "logl.etienne()", "logl_etienne()")
  logl_etienne(
    input_abundances = input_abundances, theta = theta, m = m,
    log_kda = log_kda, full = full
  )
}

#' @rdname neutr-deprecated
#' @export
optim.etienne <- function(input_abundances, init_vals = c(10, 0.5),
                          lower = c(1e-6, 1e-8), upper = NULL, ci = FALSE,
                          level = 0.95) {
  lifecycle::deprecate_soft("0.2.1", "optim.etienne()", "optim_etienne()")
  optim_etienne(
    input_abundances = input_abundances, init_vals = init_vals,
    lower = lower, upper = upper, ci = ci, level = level
  )
}

#' @rdname neutr-deprecated
#' @export
optim.ewens <- function(input_abundances) {
  lifecycle::deprecate_soft("0.2.1", "optim.ewens()", "optim_ewens()")
  optim_ewens(input_abundances = input_abundances)
}

#' @rdname neutr-deprecated
#' @export
optim.multideme <- function(input_abundance_matrix, verbose = FALSE,
                            I_max = NULL, ci = FALSE, level = 0.95) {
  lifecycle::deprecate_soft(
    "0.2.1", "optim.multideme()", "optim_multideme()"
  )
  optim_multideme(
    input_abundance_matrix = input_abundance_matrix, verbose = verbose,
    I_max = I_max, ci = ci, level = level
  )
}

#' @rdname neutr-deprecated
#' @export
optim.pitman <- function(input_abundances, init_vals = c(10, 0.1),
                         lower = c(1e-6, 1e-8), upper = NULL) {
  lifecycle::deprecate_soft("0.2.1", "optim.pitman()", "optim_pitman()")
  optim_pitman(
    input_abundances = input_abundances, init_vals = init_vals,
    lower = lower, upper = upper
  )
}
