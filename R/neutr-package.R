#' neutr: Inference and Fast Generation of Neutral Species Assemblages
#'
#' Parameter inference and fast generation for neutral species assemblage
#' models: the Ewens model, the Etienne model (a local community with dispersal
#' limitation), the multideme model and the Pitman model.
#'
#' @section Estimation:
#' Parameters are estimated by maximising the likelihood: [optim_ewens()] for
#' the Ewens model, [optim_etienne()] for the Etienne model, [optim_pitman()]
#' for the Pitman model and [optim_multideme()] for the multideme model. The
#' Ewens model is nested in both the Etienne model (\eqn{m = 1}) and the Pitman
#' model (\eqn{\sigma = 0}), and the three functions return their log-likelihood
#' on the same scale.
#'
#' @section Etienne sampling formula:
#' [logl_etienne()] evaluates the likelihood of Etienne's sampling formula and
#' [logkda()] its coefficients \eqn{K(D,A)}, computed on the log scale so that
#' large samples cannot overflow.
#'
#' @section Generation:
#' [generate_hoppe_urn0()] simulates Hoppe's urn exactly, individual by
#' individual; [generate_hoppe_urn()] and [generate_pitman_urn()] use the
#' Griffiths-Engen-McCloskey representation, whose cost does not depend on the
#' sample size, and scale to 10^12 individuals and beyond.
#'
#' @section Expected richness:
#' [kest()], [kest_gt()] and [kest_gt0()] give the expected number of species
#' (possibly above an abundance threshold).
#'
#' @author Jerome Chave
#' @keywords internal
"_PACKAGE"
