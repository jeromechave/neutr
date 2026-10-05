#' neutr: Inference and Fast Generation of Neutral Species Assemblages
#'
#' Parameter inference and fast generation for three classes of neutral species
#' assemblage models: the Ewens model, the multideme model and the Pitman model.
#'
#' @section Estimation:
#' Parameters are estimated by maximising the likelihood: [optim.ewens()] for the
#' Ewens model, [optim.pitman()] for the Pitman model and [optim.multideme()] for
#' the multideme model.
#'
#' @section Generation:
#' [generate.hoppe.urn0()] simulates Hoppe's urn exactly, individual by individual;
#' [generate.hoppe.urn()] and [generate.pitman.urn()] use the
#' Griffiths-Engen-McCloskey representation, whose cost does not depend on the
#' sample size, and scale to 10^12 individuals and beyond.
#'
#' @section Expected richness:
#' [kest()], [kest.gt()] and [kest.gt0()] give the expected number of species
#' (possibly above an abundance threshold).
#'
#' @author Jerome Chave
#' @keywords internal
"_PACKAGE"
