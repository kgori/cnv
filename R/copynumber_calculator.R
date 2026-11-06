#' Copy Number Calculator
#'
#' @description
#' The CopyNumberCalculator class calculates tumour copy number from the read depth ratio
#' (R) and the tumour purity and ploidy. It uses the following formula:
#' \eqn{N^T_i = α_i \psi^T_i + β_i (1 - \rho) / \rho}, where \cr
#' \eqn{α_i = R_i N^P_i / \psi^P_i}, \cr
#' \eqn{β_i = α_i \psi^H_i - N^H_i}, \cr
#' and \eqn{N^T_i} is the tumour copy number,
#' \eqn{\psi^T_i} is the tumour ploidy, \eqn{\rho} is the tumour purity, \eqn{N^P_i} is
#' the panel copy number, \eqn{\psi^P_i} is the panel ploidy, \eqn{\psi^H_i} is the host ploidy,
#' and \eqn{N^H_i} is the host copy number.
#'
#' Usage:
#' cnc <- CopyNumberCalculator$new(R, host_copynumber, host_ploidy, panel_copynumber, panel_ploidy) \cr
#' n <- cnc$copynumber(purity, ploidy)
#'
#' @importFrom "R6" R6Class
#' @export
CopyNumberCalculator <- R6Class(
    "CopyNumberCalculator",
    public = list(
        #' @field alpha The alpha parameter for the copy number calculation.
        #' Computed once and reused for any purity and ploidy.
        alpha = NULL,

        #' @field beta The beta parameter for the copy number calculation
        #' Computed once and reused for any purity and ploidy.
        beta = NULL,

        #' @description
        #' Constructor. Builds new CopyNumberCalculator using read depth ratio and host/panel parameters.
        #' @param R Read depth ratio of tumour mixture relative to panel of normals.
        #' @param host_copynumber The expected integer copy number of each locus for the mixed host
        #' @param host_ploidy The expected average autosomal copy number of the mixed host, e.g. 2 for diploid
        #' @param panel_copynumber The expected integer copy number of each locus for the panel of normals
        #' @param panel_ploidy The expected average copy number of the panel of normals
        initialize = function(
                R,
                host_copynumber = 2,
                host_ploidy = 2,
                panel_copynumber = 2,
                panel_ploidy = 2
        ) {
            self$alpha = R * panel_copynumber / panel_ploidy
            self$beta = self$alpha * host_ploidy - host_copynumber
        },

        #' @description
        #' Calculate tumour copy number for the given tumour purity, and tumour ploidy values.
        #' @param purity Tumour purity, between 0 and 1.
        #' @param ploidy Tumour ploidy, greater than 0.
        copynumber = function(purity, ploidy) {
            stopifnot(purity >= 0 & purity <= 1 & ploidy > 0)
            self$alpha * ploidy + self$beta * ((1 - purity) / purity)
        }
    ))

