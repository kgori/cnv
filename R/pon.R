#' Clip numeric vector so that its lowest element is >= low,
#' and its highest is <= high
clip <- function(v, low, high) {
    pmax(low, pmin(v, high))
}

#' It's log2, but has a minimum or maximum of ±31
safelog2 <- function(x) {
    clip(log2(x), -31, 31)
}

#' Panel Of Normals
#'
#' @description Build this from a matrix of normal read counts (loci in rows, samples in columns). It
#' is used to compute normalised read count ratios from external samples, usually tumours.
#' It represents the idealised copy number of the organism, with no aberration. So for a
#' diploid, the ploidy is 2, autosomal loci have a copy number of 2 and sex chromosomes
#' have a copy number of 1, 2 or 0, depending on the sex. for this reason the panel should
#' be sex-specific.
#'
#' Denoising:
#' Implements a version of GATK's logR denoising using a panel of normals.
#' This is simplified compared to the GATK procedure, as it doesn't apply any filters to the
#' data, except for an optional winsorization. As such, it is advisable to filter out rows of
#' the host matrix that have very low or high coverage, or which have row-median of zero.
#'
#' @importFrom "R6" R6Class
#' @importFrom "matrixStats" rowMedians colMedians
#' @importFrom "irlba" irlba
#' @export
PanelOfNormals = R6Class(
    "PanelOfNormals",

    public=list(

        #' @field panel_reference Row medians of the input host read count matrix, or
        #' N^P / ψ^P
        panel_reference = NULL,

        #' @field panel_copynumber The expected integer copy number of each locus for
        #' the organism represented by the panel.
        panel_copynumber = NULL,

        #' @field panel_ploidy The expected average copy number of the organism, e.g. 2
        #' for diploid
        panel_ploidy = NULL,

        #' @field U Result of Singular Value Decomposition of processed host read counts,
        #' used for denoising tumour logR.
        U = NULL,

        #' @field D Singular Values of the decomposition
        D = NULL,

        #' @description
        #' Constructor. Builds new PanelOfNormals using host read count matrix as input
        #' @param host_matrix Matrix of read counts, with samples as columns.
        #' @param winsorization_quantile Values outside the quantile range set here
        #' (q, 1-q) are clamped to the values at those quantiles.
        #' Can be set to 0 to disable.
        #' @param max_svs Maximum number of singular values to retain. The singular value
        #' decomposition is done using irlba
        initialize = function(
                host_matrix,
                panel_copynumber = NULL,
                panel_ploidy = 2,
                winsorization_quantile=0.001,
                max_svs = 10) {
            stopifnot(winsorization_quantile >= 0 & winsorization_quantile <= 1)

            # Normalised panel
            host_matrix <-  sweep(host_matrix, 2, colMedians(host_matrix), "/")

            ref <- rowMedians(host_matrix)
            ref[ref == 0] <- 1
            self$panel_reference <- ref
            self$panel_ploidy <- panel_ploidy

            if (is.null(panel_copynumber)) {
                panel_copynumber <- rep(panel_ploidy)
            }

            self$panel_copynumber <- panel_copynumber

            # Learn the denoising basis
            # row normalise
            host_matrix <- sweep(host_matrix, 1, ref, "/")
            # winsorise
            host_matrix <- apply(host_matrix, 2, function(col) clip(col,
                                                                    quantile(col, winsorization_quantile),
                                                                    quantile(col, 1 - winsorization_quantile)))
            # impute missing values
            ix=which(host_matrix==0, arr.ind=TRUE)
            host_matrix[ix] <- rowMedians(host_matrix)[ix[, "row"]]

            # renormalise
            host_matrix <- sweep(host_matrix, 2, colMedians(host_matrix), "/")

            # move to log space
            host_matrix <- apply(host_matrix, 2, safelog2)

            # decompose
            decomp <- irlba(host_matrix, max_svs)
            self$U = decomp$u
            self$D = decomp$d
        },

        #' @description
        #' Calculate standardized (not denoised) read count ratio (R) from tumour read counts
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        ratio = function(tumour_matrix) {
            stopifnot(nrow(tumour_matrix) == length(self$panel_reference))
            tumour_matrix <- sweep(tumour_matrix, 2, colMedians(tumour_matrix), "/")
            sweep(tumour_matrix, 1, self$panel_reference, "/")
        },

        #' @description
        #' Calculate standardized (not denoised) logR from tumour read counts
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        logR = function(tumour_matrix) {
            apply(self$ratio(tumour_matrix), 2, safelog2)
        },

        #' @description
        #' Calculate denoised logR from tumour read counts.
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        #' @param num_svs Number of singular values to use for denoising.
        #' @importFrom "logging" logwarn
        denoised_logR = function(tumour_matrix, num_svs) {
            stopifnot(nrow(tumour_matrix) == nrow(self$U))
            if (num_svs > ncol(self$U)) {
                logwarn("num_svs is larger than ncol(U) [%d]", ncol(self$U))
                logwarn("num_svs reduced to %d", ncol(self$U))
                num_svs <- ncol(self$U)
            }
            if (num_svs < 1) {
                logwarn("num_svs set to minimum value: 1")
                num_svs <- 1
            }

            x <- self$logR(tumour_matrix)
            u <- self$U[, seq_len(num_svs), drop = FALSE]

            x - u %*% (t(u) %*% x)
        }
    )
)
