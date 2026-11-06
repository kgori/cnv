#' Clip numeric vector so that its lowest element is >= low,
#' and its highest is <= high
clip <- function(v, low, high) {
    pmin(pmax(v, low), high)
}

#' It's log2, but has a minimum or maximum of ±31
safelog2 <- function(x) {
    clip(log2(x), -31, 31)
}

#' Validate a vector of sex labels ("F" / "M", case-insensitive)
check_sex <- function(sex, n) {
    if (length(sex) != n) {
        stop(sprintf("Input must have one entry per sample (expected %d, got %d)",
                     n, length(sex)))
    }
    sex <- toupper(as.character(sex))
    if (anyNA(sex) | !all(sex %in% c("F", "M"))) {
        stop("Input must contain only 'F' or 'M', with no NAs")
    }
    sex
}

#' Per-sample size factors: median count over autosomal / PAR loci only, so that
#' scaling does not depend on the sex of the sample.
autosomal_size_factors <- function(counts, is_auto) {
    colMedians(counts[is_auto, , drop = FALSE])
}

#' Expected copy number of every locus (rows) in each sample (columns), given
#' per-locus expected copy number for females and males and a sex label per sample.
sample_copynumber <- function(cn_female, cn_male, sex) {
    cn <- matrix(cn_female, nrow = length(cn_female), ncol = length(sex))
    cn[, sex == "M"] <- cn_male
    cn
}

#' Panel Of Normals
#'
#' @description Build this from a matrix of normal read counts (loci in rows,
#' samples in columns), possibly of mixed sex. It is used to compute normalised
#' read count ratios from external samples, usually tumours.
#'
#' The panel represents the copy number of an idealised, non-biological
#' reference genome: copy number 2 at autosomal and pseudo-autosomal loci, 2 at
#' X and 1 at Y. We use the mammalian 'X' and 'Y' sex determination system here
#' for convenience, but it is not a necessity, as all that is needed is two per-
#' locus vectors of expected copy number, for female and for male samples. The
#' panel copy number in the reference of each locus is the larger of the two.
#' Only samples expected to be at that copy number contribute to the locus
#' reference, so autosomes and PAR use all samples, non-PAR X uses females, and
#' non-PAR Y uses males. Again, these terms are from mammalian sex
#' determination, but the panel can be more general than this.
#'
#' Tumour ratios are computed against this fixed reference and do not depend on
#' the sex of the tumour or of its host. A male tumour on a female host (or any
#' other combination) is used directly: the expected ratio is 1 at autosomes,
#' 0.5 for a single X copy, 0 for no Y, and so on. Copy number is estimated as
#' ratio x panel_copynumber (up to the usual ploidy / purity scaling).
#'
#' Denoising:
#' Implements a version of GATK's logR denoising using a panel of normals.
#' This is simplified compared to the GATK procedure, as it doesn't apply any
#' filters to the data, except for an optional winsorization. As such, it is
#' advisable to filter out rows of the host matrix that have very low or high
#' coverage, or which have row-median of zero. To compute the basis, each panel
#' sample is first brought onto the reference scale by dividing by (its expected
#' copy number / panel copy number), e.g. male X is doubled. Entries where a
#' panel sample is expected to have copy number 0 (female Y) carry no
#' information and are set to logR 0 for the decomposition.
#'
#' @importFrom "R6" R6Class
#' @importFrom "matrixStats" rowMedians colMedians
#' @importFrom "irlba" irlba
#' @importFrom "logging" logwarn
#' @importFrom "stats" median quantile
#' @export
SexAgnosticPanelOfNormals = R6Class(
    "SexAgnosticPanelOfNormals",

    public=list(

        #' @field loci data.frame of locus annotation (CHROM, START, END), one row per
        #' row of the host matrix. Stored for reference; not used in calculations.
        loci = NULL,

        #' @field panel_sex Sex labels ("F"/"M") of the panel samples
        panel_sex = NULL,

        #' @field panel_reference Per-locus reference level: median of size-factor
        #' normalised counts over the reference samples for that locus. Expressed at
        #' the panel copy number.
        panel_reference = NULL,

        #' @field panel_copynumber The copy number of each locus in the idealised
        #' poly-sex panel genome (pmax of the female and male vectors): 2 for autosomes,
        #' PAR and X, 1 for Y.
        panel_copynumber = NULL,

        #' @field panel_ploidy Copy number of the autosomal / PAR loci, e.g. 2 for diploid
        panel_ploidy = NULL,

        #' @field is_auto Logical per locus: TRUE where expected copy number is the same
        #' and non-zero in both sexes (autosomes and PAR). Used for size factors.
        is_auto = NULL,

        #' @field U Result of Singular Value Decomposition of processed host read counts,
        #' used for denoising tumour logR.
        U = NULL,

        #' @field D Singular Values of the decomposition
        D = NULL,

        #' @description
        #' Constructor. Builds new PanelOfNormals using host read count matrix as input
        #' @param host_matrix Matrix of read counts, with samples as columns.
        #' @param loci data.frame or data.table with columns CHROM, START, END, one row per
        #' row of host_matrix.
        #' @param sex Vector of sex labels ("F" or "M"), one per column of host_matrix.
        #' @param copynumber_female Numeric vector, expected copy number of each locus in
        #' a female sample.
        #' @param copynumber_male Numeric vector, expected copy number of each locus in
        #' a male sample.
        #' @param winsorization_quantile Values outside the quantile range set here
        #' (q, 1-q) are clamped to the values at those quantiles.
        #' Can be set to 0 to disable.
        #' @param max_svs Maximum number of singular values to retain. The singular value
        #' decomposition is done using irlba
        #' @param min_ref_samples Warn if any locus has fewer than this many samples
        #' contributing to its reference.
        #' @param verbose Print progress updates during construction
        initialize = function(
        host_matrix,
        loci,
        sex,
        copynumber_female,
        copynumber_male,
        winsorization_quantile = 0.001,
        max_svs = 10,
        min_ref_samples = 3,
        verbose = FALSE) {
            stopifnot(winsorization_quantile >= 0 & winsorization_quantile <= 1)
            if (verbose) {
                cat(sprintf("Winsorization Quantile %.6f is valid\n", winsorization_quantile))
            }

            # validate inputs
            n_loci <- nrow(host_matrix)
            if (!is.data.frame(loci) || nrow(loci) != n_loci ||
                !all(c("CHROM", "START", "END") %in% names(loci))) {
                stop("loci must be a data.frame or data.table with columns CHROM, START, END and one row per row of host_matrix")
            }
            if (verbose) {
                cat("Loci and input matrix sizes are compatible\n")
            }

            sex <- check_sex(sex, ncol(host_matrix))
            if (verbose) {
                cat("sex labels and input matrix sizes are compatible\n")
            }

            for (cn in list(copynumber_female, copynumber_male)) {
                if (length(cn) != n_loci || anyNA(cn) || any(cn < 0)) {
                    stop("copynumber_female and copynumber_male must be non-negative, NA-free, and have one entry per row of host_matrix")
                }
            }
            if (verbose) {
                cat("Panel copynumber inputs are valid\n")
            }

            cn_panel <- pmax(copynumber_female, copynumber_male)
            self$loci <- loci
            self$panel_sex <- sex
            self$panel_copynumber <- cn_panel
            self$is_auto <- copynumber_female == copynumber_male & copynumber_female > 0
            if (!any(self$is_auto)) {
                stop("No autosomal/PAR loci (equal, non-zero copy number in both sexes) to compute size factors from")
            }
            if (verbose) {
                cat("Computing overall ploidy\n")
            }
            self$panel_ploidy <- median(cn_panel[self$is_auto])

            # reference: median over the samples at the panel copy number
            if (verbose) {
                cat("Median over samples\n")
            }
            cn <- sample_copynumber(copynumber_female, copynumber_male, sex)
            ref_valid <- (cn == cn_panel) & (cn > 0)

            n_ref <- rowSums(ref_valid)
            needed <- cn_panel > 0
            if (any(needed & n_ref == 0)) {
                stop(sprintf("%d loci have no samples of the required sex to build a reference; check the sex labels",
                             sum(needed & n_ref == 0)))
            }
            if (any(needed & n_ref < min_ref_samples)) {
                logwarn("%d loci have fewer than %d reference samples",
                        sum(needed & n_ref < min_ref_samples), min_ref_samples)
            }

            if (verbose) {
                cat("Scaling by per-sample size factors\n")
            }
            size_factors <- autosomal_size_factors(host_matrix, self$is_auto)
            scaled <- sweep(host_matrix, 2, size_factors, "/")
            scaled[!ref_valid] <- NA
            ref <- rowMedians(scaled, na.rm = TRUE)
            ref[!is.na(ref) & ref == 0] <- 1
            self$panel_reference <- ref

            # denoising
            # ratio to the panel reference, then bring each sample onto the panel scale
            # (e.g. male X x2); entries expected to be absent (female Y) are NA
            if (verbose) {
                cat("Computing the denoising basis by SVD\n")
            }
            if (verbose) {
                cat(" - Calculating ratios\n")
            }
            # Same as x <- self$ratio(host_matrix), but do it inline to avoid recomputing size_factors
            x <- sweep(host_matrix, 2, size_factors, "/")
            x <- sweep(x, 1, ref, "/")
            if (verbose) {
                cat(" - Scaling to panel copy number\n")
            }
            x <- x / (cn / cn_panel)
            if (verbose) {
                cat(" - Masking loci with expected copy number of 0 (setting to NA)\n")
            }
            x[cn == 0] <- NA

            # winsorise (NAs ignored when finding quantiles, and left as NA)
            if (verbose) {
                cat(" - Winsorization\n")
            }
            x <- apply(x, 2, function(col) clip(col,
                                                quantile(col, winsorization_quantile, na.rm = TRUE),
                                                quantile(col, 1 - winsorization_quantile, na.rm = TRUE)))

            # impute missing coverage (zeros) with the locus median
            if (verbose) {
                cat(" - Imputation of missing values\n")
            }
            ix <- which(x == 0, arr.ind = TRUE)
            x[ix] <- rowMedians(x, na.rm = TRUE)[ix[, "row"]]

            # renormalise on autosomal / PAR loci
            if (verbose) {
                cat(" - Renormalising\n")
            }
            x <- sweep(x, 2, colMedians(x[self$is_auto, , drop = FALSE]), "/")

            # move to log space; absent entries are neutral (logR 0)
            if (verbose) {
                cat(" - Moving to log space\n")
            }
            x <- safelog2(x)
            x[is.na(x)] <- 0

            stopifnot(is.matrix(x))

            # decompose
            if (verbose) {
                cat(" - Computing SVD\n")
            }
            decomp <- irlba(x, max_svs)
            self$U = decomp$u
            self$D = decomp$d
        },

        #' @description
        #' Calculate standardized (not denoised) read count ratio (R) from tumour read
        #' counts, relative to the poly-sex panel reference. Independent of tumour sex.
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        ratio = function(tumour_matrix) {
            stopifnot(nrow(tumour_matrix) == length(self$panel_reference))
            r <- sweep(tumour_matrix, 2, autosomal_size_factors(tumour_matrix, self$is_auto), "/")
            sweep(r, 1, self$panel_reference, "/")
        },

        #' @description
        #' Calculate standardized (not denoised) logR from tumour read counts
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        logR = function(tumour_matrix) {
            safelog2(self$ratio(tumour_matrix))
        },

        #' @description
        #' Calculate denoised logR from tumour read counts. Any NA entries are left out
        #' of the projection and returned as NA.
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        #' @param num_svs Number of singular values to use for denoising.
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
            missing <- is.na(x)
            x[missing] <- 0
            u <- self$U[, seq_len(num_svs), drop = FALSE]

            d <- x - u %*% (t(u) %*% x)
            d[missing] <- NA
            d
        }
    )
)


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
