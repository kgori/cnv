objective_maker <- function(cncalculator, segmentation, lambda) {
    ## Prefix sums (length = nrow(cncalculator$alpha) + 1)
    alpha_ps <- c(0, cumsum(cncalculator$alpha))
    beta_ps  <- c(0, cumsum(cncalculator$beta))

    ## Segment statistics
    segment_alpha <- (
        alpha_ps[segmentation$ends + 1] -
            alpha_ps[segmentation$starts + 1 - 1]
    ) / segmentation$lengths

    segment_beta <- (
        beta_ps[segmentation$ends + 1] -
            beta_ps[segmentation$starts + 1 - 1]
    ) / segmentation$lengths

    segment_lengths <- segmentation$lengths

    function(purity, ploidy) {
        c <- (1 - purity) / purity

        means <-
            segment_alpha * ploidy +
            segment_beta * c

        cn_int <- round(means)
        integer_loss <- sum(segment_lengths * (means - cn_int)^2)
        nonneg_loss <- sum(segment_lengths * (means - pmax(0, cn_int))^2)

        loss <- (1 - lambda) * integer_loss + lambda * nonneg_loss
        loss
    }
}

#' @export
find_local_minima <- function(grid) {
    purities <- sort(unique(grid$purity))
    ploidies <- sort(unique(grid$ploidy))
    nr <- length(purities)
    nc <- length(ploidies)
    m <- matrix(grid$result, nr, nc)
    minima <- matrix(FALSE, nr, nc)
    for (i in 1+seq_len(nr - 2)) {
        for (j in 1+seq_len(nc - 2)) {
            centre <- m[i, j]
            neighbours <- m[(i-1):(i+1), (j-1):(j+1)]
            minima[i, j] <- all(centre <= neighbours)
        }
    }
    min_ix <- which(minima, arr.ind = TRUE)
    mins <- m[min_ix]
    result <- cbind(
        purities[min_ix[, 1]],
        ploidies[min_ix[, 2]],
        mins
    )
    colnames(result) <- c("purity", "ploidy", "result")
    result[order(result[, 3]), ]
}

#' @export
grid_search <- function(
        copynumber_calculator,
        segmentation,
        purity_range = c(0.1, 1.0),
        ploidy_range = c(1, 8),
        purity_step = 0.01,
        ploidy_step = 0.05,
        nonneg_pen_lambda = 0.1) {
    stopifnot(all(purity_range >= 0 & purity_range <= 1))
    stopifnot(all(ploidy_range > 0))
    stopifnot(purity_step > 0 & ploidy_step > 0)

    grid <- expand.grid(
        purity = seq(purity_range[1], purity_range[2], purity_step),
        ploidy = seq(ploidy_range[1], ploidy_range[2], ploidy_step))

    obj <- objective_maker(copynumber_calculator, segmentation, nonneg_pen_lambda)

    grid$result <- mapply(
        obj,
        grid$purity,
        grid$ploidy
    )
    grid
}

#' @importFrom "dfoptim" nmkb
#' @export
neldermead <- function(
        startval,
        copynumber_calculator,
        segmentation,
        lower_bound = c(0, 1),
        upper_bound = c(1, 8),
        nonneg_pen_lambda = 0.1
) {
    # Objective function to optimise (1 arg version, takes x = c(pur, plo))
    fn <- objective_maker(copynumber_calculator, segmentation, nonneg_pen_lambda)
    obj <- function(x) fn(x[1], x[2])
    nmkb(startval, obj, lower = lower_bound, upper = upper_bound)
}

#' @importFrom "ggplot2" ggplot aes geom_tile geom_point scale_fill_viridis_c coord_fixed theme_bw
#' @export
plot_grid <- function(grid) {
    foundbest <- as.data.frame(grid[which.min(grid$result),])
    pur_width  <- min(diff(sort(unique(grid$purity))))
    plo_height <- min(diff(sort(unique(grid$ploidy))))
    minima <- find_local_minima(grid)
    
    ggplot(grid, aes(purity, ploidy, fill = -log10(result+1))) +
        geom_tile(width = pur_width, height = plo_height) +
        geom_point(data = minima[seq_len(min(10, nrow(minima))), ],
                   aes(purity, ploidy, size = -log10(result+1)),
                   inherit.aes = FALSE,
                   shape = 3,
                   colour = "black") +
        geom_point(data = foundbest,
                   aes(purity, ploidy),
                   inherit.aes = FALSE,
                   shape = 4,
                   colour = "blue",
                   size = 3) +
        scale_fill_viridis_c() +
        coord_fixed() +
        theme_bw()
}

#' @importFrom "R6" R6Class
#' @export
GridSearch <- R6Class(
    "GridSearch",
    public = list(
        #' @field copynumber_calculator A CopyNumberCalculator object.
        copynumber_calculator = NULL,
        #' @field segmentation A segmentation object, used when calculating the penalty function
        #' over the grid.
        segmentation = NULL,

        #' @field nonneg_pen_lambda Penalty for non-integer copy number solutions.
        nonneg_pen_lambda = NULL,

        #' @description
        #' Constructor. Builds new GridSearch object using a CopyNumberCalculator and a segmentation.
        #' The penalty calculated at each parameter combination is the sum of the squared differences between the calculated copy number and the nearest integer copy number, weighted by the segment lengths.
        #' We can additionally penalise negative integer solutions, by summing the squared differences between the calculated copy number and the nearest non-negative integer copy number, weighted by the segment lengths.
        #' The penalty term is a weighted sum of the integer penalty and the non-negative penalty, with the weighting factor given by nonneg_pen_lambda.
        #' @param copynumber_calculator A CopyNumberCalculator object.
        #' @param segmentation A segmentation object, as returned by segment_logR().
        #' @param purity_range Range of tumour purity values to search over.
        #' @param ploidy_range Range of tumour ploidy values to search over.
        #' @param purity_step Step size for tumour purity values.
        #' @param ploidy_step Step size for tumour ploidy values.
        #' @param nonneg_pen_lambda Penalty for non-integer copy number solutions. A weighting factor for the non-negative penalty term in the objective function.
        #' when lambda is 0, the objective function is purely the integer penalty, and when lambda is 1, the objective function is purely the non-negative penalty.
        initialize = function(
                copynumber_calculator,
                segmentation,
                nonneg_pen_lambda = 0
                ) {
            stopifnot(nonneg_pen_lambda >= 0 & nonneg_pen_lambda <= 1)
            self$copynumber_calculator <- copynumber_calculator
            self$segmentation <- segmentation
            self$nonneg_pen_lambda <- nonneg_pen_lambda
        },    

        #' @description
        #' Perform a grid search over the specified ranges of tumour purity and ploidy values.
        #' @param purity_range Range of tumour purity values to search over.
        #' @param ploidy_range Range of tumour ploidy values to search over.
        #' @param purity_step Step size for tumour purity values.
        #' @param ploidy_step Step size for tumour ploidy values.
        grid_search = function(
                purity_range = c(0.1, 1.0),
                ploidy_range = c(1, 8),
                purity_step = 0.01,
                ploidy_step = 0.05,
        ) {
            grid_search(
                self$copynumber_calculator,
                self$segmentation,
                purity_range,
                ploidy_range,
                purity_step,
                ploidy_step,
                self$nonneg_pen_lambda
            )
        },

        #' @description
        #' Optimise the objective function around a given starting point using the Nelder-Mead algorithm
        #' with box constraints.
        #' @param startval A numeric vector of length 2, containing the starting values for purity and ploidy.
        #' @param lower_bound A numeric vector of length 2, containing the lower bounds for purity and ploidy.
        #' @param upper_bound A numeric vector of length 2, containing the upper bounds for purity and ploidy.
        optimise_around_point = function(
                startval,
                lower_bound = c(0, 1),
                upper_bound = c(1, 8),
        ) {
            neldermead(
                startval,
                self$copynumber_calculator,
                self$segmentation,
                lower_bound,
                upper_bound,
                self$nonneg_pen_lambda
            )
        }

        #' @description
        #' Find local minima in the grid search results.
        #' @param grid A data frame containing the grid search results, as returned by grid_search
        find_local_minima = function(grid) {
            find_local_minima(grid)
        }
    )
)
