#' @export
tumour_copy_number <- function(
        R,
        purity,
        tumour_ploidy,
        host_copynumber = 2,
        host_ploidy = 2,
        panel_copynumber = 2,
        panel_ploidy = 2) {
    stopifnot(
        purity >= 0,
        purity <= 1,
        tumour_ploidy > 0,
        host_ploidy > 0,
        panel_ploidy > 0,
        host_copynumber >= 0,
        panel_copynumber >= 0
    )
    rnp <- R * panel_copynumber / panel_ploidy
    rnp * tumour_ploidy + (1 - purity) / purity *
        (rnp * host_ploidy - host_copynumber)
}

#' @export
prob_read_came_from_tumour <- function(
        R,
        purity,
        tumour_ploidy,
        host_copynumber = 2,
        host_ploidy = 2,
        panel_copynumber = 2,
        panel_ploidy = 2) {
    stopifnot(
        purity >= 0,
        purity <= 1,
        tumour_ploidy > 0,
        host_ploidy > 0,
        panel_ploidy > 0,
        host_copynumber >= 0,
        panel_copynumber >= 0
    )
    alpha <- R * panel_copynumber / panel_ploidy
    mixture_ploidy <-
        purity * tumour_ploidy +
        (1 - purity) * host_ploidy
    mixture_copynumber <- alpha * mixture_ploidy

    prob <- 1 - (1 - purity) * host_copynumber / mixture_copynumber
    pmax(0, pmin(1, prob))
}

#' @export
prob_read_came_from_host <- function(
        R,
        purity,
        tumour_ploidy,
        host_copynumber = 2,
        host_ploidy = 2,
        panel_copynumber = 2,
        panel_ploidy = 2) {
    1 - prob_read_came_from_tumour(R, purity,
                                   tumour_ploidy, host_copynumber, host_ploidy,
                                   panel_copynumber, panel_ploidy)
}
