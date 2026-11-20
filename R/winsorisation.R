#' Copied from ASCAT
medianFilter <- function (x, k)
{
    n <- length(x)
    filtWidth <- 2 * k + 1
    if (filtWidth > n) {
        if (n == 0) {
            filtWidth <- 1
        }
        else if (n%%2 == 0) {
            filtWidth <- n - 1
        }
        else {
            filtWidth <- n
        }
    }
    runMedian <- runmed(x, k = filtWidth, endrule = "median")
    return(runMedian)
}

#' Copied from ASCAT
madWins <- function (x, tau, k)
{
    xhat <- medianFilter(x, k)
    d <- x - xhat
    SD <- mad(d)
    z <- tau * SD
    xwin <- xhat + psi(d, z)
    outliers <- rep(0, length(x))
    outliers[x > xwin] <- 1
    outliers[x < xwin] <- -1
    return(list(ywin = xwin, sdev = SD, outliers = outliers))
}

#' Copied from ASCAT
psi <- function(x, z) {
    xwin <- x
    xwin[x < -z] <- -z
    xwin[x > z] <- z
    return(xwin)
}

#' Winsorizes an input vector. Winsorizing means to clip outlier values to a level closer to the central tendency of the distribution.
#' @param tau Multiplier of the local standard deviation: smaller values result in a stronger clipping effect
#' @param k The local window size
#' @export
winsorize <- function(x, tau=2.5, k=5) {
    madWins(x, tau, k)$ywin
}
