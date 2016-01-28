#' Single-Point One-Dimensional Gausssian Kernel Density Estimate
#' @description Computes the one-dimensional guassian KDE for a single point.
#' @param p A real number
#' @param d A vector of observations
#' @param bw Bandwidth
#' @return A numeric density estimate
KDE_height <- function(p, d, bw){
  sum(dnorm(0, mean = (p-d)/bw))/(bw * length(d))
}

#' General One-Dimensional Gausssian Kernel Density Estimate
#' @description Computes the one-dimensional guassian KDE for an entire region.
#' @param d A vector of observations.
#' @param bw Bandwidth. Default is the unbiased cross-validation estimate.
#' @param dmin The lowest density you'd like to estimate. Default somewhat
#'   arbitrary. Used two-standard deviations from the relevant value.
#' @param dmax The highest quantity you'd like to estimate. Default somewhat
#'   arbitrary. Used two-standard deviations from the relevant value.
#' @param accuracy Numer of bins you'd like. Default somewhat arbitrary.
#' @return A matrix with locations and heights.
#' @examples
#' x <- rnorm(100)
#' a <- KDE_1(x)
#' plot(a, type ='l')
#'
#' y <- rgamma(30, 2, 4)
#' KDE_1(x, 0.25)
#' b <- KDE_1(y, 5)
#' c <- KDE_1(y)
#' plot(b, type = 'l')
#' lines(c, type = 'l', col = 'red')
KDE_1 <- function(d, bw = my_ucv(d)$par, dmin = min(d) - bw * abs(pnorm(0.025, min(d)/bw)), dmax = max(d) + bw * abs(pnorm(0.975, max(d)/bw)), accuracy = 1000){
  bins <- seq(dmin, dmax, length.out = accuracy)
  heights <- sapply(bins, KDE_height, d = d, bw = bw)
  return(cbind(bins, heights))
}
