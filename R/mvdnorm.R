#' Bivariate Guassian Density Computation Compute the density of a point under a
#' bivariate gaussian model
#' @param p A (x,y) vector of real numbers, i.e. a point.
#' @param mu A vector indicating the mean.
#' @param sigma Covariance Matrix
#' @return The density at point p of a bivariate gaussian with mean mu and
#'   covariance matrix sigma.
twod_dnorm <- function(p, mu, sigma){
 (2*pi)^{-1} * det(sigma)^{-1/2} * exp(-1/2 * t(p-mu) %*% solve(sigma) %*% (p-mu))
}
