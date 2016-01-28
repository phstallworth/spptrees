#' Oversmoothed bandwidth Estimate for One-Dimensional Gaussian KDE
#'
#' @param d A vector
#' @return Computes the oversmoothed bandwidth estimate for the one-dimensional
#'   Guassian kernel density estimate.
my_oversmooth <- function(d) 1.06*sd(d)*length(d)^(-1/5)

#' Unbiased Cross-Validation for One-Dimensional Gaussian Kernel Bandwidth
#' Selection.
#'
#' @param d A vector.
#' @param lower A real number.
#' @param upper A real number.
#' @return Computes the unbiased cross-validation bandwidth of a one-dimensional
#'   gaussian kernel density estimator.
my_ucv <- function(d, lower = 0.0025, upper= 2){
  optim(1, function(h, d) hucv(h, d), d = d, lower = lower, upper = upper)
}

#' Biased Cross-Validation for One-Dimensional Gaussian Kernel Bandwidth
#' Selection.
#'
#' @param d A vector.
#' @param lower A real number.
#' @param upper A real number.
#' @return Computes the biased cross-validation bandwidth of a one-dimensional
#'   gaussian kernel density estimator.
my_bcv <- function(d, lower = 0.0025, upper= my_oversmooth(d)){
  optim(1, hbcv(h, d), d = d, lower = lower, upper = upper)
}

delta <- function(x, y, h) (x-y)/h

hucv <- function(h, d){
  c <- 0
  n <- length(d)
  for(i in 2:length(d)){
    j = 1
    while(j != i){
      c <- c + exp(-delta(d[i], d[j], h)^2/4) - sqrt(8)*exp(-delta(d[i], d[j], h)^2/2)
      j <- j + 1
    }
  }
  return(1/(2*n * h * sqrt(pi)) + 1/(n^2*h*sqrt(pi)) * c)
}

hbcv <- function(h, d){
  c <- 0
  n <- length(d)
  for(i in 2:length(d)){
    j = 1
    while(j != 1){
      c <- c + (delta(d[i], d[j], h)^4 - 12 * delta(d[i], d[j], h)^2 + 12)*exp^(-delta(d[i],d[j],h)/4)
    }
    return(1/(2*n*h*sqrt(pi))+1/(64*n^2*h*sqrt(pi))*c)
  }
}






