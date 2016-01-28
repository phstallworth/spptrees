#' Calculate non-parametric intensity of an IPP at a point
#'
#' @param pt_x An x value
#' @param pt_y A y value
#' @param my_set A set containing information on validation set(needs to be removed),
#'  lattice points, 1-st tree neighbor, and distance from lattice to tree neighbor
#' @param bw The desired bw.
#' @param my_region The portion of the space already explored
#' @param my_Area The area of the explored space
#' @return A non-parametric estimate of the intensity at the point (pt_x, pt_y)
lambda_hat <- function(pt_x, pt_y, my_set, bw, my_region, my_Area){
  pt <- c(pt_x, pt_y)
  num_lambda_hat(pt, my_set, bw)/lambda_A_checked(pt, my_region, my_Area, bw)
}

#' Calculates the denominator for lamda_hat
lambda_A_checked <- function(pt, my_region, my_Area, bw){
  mu <- 0
  sigma <- matrix(c(bw, 0, 0, bw), nrow = 2)
  for(i in 1:nrow(my_region)){
    mu <- (mu * (i - 1) + twod_dnorm(my_region[i,], pt, sigma))/i
  }
  mu * my_Area
}

#' Calculates the numerator of the non-parametric intensity estimate
num_lambda_hat <- function(pt, my_set, bw){
  sum <- 0
  covMatrix <- matrix(c(bw, 0, 0, bw), nrow = 2)
  for(i in 1:nrow(my_set)){
    sum <- sum + twod_dnorm(my_set[i, 4:5], pt, covMatrix)
  }
  sum
}

#' Drop points in the x-y plane
drop_points <- function(npoints, xmax, xmin, ymax, ymin){
  cbind(runif(npoints, xmin, xmax), runif(npoints, xmin, xmax))
}

#' Discover whether points are in a given region
#'
#' @param pts A two-column matrix with x-values in column 1 and y-values in column 2.
#' @param region A set containing information on validation set(needs to be removed),
#'  lattice points, 1-st tree neighbor, and distance from lattice to tree neighbor
in_Area <- function(pts, region){
  v <- rep(NA, nrow(pts))
  for(i in 1:length(v)){
    for(j in 1:nrow(region)){
      if(d2(pts[i, 1:2], region[j, 2:3]) < region[j, 6]){
        v[i] <- T
      }
    }
  }
  return(v)
}



