#' Create a vector with validation set information
#'
#' @param sets_number How many validation sets.
#' @param total Total number of items.
val_sets <- function(sets_number, total){
  sample(1:total, size = total, replace = FALSE) %/% sets_number + 1
}

#' Compute the distance between two points
#'
#' @param a,b 2-dimeninsional vectors
d2 <- function(a, b) sqrt((a[1]-b[1])^2 + (a[2]- b[2])^2)

#' Find the data point nearest to (x,y)
#'
#' @param x,y Real numbers.
#' @param d A two column  matrix with each row a point
nearest_point <- function(x, y, d){
  my_min <- c(NA, NA, Inf)
  for(i in 1:nrow(d)){
    current_distance <- d2(c(d[i,1], d[i,2]), c(x, y))
    current_min <- my_min[3]
    if(current_distance < current_min){
      my_min <- c(d[i,1], d[i,2], current_distance)
    }
  }
  my_min
}

#' Find the duplicate rows in a "final" matrix
#'
#' @param my_matrix A 6 row matrix with the first row validation set info,
#'    the second and third row lattice points, the 4th and 5th row
#'    nearest tree neighbor, and the 6th distance
#' @return A vector with the row number of each non-duplicate entry
find_duplicates <- function(my_matrix){
  for_keeps <- 1:nrow(my_matrix)
  for(i in 2:nrow(my_matrix)){
    for(j in 1:(i-1)){
      if(my_matrix[j, 4] == my_matrix[i, 4] & my_matrix[j, 5] == my_matrix[i, 5]){
        for_keeps[i] <- NA
      }
    }
  }
  for_keeps[!is.na(for_keeps)]
}

#' Compute the exponentiated part of the bandwidth optimizer
exp_value <- function(validationRegion, validationRegionArea, checkedRegion, checkedArea, checkset, bw){
  mu <- 0
  for(i in 1:nrow(validationRegion)){
    mu <- (mu * (i - 1) + lambda_hat(validationRegion[i, 1], validationRegion[i, 2], checkset, bw, checkedRegion, checkedArea)) / i
  }
  mu * validationRegionArea
}

#' Compute the product part of the bandwidth optimizer
prod_part <- function(valid_set, checkedRegion, checkedArea, checkset, bw){
  prod <- 1
  for(i in 1:nrow(valid_set)){
    prod <- prod * lambda_hat(valid_set[i, 4], valid_set[i,5], checkset, bw, checkedRegion, checkedArea)
  }
  prod
}


#' Compute the full bandwidth optimizer for a single validation set
final_comp <- function(validSet, checkedSet, validationRegion, validationArea, checkedRegion, checkedArea, bw){
  exp(-exp_value(validationRegion, validationArea, checkedRegion, checkedArea, checkedSet, bw)) * prod_part(validSet, checkedRegion, checkedArea, checkedSet, bw)
}

#' Perform the whole validation set bandwidth probability computation from start to finish
prob_O_Si <- function(val_num, d, pts, bw, xmax, xmin, ymax, ymin){
  #Break up into checksets and valset
  checkset <- d[d[,1] != val_num,]
  valset <- d[d[,1] == val_num,]
  #Learn Regions
  checked_Region <- inCheckedRegion(pts, checkset)
  validation_Region <- inValidationRegion(pts, valset, checkset)
  check_region_area <- nrow(checked_Region)/nrow(pts) * (xmax - xmin) * (ymax - ymin)
  validation_region_area <- nrow(validation_Region)/nrow(pts) * (xmax - xmin) * (ymax - ymin)
  #Remove duplicates
  checkset <- checkset[find_duplicates(checkset),]
  valset <- valset[find_duplicates(valset),]
  final_comp(valset, checkset, validation_Region, validation_region_area, checked_Region, check_region_area, bw)
}

#' Compute the score for a specific bandwidthfunction
bw_level <- function(bw, nvalsets, d, pts, xmax, xmin, ymax, ymin){
  prod(sapply(seq(1:nvalsets), prob_O_Si, d, pts, bw, xmax, xmin, ymax, ymin))
}

#' Find the optimal bandwidth.
#'  @description Currently to slow to work well
opt_bw <- function(npts, nvalsets, d, xmax, xmin, ymax, ymin, min_guess, max_guess, step_value){
  mypts <- drop_points(npts, xmax = xmax, xmin = xmin, ymax = ymax, ymin = ymin)
  v <- sapply(seq(min_guess, max_guess, by = step_value), bw_level, nvalsets, d, mypts, xmax, xmin, ymax, ymin)
  min_guess + step_value * which.max(v)
}

#' Determine which points are in the checked region
inCheckedRegion<- function(pts, checkset){
  pts[!is.na(in_Area(pts, checkset)),]
}

#' Determine which points are in the validation region
inValidationRegion <- function(pts, valset, checkset){
  pts[!is.na(in_Area(pts, valset)) & is.na(in_Area(pts, checkset)),]
}


