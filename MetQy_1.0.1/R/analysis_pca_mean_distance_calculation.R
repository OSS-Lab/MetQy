#' Given a set of coordinates, calculate the mean distance between all points.
#'
#' @param MATRIX - matrix with the rows refering to the points with N columns containing the coordinates (and therefore N dimensions).
#'
#' @details The mean distance of \eqn{p} points is calculated by the sum of the individual Euclidean distances
#' divided by the total numer of distances (given by \eqn{p*(p-1)/2}).
#'
#' @return The mean distance (\code{numeric}).
#' @export

############################################################################################################################################

analysis_pca_mean_distance_calculation <- function(MATRIX){

  stopifnot(ncol(MATRIX)>1)
  stopifnot(nrow(MATRIX)>1) # minimum two points

  p <- nrow(MATRIX)
  N <- sum(1:(p-1)) # equivalent to p*(p-1)/2
  distances <- numeric(length = N)
  index <- 0

  for(ii in 1:(p-1)){
    for(jj in (ii+1):p){
      index <- index + 1
      distances[index] <- sqrt(sum((MATRIX[jj,]-MATRIX[ii,])^2))
    }
  }
  mean_distance <- sum(distances)/N
  return(mean_distance)
}
