
#' Internal function to compute the distance of point b to segment ac
#'
#' @param a a vector containing the xyz coordinates of a point
#' @param b a vector containing the xyz coordinates of a point
#' @param c a vector containing the xyz coordinates of a point
#'
#' @return the distance of point b to segment ac
#' @export
#'
#' @keywords internal


dist2line <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  v3 <- cross_prod_3d(v1,v2)
  area <- sqrt(sum(v3*v3))/2
  d <- 2*area/sqrt(sum(v1*v1))
  return(d)
}

