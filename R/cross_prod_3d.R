
#' Internal function that performs a cross product for the dist2line function
#'
#' @param v1 vector 1
#' @param v2 vector 2
#'
#' @export
#'
#' @keywords internal

cross_prod_3d <- function(v1,v2){
  v3 <- vector()
  v3[1] <- v1[2]*v2[3]-v1[3]*v2[2]
  v3[2] <- v1[3]*v2[1]-v1[1]*v2[3]
  v3[3] <- v1[1]*v2[2]-v1[2]*v2[1]
  return(v3)
}
