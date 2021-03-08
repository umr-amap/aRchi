#' Calculate the zenith angle from xyz coordinates
#'
#' Calculate a zenith angle between two segments from 3d (i.e x,y,z) coordinates
#'
#' @param o 3d coordinates of the common point of the two segments
#' @param a 3d coordinates of the other point of segment a
#' @param b 3d coordinates of the other point of segment b
#' @return  The angle in degree
#' @include aRchiClass.R
#' @examples
#' origin=c(0,0,0)
#' a=c(0,0,1)
#' b=c(1,0,0)
#'
#' angle3d(o=origin,a=a,b=b)
#'

angle3d=function (o,a,b) {

  xo=o[1]
  yo=o[2]
  zo=o[3]

  xa=a[1]
  ya=a[2]
  za=a[3]

  xb=b[1]
  yb=b[2]
  zb=b[3]

  v_OA = c((xa-xo),(ya-yo),(za-zo))
  v_OB = c((xb-xo),(yb-yo),(zb-zo))

  n_OA = sqrt(v_OA[1]^2+v_OA[2]^2+v_OA[3]^2)
  n_OB = sqrt(v_OB[1]^2+v_OB[2]^2+v_OB[3]^2)

  scal_OAOB= (v_OA[1]*v_OB[1])+(v_OA[2]*v_OB[2])+(v_OA[3]*v_OB[3])

  theta=acos(scal_OAOB/(n_OA*n_OB))
  return(theta=theta)}
