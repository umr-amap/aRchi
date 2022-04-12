

#' Transform a QSM into a triangular or quadrangular mesh
#'
#' @param qsm a qsm.
#' @param tmesh logical. If TRUE a triangular mesh is returned, otherwise a quadrangular mesh is returned.
#' @param sides integer. The number of faces that compose a cylinder.
#'
#' @return a 3D mesh.
#' @export
#'
#' @keywords internal

QSM2mesh <- function(qsm,tmesh = FALSE,sides = 16){

  # to pass CRAN checks
  startX=startY=startZ=endX=endY=endZ=radius_cyl=.=NULL

  # Build a matrix with two lines by cylinder
  mat = matrix(ncol=4,nrow = 2*nrow(qsm))
  mat[seq(1,nrow(mat)-1,2),1:4] = as.matrix(qsm[,.(startX,startY,startZ,radius_cyl)])
  mat[seq(2,nrow(mat),2),1:4] = as.matrix(qsm[,.(endX,endY,endZ,radius_cyl)])

  # build empty matrix for vertices and indices
  vb = matrix(ncol = nrow(qsm)*sides*2,nrow=4)
  ib = matrix(ncol = nrow(qsm)*sides,nrow=4)

  # for each cylinder store the vertices and indices
  n_vb = 1 ; n_ib = 1 ; iter = 1
  for(i in seq(1,nrow(mat-1),2)){

    cyl = rgl::cylinder3d(mat[i:(i+1),1:3],radius= mat[i,4],sides=sides)

    vb[1:4,n_vb:(n_vb+(sides*2-1))] = cyl$vb
    ib[1:4,n_ib:(n_ib+(sides-1))] = cyl$ib + (iter-1)*sides*2
    n_vb = n_vb+sides*2
    n_ib = n_ib+sides
    iter = iter+1
  }

  # build the mesh
  if(tmesh){
    # if the desired output is a triangular mesh
    it = matrix(ncol= 2*ncol(ib),nrow = 3)
    it[1:3,seq(1,ncol(it)-1,2)] = ib[1:3,]
    it[1:3,seq(2,ncol(it),2)] = ib[c(1,3,4),]

    return(rgl::tmesh3d(vertices = vb, indices = it))
  }else{
    # if the desired output is a quadrangular mesh
    return(rgl::qmesh3d(vertices = vb, indices = ib))
  }
}





