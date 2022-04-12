#' Generate a QSM in treeQSM format
#' @docType methods
#' @rdname aRchi2treeQSM
#' @description a data.table with treeQSM format (v2.3.3)
#' @param aRchi the object of class aRchi to write
#' @include aRchiClass.R


setGeneric("aRchi2treeQSM",
           function(aRchi){standardGeneric("aRchi2treeQSM")}
)

#' @rdname aRchi2treeQSM
#' @export

setMethod("aRchi2treeQSM",
          signature(aRchi="aRchi"),
          function(aRchi){

        radius_cyl=cyl_ID=length=startX=startY=startZ=axisX=axisY=axisZ=parent_ID=extension_ID=added=radius_cyl=axis_ID=branching_order=position_in_branch=NULL

        if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")

  aRchi_QSM=aRchi@QSM
  aRchi_QSM$axisX=-(aRchi_QSM$startX-aRchi_QSM$endX)/aRchi_QSM$length
  aRchi_QSM$axisY=-(aRchi_QSM$startY-aRchi_QSM$endY)/aRchi_QSM$length
  aRchi_QSM$axisZ=-(aRchi_QSM$startZ-aRchi_QSM$endZ)/aRchi_QSM$length
  aRchi_QSM$added=0
  aRchi_QSM$position_in_branch=0
  setorder(aRchi_QSM,cols=cyl_ID)
  aRchi_QSM=aRchi_QSM[,c("radius_cyl","length","startX","startY", "startZ","axisX","axisY","axisZ","parent_ID","extension_ID","added","radius_cyl","axis_ID","branching_order","position_in_branch")]
  return(aRchi_QSM)

}
)
