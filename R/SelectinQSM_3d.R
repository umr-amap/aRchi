#' SelectinQSM_3d
#'
#' @export
#' @docType methods
#' @rdname SelectinQSM_3d
#' @description Select interactively a sub-part of a QSM (cylinder, segment, node, axis, branch, subtree) in a 3d device and return its characteristics.
#' @param aRchi An object of class aRchi
#' @param skeleton logical. Display the skeleton only. Default is \code{TRUE}. Faster than displaying the QSM with the fleshed cylinders.
#' @param level character. \code{cylinder} (default), \code{segment}, \code{node}, \code{axis}, \code{branch}, \code{subtree}.
#' @return a data.table with the cylinders characteristics at the requested level (i.e sub-part of the original QSM).
#' @details
#'
#' The selection is performed in two times: i) Identifying the zone of interest in the 3d device and zoom into it if needed. When identified, the user has to hit enter in the R console. At this point, it is impossible to rotate the displayed QSM anymore as the left button of the mouse is used for the selection. However translation are still possible with the right button. ii) Draw a rectangle with the left button of the mouse in the zone of interest.
#'
#' Some details about the level of organization are given below.
#'
#'  "cylinder": return characteristics for the cylinders selected only
#'
#'  "segment:" return characteristics for the cylinders of the segments selected
#'
#'  "node": return the characteristics for the cylinders of the node selected. For a specific node, select the mother.
#'
#'  "axis": return the characteristics for the cylinders of the axis selected. An axis is a continuous succession of cylinder having a same branching order value.
#'
#'  "branch": return the characteristics for the cylinders of the branch selected. A branch is similar to an axis but regroup also everything that is upstream the axis (i.e all that the axis carries)
#'
#'  "Subtree": return the characteristics for the cylinders of the subtree selected. A subtree is similar to a branch but starting from the cylinder selected and not from the point of insertion of the selected axis. In other word, when the user draw a rectangle on a cylinder, the subtree selection return all that the cylinder carries. If several cylinders are selected, the subtree selection return all that the most downstream cylinder carries.
#' @examples
#' # Read an aRchi file with at least a QSM
#' if(interactive()){
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Select a branch
#' SelectinQSM_3d(Tree1_aRchi,level="branch")
#' # Same with the fleshed cylinder and keep the branch QSM in an object
#' My_branch=SelectinQSM_3d(Tree1_aRchi,level="branch",skeleton=FALSE)
#' My_branch
#' # Compute the moment of force
#' Tree1_aRchi=Compute_Mf(Tree1_aRchi,WoodDensity=550)
#' #Select a cylinder to return the moment of force at his position
#' SelectinQSM_3d(Tree1_aRchi,skeleton=FALSE)
#'}
#' @include aRchiClass.R
setGeneric("SelectinQSM_3d",
           function(aRchi,skeleton=TRUE,level="cylinder"){standardGeneric("SelectinQSM_3d")}
)

#' @rdname SelectinQSM_3d
#' @export

setMethod("SelectinQSM_3d",
          signature = "aRchi",
          function(aRchi,skeleton,level){
            axis_ID=cyl_ID=startX=startY=startZ=endX=endY=endZ=segment_ID=node_ID=parent_ID=radius=axis_ID=ID_Path=NULL

            QSM=aRchi@QSM
            if(level%in%c("branch","subtree")){
              Paths=aRchi@Paths
              if(is.null(Paths)){stop("Paths are needed for branch or subtree level. Please run Make_Path() on your object aRchi before selecting a branch or subtree in 3d.")}
            }
            if (skeleton==FALSE){
              dat_plot=plyr::alply(QSM,1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
            }

            if (skeleton==TRUE){
              # build the data for segment ploting
              dat_plot = data.frame(matrix(ncol=8,nrow=2*nrow(QSM)))
              dat_plot[seq(1,nrow(dat_plot)-1,2),] = QSM[,c(startX,startY,startZ,cyl_ID,segment_ID,node_ID,axis_ID,parent_ID)]
              dat_plot[seq(2,nrow(dat_plot),2),] = QSM[,c(endX,endY,endZ,cyl_ID,segment_ID,node_ID,axis_ID,parent_ID)]
              dat_plot = data.table::data.table(X = dat_plot[,1],Y=dat_plot[,2],Z=dat_plot[,3],cyl_ID=dat_plot[,4],segment_ID=dat_plot[,5],node_ID=dat_plot[,6],axis_ID=dat_plot[,7],parent_ID=dat_plot[,8])
            }





            if (interactive()) {
              pc=QSM[startX==min(startX)|startX==max(endX)|startY==min(startY)|startY==max(endY)|startZ==min(startZ)|startZ==max(endZ),1:3]
              names(pc)=c("X","Y","Z")
              pc = pkgcond::suppress_messages( lidR::LAS(pc)) # pkgcond::supress_messages removes messages from the LAS building
              lidR::plot(pc,bg="black",colorPalette="black",size=0,clear_artifacts=FALSE,axis=T)
              ifelse(skeleton,rgl::segments3d(dat_plot,lwd=3,col="white",add=TRUE), rgl::shapelist3d(dat_plot,color="white",alpha=1,add=TRUE,lit=TRUE))

              valid=gtools::ask("Find and Zoom into the zone of interest, then, hit Enter:\n")
              cat("Now select the zone of interest by drawing a rectangle.\n")
              f <- rgl::select3d()
              if (!is.null(f)) {
                keep <- which(f(QSM[,1:3]))
                if(length(keep)==0){
                  QSM$axisX=QSM$endX-QSM$startX
                  QSM$axisY=QSM$endY-QSM$startY
                  QSM$axisZ=QSM$endZ-QSM$startZ
                  incr=0
                  for (i in 1:100) {

                    incr=incr+0.05
                    QSM$x_1cm=QSM$startX+(incr*QSM$axisX)
                    QSM$y_1cm=QSM$startY+(incr*QSM$axisY)
                    QSM$z_1cm=QSM$startZ+(incr*QSM$axisZ)

                    keep <-which(f(QSM[,c("x_1cm","y_1cm","z_1cm")]))
                    if(length(keep)!=0){break}
                  }
                }

                if (skeleton==TRUE){
                  if(level=="cylinder"){
                  ls_keep=dat_plot[cyl_ID%in%QSM[keep]$cyl_ID]
                  ls_noKeep=dat_plot[!cyl_ID%in%QSM[keep]$cyl_ID]
                  }
                  if(level=="segment"){
                  ls_keep=dat_plot[segment_ID%in%QSM[keep]$segment_ID]
                  ls_noKeep=dat_plot[!segment_ID%in%QSM[keep]$segment_ID]
                  }
                  if(level=="node"){
                  node_cyl=rbind(QSM[segment_ID%in%unique(QSM[keep]$segment_ID)],QSM[node_ID%in%unique(QSM[keep]$segment_ID)])
                  ls_keep=dat_plot[cyl_ID%in%node_cyl$cyl_ID]
                  ls_noKeep=dat_plot[!cyl_ID%in%ls_keep$cyl_ID]
                  }
                  if(level=="axis"){
                    ls_keep=dat_plot[axis_ID%in%QSM[keep]$axis_ID ]
                    ls_noKeep=dat_plot[!axis_ID %in%QSM[keep]$axis_ID ]
                  }
                  if(level=="branch"){
                    Cyl_ID_start=min(QSM[axis_ID%in%unique(QSM[keep]$axis_ID)]$cyl_ID)
                    branch_cyl=QSM[cyl_ID%in%unique(Paths[ID_Path%in%Paths[cyl_ID==Cyl_ID_start]$ID_Path&cyl_ID>=Cyl_ID_start]$cyl_ID)]
                    ls_keep=dat_plot[cyl_ID%in%branch_cyl$cyl_ID]
                    ls_noKeep=dplyr::anti_join(dat_plot,ls_keep, by = c("X", "Y", "Z", "cyl_ID", "segment_ID", "node_ID", "axis_ID", "parent_ID"))
                  }
                  if(level=="subtree"){
                    Cyl_ID_start=min(unique(QSM[keep]$cyl_ID))
                    subtree_cyl=QSM[cyl_ID%in%unique(Paths[ID_Path%in%Paths[cyl_ID==Cyl_ID_start]$ID_Path&cyl_ID>=Cyl_ID_start]$cyl_ID)]
                    ls_keep=dat_plot[cyl_ID%in%subtree_cyl$cyl_ID]
                    ls_noKeep=dplyr::anti_join(dat_plot,ls_keep, by = c("X", "Y", "Z", "cyl_ID", "segment_ID", "node_ID", "axis_ID", "parent_ID"))
                  }
                  rgl::clear3d()
                  rgl::segments3d(ls_keep,lwd=3,col="red",add=TRUE)
                  rgl::segments3d(ls_noKeep,lwd=3,col="white",add=TRUE)


                }
                if (skeleton==FALSE){
                  if(level=="cylinder"){
                    cyl=QSM[keep]
                    ls_keep=dat_plot[which(QSM$cyl_ID%in%cyl$cyl_ID)]
                    ls_noKeep=dat_plot[-which(QSM$cyl_ID%in%cyl$cyl_ID)]
                  }
                  if(level=="segment"){
                    segment_cyl=QSM[segment_ID%in%unique(QSM[keep]$segment_ID)]
                    ls_keep=dat_plot[which(QSM$cyl_ID%in%segment_cyl$cyl_ID)]
                    ls_noKeep=dat_plot[-which(QSM$cyl_ID%in%segment_cyl$cyl_ID)]
                  }
                  if(level=="node"){
                    node_cyl=rbind(QSM[segment_ID%in%unique(QSM[keep]$segment_ID)],QSM[node_ID%in%unique(QSM[keep]$segment_ID)])
                    ls_keep=dat_plot[which(QSM$cyl_ID%in%node_cyl$cyl_ID)]
                    ls_noKeep=dat_plot[-which(QSM$cyl_ID%in%node_cyl$cyl_ID)]
                  }
                  if(level=="axis"){

                    axis_cyl=QSM[axis_ID%in%unique(QSM[keep]$axis_ID)]
                    ls_keep=dat_plot[which(QSM$cyl_ID%in%axis_cyl$cyl_ID)]
                    ls_noKeep=dat_plot[-which(QSM$cyl_ID%in%axis_cyl$cyl_ID)]
                  }
                  if(level=="branch"){
                    Cyl_ID_start=min(QSM[axis_ID%in%unique(QSM[keep]$axis_ID)]$cyl_ID)
                    branch_cyl=QSM[cyl_ID%in%unique(Paths[ID_Path%in%Paths[cyl_ID==Cyl_ID_start]$ID_Path&cyl_ID>=Cyl_ID_start]$cyl_ID)]
                    ls_keep=dat_plot[which(QSM$cyl_ID%in%branch_cyl$cyl_ID)]
                    ls_noKeep=dat_plot[-which(QSM$cyl_ID%in%branch_cyl$cyl_ID)]
                  }
                  if(level=="subtree"){
                    Cyl_ID_start=min(unique(QSM[keep]$cyl_ID))
                    subtree_cyl=QSM[cyl_ID%in%unique(Paths[ID_Path%in%Paths[cyl_ID==Cyl_ID_start]$ID_Path&cyl_ID>=Cyl_ID_start]$cyl_ID)]
                    ls_keep=dat_plot[which(QSM$cyl_ID%in%subtree_cyl$cyl_ID)]
                    ls_noKeep=dat_plot[-which(QSM$cyl_ID%in%subtree_cyl$cyl_ID)]
                    }

                  rgl::clear3d()
                  rgl::shapelist3d(ls_keep,color="red",alpha=1,add=TRUE,lit=TRUE)

                  if(length(ls_noKeep)!=0){
                  rgl::shapelist3d(ls_noKeep,color="white",alpha=1,add=TRUE,lit=TRUE)}


                }

              }

            }
            if(level=="cylinder"){return(QSM[keep])}
            if(level=="segment"){return(QSM[segment_ID%in%unique(QSM[keep]$segment_ID)])}
            if(level=="node"){return(node_cyl)}
            if(level=="axis"){return(QSM[axis_ID%in%unique(QSM[keep]$axis_ID)])}
            if(level=="branch"){return(branch_cyl)}
            if(level=="subtree"){return(subtree_cyl)}
          }
)


