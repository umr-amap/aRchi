#' Plot an object of class aRchi
#'
#' Plot an object of class aRchi in a 3d device. The QSM can be plotted according to different level of organization and the point cloud can be displayed if available.
#'
#' @export
#' @docType methods
#' @rdname plot
#' @description Plot an object of class aRchi.
#' @param x An aRchi object
#' @param y Unused (inherited from R base)
#' @param skeleton logical (Default is \code{TRUE}). Display the skeleton only (i.e segments). Faster than displaying the whole QSM with the fleshed cylinders.
#' @param leaves logical (Default is \code{FALSE}). Display the leaves ?
#' @param color The color of the cylinders. Can be either a single color or a level of organization:
#' "branching_order" for branching branching_order, "cylinder" to coloryze each cylinder independently, "segment" to coloryze the branch segments, "axis" to coloryze the axis, "A0" to colorize only the main axis from \code{\link{Compute_A0}} function
#' @param lwd line width of the skeleton
#' @param transparency The transparency of the cylinders
#' @param bg The background color
#' @param show_point_cloud logical (Default = \code{FALSE}). Display the point cloud ?
#' @param lit logical (Default is \code{TRUE}). Specify if lighting calculation should take place on the geometry. Only applies if \code{skeleton = FALSE}.
#' @include aRchiClass.R
#' @examples
#' \donttest{
#' # Read an aRchi file with at least a QSM
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Plot the QSM by coloring the branching order
#' plot(Tree1_aRchi,color="branching_order")
#' # Same with the fleshed cylinder and the point cloud
#' plot(Tree1_aRchi,color="branching_order",skeleton=FALSE,show_point_cloud=TRUE)
#'}
#'
#'
setMethod("plot",
          "aRchi",
          function(x,y,transparency=1,color = "white",bg="black",lwd = 3,show_point_cloud = FALSE,skeleton=TRUE,leaves=FALSE,lit=TRUE){

            pc_col=startX=endX=endY=startY=startZ=endZ=segment_ID=node_ID=parent_ID=radius=axis_ID=ID_Path=.=NULL


            # make sure the data is an aRchi file and contains a QSM
            if(class(x) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(x@QSM)) stop("This aRchi file does not contains a QSM")
            if(skeleton){
              skel = x@QSM

              if(color %in% c("cylinder","segment","axis","reconstructed","node","annual_shoots","physiological_age","branching_order","A0")){
                if (color=="cylinder") col = rep(skel$cyl_ID,each=2)
                if (color=="axis") col = rep(skel$axis_ID+1,each=2)
                if (color=="reconstructed") col = c("white","red")[rep(skel$reconstructed,each=2)]
                if (color=="segment") col = rep(skel$cyl_ID,each=2)
                if (color=="node") col = rep(skel$node_ID,each=2)+1
                if (color=="annual_shoots") col = rep(skel$AS,each=2)
                if (color=="branching_order") col = rep(skel$branching_order+2,each=2)
                if (color=="physiological_age") col = rep(skel$PA,each=2)
                if (color=="A0"){col = rep(skel$A0,each=2)
                col=ifelse(col==1,"white","red")}
              }else{
                col = color
              }

              # build the data for segment ploting
              dat_plot = data.frame(matrix(ncol=3,nrow=2*nrow(skel)))
              dat_plot[seq(1,nrow(dat_plot)-1,2),] = skel[,c(startX,startY,startZ)]
              dat_plot[seq(2,nrow(dat_plot),2),] = skel[,c(endX,endY,endZ)]
              dat_plot = data.table::data.table(X = dat_plot[,1],Y=dat_plot[,2],Z=dat_plot[,3])

              if(show_point_cloud){
                if(is.null(x@pointcloud)){warning("There is no point cloud to plot")}
                lidR::plot(x@pointcloud,bg=bg,clear_artifacts=FALSE,pal = height.colors,axis=T)
              }else{
                # translate the data for ploting with lidR tools
                #  dat_plot[,':='(X = X-min(X),Y = Y - min(Y),Z = Z-min(Z))]
                # empty window

                ########## modification to plot correct axes
                #pc=skel[startX==min(startX)|startX==max(endX)|startY==min(startY)|startY==max(endY)|startZ==min(startZ)|startZ==max(endZ),1:3]
                pc = skel[,.(startX,startY,startZ)]

                names(pc)=c("X","Y","Z")
                pc = pkgcond::suppress_messages( lidR::LAS(pc)) # pkgcond::supress_messages removes messages from the LAS building
                lidR::plot(pc,bg=bg,clear_artifacts=FALSE,axis=T,alpha = 0)
              }

              rgl::segments3d(dat_plot,lwd=lwd,col=col,add=TRUE)

            }
            if(skeleton==FALSE){
            QSM=x@QSM
            # particular color values
            if(color %in% c("branching_order","annual_shoots","axis","cylinder","segment","node","A0")){
              if (color=="branching_order"){col=QSM$branching_order+2}
              if (color=="annual_shoots") col = QSM$AS
              if (color=="cylinder"){col=QSM$cyl_ID}
              if (color=="axis") {col = QSM$axis_ID+1}
              if (color=="segment"){col=QSM$segment_ID+1}
              if (color=="node"){col=QSM$node+1}
              if (color=="A0"){col=ifelse(QSM$A0==1,"white","red")}
            }else{
              col = color
            }

            #ls_cyl=plyr::alply(QSM,1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
            mesh = aRchi::QSM2mesh(QSM)

            if(show_point_cloud){
              if(is.null(x@pointcloud)){warning("There is no point cloud to plot")}
              lidR::plot(x@pointcloud,bg=bg,clear_artifacts=FALSE,axis=T) #plot the point cloud
            }else{
              # empty window. to check
              pc=QSM[startX==min(startX)|startX==max(startX)|startY==min(startY)|startY==max(startY)|startZ==min(startZ)|startZ==max(startZ),1:3]
              names(pc)=c("X","Y","Z")
              pc = pkgcond::suppress_messages( lidR::LAS(pc)) # pkgcond::supress_messages removes messages from the LAS building
              lidR::plot(pc,bg=bg,size=0,clear_artifacts=FALSE,axis=T,col=bg)
            }
            #rgl::shapelist3d(ls_cyl,color=col,alpha=transparency,add=TRUE,lit=lit) # plot the list
            rgl::shade3d(mesh,col=rep(col,each=32))
            }
            if(leaves == TRUE){
              if(is.null(x@leaves$mesh)){stop("There is no point leaves to display")}
              rgl::shade3d(x@leaves$mesh,col = "forestgreen",add=T)
            }
          }
)





