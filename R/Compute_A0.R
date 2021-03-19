#' Find the principal axis of a tree
#'
#' Find the principal axis of a tree (i.e A0) and add a column to the QSM of an aRchi object. This is alternative method to the default branch order proposed in a QSM for the principal axis only.
#'
#' @export
#' @docType methods
#' @rdname Compute_A0
#' @param aRchi an object of class aRchi with at least a QSM and the Paths table
#' @param plotresult logical (default = FALSE). Show the results in a 3d plot if TRUE
#' @return The aRchi object with the QSM having a new column A0.
#' @details
#'
#' The method used to find the principal axis consist in finding the highest vertical path with consecutive segments of similar diameters and orientations. An index called A0 of the probability of being the principal axis is thus computed for each path of the tree and the path with the highest value is considered as the principal axis (see Martin-Ducup et al. 2020 for more information).
#'
#' A0 ranges between 0 and 4 with 0 indicating a path with a low probability of being the principal axis and 4 indicating a high probability of being the principal axis, i.e. the highest vertical path with consecutive segments of similar diameters and orientations.The path with the maximum A0 value was selected as the principal axis.
#'
#' The new column A0 of the QSM slot take the value 2 if the cylinder is part of the principal axis or 1 if not.
#'
#' #'@references
#'
#' 	Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' @seealso \code{\link{DAI}} to compute the dominance of a principal axis index that uses the A0 index.
#' @include aRchiClass.R
#'
#' @examples
#' \dontrun{
#' # Read an aRchi file with at least the QSM and the paths table
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' Tree1_aRchi=Compute_A0(Tree1_aRchi,plotresult=TRUE)
#'}
setGeneric("Compute_A0",
function(aRchi,plotresult=FALSE){standardGeneric("Compute_A0")}
)


#' @rdname Compute_A0
#' @export
#'
setMethod("Compute_A0",
          signature = "aRchi",

          function(aRchi,plotresult){
            segment_ID=ID_Path=V1=cyl_ID=A0=H_segment_rel=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(aRchi@Paths)) stop("The archi file does not contains Paths")

            QSM=aRchi@QSM
            Paths=aRchi@Paths
            # Make a table with information at the segment level (i.e length, radius of the first cylinder, radius of the last cylinder of the parent, length of the parent, angle and angle of the parent). This table is needed to identify A0
            segment_table=plyr::ddply(QSM,("segment_ID "),function(x){
              segment= x
              segment_P=QSM[segment_ID==unique(x$node_ID)]
              cbind(
                segment_ID=unique(segment$segment_ID),
                Angle=circular::deg(angle3d(as.numeric(segment[1,c("startX","startY","startZ")]),as.numeric(c(segment[1,c("startX","startY")],segment[1,"startZ"]+1)),as.numeric(segment[nrow(segment),c("endX","endY","endZ")]))),
                Angle_P=circular::deg(angle3d(as.numeric(segment_P[1,c("startX","startY","startZ")]),as.numeric(c(segment_P[1,c("startX","startY")],segment_P[1,"startZ"]+1)),as.numeric(segment_P[nrow(segment_P),c("endX","endY","endZ")]))),
                radius=segment[1,"radius_cyl"],
                radius_parent=segment_P[nrow(segment_P),"radius_cyl"],
                length=sum(segment$length),
                length_parent=sum(segment_P$length),
                z_min=segment[1,"startZ"],
                z_max=segment[nrow(segment),"endZ"])},.progress = "text")
            names(segment_table)[5]="radius_parent"

            # Merge segment table and paths
            Identification_trunk_table=merge(segment_table,unique(Paths[,c("segment_ID","node_ID","ID_Path")]),by="segment_ID")

            # Compute angle asymetry between a segment and its parent
            Identification_trunk_table=plyr::adply(Identification_trunk_table,1,function(x){1-(abs(x$Angle_P-x$Angle)/max(x$Angle_P,x$Angle))})
            names(Identification_trunk_table)[12]="AS_angle"
            # Compute radius asymmetry between a segment and its parent
            Identification_trunk_table=plyr::adply(Identification_trunk_table,1,function(x){AS_diam=1-(abs(x$radius_parent-x$radius )/max(x$radius_parent,x$radius))})
            names(Identification_trunk_table)[13]="AS_diam"
            Identification_trunk_table$vert_angle=1-(Identification_trunk_table$Angle/90)

            # Minimum Height
            H_min=min(segment_table$z_min)
            # Segment Height
            Identification_trunk_table$H_segment_rel=Identification_trunk_table$z_max-H_min
            #Path height
            H_path=plyr::ddply(Identification_trunk_table,("ID_Path"),plyr::summarise,H_path=max(H_segment_rel))
            # Merge path Height to tablle with segment and path ID
            Identification_trunk_table=merge(Identification_trunk_table,H_path,by="ID_Path")
            #Highest path
            max_H_path=max(Identification_trunk_table$H_path)
            # relative path height
            Identification_trunk_table$H_path_rel=Identification_trunk_table$H_path/max_H_path
            # A0 index (added) per path
            Index_A1_add=data.table::data.table(plyr::ddply(Identification_trunk_table,("ID_Path"),function(x){mean(x$AS_angle+x$AS_diam+x$vert_angle+x$H_path_rel,na.rm=T)}))
            # A0 index (multiplied) per path
            Index_A1_mult=data.table::data.table(plyr::ddply(Identification_trunk_table,("ID_Path"),function(x){mean(x$AS_angle*x$AS_diam*x$vert_angle*x$H_path_rel,na.rm=T)}))

            # A0 cylinder ID
            A1_cyl_ID=Paths[ID_Path==Index_A1_add[V1==max(V1)]$ID_Path]$cyl_ID

            QSM$A0=1
            QSM[cyl_ID%in%A1_cyl_ID]$A0=2


            if(plotresult){

              pc = pkgcond::suppress_messages( lidR::LAS(data.frame(X=mean(QSM$startX),Y=mean(QSM$startY),Z=mean(QSM$startZ)))) # pkgcond::supress_messages removes messages from the LAS building
              lidR::plot(pc,bg="black",colorPalette="black",size=0,clear_artifacts=F)

              ls_cyl=plyr::alply(QSM[A0==2],1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
              rgl::shapelist3d(ls_cyl,color=2,alpha=1,add=T,lit=T) # plot the list

              ls_cyl=plyr::alply(QSM[A0==1],1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
              rgl::shapelist3d(ls_cyl,color="white",alpha=1,add=T,lit=T) # plot the list
              rgl::bbox3d(color="white")

            }

            aRchi@QSM=QSM
            aRchi@operations$Compute_A0=c(Computed="YES")
            return(aRchi)



          }
)


