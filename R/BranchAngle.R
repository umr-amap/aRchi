#' Estimation of the tree branch angle from an aRchi file
#'
#' Estimate the branch angle of a QSM. Two methods are possible (see method argument)
#'
#' @export
#' @docType methods
#' @rdname BranchAngle
#' @param aRchi an object of class aRchi with at least a QSM and a path table
#' @param method character. \code{SegmentAngle} or code{King98}
#' @param level character. The level at which the branch angle is computed. \code{Tree} for tree level; \code{branching_order} for branch order level; \code{Axis} one angle value per axis.
#' @param A0 logical (default = FALSE). If TRUE the main axis to remove from the calculation is re-estimated using the \code{\link{Compute_A0}} function. If false the default branch order 0 is kept.
#' @return a numeric or data.table. The branch angle in degree at the selected level. with 0 a perfectly vertical branch angle, 90 a perfectly horizontal branch angle and >90 a downward branch angle
#' @details
#'
#' The method "SegmentAngle" compute the angle by considering the first and the last cylinder or each segment, mean is then used for the level of organization selected.
#'
#' The method "King98" compute the angle by considering the first and the last cylinder of each axis mean is then used for the level of organization selected.
#'
#' The main axis is always removed.

#'@references
#'
#' 	Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' @include aRchiClass.R
#' @examples
#' \donttest{
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Compute the branch angle at various level
#' BranchAngle(Tree1_aRchi,method="SegmentAngle")
#' BranchAngle(Tree1_aRchi,level="branching_order",method="SegmentAngle",A0=TRUE)
#'}
setGeneric("BranchAngle",
function(aRchi,method=NULL,A0=FALSE,level="Tree"){standardGeneric("BranchAngle")}
)


#' @rdname BranchAngle
#' @export
setMethod("BranchAngle",
          signature = "aRchi",

          function(aRchi,method,A0,level){
            branching_order=cyl_ID=segment_ID=Angle=axis_ID=branch_angle=.=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(aRchi@Paths)) stop("The archi file does not contains Paths")
            if(is.null(method)) stop("Please choose a method, either SegmentAngle or King98")
            if(method %in% c("SegmentAngle","King98")==FALSE) stop("Incorrect method argument")
            if(level %in% c("Tree","branching_order","Branch")==FALSE) stop("Incorrect level argument")
            if(A0==TRUE){
              aRchi=Compute_A0(aRchi,plotresult=FALSE)
              QSM=aRchi@QSM
              cyl_ID_A1=dplyr::anti_join(QSM[A0==2],QSM[branching_order==0],by="cyl_ID")$cyl_ID
              QSM[cyl_ID%in%cyl_ID_A1]$branching_order=0
              cyl_ID_A1F=dplyr::anti_join(QSM[branching_order==0],QSM[A0==2],by="cyl_ID")$cyl_ID
              QSM[cyl_ID%in%cyl_ID_A1F]$branching_order=1
              aRchi@QSM=QSM
            }
            QSM=aRchi@QSM
            Paths=aRchi@Paths

            QSM_branch_angle=QSM[branching_order!=0]
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
                z_max=segment[nrow(segment),"endZ"])})
            names(segment_table)[5]="radius_parent"
            # Add the branch order to segment_table
            segment_table=data.table::data.table(merge(segment_table,unique(QSM_branch_angle[,c("axis_ID","branching_order","segment_ID")]),by="segment_ID"))

            if(method=="SegmentAngle"){
              if(level=="Tree"){return(mean(segment_table$Angle))}
              if(level=="branching_order"){return(segment_table[,.(branch_angle=mean(Angle)),by=branching_order])}
              if(level=="Branch"){
                Branch_angle_table=segment_table[,.(branch_angle=mean(Angle)),by=axis_ID]
                Branch_angle_table=data.table::setorder(Branch_angle_table,axis_ID)
                return(Branch_angle_table)
              }
            }
            if(method=="King98"){
              Branch_angle_table=data.table::data.table(plyr::ddply(QSM_branch_angle,("axis_ID"),function(x){ cbind(branch_angle=circular::deg(angle3d(as.numeric(x[1,c("startX","startY","startZ")]),as.numeric(c(x[1,c("startX","startY")],x[1,"startZ"]+1)),as.numeric(x[nrow(x),c("endX","endY","endZ")]))),branching_order=x$branching_order[1])}))
              Branch_angle_table=data.table::setorder(Branch_angle_table,axis_ID)
              if(level=="Branch"){return(Branch_angle_table[,1:2])}
              if(level=="Tree"){return(branch_angle=mean(Branch_angle_table$branch_angle))}
              if(level=="branching_order"){return(Branch_angle_table[,.(branch_angle=mean(branch_angle)),by=branching_order])}
            }

          }
)




