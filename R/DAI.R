#' Estimate the index of dominance of the principal axis: DAI.
#'
#' Estimate the index of dominance of the principal axis (DAI) from an aRchi object.
#'
#' @export
#' @docType methods
#' @rdname DAI
#' @param aRchi an object of class aRchi with at least the QSM and the Paths table
#' @return Numeric. The value of DAI.
#' @details
#'
#' The idea of DAI is to disentangle architectures of trees with a strong apical dominance, i.e. with a central main stem growing more strongly than other side axes (such as in most conifers, for instance), from those having a spread out branching pattern with similar axes and no obvious main stem.
#'
#' The higher the index the more dominant the principal axis. DAI is computed based on the indices A0 which is an index of the probability of being the principal axis for each path of the tree (see function [Compute_A0()]).
#'
#' DAI is thus computed using the following formula:
#' DAI = equation with n_p the number of path, A0_i the A0 value of the path i and maxA0 the maximum A0 value (i.e the principal axis A0 value).
#'
#' @references
#'
#' Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' @seealso [Compute_A0())] to identify the principal axis (i.e max(A0)).
#'
#' @include aRchiClass.R
#'
#' @examples
#' # Read an aRchi file with at least the QSM and the paths table
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' DAI(Tree1_aRchi)
#'
setGeneric("DAI",
           function(aRchi){standardGeneric("DAI")}
)

setMethod("DAI",
          signature = "aRchi",
          function(aRchi){
            H_segment_rel=segment_ID=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(aRchi@Paths)) stop("The archi file does not contains Paths")

            QSM=aRchi@QSM
            Paths=aRchi@Paths
            # Make a table with information at the segment level (i.e length, radius of the first cylinder, radius of the last cylinder of the parent, length of the parent, angle and angle of the parent). This table is needed to identify A1
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
            # A1 index (added) per path
            Index_A1_add=data.table::data.table(plyr::ddply(Identification_trunk_table,("ID_Path"),function(x){mean(x$AS_angle+x$AS_diam+x$vert_angle+x$H_path_rel,na.rm=T)}))
            # A1 index (multiplied) per path
            Index_A1_mult=data.table::data.table(plyr::ddply(Identification_trunk_table,("ID_Path"),function(x){mean(x$AS_angle*x$AS_diam*x$vert_angle*x$H_path_rel,na.rm=T)}))

            # Index of dominacne of the principal axis.
            DA=max(Index_A1_add$V1)-mean(Index_A1_add$V1)
            return(DA)

          }
)

