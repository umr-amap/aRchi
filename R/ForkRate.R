#' Compute the fork rate of a tree
#'
#' Compute the fork rate from an aRchi object
#'
#' @export
#' @docType methods
#' @rdname ForkRate
#' @param aRchi a file of class aRchi with at least a QSM and a Path table.
#' @return a vector with two numeric value. The number of Fork and the fork rate.
#' @details
#'
#' The fork rate is the mean number of forks per meter of tree height. This metric is computed by browsing tree QSM from the base to the top. Each time a ramification point is encountered it is evaluated as a fork if at least one daughter had a radius not less than 75% of the diameter of the largest daughter. This threshold was chosen in order to exclude non-perennial structures, such as traumatic or epicormic shoots, or branches that will not last on the tree. If the ramification point is a fork, all retained daughter branches are browsed through until the next ramification point and further until the path end. If a daughter is rejected, it is removed as well as all the paths passing through it.
#'
#'@references
#'
#' 	Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' @include aRchiClass.R
#' @seealso [Clean_QSM())] to clean a QSM based on a threshold of percentage of largest daughter's diameter.
#' @examples
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Compute the fork rate of Tree1
#' ForkRate(Tree1_aRchi)
setGeneric("ForkRate",
           function(aRchi){standardGeneric("ForkRate")}
)

setMethod("ForkRate",
          signature = "aRchi",
          function(aRchi){
            radius=percent_diam=segment_ID=axis_ID=V1=NULL

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
                z_max=segment[nrow(segment),"endZ"],
                node_ID=x$node_ID[1])})
            names(segment_table)[5]="radius_parent"

            # Search for the maximum radius between the daughters of a same node
            node_max_radius_table=plyr::ddply(segment_table,("node_ID"),dplyr::summarise,node_max_radius=max(radius))
            # Add it to the segment table
            segment_table=data.table::data.table(merge(segment_table,node_max_radius_table,by="node_ID"))
            # Compute daughter radius/ max(daughter radius)
            segment_table$percent_diam=segment_table$radius/segment_table$node_max_radius
            # Retrieve Paths table with all the variable of Segment table
            Perenial_structure_table=merge(segment_table,unique(Paths[,c("segment_ID","ID_Path")]),by="segment_ID")
            # Keep only paths for which  daughter radius/ max(daughter radius) > to the threshold
            Perenial_structure_table=data.table::data.table(dplyr::anti_join(Perenial_structure_table,unique(Perenial_structure_table[percent_diam<0.75,"ID_Path"]),by="ID_Path"))
            Perenial_structure_table=unique(Perenial_structure_table[,-"ID_Path"])


            cleanedQSM=QSM[segment_ID%in%Perenial_structure_table$segment_ID]


            # Rectification of the segment and node ID following the suppression of certain subtree
            # How many daughters per node ?
            nb_br_node=cleanedQSM[,(length(unique(segment_ID))),by="node_ID"]
            # How many fork in the QSM
            N_Fork=nrow(nb_br_node[V1>1])
            Forkrate=N_Fork/(max(QSM$startZ)-min(QSM$endZ))
            return(c(N_Fork=N_Fork,Forkrate=Forkrate))
          }
)


