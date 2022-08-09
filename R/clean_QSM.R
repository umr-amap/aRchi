#' Cleans a QSM
#'
#' Cleans the QSM in an object of class \code{aRchi} by removing branches that have a disproportionate lower radius than their siblings.
#'
#' @export
#' @docType methods
#' @rdname Clean_QSM
#' @param aRchi an object of class \code{aRchi} with at least a QSM and a Paths table.
#' @param threshold numeric. The proportion of the largest daughter diameter (between 0 and 1) under which a branch is removed.
#' @param plotresult logical (default = FALSE). Show the results in a 3d plot if \code{TRUE}
#' @return An object of class \code{aRchi} with the cleaned QSM.
#' @details This cleaning is done by browsing the tree QSM from the base to the top. Each time a ramification point is encountered a daughter branch is removed if its radius is lower than a selected (i.e \code{threshold}) proportion of radius of the largest daughter. This allows removing small branches on large branches that can be for example traumatic or epicormic shoots or false branches due to noise in QSM. In \code{\link{ForkRate}} function the same approach is used with a \code{threshold} of 75\% (i.e 0.75) to count the number of fork and compute the fork rate.
#' @seealso \code{\link{ForkRate}} to compute the fork rate; \code{\link{Truncate_QSM}} to truncate a QSM at a specific diameter threshold
#' @include aRchiClass.R
#' @examples
#' \donttest{
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Clean the QSM: threshold of 0.5
#' Cleaned_Tree1_aRchi=Clean_QSM(Tree1_aRchi,threshold = 0.5,plotresult = TRUE)
#' # show the cleaned QSM data.table
#' get_QSM(Cleaned_Tree1_aRchi)
#'}
setGeneric("Clean_QSM",
           function(aRchi,threshold=NULL,plotresult=FALSE){standardGeneric("Clean_QSM")}
)

#' @rdname Clean_QSM
#' @export

setMethod("Clean_QSM",
          signature = "aRchi",
          function(aRchi,threshold,plotresult){
            radius=percent_diam=segment_ID=axis_ID=V1=node_ID=startX=startY=startZ=endX=endY=endZ=NULL

            if(is.null(threshold)) stop("Please provide a threshold")
            if(inherits(aRchi,"aRchi")==F) stop("The provided data is not of class aRchi")
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
                node_ID=x$node_ID[1])},.progress = "text")
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
            Perenial_structure_table=data.table::data.table(dplyr::anti_join(Perenial_structure_table,unique(Perenial_structure_table[percent_diam<threshold,"ID_Path"]),by="ID_Path"))
            Perenial_structure_table=unique(Perenial_structure_table[,-"ID_Path"])

            cleanedQSM=QSM[segment_ID%in%Perenial_structure_table$segment_ID]


            # Rectification of the segment and node ID following the suppression of certain subtree
            # How many daughters per node ?
            nb_br_node=cleanedQSM[,(length(unique(segment_ID))),by="node_ID"]

            # Node with only one daughters are not node anymore. We need to attribute the the segment ID and the node ID of these daughters who lost her sister. Except for the first node (node_ID = 0) which has no mother neither daughter..

            one_br_node=nb_br_node[V1==1][-1,]
            for (i in one_br_node$node_ID) {

              new_segm_ID= unique(cleanedQSM[node_ID==i]$segment_ID)
              cleanedQSM[node_ID==i]$axis_ID=unique(cleanedQSM[segment_ID==i]$axis_ID)
              cleanedQSM[node_ID==i]$branching_order=unique(cleanedQSM[segment_ID==i]$branching_order)
              cleanedQSM[node_ID==i]$node_ID=unique(cleanedQSM[segment_ID==i]$node_ID)
              cleanedQSM[segment_ID==i]$segment_ID=new_segm_ID



            }
            # Remove sub_treebiomass Mf and Mf_r from QSM because they are obsolete with the new QSM
            if(is.null(cleanedQSM$Mf)==FALSE){cleanedQSM=cleanedQSM[,-c("sub_tree_biomass","Mf","Mf_r")]}

            aRchi@QSM=cleanedQSM
            aRchi=Make_Path(aRchi)
            message("\nPaths table has been re-estimated according to the new cleaned QSM")
            if(is.null(aRchi@Nodes)==FALSE){
              Make_Node(aRchi)

              message("\nNodes table has been re-estimated according to the new cleaned QSM")
            }


            if(plotresult){
              pc=QSM[startX==min(startX)|startX==max(endX)|startY==min(startY)|startY==max(endY)|startZ==min(startZ)|startZ==max(endZ),1:3]
              names(pc)=c("X","Y","Z")
              pc = pkgcond::suppress_messages( lidR::LAS(pc)) # pkgcond::supress_messages removes messages from the LAS building
              lidR::plot(pc,bg="black",size=0,clear_artifacts=FALSE,axis=T)

              ls_cyl=plyr::alply(cleanedQSM,1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
              rgl::shapelist3d(ls_cyl,color=2,alpha=1,add=TRUE,lit=TRUE) # plot the list

              rest_of_QSM=QSM[!segment_ID%in%Perenial_structure_table$segment_ID]
              if(nrow(rest_of_QSM)==0){
                message("\nNo branches to remove")
                return(aRchi)
              }
              ls_cyl=plyr::alply(rest_of_QSM,1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
              rgl::shapelist3d(ls_cyl,color="white",alpha=1,add=TRUE,lit=TRUE) # plot the list

            }
            aRchi@operations$Clean_QSM = c(threshold = threshold)
            return(aRchi)
          }
)


