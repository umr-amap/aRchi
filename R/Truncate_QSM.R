#' Truncate a QSM
#'
#' Truncate a QSM at a radius threshold
#'
#' @export
#' @docType methods
#' @rdname Truncate_QSM
#' @param aRchi an object of class aRchi with at least a QSM and a Paths table.
#' @param threshold numeric. The radius threshold in meter.
#' @param Keepdaughters logical (default = \code{FALSE}). Keep the daughters of the last segment retained even if they are lower than the threshold.
#' @param plotresult logical (default = \code{FALSE}). Show the results in a 3d plot if \code{TRUE}
#' @return An aRchi file with the QSM truncated
#' @details The threshold is applied to a whole segments. In other word, if a segment has at least one cylinder lower than the threshold it is removed as well as everything upstream (except the direct daughters if \code{Keepdaughters=TRUE}).
#' @seealso \code{\link{Clean_QSM}} to clean a QSM of an object aRchi.
#' @examples
#' \donttest{
#' # Read an aRchifile with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Truncate the QSM: 5cm radius threshold
#' Truncated_Tree1_aRchi=Truncate_QSM(Tree1_aRchi,plotresult = TRUE,threshold = 0.05)
#' }
#' @include aRchiClass.R
setGeneric("Truncate_QSM",
           function(aRchi,threshold=NULL,Keepdaughters=FALSE,plotresult=FALSE){standardGeneric("Truncate_QSM")}
)

#' @rdname Truncate_QSM
#' @export

setMethod("Truncate_QSM",
          signature = "aRchi",
          function(aRchi,threshold,Keepdaughters,plotresult){
            ID_Path=radius_cyl=segment_ID=cyl_ID=node_ID=axis_ID=V1=startX=startY=startZ=endX=endY=endZ=NULL

            if(is.null(threshold)) stop("Please provide a threshold")
            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(aRchi@Paths)) stop("The archi file does not contains Paths")


            QSM=aRchi@QSM
            Paths=aRchi@Paths

            Paths_subSample_tree=NULL

            pb <- progress::progress_bar$new(total =length(unique(Paths$ID_Path)),width = 60,
                                             format = " Browsing paths [:bar] :percent",clear=FALSE)
            for (j in unique(Paths$ID_Path)) {

              # Path j
              pb$tick()
              Path_j=Paths[ID_Path==j]

              # Find the list of cylinders ID with a radius lower to the threshold
              cyl_ID_threshold=Path_j[radius_cyl<(threshold)]$cyl_ID
              if(length(cyl_ID_threshold)==0) {
                sub_path_j=Path_j
                Paths_subSample_tree=rbind(Paths_subSample_tree,sub_path_j)
                next
              }

              # Take the more upstream cylinder
              cyl_ID_min=min(cyl_ID_threshold)

              # Find the first segment_ID with a radius lower than the threshold in path j
              Internode_Diam_lim=Path_j[segment_ID==Path_j[cyl_ID==cyl_ID_min]$segment_ID[1]]

              # Select the segments having more upstream than Internode_Diam_lim.
              sub_path_j= Path_j[segment_ID<unique(Internode_Diam_lim$segment_ID)]


              # Find the daugthers of the kept segment in the path j. If keepdaughters= TRUE
              if(Keepdaughters){
                sub_path_j_Daughters=NULL
                for (k in unique(sub_path_j$segment_ID)){

                  sub_path_j_k=sub_path_j[segment_ID==k]
                  Daugthers_internode=QSM[node_ID ==tail(QSM[segment_ID==unique(sub_path_j_k$segment_ID)]$cyl_ID,n=1)]
                  Daugthers_internode=cbind(Daugthers_internode,ID_Path=rep(j,nrow(Daugthers_internode)))

                  sub_path_j_k=rbind(Daugthers_internode,sub_path_j_k)

                  sub_path_j_Daughters=rbind(sub_path_j_Daughters, sub_path_j_k)

                }
                sub_path_j=rbind(sub_path_j_Daughters,sub_path_j)
              }

              # rbind the retained path.
              Paths_subSample_tree=rbind(Paths_subSample_tree,sub_path_j)

            }
            # Take the retained cyl_ID
            cyl_ID_subSample_tree=unique(Paths_subSample_tree$cyl_ID)

            # Find them in the QSM to have the truncated QSM !
            TruncatedQSM=QSM[cyl_ID%in%sort(cyl_ID_subSample_tree)]

            # Find the cylinder ID of the new apex cylinder to attribute them an extension ID of 0. This is needed for path recalculation.
            Cyl_ID_term=unique(plyr::ddply(Paths_subSample_tree,("ID_Path"),function(x){x[which(x$cyl_ID==max(x$cyl_ID)),]})$cyl_ID)
            Cyl_ID_term=unique(sapply(Cyl_ID_term,function(x){max(Paths_subSample_tree[ID_Path%in%Paths_subSample_tree[cyl_ID==x]$ID_Path]$cyl_ID)}))
            TruncatedQSM[cyl_ID%in%Cyl_ID_term]$extension_ID=0

            # Rectification of the segment and node ID following the suppression of certain subtree
            # How many daughters per node ?
            nb_br_node=TruncatedQSM[,(length(unique(segment_ID))),by="node_ID"]

            # Node with only one daughters are not node anymore. We need to attribute the segment ID and the node ID of these daughters who lost her sister. Except for the first node (node_ID = 0) which has no mother neither daughter..

            one_br_node=nb_br_node[V1==1][-1,]
            for (i in one_br_node$node_ID) {

              new_segm_ID= unique(TruncatedQSM[node_ID==i]$segment_ID)
              TruncatedQSM[node_ID==i]$axis_ID=unique(TruncatedQSM[segment_ID==i]$axis_ID)
              TruncatedQSM[node_ID==i]$branching_order=unique(TruncatedQSM[segment_ID==i]$branching_order)
              TruncatedQSM[node_ID==i]$node_ID=unique(TruncatedQSM[segment_ID==i]$node_ID)
              TruncatedQSM[segment_ID==i]$segment_ID=new_segm_ID
            }
            if(is.null(aRchi@QSM$Mf)==FALSE){TruncatedQSM=TruncatedQSM[,-c("sub_tree_biomass","Mf","Mf_r")]}

            aRchi@QSM=TruncatedQSM
            if(length(unique(aRchi@QSM$branching_order)>1)){
            aRchi=Make_Path(aRchi)
            message("\nPaths table has been re-estimated according to the new truncated QSM")}
            if(length(unique(aRchi@QSM$branching_order)==1)){aRchi@Paths=NULL}
            if(is.null(aRchi@Nodes)==FALSE){
              aRchi=Make_Node(aRchi)

              message("\nNodes table has been re-estimated according to the new truncated QSM")
            }

            if(plotresult){

              pc=QSM[startX==min(startX)|startX==max(endX)|startY==min(startY)|startY==max(endY)|startZ==min(startZ)|startZ==max(endZ),1:3]
              names(pc)=c("X","Y","Z")
              pc = pkgcond::suppress_messages( lidR::LAS(pc)) # pkgcond::supress_messages removes messages from the LAS building

              lidR::plot(pc,bg="black",size=0,clear_artifacts=FALSE,axis=T)

              ls_cyl=plyr::alply(TruncatedQSM,1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
              rgl::shapelist3d(ls_cyl,color=2,alpha=1,add=TRUE,lit=TRUE) # plot the list

              rest_of_QSM=QSM[!cyl_ID%in%cyl_ID_subSample_tree]
              ls_cyl=plyr::alply(rest_of_QSM,1,function(x){rgl::cylinder3d(rbind(as.matrix(x[,c("startX","startY","startZ")]),as.matrix(x[,c("endX","endY","endZ")])),radius= x[,"radius_cyl"][[1]],sides=8,closed=-2)}) # a list of cylinder
              rgl::shapelist3d(ls_cyl,color="white",alpha=1,add=TRUE,lit=TRUE) # plot the list

            }
            aRchi@operations$Truncate_QSM=c(threshold=threshold,Keepdaughters=Keepdaughters)

            return(aRchi)

          }
)


