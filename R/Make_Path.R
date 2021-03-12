#' Make the path of a QSM in an aRchi object
#'
#' Identify and record the paths of a QSM in an object of class aRchi. The path are needed for several tree metrics estimation of the aRchi package.
#'
#' @export
#' @docType methods
#' @rdname Make_Path
#' @param aRchi a file of class aRchi
#' @return The aRchi object with the table of paths
#' @details
#'
#' A path is a continuous succession of cylinders from a terminal segment (i.e a branch tip) to the trunk base. Thus, there is as many path as terminal segments in a QSM.
#'
#' This function fill the slot Paths of an object of class aRchi with a data.table of the paths. This data.table contains the same variables as a classic QSM plus an ID_path column. This table is thus larger than the QSM table as each cylinder is repeated as many times as it appears in a path. For example, the first cylinder of the QSM (i.e the beginning of the trunk) is repeated N path times, with N path the number of path of the QSM.
#'
#' Many function of the aRchi packages request a path table (check see also).
#' @seealso \code{\link{BranchAngle}};\code{\link{Truncate_QSM}}; \code{\link{Clean_QSM}}; \code{\link{ForkRate}}; \code{\link{PathFraction}}
#'
#' @include aRchiClass.R
#' @examples
#' # Read a QSM file
#' file=system.file("extdata","Tree_1_TreeQSM.txt",package = "aRchi")
#' QSM=read_QSM(file,model="treeQSM")
#' # Build an object of class aRchi
#' Tree1_aRchi<-build_aRchi(QSM=QSM)
#' Tree1_aRchi
#' # Make the path table
#' Tree1_aRchi<-Make_Path(Tree1_aRchi)
#' Tree1_aRchi
#' PathFraction(Tree1_aRchi)


setGeneric("Make_Path",
           function(aRchi){standardGeneric("Make_Path")}
)

#' @rdname Make_Path
#' @export


setMethod("Make_Path",
          signature = "aRchi",

          function(aRchi){
            extension_ID=segment_ID=node_ID=cyl_ID=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            QSM=aRchi@QSM
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            # Terminal segments
            segment_term=QSM[extension_ID==0]$segment_ID
            if(length(unique(segment_term))==1) stop("The QSM contains only a trunk")


            # Paths. A data.table
            Path_complet=NULL

            #List of data.table. Each element of the list has a name corresponding to a segment and contains all what the segment bears.
            Mf_var=NULL
            pb <- progress::progress_bar$new(total =length(unique(segment_term)),width = 60,
                                             format = " Searching for path [:bar] :percent",clear=FALSE)
            for (j in unique(segment_term)) {
              pb$tick()

              # Terminal segment j
              internode_Mf_var=list(QSM[segment_ID==j])

              # We attribute the ID of the segment the the element of the list. Here it is the first segment. it has no child !
              names(internode_Mf_var)=j

              # Segment_ID of the bearer
              num_internode_mere=unique(QSM[segment_ID==j]$node_ID)

              # Segment_ID of the bearer again. We need it two times with different names
              num_internode_mere_tot=num_internode_mere

              # We start the construction of Mf_var with the terminal segment. The first element of the list
              Mf_var=c(Mf_var,internode_Mf_var)

              # When the bearer of the terminal segment is directly the trunk, we add the trunk to the segment and the path is completed -> go to next path (j)
              if(unique(QSM[segment_ID==num_internode_mere]$node_ID)==0){
                last_step_path=cbind(internode_Mf_var[[1]],ID_Path=rep(j,nrow(internode_Mf_var[[1]])))
                last_step_path=rbind(last_step_path,cbind(QSM[node_ID==0],ID_Path=sort(rep(unique(last_step_path$ID_Path),nrow(QSM[node_ID ==0])))))
                Path_complet= rbind(Path_complet,last_step_path)
                next()
              }

              # We continue the construction of Mf_var with the bearer segment, the second element of the list
              internode_Mf_var=list(rbind(Mf_var[[length(Mf_var)]],QSM[segment_ID==num_internode_mere]))

              # We assign the ID to the element of the list.
              names(internode_Mf_var)=num_internode_mere

              # Mf_var has its second element
              Mf_var=c(Mf_var,internode_Mf_var)

              # Same approach as above through a loop for all the following segment until the trunk
              test_internode_1="start"
              while (test_internode_1!=0) {
                seg_mere=unique(QSM[segment_ID==num_internode_mere]$node_ID)
                test_internode_1=QSM[cyl_ID ==seg_mere]$node_ID
                num_internode_mere=unique(QSM[cyl_ID==seg_mere]$segment_ID)
                num_internode_mere_tot=c(num_internode_mere_tot,num_internode_mere)
                internode_Mf_var=list(rbind(Mf_var[[length(Mf_var)]],QSM[segment_ID==num_internode_mere]))
                names(internode_Mf_var)=num_internode_mere
                Mf_var=c(Mf_var,internode_Mf_var)
              }

              # Add a column with the path name
              last_step_path=cbind(internode_Mf_var[[1]],ID_Path=rep(j,nrow(internode_Mf_var[[1]])))

              # Add the trunk
              last_step_path=rbind(last_step_path,cbind(QSM[node_ID==0],ID_Path=sort(rep(unique(last_step_path$ID_Path),nrow(QSM[node_ID ==0])))))

              # Add path j to Path_complet
              Path_complet= rbind(Path_complet,last_step_path)


            }
            aRchi@Paths=Path_complet
            aRchi@operations$Make_Path=c(Computed="YES")
            return(aRchi)
          }
)




