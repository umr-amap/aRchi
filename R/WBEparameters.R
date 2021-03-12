#' Compute WBE parameters
#'
#'   Compute from an object of class aRchi the West, Brown and Enquist (WBE) scaling exponent (alpha and beta) and the subsequent estimated metabolic rate (theta) at node, branch, branch order or tree level.
#'
#' @export
#' @docType methods
#' @rdname WBEparameters
#' @param aRchi an object of class aRchi with at least a QSM and the Nodes table (see function \code{\link{Make_Node}})
#' @param level characters. At which level R_ratio has to be estimated. \code{Node} for node level, \code{Axis} for axis level, \code{branching_order} for branch order level and \code{Tree} for tree level (default).
#' @param position At which position from the node WBE parameters have to be estimated. Either a numeric or a character. Use a numeric multiple of ten to select the distance from the node in cm where R ratio has to be estimated (e.g 10 for 10cm from the node). Use the \% sign after a multiple of ten to select the distance from the node in percentage of the length of the parent and daughters segments (e.g 50\% for an estimation at mid-length of the segments). Note that 0 is accepted and correspond to the closest position from the node.
#' @return Data.table of the summary (median and mean) of WBE parameters at the selected level. The median should be used according to Bentley et al. 2013.
#' @include aRchiClass.R
#' @details
#' Details for WBE parameters calculation are given in the details part of function \code{\link{Make_Node}}.
#' @references
#'
#' 	Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' 	Lau, A. et al. Estimating architecture-based metabolic scaling exponents of tropical trees using terrestrial LiDAR and 3D modelling. Forest Ecology and Management 439, 132–145 (2019).
#'
#' 	Bentley, L. P. et al. An empirical assessment of tree branching networks and implications for plant allometric scaling models. Ecology Letters 16, 1069–1078 (2013).
#'

#' @seealso \code{\link{Make_Node}} for node metrics estimation; \code{\link{LeonardoRatio}} to estimates Leonardo Da Vinci's ratio at different level.
#' @include aRchiClass.R
#' @examples
#' # Read an aRchifile with a QSM and node tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # WBE parameters at the branching order level estimated at midlength of the segments
#' LeonardoRatio(Tree1_aRchi,level="Tree", position="50%")
#' # WBE parameters at the tree level estimated at 10 cm from the node
#' LeonardoRatio(Tree1_aRchi)
#'
setGeneric("WBEparameters",
           function(aRchi,level="Tree",position=10){standardGeneric("WBEparameters")}
)

#' @rdname WBEparameters
#' @export
setMethod("WBEparameters",
          signature = "aRchi",
          function(aRchi, level, position){
            node_ID=pos_parent=pos_daugthers=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(aRchi@Nodes)) stop("The archi file does not contains node tables. Please use the MakeNode function")
            test_position=as.numeric(stringr::str_extract_all(position, "\\d+"))
            if((test_position/10)%%1!=0) stop("position should be a multiple of 10 (e.g 10) or a multiple of 10 and a % sign (e.g 50%)")
            if(level %in% c("Tree","branching_order","Node","Axis")==F) stop("Incorrect level argument")

            if(is.numeric(position)){
              NodeTable=aRchi@Nodes$Absolute_positions
            }else{
              NodeTable=aRchi@Nodes$Relative_positions
            }
            position=as.numeric(stringr::str_extract_all(position, "\\d+"))/100

            NodeTable_pos=NodeTable[pos_parent==position&pos_daugthers==position]

            if(level=="Node"){
              return(NodeTable_pos[,c("node_ID","pos_parent","pos_daugthers","radius_parent","R_ratio" )])
            }

            if(level=="Axis"){
              summary_table=plyr::ddply(aRchi@QSM[node_ID%in%NodeTable_pos$node_ID],("axis_ID"),function(x){cbind(Median_alpha=median(NodeTable_pos[node_ID%in%unique(x$node_ID )]$alpha),
                                                                                                                    Mean_alpha=mean(NodeTable_pos[node_ID%in%unique(x$node_ID )]$alpha),
                                                                                                                    Median_beta=median(NodeTable_pos[node_ID%in%unique(x$node_ID )]$beta),
                                                                                                                    Mean_beta=mean(NodeTable_pos[node_ID%in%unique(x$node_ID )]$beta),
                                                                                                                    Median_theta=median(NodeTable_pos[node_ID%in%unique(x$node_ID )]$theta),
                                                                                                                    Mean_theta=mean(NodeTable_pos[node_ID%in%unique(x$node_ID)]$theta),
                                                                                                                    N_node=length(unique(x$node_ID)))

              })

              return(data.table::data.table(summary_table))
            }


            if(level=="branching_order"){
              summary_table=plyr::ddply(aRchi@QSM[node_ID%in%NodeTable_pos$node_ID],("branching_order"),function(x){cbind(Median_alpha=median(NodeTable_pos[node_ID%in%unique(x$node_ID )]$alpha),
                                                                                                                      Mean_alpha=mean(NodeTable_pos[node_ID%in%unique(x$node_ID )]$alpha),
                                                                                                                      Median_beta=median(NodeTable_pos[node_ID%in%unique(x$node_ID )]$beta),
                                                                                                                      Mean_beta=mean(NodeTable_pos[node_ID%in%unique(x$node_ID )]$beta),
                                                                                                                      Median_theta=median(NodeTable_pos[node_ID%in%unique(x$node_ID )]$theta),
                                                                                                                      Mean_theta=mean(NodeTable_pos[node_ID%in%unique(x$node_ID )]$theta),
                                                                                                                      N_node=length(unique(x$node_ID)))
              })
              return(data.table::data.table(summary_table))
            }

            if(level=="Tree"){

              summary_table=cbind(Median_alpha=median(NodeTable_pos$alpha),
                                  Mean_alpha=mean(NodeTable_pos$alpha),
                                  Median_beta=median(NodeTable_pos$beta),
                                  Mean_beta=mean(NodeTable_pos$beta),
                                  Median_theta=median(NodeTable_pos$theta),
                                  Mean_theta=mean(NodeTable_pos$theta),
                                  N_node=length(unique(NodeTable_pos$node_ID)))
              return(summary_table)
            }
          }
)



