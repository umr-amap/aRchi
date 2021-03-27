#' Compute Leonardo's ratio
#'
#' Compute from an object of class aRchi the Leonardo's ratio (i.e R_ratio) at node, axis, branch order or tree level.
#'
#' @export
#' @docType methods
#' @rdname LeonardoRatio
#' @param aRchi an object of class aRchi with at least a QSM and the Nodes table (see function \code{\link{Make_Node}})
#' @param level characters. At which level R_ratio has to be estimated. \code{Node} for node level, \code{Axis} for the axis level, \code{branching_orde} for branch order level and \code{Tree} for tree level (default).
#' @param position At which position from the node R_ratio had to be estimated. Either a numeric or a character. Use a numeric multiple of ten to select the distance from the node in cm where R ratio has to be estimated (e.g 10 for 10cm from the node). Use the \% sign after a multiple of ten to select the distance from the node in percentage of the length of the parent and daughters segments (e.g 50\% for an estimation at mid-length of the segments). Note that 0 is accepted and correspond to the closest position from the node.
#' @return Data.table of the summary of R_ratio for the selected level.
#' @details
#' Details for Leonardo Da Vinci's ratio calculation are given in the details part of function \code{\link{Make_Node}}.
#' @include aRchiClass.R
#' @seealso \code{\link{Make_Node}} for node metrics estimation; \code{\link{WBEparameters}} to estimates WBE parameters at different level;
#' @include aRchiClass.R
#' @examples
#' # Read an aRchifile with a QSM and node tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Leonardo'ratio at the branching order level estimated at midlength of the segments
#' LeonardoRatio(Tree1_aRchi,level="Tree", position="50%")
#' # Leonardo'ratio at the node level estimated at 30 cm from the node
#' LeonardoRatio(Tree1_aRchi,level="Node", position=30)
setGeneric("LeonardoRatio",
           function(aRchi,level="Tree", position=10){standardGeneric("LeonardoRatio")}
)

#' @rdname LeonardoRatio
#' @export

setMethod("LeonardoRatio",
          signature = "aRchi",
          function(aRchi, level, position){
            pos_parent=pos_daugthers=node_ID=NULL
            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(aRchi@Nodes)) stop("The archi file does not contains node tables. Please use the MakeNode function")
            test_position=as.numeric(stringr::str_extract_all(position, "\\d+"))
            if((test_position/10)%%1!=0) stop("position should be a multiple of 10 (e.g 10) or a multiple of 10 and a % sign (e.g 50%)")
            if(level %in% c("Tree","branching_order","Node","Axis")==FALSE) stop("Incorrect level argument")
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
              summary_table=plyr::ddply(aRchi@QSM[node_ID%in%NodeTable_pos$node_ID],("axis_ID"),function(x){cbind(TRUE(as.matrix(summary(NodeTable_pos[node_ID%in%unique(x$node_ID)]$R_ratio))),N_node=length(unique(x$node_ID)))})
              return(data.table::data.table(summary_table))
            }


            if(level=="branching_order"){
              summary_table=plyr::ddply(aRchi@QSM[node_ID%in%NodeTable_pos$node_ID],("branching_order"),function(x){cbind(TRUE(as.matrix(summary(NodeTable_pos[node_ID%in%unique(x$node_ID)]$R_ratio))),N_node=length(unique(x$node_ID)))})
              return(data.table::data.table(summary_table))
            }

            if(level=="Tree"){
              return(summary(NodeTable_pos$R_ratio))
            }
          }
)



