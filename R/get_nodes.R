#' Get the nodes from an object of class aRchi
#'
#' @docType methods
#' @rdname get_nodes
#' @description Get the nodes from an object of class aRchi
#' @param aRchi The object of class aRchi
#' @include aRchiClass.R
#' @seealso \code{\link{get_QSM}}; \code{\link{get_pointcloud}}; \code{\link{get_paths}}
#' @examples
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' # get the nodes (a list of two data.table)
#' get_nodes(Tree1_aRchi)
#'
setGeneric("get_nodes",
           function(aRchi){standardGeneric("get_nodes")}
)

#' @rdname get_nodes
#' @export
setMethod("get_nodes",
          signature = "aRchi",
          function(aRchi){
            return(aRchi@Nodes)
          }
)

