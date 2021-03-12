#' Get the QSM from an object of class aRchi
#'
#' @docType methods
#' @rdname get_QSM
#' @description Get the QSM from an object of class aRchi
#' @param aRchi The object of class aRchi
#' @include aRchiClass.R
#' @seealso \code{\link{get_pointcloud}}; \code{\link{get_paths}}; \code{\link{get_nodes}}
#' @examples
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' # show the QSM data.table
#' get_QSM(Tree1_aRchi)
#'
setGeneric("get_QSM",
           function(aRchi){standardGeneric("get_QSM")}
)

#' @rdname get_QSM
#' @export

setMethod("get_QSM",
          signature = "aRchi",
          function(aRchi){
            return(aRchi@QSM)
          }
)

