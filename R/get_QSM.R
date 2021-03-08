#' Get the QSM from an object of class aRchi
#'
#' @docType methods
#' @rdname get_QSM
#' @description Get the QSM from an object of class aRchi
#' @param aRchi The object of class aRchi
#' @include aRchiClass.R
#' @seealso [get_pointcloud()]; [get_paths()]; [get_nodes()];
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

setMethod("get_QSM",
          signature = "aRchi",
          function(aRchi){
            return(aRchi@QSM)
          }
)

