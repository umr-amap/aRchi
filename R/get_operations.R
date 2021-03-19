#' Get the nodes from an object of class aRchi
#'
#' @docType methods
#' @rdname get_operations
#' @description Show a list with all the operations (and their parameters) that have been performed on an object of class aRchi
#' @param aRchi The object of class aRchi
#' @include aRchiClass.R
#' @seealso \code{\link{get_QSM}}; \code{\link{get_pointcloud}}; \code{\link{get_paths}}
#' @examples
#' \dontrun{
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#'# Making some operations
#' Tree1_aRchi<-Make_Path(Tree1_aRchi)
#' Tree1_aRchi=Compute_A0(Tree1_aRchi)
#' Tree1_aRchi=Compute_Mf(Tree1_aRchi,WoodDensity = 550)
#' Tree1_aRchi=Clean_QSM(Tree1_aRchi,threshold = 0.5)
#'
#' # Show the oprations and their parameters
#' get_operations(Tree1_aRchi)
#'}
setGeneric("get_operations",
           function(aRchi){standardGeneric("get_operations")}
)
#' @rdname get_operations
#' @export
setMethod("get_operations",
          signature = "aRchi",
          function(aRchi){
            return(aRchi@operations)
          }
)

