#' Get the pointcloud from an object of class aRchi
#'
#' @docType methods
#' @rdname get_pointcloud
#' @description Get the pointcloud from an object of class aRchi
#' @param aRchi The object of class aRchi
#' @include aRchiClass.R
#' @seealso \code{\link{get_QSM}}; \code{\link{get_paths}}; \code{\link{get_nodes}}
#' @examples
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' # get the pointcloud
#' get_pointcloud(Tree1_aRchi)
#'
setGeneric("get_pointcloud",
           function(aRchi){standardGeneric("get_pointcloud")}
)

#' @rdname get_pointcloud
#' @export


setMethod("get_pointcloud",
          signature = "aRchi",
          function(aRchi){
            if (utils::packageVersion("lidR") < "4.0.0")
              return(aRchi@pointcloud)
            else
              return(suppressWarnings(lidR::LAS(aRchi@pointcloud@data,
                                                aRchi@pointcloud@header@PHB,
                                                proj4string = aRchi@pointcloud@proj4string,
                                                check = FALSE)))
          }
)

