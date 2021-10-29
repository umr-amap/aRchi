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
            # Transparent backward compatibility with lidR 3 and lidR 4
            if (is(aRchi@pointcloud, "LAS"))
            {
              if (!.hasSlot(aRchi@pointcloud, "crs"))
              {
                pts <- aRchi@pointcloud@data
                phb <- aRchi@pointcloud@header@PHB
                crs <- aRchi@pointcloud@proj4string
                return(suppressWarnings(lidR::LAS(pts, phb, proj4string = crs, check = FALSE)))
              }
            }

            return(aRchi@pointcloud)
          }
)

