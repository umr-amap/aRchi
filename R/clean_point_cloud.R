#' Filter noise from a point cloud
#'
#' @description Uses the \code{\link[VoxR]{filter_noise}} function to filter noise from the point cloud.
#'
#' @param aRchi An aRchi object containing a point cloud
#' @param k numeric. The number of nearest neighbours to use.
#' @param sigma numeric. The multiplier of standard deviation to consider a point as noise.
#'
#' @return The aRchi file with a clean point cloud.
#' @export
#'
#' @examples
#' # import aRchi file
#' aRchi=system.file("extdata","Tree_2.aRchi",package = "aRchi")
#' aRchi = aRchi::read_aRchi(aRchi)
#'
#' # clean point cloud
#' aRchi = aRchi::clean_point_cloud(aRchi)
#'

setGeneric("clean_point_cloud",
           function(aRchi, k = 5, sigma = 1.5){standardGeneric("clean_point_cloud")}
)

#' @rdname clean_point_cloud
#' @export

setMethod("clean_point_cloud",
          signature = "aRchi",
          function(aRchi, k = 5, sigma = 1.5){

            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@pointcloud)) stop("This aRchi file does not contains a point cloud")

            aRchi@pointcloud@data = VoxR::filter_noise(aRchi@pointcloud@data,k=k,sigma=sigma,
                                                       store_noise=FALSE,message=FALSE)

            aRchi@operations$clean_point_cloud = list(k=k,sigma=sigma)

            return(aRchi)
          }
)
