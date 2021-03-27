#' Add a point cloud to an object of class aRchi
#'
#' @docType methods
#' @rdname add_pointcloud
#' @description Add a point cloud to an object of class aRchi
#' @param aRchi The object of class aRchi
#' @param point_cloud The point cloud. Either a las or a data.frame
#' @include aRchiClass.R

setGeneric("add_pointcloud",
           function(aRchi,point_cloud){standardGeneric("add_pointcloud")}
)
#' @rdname add_pointcloud
#' @export


setMethod("add_pointcloud",
          signature = "aRchi",
          function(aRchi,point_cloud){
            ##########################
            # import the point cloud #
            ##########################
              if(class(point_cloud)[1] == "LAS"){
                aRchi@pointcloud = point_cloud
              }else{
                if(!is.data.frame(point_cloud)) stop("point_cloud must be a data.frame, data.table or LAS")
                if(ncol(point_cloud)<3) stop("point_cloud must have at least 3 columns")
                if(ncol(point_cloud)>3) warning("the point cloud contains more than three columns, three first used")

                # ensure the point cloud is a data.table with three columns and with the right columns names
                point_cloud = data.table::data.table(point_cloud[,1:3])
                data.table::setnames(point_cloud,c("X","Y","Z"))

                # build the LAS file and fill the aRchi slot
                aRchi@pointcloud = pkgcond::suppress_messages(lidR::LAS(point_cloud))
              }
            return(aRchi)
          }
)

