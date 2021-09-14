#' Add the radius to a skeleton based on point distance to the skeleton
#'
#' @param aRchi an object of class aRchi containing at least a point cloud and a QSM.
#' @param sec_length numeric. The length of the section to compute the radius. See details.
#' @param by_axis logical. If \code{TRUE} the radius is calculated for each section within each axis, if \code{FALSE}
#'                the radius is calculated within section for all axes simultaneously. Default = TRUE.
#' @param method character. The axis radius can be either the mean (\code{method = "mean"}) or the
#'               median (\code{method = "median"}) distance to the skeleton. Default = "median".
#'
#' @details The point distance to skeleton is likely to vary considerably over short
#'          distances (e.g. at a branching point). Therefore, the radius is computed
#'          by averaging the point distance to the skeleton over sections of user defined
#'          length.
#' @return The aRchi file with an updated radius in the QSM slot.
#' @export
#'
#' @examples
#' \donttest{
#' # import a point cloud
#' tls=system.file("extdata","Tree_2_point_cloud.las",package = "aRchi")
#' tls = lidR::readLAS(tls)
#'
#' aRchi = aRchi::build_aRchi()
#' aRchi = aRchi::add_pointcloud(aRchi,point_cloud = tls)
#'
#' # build a skeleton from the point cloud
#' aRchi = aRchi::skeletonize_pc(aRchi)
#'
#' # smooth the skeleton
#' aRchi = aRchi::smooth_skeleton(aRchi)
#'
#' # add the diameter to the skeleton
#' aRchi = aRchi::add_radius(aRchi)
#'
#' # plot the QSM
#' plot(aRchi,skeleton = FALSE,show_point_cloud = FALSE)
#' }

setGeneric("add_radius",
           function(aRchi,sec_length = 0.5,by_axis = TRUE,method = "median"){standardGeneric("add_radius")}
)

#' @rdname add_radius
#' @export

setMethod("add_radius",
          signature = "aRchi",
          function(aRchi,sec_length = 0.5,by_axis = TRUE,method = "median"){

            # to pass CRAN checks
            .=X=Y=Z=axis_ID=cyl_ID=dist=dist_tip=endX=endY=endZ=radius_cyl=radius=sec=seg=startX=startY=startZ=NULL

            if(!is.numeric(sec_length)) stop("sec_length must be numeric")

            pc = aRchi@pointcloud@data[,.(X,Y,Z)]
            skel = aRchi@QSM

            #- segment distance from the axis tip
            skel[,dist_tip := rev(cumsum(length)),by=axis_ID]

            #### to which segment belongs each point ?
            # point distance to the two nearest nodes
            nearest = FNN::knnx.index(data = skel[,4:6],query = pc[,1:3], algorithm = "kd_tree", k = 2)
            # the coordinates of the two nearest nodes
            near_coord = cbind(skel[nearest[,1],4:6],skel[nearest[,2],4:6])
            rm(nearest)
            data.table::setnames(near_coord,c("X1","Y1","Z1","X2","Y2","Z2"))
            near_coord[,index := 1:nrow(near_coord)]

            # as a node can be either the start or the end of the segment
            # first possibility
            data.table::setkeyv(skel,c("startX","startY","startZ","endX","endY","endZ"))
            data.table::setkeyv(near_coord,c("X1","Y1","Z1","X2","Y2","Z2"))
            near_coord[skel,seg := cyl_ID]
            # second possibility
            data.table::setkeyv(near_coord,c("X2","Y2","Z2","X1","Y1","Z1"))
            near_coord[skel,seg := cyl_ID]

            # add the seg ID to each point
            pc[near_coord$index,seg := near_coord$seg]
            rm(near_coord)

            pc[is.na(seg),seg := skel$cyl_ID[FNN::knnx.index(data = skel[,4:6],query = pc[is.na(seg),1:3], algorithm = "kd_tree", k = 1)]]


            # add distance from axis tip, axis_ID and cylinder start and end coordinates for each point
            pc[,':='(
              axis_ID = skel$axis_ID[skel$cyl_ID == seg],
              dist_tip = skel$dist_tip[skel$cyl_ID == seg],
              startX = skel$startX[skel$cyl_ID == seg],
              startY = skel$startY[skel$cyl_ID == seg],
              startZ = skel$startZ[skel$cyl_ID == seg],
              endX = skel$endX[skel$cyl_ID == seg],
              endY = skel$endY[skel$cyl_ID == seg],
              endZ = skel$endZ[skel$cyl_ID == seg]
            ),by=seg]


            # compute point distance to skeleton
            pc[,dist := dist2line(
              a=c(X,Y,Z),
              b=c(endX,endY,endZ),
              c=c(startX,startY,startZ)),
              by=row.names(pc)]

            # remove useless variable
            pc[,':='(startX = NULL,startY = NULL,startZ = NULL,
                     endX = NULL,endY = NULL,endZ = NULL)]


            pc[,sec:=round((dist_tip-(sec_length/2))/sec_length)*sec_length] # sections


            # segment diameter is the average distance of the closest points
            if(by_axis){
              if(method == "median") pc[,radius:= median(dist),by=.(sec,axis_ID)]
              if(method == "mean") pc[,radius:= mean(dist),by=.(sec,axis_ID)]
            }else{
              if(method == "median") pc[,radius:= median(dist),by=sec]
              if(method == "mean") pc[,radius:= mean(dist),by=sec]
            }

            data.table::setkey(skel,"cyl_ID")
            data.table::setkey(pc,"seg")

            skel[pc,radius_cyl := radius]

            skel[cyl_ID == 1, radius_cyl := skel$radius_cyl[skel$cyl_ID == 2]]

            ##### correct diameter
            skel[,sec:=round((dist_tip-(sec_length/2))/sec_length)*sec_length] # sections
            skel[,radius_cyl := mean(radius_cyl,na.rm=TRUE),by=.(sec)] # mean dist by section and axis

            # segment diameter is the average distance of the closest points
            if(by_axis){
              for(j in unique(skel$axis_ID)){
                sections = sort(unique(skel$sec[skel$axis_ID==j]))
                if(length(sections) > 2){
                  for(i in 2:(length(sections)-1)){
                    if(unique(skel$radius_cyl[skel$sec==sections[i] & skel$axis_ID==j]) < unique(skel$radius_cyl[skel$sec==sections[i-1]& skel$axis_ID==j])){
                      skel[sec == sections[i] & skel$axis_ID==j,radius_cyl := (unique(skel$radius_cyl[skel$sec==sections[i-1]& skel$axis_ID==j])+unique(skel$radius_cyl[skel$sec==sections[i+1]& skel$axis_ID==j]))/2]
                    }
                  }
                }
              }
            }else{
              sections = sort(unique(skel$sec))
              for(i in 2:(length(sections)-1)){
                if(unique(skel$radius_cyl[skel$sec==sections[i]]) < unique(skel$radius_cyl[skel$sec==sections[i-1]])){
                  skel[sec == sections[i],radius_cyl := (unique(skel$radius_cyl[skel$sec==sections[i-1]])+unique(skel$radius_cyl[skel$sec==sections[i+1]]))/2]
                }
              }
            }


            skel[,':='(dist_tip = NULL,sec = NULL,volume = length*pi*radius_cyl^2)]

            aRchi@QSM = skel

            aRchi@operations$add_radius = list(sec_length = sec_length)

            return(aRchi)
          }
)

