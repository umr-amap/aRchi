#' Smooth a tree skeleton
#'
#' @param aRchi a file of class aRchi containing at least a skeleton
#' @param niter integer. number of iterations to perform
#' @param th numeric. The distance threshold to correct the segments tips position
#'
#' @return a file of class aRchi with a smoothed skeleton
#' @export
#'
#' @examples
#' # import aRchi file
#' aRchi=system.file("extdata","Tree_2.aRchi",package = "aRchi")
#' aRchi = aRchi::read_aRchi(aRchi)
#'
#' # plot unsmoothed skeleton
#' plot(aRchi)
#'
#' # smooth skeleton
#' aRchi = aRchi::smooth_skeleton(aRchi)
#'
#' # plot smoothed skeleton
#' plot(aRchi)


setGeneric("smooth_skeleton",
           function(aRchi,niter = 1,th = 0){standardGeneric("smooth_skeleton")}
)

#' @rdname smooth_skeleton
#' @export

setMethod("smooth_skeleton",
          signature = "aRchi",
          function(aRchi,niter = 1,th = 0){

            # to pass CRAN check
            startX=startY=startZ=endX=endY=endZ=parent_ID=.=NULL

            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")

            data = aRchi@QSM

            if(!is.numeric(th)) stop("th must be numeric")
            if(!is.numeric(niter)) stop("niter must be numeric")

            for(i in 1:niter){
              #- progress
              svMisc::progress(i,max.value = niter)

              # for each axis
              for(ax in unique(data$axis_ID)){
                # rows corresponding to the target axis
                axis = which(data$axis_ID == ax)
                # correction only applies if the axis is composed of at least 2 segments
                if(length(axis)>=2){
                  for(j in 2:length(axis)){
                    # computes the distance of b relative to the line ac
                    d = aRchi::dist2line(b=data[axis[j-1],c(startX,startY,startZ)],
                                         c=data[axis[j-1],c(endX,endY,endZ)],
                                         a=data[axis[j],c(endX,endY,endZ)])
                    # if the distance is greater than the threshold -> apply the correction
                    if(d > th){
                      # new coordinates are the mean between the extremities of the segments
                      new_coord = c( (data$startX[axis[j-1]]+data$endX[axis[j]])/2,
                                     (data$startY[axis[j-1]]+data$endY[axis[j]])/2,
                                     (data$startZ[axis[j-1]]+data$endZ[axis[j]])/2 )
                      # new coordinates assigned to the end of the previous segment and
                      # to the start of the curent segment and child segments of j-1
                      data[axis[j-1],':='(endX = new_coord[1],
                                          endY = new_coord[2],
                                          endZ = new_coord[3])]
                      data[c(axis[j],which(parent_ID == data$cyl_ID[axis[j-1]])),
                           ':='(startX = new_coord[1],
                                startY = new_coord[2],
                                startZ = new_coord[3])]
                    }
                  }
                }
              }
            }

            data[,length := sqrt( (startX - endX)^2 +
                                    (startY - endY)^2 +
                                    (startZ - endZ)^2)]

            aRchi@QSM = data

            aRchi@operations$"Smoothing" = c(niter = niter,th = th)

            return(aRchi)
          }
)






