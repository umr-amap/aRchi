#' Add the physiological age of an annual shoots based on their length
#'
#' @param aRchi an object of class aRchi containing at least a QSM.
#' @param th numeric. The length thresholds used to segment annual shoots into physiological ages.
#' @param correct_PA logical. Should a correction of the physiological age be performed ?
#' This correction is based on the assumption that the child annual shoot can not be of lower
#' order than its parent annual shoot. The parent annual shoot physiological age is modified accordingly.
#'
#' @return The aRchi file with a physiological age field added to the QSM slot
#'
#' @references Lecigne, B., Delagrange, S., & Taugourdeau, O. (2021). Annual Shoot Segmentation and
#' Physiological Age Classification from TLS Data in Trees with Acrotonic Growth. Forests, 12(4), 391.
#' https://doi.org/10.3390/f12040391
#'
#' @export
#'
#' @examples
#' \donttest{
#' # import aRchi file
#' aRchi=system.file("extdata","Tree_2.aRchi",package = "aRchi")
#' aRchi = aRchi::read_aRchi(aRchi)
#'
#' # smooth skeleton
#' aRchi = smooth_skeleton(aRchi)
#'
#' # segment annual shoots
#' aRchi = aRchi::segment_annual_shoots(aRchi,tree_age = 13)
#'
#' # add physiological ages
#' aRchi = aRchi::add_physiological_ages(aRchi)
#'
#' plot(aRchi,color="physiological_age",bg = "white")
#' }

setGeneric("add_physiological_ages",
           function(aRchi,th = c(0.1,0.2,0.5), correct_PA = TRUE){standardGeneric("add_physiological_ages")}
)
#' @rdname add_physiological_ages
#' @export

setMethod("add_physiological_ages",
          signature = "aRchi",
          function(aRchi,th = c(0.1,0.2,0.5), correct_PA = TRUE){

            # to pass CRAN check
            .=AS=PA=axis_ID=endX=endY=endZ=physiological_age=startX=startY=startZ=NULL

            if(!is.numeric(th)) stop("th must be numeric")

            if(inherits(aRchi,"aRchi")==F) stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("This aRchi file does not contains a QSM")
            if(is.null(aRchi@QSM$AS)) stop("Annual shoots were not segmented in this QSM.")

            # the skeleton
            skel = aRchi@QSM

            # cylinder length
            skel[,length := sqrt( (startX-endX)^2 + (startY-endY)^2 + (startZ-endZ)^2 )]

            # length of each annual shoot
            AS_length = skel[,sum(length),by = .(axis_ID,AS)]
            data.table::setnames(AS_length,c("axis_ID","AS","length"))

            ## compute physiological based on annual shoot length
            # add 0 to the thresholds
            th = c(0,sort(th))
            for(i in 1:(length(th))){
              if(i == (length(th))){
                AS_length[length >= th[i],physiological_age := i]
              }else{
                AS_length[length >= th[i] & length < th[i+1],physiological_age := i]
              }
            }

            # the longest AS are the smallest physiological age
            AS_length[,physiological_age := (length(th)+1)-physiological_age ]

            # add physiological ages to the skeleton
            data.table::setkeyv(AS_length,c("axis_ID","AS"))
            data.table::setkeyv(skel,c("axis_ID","AS"))

            skel[AS_length,PA := physiological_age]

            if(correct_PA){
              for(i in sort(unique(skel$axis_ID),decreasing = T)){

                first = TRUE
                for(j in sort(unique(skel$AS[skel$axis_ID == i]))){
                  if(!first){
                    childs_PA = skel$PA[which(skel$parent_ID %in% skel$cyl_ID[which(skel$axis_ID == i & skel$AS == j)] & skel$AS != j)]
                    if(length(childs_PA) >= 1){
                      if(min(childs_PA) < AS_length[axis_ID == i & AS == j, unique(physiological_age)]){
                        skel[axis_ID == i & AS == j, PA := min(childs_PA)]
                      }
                    }
                  }
                  first = FALSE
                }
              }
            }

            aRchi@QSM = skel

            aRchi@operations$add_physiological_ages = list(th = th, correct_PA = correct_PA)

            return(aRchi)
          }
)
