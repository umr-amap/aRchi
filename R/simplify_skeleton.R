
#' Simplify a skeleton by removing unnecessary cylinders
#'
#' @param aRchi an object of class aRchi containing at least a QSM.
#' @param seg_length numeric. The target maximal cylinder length.
#'
#' @return the simplified QSM.
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
#' # simplyfy the skeleton
#' aRchi = aRchi::simplify_skeleton(aRchi,seg_length = 0.05)
#'
#' # plot the simplifyed skeleton
#' plot(aRchi,show_point_cloud = TRUE)
#' }

setGeneric("simplify_skeleton",
           function(aRchi,seg_length = 0.1){standardGeneric("simplify_skeleton")}
)

#' @rdname simplify_skeleton
#' @export

setMethod("simplify_skeleton",
          signature = "aRchi",
          function(aRchi,seg_length = 0.1){

            # to pass CRAN check
            .=axis_ID=cyl_ID=endX=endY=endZ=parent_ID=startX=startY=startZ=NULL

            skel = aRchi@QSM

            ###########################
            #- simplify the skeleton -#
            ###########################

            # create the output data depending on if a radius exists or not
            if(!is.null(skel$radius)){
              skel_new = data.table::data.table(matrix(ncol=12))
              data.table::setnames(skel_new,c("startX","startY","startZ","endX","endY","endZ","cyl_ID","parent_ID","length","axis_ID","branching_order","radius"))
              rad_store = c()
              rad = TRUE
            }else{
              skel_new = data.table::data.table(matrix(ncol=11))
              data.table::setnames(skel_new,c("startX","startY","startZ","endX","endY","endZ","cyl_ID","parent_ID","length","axis_ID","branching_order"))
              rad = FALSE
            }

            # the axes index
            axes = sort(unique(skel$axis_ID))

            # progress bar
            pb <- progress::progress_bar$new(total = nrow(skel),width = 60,
                                             format = " Simplifying skeleton [:bar] :percent",clear=FALSE)
            for(i in axes){
              # for each axis go from the base to the tip
              axis = skel[axis_ID == i]

              # the curent cumulasted cylinder length
              cur_L = 0

              # the coordinates of the cylinder starting point
              start = axis[1,.(startX,startY,startZ)]
              for(j in 1:nrow(axis)){
                pb$tick() # progress bar

                # update the curent length
                cur_L = cur_L + axis$length[j]

                # if the radius exist, store it to attibute it to the new cylinder
                if(rad) rad_store = c(rad_store,axis$radius[j])

                # compute the number of child segments
                n_child = nrow(skel[parent_ID == axis$cyl_ID[j]])


                if(cur_L >= seg_length | n_child == 0 | n_child > 1){
                  # if the curent cumulated length is target length or that there is a child axis
                  # store the coordinates of the segment tip
                  end = axis[j,.(endX,endY,endZ)]

                  # add the segments to the new segment table
                  if(rad){
                    skel_new = data.table::rbindlist(list(skel_new,data.table::data.table(start[,1:3],end[,1:3],NA,NA,NA,i,axis$branching_order[j],mean(rad_store))),use.names = F)
                    rad_store = c()
                  }else{
                    skel_new = data.table::rbindlist(list(skel_new,data.table::data.table(start[,1:3],end[,1:3],NA,NA,NA,i,axis$branching_order[j])),use.names = F)
                  }

                  # define the starting point of the next segment
                  start = end

                  # after a new segment was created, the curent length is back to 0
                  cur_L = 0
                }
              }
            }

            # remove the first row that only contains NAs
            skel_new = skel_new[-1]

            ######################
            #- compute topology -#
            ######################

            # create segment ID
            skel_new[,cyl_ID := 1:nrow(skel_new)]

            # assign bearer ID
            skel_new[,parent_ID := 0]
            skel_new[,parent_ID := FNN::knnx.index(data = skel_new[,4:6],query = skel_new[,1:3], algorithm = "kd_tree", k = 1)]
            skel_new[cyl_ID == parent_ID,parent_ID := 0]

            # segments length
            skel_new[,length := sqrt( (startX - endX)^2 + (startY - endY)^2 + (startZ - endZ)^2)]

            aRchi@QSM = skel_new

            aRchi@operations$simplify_skeleton = list(seg_length = seg_length)

            return(aRchi)
          }
)


