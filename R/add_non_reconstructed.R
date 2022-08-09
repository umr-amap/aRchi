
#' Detect and add the axes that were not reconstructed by the QSM method
#'
#' @param aRchi a file of class aRchi.
#' @param max_dist numeric. The maximum distance of a point to the skeleton to be considered in the computation.
#' @param sec_length numeric. The length of an axis section to compute the average distance used in the computation of the point
#'                   relative distance to skeleton.
#' @param method character. Defines the method to use to identify the non reconstructed portions of the point cloud, see details.
#' @param th numeric. The threshold to use to identify the non reconstructed portions of the point cloud. See details.
#' @param d_clust numeric. The distance to use for axes clustering from the identifyed non reconstructed portions of the point cloud.
#'
#' @details  This function detects the non reconstructed axes by analyzing the point cloud distance to the skeleton.
#'           In a first step, the distance of each point to the skeleton is computed. At this step, the
#'           points that are too far from the skeleton can be removed with the parameter \code{max_dist}.
#'           The relative distance of the points is then computed by dividing their respective distance by
#'           the average distance of all the points within axes sections of a given length (\code{sec_length}). The
#'           points that are far from the skeleton (in terms of relative distance) are then considered as
#'           being part of a non reconstructed axis (NR). To do so two parameters are available. First, the
#'           \code{method} parameter defines which method should be used to identify the NR points and \code{th}
#'           sets the threshold to use. If \code{method = "absolute"} then \code{th} is a cut-off distance
#'           for the relative distance (1 being the average distance). If \code{method = "statistical"} then
#'           \code{th} is a multiplier of the standard deviation so that the points further than th*sd (sd = the
#'           standard deviation of the distribution of relative distances) are considered as NR points. NR axes are then
#'           segmented by clustering the NR points based on distance. To do so, the \code{d_clust} parameter
#'           sets the minimal distance between two points to be considered as part of two different axes. The
#'           NR axes are then reconstructed as a single segment.
#'
#' @note The tree topology is fully recomputed. Therefore any existing topology will be lost. It is thus recommended
#'       to use this function at the beginning of the processing pipeline.
#'
#' @references Lecigne, B., Delagrange, S., Lauri, P. É., & Messier, C. (2022). Trimming influences tree light interception
#'             and space exploration: contrasted responses of two cultivars of Fraxinus pennsylvanica at various scales of
#'             their architecture. Trees, 1-17. https://doi.org/10.1007/s00468-022-02273-5
#'
#' @return the aRchi file with reconstructed axes added to the skeleton with a field "reconstructed".
#'
#' @export
#'
#' @examples
#' \donttest{
#' # import aRchi file
#' aRchi=system.file("extdata","Tree_2.aRchi",package = "aRchi")
#' aRchi = aRchi::read_aRchi(aRchi)
#'
#' # smooth the skeleton
#' aRchi = aRchi::smooth_skeleton(aRchi,th = 0.01)
#'
#' # clean point cloud
#' aRchi = aRchi::clean_point_cloud(aRchi)
#'
#' # add missing axes
#' aRchi = aRchi::add_non_reconstructed(aRchi,th = 3)
#'
#' # plot the skeleton with reconstructed axes in red
#' plot(aRchi,color="reconstructed",show_point_cloud = TRUE)
#' }

setGeneric("add_non_reconstructed",
           function(aRchi,max_dist = 99999,sec_length = 0.2,method = "statistical",th = 2.5,d_clust = 0.02){standardGeneric("add_non_reconstructed")}
)
#' @rdname add_non_reconstructed
#' @export

setMethod("add_non_reconstructed",
          signature = "aRchi",

          function(aRchi,max_dist = 99999,sec_length = 0.2,method = "statistical",th = 2.5,d_clust = 0.02){

            # to pass CRAN checks
            X=Y=Z=startX=startY=startZ=endX=endY=endZ=axis_ID=seg=cyl_ID=dist=sec=mean_dist=.=rel_dist=cluster=tip=parent_ID=
              branching_order=dist_tip=section=bearer_ID=seg_ID=bear_length=segment_ID=node_ID=reconstructed=NULL

            if(inherits(aRchi,"aRchi")==F) stop("The provided data is not of class aRchi")

            # checking arguments
            if(is.null(aRchi@QSM)) stop("The archi file does not contains skeleton")
            if(is.null(aRchi@pointcloud)) stop("The archi file does not contains a pointcloud")
            if(!is.numeric(max_dist)) stop("max_dist must be numeric")
            if(!is.numeric(sec_length)) stop("sec_length must be numeric")
            if(!is.numeric(th)) stop("th must be numeric")
            if(!is.numeric(d_clust)) stop("d_clust must be numeric")
            if(!is.character(method)) stop("method must be a character string")
            if(!method %in% c("statistical","absolute")) stop("method must be either statistical or absolute")

            skel = aRchi@QSM
            pc = aRchi@pointcloud@data
            pc[,index := 1:nrow(pc)] # add index to pc for matching

            #-
            if(!is.null(skel$reconstructed)) skel = skel[reconstructed == 1]

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
            pc = pc[!is.na(seg)] #- ambiguous points are removed

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
            pc[,dist := aRchi::dist2line(
              a=c(X,Y,Z),
              b=c(endX,endY,endZ),
              c=c(startX,startY,startZ)),
              by=row.names(pc)]

            # remove useless variable
            pc[,':='(startX = NULL,startY = NULL,startZ = NULL,
                     endX = NULL,endY = NULL,endZ = NULL)]

            #- remove points that are too far from the skeleton
            pc = pc[dist < max_dist]

            # compute the relative distance to skeleton by section
            pc[,sec:=round((dist_tip-(sec_length/2))/sec_length)*sec_length] # sections
            pc[,mean_dist := mean(dist),by=.(axis_ID,sec)] # mean dist by section and axis
            pc[,rel_dist := dist/mean_dist] # relative distance

            # identify non reconstructed
            if(method == "absolute"){
              pc[,NR := 1]
              pc[rel_dist > th, NR:=2]
            }
            if(method == "statistical"){
              pc[rel_dist > mean(rel_dist)+th*stats::sd(rel_dist), NR:=2]
            }

            ####### non reconstructed axes reconstruction
            # select non reconstructed axes
            NR = pc[NR == 2,.(X,Y,Z,seg,dist)]

            # non reconstructed axes clustering
            NR[,cluster := stats::cutree(fastcluster::hclust(stats::dist(NR[,1:3]), method = "single"), h = d_clust)]

            # identify the root and tip of the cluster
            NR[,tip := dist == max(dist),by=cluster] # tip is the farest
            NR[,seg := min(seg),by=cluster]
            # keep root and tip
            NR = NR[tip == TRUE]
            #NR = NR[order(cluster,seg,tip)] # sort

            # reconstruct NR as single segments
            NR[, ':='(
              startX = skel$startX[skel$cyl_ID == seg],
              startY = skel$startY[skel$cyl_ID == seg],
              startZ = skel$startZ[skel$cyl_ID == seg]
            ),by=cluster]

            NR[, ':='(
              endX = NR$X,
              endY = NR$Y,
              endZ = NR$Z
            )]

            NR[,cyl_ID := 1:nrow(NR)]
            NR[,cyl_ID := cyl_ID+max(skel$cyl_ID)]

            skel = data.table::rbindlist(list(skel[order(cyl_ID),.(startX,startY,startZ,endX,endY,endZ,cyl_ID)],NR[,.(startX,startY,startZ,endX,endY,endZ,cyl_ID)]),use.names = TRUE)
            data.table::setnames(skel,old="cyl_ID",new="seg_ID")

            ########################
            ##- COMPUTE TOPOLOGY -##
            ########################

            # assign bearer ID
            skel[,bearer_ID := 0]
            skel[,bearer_ID := FNN::knnx.index(data = skel[,4:6],query = skel[,1:3], algorithm = "kd_tree", k = 1)]
            skel[seg_ID == bearer_ID,bearer_ID := 0]

            # segments length
            skel[,length := sqrt( (startX - endX)^2 + (startY - endY)^2 + (startZ - endZ)^2)]

            # CLASS SEGMENTS INTO AXES
            ## compute total length beared by each segment
            skel[, bear_length := 0]

            pb <- progress::progress_bar$new(total = nrow(skel)*2,width = 60,
                                             format = "Computing topology  [:bar] :percent",clear=FALSE)
            for(s in rev(sort(skel$seg_ID))){
              pb$tick()
              childs = skel[seg_ID == s | bearer_ID == s]
              skel[seg_ID == s, bear_length := length+sum(childs$bear_length)]
            }

            ## the axis follows the longest bear_length
            skel[, axis_ID := 0]

            cur_seg = skel[bearer_ID == 0] # start with the trunk base
            cur_ID = 1 # curent axis ID
            cur_sec = 1 # curent section (for segment computation)

            skel[bearer_ID == 0, axis_ID := cur_ID]

            queue = c()
            while(min(skel$axis_ID)==0){
              pb$tick()
              childs = skel[seg_ID == cur_seg$seg_ID,section := cur_sec]
              childs = skel[bearer_ID == cur_seg$seg_ID]
              if(nrow(childs) >= 1){
                # if only one child -> it's in the same axis
                if(nrow(childs) == 1){
                  skel[seg_ID == childs$seg_ID, axis_ID := cur_ID]
                  cur_seg = childs
                }
                # if more than one child -> the one that bare longest structure belongs
                # to the same axis, the other ones goes to the queue
                if(nrow(childs) > 1){
                  skel[seg_ID == childs$seg_ID[which.max(childs$bear_length)], axis_ID := cur_ID]
                  cur_seg = childs[which.max(childs$bear_length),]
                  queue = c(queue , childs$seg_ID[-which.max(childs$bear_length)])

                  cur_sec = cur_sec + 1
                }
              }else{
                # if there is no child -> pick the first one in the queue and increment cur_ID
                cur_ID = cur_ID + 1 # increment cur_ID
                cur_seg = skel[seg_ID == queue[1]] # select curent segment
                skel[seg_ID == cur_seg$seg_ID, axis_ID := cur_ID]

                queue = queue[-1]

                cur_sec = cur_sec + 1
              }
            }

            # ADD BRANCHING ORDER
            cur_BO = skel[skel$axis_ID == 1] # axes of branching order 1
            skel[axis_ID == 1,branching_order := 1]
            BO = 2 # first branching order to detect
            while(nrow(cur_BO)>0){
              # find all child axes of the axes of the curent BO
              child_axes=skel[bearer_ID %in% cur_BO$seg_ID &
                                !(axis_ID %in% unique(cur_BO$axis_ID)),c(axis_ID)]
              # add the new BO to the child axes
              skel[axis_ID %in% child_axes,branching_order := BO]
              # select the child axes for the next round
              cur_BO = skel[skel$axis_ID %in% child_axes]
              BO = BO+1
            }

            skel[,reconstructed := 2L]
            skel[1:nrow(aRchi@QSM),reconstructed := 1L]

            # produce output data
            aRchi@QSM = data.table::data.table(
              startX = skel$startX,
              startY = skel$startY,
              startZ = skel$startZ,
              endX = skel$endX,
              endY = skel$endY,
              endZ = skel$endZ,
              cyl_ID = skel$seg_ID,
              parent_ID = skel$bearer_ID,
              radius = 0,
              length = skel$length,
              volume = 0,
              axis_ID = skel$axis_ID,
              segment_ID = skel$section,
              node_ID = 0,
              branching_order = skel$branching_order,
              reconstructed = skel$reconstructed
            )


            rm(skel)

            aRchi@QSM[,segment_ID := max(cyl_ID),by=segment_ID] # right segment_ID

            aRchi@QSM[,node_ID := aRchi@QSM$segment_ID[which(
              aRchi@QSM$cyl_ID %in% parent_ID & aRchi@QSM$segment_ID != segment_ID
            )],by = segment_ID]


            aRchi@pointcloud@data[pc$index,':='(dist_to_skel = pc$dist,rel_dist_to_skel = pc$rel_dist)]

            aRchi@operations$add_non_reconstructed = c(max_dist = max_dist,sec_length = sec_length,
                                                       method=method,th=th,d_clust=d_clust)
            return(aRchi)
          }


)

