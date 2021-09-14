
#' Build the skeleton of a tree point cloud
#'
#' @param aRchi an object of class aRchi containing at least a point cloud
#' @param D numeric. The distance of research for point neighborhood. Sets the
#' layer tickness. See description for details. Default = 0.02
#' @param progressive logical. Should the clustering distance be progressive ?
#' See description for details. Default = TRUE.
#' @param cl_dist numeric. The clustering distance. If \code{pregressive = FALSE}
#' sets the clustering distance for all the point cloud. If \code{pregressive = TRUE}
#' sets the minimum clustering distance to be used. See description for details.
#' Default = 0.02.
#' @param max_d The maximum searching distance for skeleton building. See
#' description for details.Default = 0.05.
#' @param K The number of points for neighboor searching. See description for
#' details. Default = 100.
#'
#' @description The skeletonization algorithm works in four steps. At STEP 1 the
#' point cloud is devided in layers of regular tickness (defined by parameter
#' \code{D}). To do so, the tree base is fisrt defined as first layer and all the points
#' of the point cloud that are within \code{D} of any points of the layer are defined
#' as the next layer. This process continues until that no more points are found.
#' If there are some remaining points (that are further than \code{D} to any
#' classifyed points),a new layer is defined as the point that is the closer of
#' the already classified point within the layer produced at earliest iteration.
#' This continues until no more points remain unclassified. at STEP 2, the layers
#' are devided into clusters based on point distance: two point that are further
#' than a given distance are considered as being part of two different objects.
#' Two possibilities are available to define the clustering distance. First, it
#' remains constant and is defined by \code{cl_dist}. Second, the distance is
#' defined as the average distance of the points of the previous layer to the
#' center of their corresponding layer (this is achieved by setting
#' \code{progressive = TRUE}). In this case, \code{cl_dist} defines the minimal
#' clustering distance that can be used. This option helps to adapt the clustering
#' distance to the size of the actual obtects (i.e. branch sections) to cluster.
#' By default the first layer is assumed as being part of a single cluster. At
#' STEP 3, the cluster centers are combined to build the skeleton. To do so, an
#' iterative and hierarchical process starting with a starting point defined as
#' the cluster center with the smallest Z value is used. A cluster is automatically
#' conected to its nearest neighbor located within \code{max_d}. If no neighbour is
#' detected (because the cluster is a branch tip or due to occlusion), the process
#' restarts with the cluster located in the earliest layer that is the closer to
#' an already used cluster as new starting point. This process ends once all cluster
#' centers are connected.
#'
#'
#' @return an aRchi file containing the original point cloud and the corresponding skeleton
#' @export
#'
#' @examples
#' \donttest{
#' # import a point cloud
#' tls=system.file("extdata","Tree_2_point_cloud.las",package = "aRchi")
#' tls = lidR::readLAS(tls)
#'
#' # build an empty aRchi file and add the point cloud
#' aRchi = aRchi::build_aRchi()
#' aRchi = aRchi::add_pointcloud(aRchi,point_cloud = tls)
#'
#' # plot the point cloud
#' plot(aRchi@pointcloud)
#'
#' # build a skeleton from the point cloud
#' aRchi = skeletonize_pc(aRchi)
#'
#' # smooth the skeleton
#' aRchi = smooth_skeleton(aRchi)
#'
#' # plot the skeleton
#' plot(aRchi,show_point_cloud = TRUE)
#' }

setGeneric("skeletonize_pc",
           function(aRchi,D = 0.03,progressive = TRUE,cl_dist = 0.02,max_d = 0.05,K=100){standardGeneric("skeletonize_pc")}
)

#' @rdname skeletonize_pc
#' @export

setMethod("skeletonize_pc",
          signature = "aRchi",
          function(aRchi,D = 0.03,progressive = TRUE,cl_dist = 0.02,max_d = 0.05,K=100){

            # to pass CRAN check
            .=ID=X=Y=Z=axis_ID=bear_length=bearer_ID=branching_order=cluster=cyl_ID=dist=endX=endY=endZ=iter=node_ID=parent_ID=radius=
              section=seg_ID=segment_ID=startX=startY=startZ=extension_ID=NULL

            data = aRchi@pointcloud@data

            ################################
            #- iteratively compute layers -#
            ################################

            #- add required columns
            data[,':='(ID=1:nrow(data),iter = -1)]

            #- first layer
            lay = data[Z <= min(Z)+0.1]

            #- the points in the first layer belong to iteration 1
            data[lay$ID,iter:=1]

            #- create the data for searching
            dat2 = data # duplicate data
            dat2 = dat2[-lay$ID] # remove the first layer

            # cat("Computing layers \n")
            i=2 # iterations

            pb <- progress::progress_bar$new(total = nrow(data),width = 60,
                                             format = " (1/4) Computing layers    [:bar] :percent",clear=FALSE)
            while(nrow(dat2)>0){

              # adjust K if there are not enough remaining points
              if(nrow(dat2)< K){k = nrow(dat2)}else{k = K}

              # compute the distance maximum distance of each point in the data to the
              # points in the layer
              dat2[,dist := Rfast::rowMaxs(FNN::knnx.dist(data = lay[,1:3],
                                                          query = dat2[, 1:3],
                                                          algorithm = "kd_tree",
                                                          k = 1),value=TRUE)]

              # keep the points that are within the layer tickness to produce the next layer
              lay = dat2[dist <= D]

              # if there is no points in the new layer, the point from the data the closest
              # to an already classified point is used a new root
              if(nrow(lay)==0) lay = dat2[which.min(FNN::knnx.dist(data = data[iter != -1,1:3],
                                                                   query = dat2[, 1:3],
                                                                   algorithm = "kd_tree",
                                                                   k = 1))]

              # the curent iteration is added to the points in the curent layer
              data[lay$ID,iter := i]

              # the points of the curent layer are removed from the search data
              dat2 = dat2[!ID %in% lay$ID,]
              i=i+1
              pb$update(1-(nrow(dat2)/nrow(data)))
            }

            #######################################################
            #- clustering non connected components in each layer -#
            #######################################################
            if(progressive == T){ # if the clustering distance vary

              first = TRUE # first iteration ?

              #pb=txtProgressBar(min=1,max=length(unique(data$iter)),width = 33,style=3,label = "Clustering layers \n")
              pb <- progress::progress_bar$new(total = length(unique(data$iter)),width = 60,
                                               format = " (2/4) Clustering layers   [:bar] :percent",clear=FALSE)
              for(i in sort(unique(data$iter))){
                #setTxtProgressBar(pb,i)
                pb$tick()

                if(first){
                  # at the first iteration the clustering distance (cl_dist) is the radius of
                  # the first layer

                  # only one cluster at iteration 1 (the tree base)
                  data[iter == i,cluster := 1]
                  # compute the radius
                  data[iter == i,radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  )]
                  # the average radius is the first clustering distance
                  cl_d = mean(data$radius[data$iter == i])
                  first = FALSE
                }else{
                  # at subsequent iterations the clustering distance is the average radius
                  # of clustered objects in the curent layer

                  # keep points in the curent layer
                  in_iter = which(data$iter==i)

                  # if there are more than one point in the layer -> cluster it
                  if(length(in_iter)>=2){

                    # clustering
                    data[in_iter,cluster := stats::cutree(fastcluster::hclust(stats::dist(data[in_iter,1:3]), method = "single"), h = cl_d)]

                    # compute the radius for each cluster
                    data[in_iter,radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  ),by = cluster]

                    # the new clustering distance is the average radius of all clusters
                    cl_d = mean(data$radius[data$iter == i])/2

                    # if the computed clustering distance is too small
                    # -> replace by the user defined minimum distance
                    if(cl_d < cl_dist) cl_d = cl_dist
                  }else{
                    # if there is only one point, there is only one cluster
                    data[in_iter,cluster := 1]
                  }
                }
              }
            }else{ # if the clustering distance remains constant
              pb <- progress::progress_bar$new(total = length(unique(data$iter)),width = 60,
                                               format = " (2/4) Clustering layers  [:bar] :percent",clear=FALSE)
              for(i in sort(unique(data$iter))){
                pb$tick()

                # keep points in the curent layer
                in_iter = which(data$iter==i)
                if(length(in_iter)>=2){
                  # clustering
                  data[in_iter,cluster := stats::cutree(fastcluster::hclust(stats::dist(data[in_iter,1:3]), method = "single"), h = cl_dist)]

                  # compute radius by layer
                  data[in_iter,radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  ),by = cluster]
                }else{
                  # if there is only one point, there is only one cluster
                  data[in_iter,cluster := 1]
                }
              }
            }

            ########################
            #- Build the skeleton -#
            ########################
            # keep only the center of each cluster
            cl = data[,.(X = mean(X),Y=mean(Y), Z=mean(Z), iter = iter),by = .(iter,cluster)]
            cl = cl[,3:6]

            # add required columns
            cl[,':='(ID = 1:nrow(cl), done = 0)]

            # create the searching data
            cl2 = cl

            # the root is the luster center with the lower Z coordinate
            root = cl2[Z == min(Z)]

            # note the root point as being done
            cl[root$ID, done := 1]

            # remove the root from the search data
            cl2 = cl2[-root$ID]

            # create a data.table to store the skeleton
            skel = data.table::data.table(matrix(ncol=6))
            data.table::setnames(skel,c("startX","startY","startZ","endX","endY","endZ"))


            pb <- progress::progress_bar$new(total = nrow(data),width = 60,
                                             format = " (3/4) Building skeleton   [:bar] :percent",clear=FALSE)
            while(nrow(cl2)>0){

              # keep the coordinates of the root as it will be the start of the segment
              start=root[,1:3]

              # keep the cluster centers that are close to the root
              # distance
              cl2[,dist:=sqrt((X-root$X)^2 + (Y-root$Y)^2 + (Z-root$Z)^2  )]

              # the new root is the cluster center the closer to the root that:
              #     - falls within the maximum search distance
              #     - is located in a layer produced after the root's layer
              root = cl2[dist == min(dist) & dist <= max_d & iter > root$iter]

              if(nrow(root)>0){
                # if the new root is selected

                # note the cluster as done
                cl[root$ID, done := 1]

                # remove the cluster from the search data
                cl2 = cl2[ID != root$ID]

                # add the segment start and tip (start and root respectively) to the skeleton data
                skel = data.table::rbindlist(list(skel,data.table::data.table(start[,1:3],root[,1:3])),use.names = F)
              }else{
                # if there is no root selected -> the new root is within the smaller remaining layer
                # the closer to an already classifyed cluster center

                # keep the points in the smaller layer
                root = cl2[which.min(iter)]

                # find the cluster center the closest to any already classifyed cluster center
                root = root[which.min(FNN::knnx.dist(data = cl[done == 1,1:3], query = root[, 1:3], algorithm = "kd_tree",k = 1))]

                # fing the start point of the segment as the closest already classifyed cluster center
                d = sqrt((cl$X[cl$done==1]-root$X)^2 + (cl$Y[cl$done==1]-root$Y)^2 + (cl$Z[cl$done==1]-root$Z)^2  )
                done = cl[done == 1]
                start = done[ d == min(d) ]

                # mark the root's cluster center as done
                cl[root$ID, done := 1]

                # remove the root from the search data
                cl2 = cl2[ID != root$ID]

                # add the new segment to the skeleton data set
                skel = data.table::rbindlist(list(skel,data.table::data.table(start[,1:3],root[,1:3])),use.names = F)
              }
              pb$update(1-(nrow(cl2)/nrow(cl)))
            }
            #- remove the first row that only contains NAs
            skel = skel[-1,]

            ######################
            #- compute topology -#
            ######################

            # create segment ID
            skel[,seg_ID := 1:nrow(skel)]

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
                                             format = " (4/4) Computing topology  [:bar] :percent",clear=FALSE)
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
              extension_ID = 0,
              radius_cyl = 0,
              length = skel$length,
              volume = 0,
              axis_ID = skel$axis_ID,
              segment_ID = skel$section,
              node_ID = 0,
              branching_order = skel$branching_order
            )


            rm(skel)

            # add segment ID
            aRchi@QSM[,segment_ID := max(cyl_ID),by=segment_ID] # right segment_ID

            # add node ID
            aRchi@QSM[,node_ID := aRchi@QSM$segment_ID[which(
              aRchi@QSM$cyl_ID %in% parent_ID & aRchi@QSM$segment_ID != segment_ID
            )],by = segment_ID]

            # add extension ID
            aRchi@QSM[,extension_ID := FNN::knnx.index(data = aRchi@QSM[,1:3],query = aRchi@QSM[,4:6], algorithm = "kd_tree", k = 1)]
            aRchi@QSM[cyl_ID == extension_ID,extension_ID := 0] # the ones that is its own child is a tip
            aRchi@QSM[! cyl_ID %in% parent_ID, extension_ID := 0]

            # keep the operation in memory
            aRchi@operations$skeletonize = list(D=D,progressive=progressive,cl_dist=cl_dist,max_d=max_d,K=K)

            return(aRchi)
          }
)


