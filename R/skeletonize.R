
#' Build the skeleton of a tree point cloud
#'
#' @param aRchi an object of class aRchi containing at least a point cloud.
#' @param D numeric. The distance of research for point neighborhood. Sets the
#' layer thickness. See description for details.
#' @param progressive logical. Should the clustering distance be progressive ?
#' See description for details.
#' @param cl_dist numeric. The clustering distance. If \code{pregressive = FALSE}
#' sets the clustering distance for all the point cloud. If \code{pregressive = TRUE}
#' sets the minimum clustering distance to be used. See description for details.
#' @param max_d The maximum searching distance for skeleton building. See
#' description for details.
#'
#' @details The skeletonization algorithm works in four steps. At STEP 1 the
#' point cloud is divided in layers of regular thickness (defined by parameter
#' \code{D}). To do so, the tree base is first defined as the first layer. Then, all the points
#' that are within \code{D} of any points of the layer \emph{N}  belong to the layer
#' \emph{N+1}. This process continues until no more points are found within \code{D}.
#' If there are some remaining points (that are further than \code{D} to any
#' already classified points), a new layer is defined as the point that is the closest of
#' the already classified point. This continues until no more points remain unclassified.
#' At STEP 2, the layers are divided into clusters based on point distance: two point that
#' are further than a given distance are considered as being part of two different objects.
#' Two possibilities are available to define the clustering distance. First, it
#' remains constant and is defined by \code{cl_dist}. Second, the distance is
#' defined as the average distance of the points of the previous layer to the
#' center of their corresponding cluster (this is achieved by setting
#' \code{progressive = TRUE}). In this latter case, \code{cl_dist} defines the minimal
#' clustering distance that can be used. This option helps to adapt the clustering
#' distance to the size of the actual objects (i.e. branch sections) that have to be to clustered.
#' By default the first layer is assumed as being part of a single cluster. At
#' STEP 3, the cluster centers are used as node of the skeleton and are combined
#' to iteratively build the skeleton. At first iteration, the node with the
#' smallest Z value is used as the root node. At subsequent iterations, the node
#' located within \code{max_d} and that is the closest to the root node
#' is integrated into the skeleton and is selected as new root node. This process
#' continues until no node is found within \code{max_d} to the root node
#' (either because the cluster is a branch tip or because there is a gap in the point cloud).
#' In this case, the node that belong to the layer that was produced at earliest
#' iteration of STEP 1 and that is the closer to a skeleton node is selected as new root node.
#' This process ends once all nodes are integrated into the skeleton. At STEP 4, a basic
#' tree topology is computed. First, an unique ID is assigned to all segments of the skeleton
#' (i.e. the segment that link two nodes) and its parent segment ID is retrieved.
#' The segments are then classified into axes. An axis being defined as a continuous set of
#' segments that always follow the segment that support the longest portion of the skeleton.
#' Axes are then partitioned into axes segments defined as the portion of an axis
#' located between two branching point or between a branching point and an extremity of
#' the axis. The axis branching order is then computed as the number of branching point
#' observed between the tree base and the axis first segment + 1.
#'
#'
#' @return an aRchi file containing the original point cloud and the corresponding skeleton
#' @export
#'
#' @examples
#' \donttest{
#' #################################################################
#' # Example with a small high quality point cloud : using default #
#' # default parameters to detect fine architectural details       #
#' #################################################################
#'# import a point cloud
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
#'
#'
#' ##############################################################
#' # Example with a large point cloud with a lot of occlusion : #
#' # parameters selected for speed                              #
#' ##############################################################
#' # import a point cloud
#' tls=system.file("extdata","Tree_1_point_cloud.las",package = "aRchi")
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
#' aRchi = skeletonize_pc(aRchi, D = 0.5, cl_dist = 0.2, max_d = 1)
#'
#' # smooth the skeleton
#' aRchi = smooth_skeleton(aRchi)
#'
#' # plot the skeleton
#' plot(aRchi,show_point_cloud = TRUE)
#' }
setGeneric("skeletonize_pc", function(aRchi,D = 0.03,progressive = TRUE,cl_dist = 0.02,max_d = 0.05)
{
  standardGeneric("skeletonize_pc")
})

#' @rdname skeletonize_pc
#' @export
setMethod("skeletonize_pc", signature = "aRchi", function(aRchi,D = 0.03,progressive = TRUE,cl_dist = 0.02,max_d = 0.05)
{
  # to pass CRAN check
  .=ID=X=Y=Z=axis_ID=bear_length=bearer_ID=branching_order=cluster=cyl_ID=dist=endX=endY=endZ=iter=node_ID=parent_ID=radius=section=seg_ID=segment_ID=startX=startY=startZ=extension_ID=NULL

  data = aRchi@pointcloud@data

  # ==================================
  # Step 1. iteratively compute layers
  # ==================================

  cat("(1/4) Computing layers\n")
  data = cpp_compute_layers(as.matrix(data), D)
  data.table::setDT(data)

  # =========================================================
  # Step 2. clustering non connected components in each layer
  # =========================================================

  # if the clustering distance vary
  progressive = T
  if (progressive == T)
  {

    first = TRUE # first iteration ?

    cat("(2/4) Clustering layers\n")
    #pb <- progress::progress_bar$new(total = length(unique(data$iter)),width = 60, format = " (2/4) Clustering layers   [:bar] :percent", clear=FALSE)

    for(i in sort(unique(data$iter)))
    {
     # pb$tick()

      if(first)
      {
        # at the first iteration the clustering distance (cl_dist) is the radius of
        # the first layer

        # only one cluster at iteration 1 (the tree base)
        data[iter == i,cluster := 1]
        # compute the radius
        data[iter == i,radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  )]
        # the average radius is the first clustering distance
        cl_d = mean(data$radius[data$iter == i])
        first = FALSE
      }
      else
      {
        # at subsequent iterations the clustering distance is the average radius
        # of clustered objects in the curent layer

        # keep points in the curent layer
        in_iter = which(data$iter==i)

        # if there are more than one point in the layer -> cluster it
        if(length(in_iter)>=2)
        {
          # clustering
          data[in_iter,cluster := stats::cutree(fastcluster::hclust(stats::dist(data[in_iter,1:3]), method = "single"), h = cl_d)]

          # compute the radius for each cluster
          data[in_iter,radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  ),by = cluster]

          # the new clustering distance is the average radius of all clusters
          cl_d = mean(data$radius[data$iter == i])/2

          # if the computed clustering distance is too small
          # -> replace by the user defined minimum distance
          if(cl_d < cl_dist) cl_d = cl_dist
        }
        else
        {
          # if there is only one point, there is only one cluster
          data[in_iter,cluster := 1]
        }
      }
    }
  }
  else # if the clustering distance remains constant
  {
    pb <- progress::progress_bar$new(total = length(unique(data$iter)),width = 60, format = " (2/4) Clustering layers  [:bar] :percent", clear=FALSE)

    for(i in sort(unique(data$iter)))
    {
      pb$tick()

      # keep points in the curent layer
      in_iter = which(data$iter==i)
      if(length(in_iter)>=2)
      {
        # clustering
        data[in_iter,cluster := stats::cutree(fastcluster::hclust(stats::dist(data[in_iter,1:3]), method = "single"), h = cl_dist)]

        # compute radius by layer
        data[in_iter,radius := sqrt((X-mean(X))^2 + (Y-mean(Y))^2 + (Z-mean(Z))^2  ),by = cluster]
      }
      else
        {
        # if there is only one point, there is only one cluster
        data[in_iter,cluster := 1]
      }
    }
  }

  # ==========================
  # Step 3. Build the skeleton
  # ==========================

  cat("(3/4) Building skeleton\n")
  skel = cpp_build_skeleton(data, max_d)

  # segments length
  dx = skel$startX - skel$endX
  dy = skel$startY - skel$endY
  dz = skel$startZ - skel$endZ
  skel$length = sqrt(dx^2 + dy^2 + dz^2)
  skel = skel[skel$length > 0,]

  # ========================
  # Step 4. compute topology
  # ========================

  cat("(4/4) Computing topology\n")

  # create segment ID
  data.table::setDT(skel)
  skel[, seg_ID := 1:nrow(skel)]

  # assign bearer ID
  skel[, bearer_ID := 0]
  skel[, bearer_ID := FNN::knnx.index(data = skel[,4:6], query = skel[,1:3], algorithm = "kd_tree", k = 1)]
  skel[seg_ID == bearer_ID, bearer_ID := 0]

  # CLASS SEGMENTS INTO AXES
  ## compute total length beared by each segment
  skel[, bear_length := 0]

  res = cpp_compute_topology(skel)
  skel = cbind(skel, res)
  data.table::setDT(skel)

  # ADD BRANCHING ORDER
  cur_BO = skel[skel$axis_ID == 1] # axes of branching order 1
  skel[axis_ID == 1,branching_order := 1]
  BO = 2 # first branching order to detect
  while(nrow(cur_BO)>0)
  {
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

  # add segment ID
  aRchi@QSM[, segment_ID := max(cyl_ID), by=segment_ID] # right segment_ID

  # add node ID
  aRchi@QSM[,node_ID := aRchi@QSM$segment_ID[which(
    aRchi@QSM$cyl_ID %in% parent_ID & aRchi@QSM$segment_ID != segment_ID
  )],by = segment_ID]

  # add extension ID
  aRchi@QSM[,extension_ID := FNN::knnx.index(data = aRchi@QSM[,1:3],query = aRchi@QSM[,4:6], algorithm = "kd_tree", k = 1)]
  aRchi@QSM[cyl_ID == extension_ID,extension_ID := 0] # the ones that is its own child is a tip
  aRchi@QSM[! cyl_ID %in% parent_ID, extension_ID := 0]

  # keep the operation in memory
  aRchi@operations$skeletonize = list(D=D,progressive=progressive,cl_dist=cl_dist,max_d=max_d)

  return(aRchi)
})


