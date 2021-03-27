#' Build a an object of class aRchi
#'
#' @description Build an object of class aRchi
#' @param QSM A data.table obtained from \code{\link{read_QSM}} function
#' @param point_cloud A point cloud. Either a LAS or a data.table with at least three columns with 3d coordinates (i.e X,Y,Z)
#' @param keep_original logical (Default = FALSE). Should the original branching order and axis be kept ? Otherwise, it is re-estimated.
#' @include aRchiClass.R
#' @seealso \code{\link{aRchi}}; \code{\link{write_aRchi}}; \code{\link{read_aRchi}}
#' @examples
#' file_QSM=system.file("extdata","Tree_1_TreeQSM.txt",package = "aRchi")
#' file_pc=system.file("extdata","Tree_1_point_cloud.las",package = "aRchi")
#' QSM=read_QSM(file_QSM,model="treeQSM")
#' pc=lidR::readLAS(file_pc)
#' # Make an object of class aRchi
#' Tree1_aRchi=build_aRchi(QSM=QSM,point_cloud=pc)


build_aRchi=function(QSM,point_cloud,keep_original = FALSE){
  branching_order=axis_ID=cyl_ID=parent_ID=extension_ID=startX=endX=endY=startY=startZ=endZ=bear_length=segment_ID=node_ID=parentID=radius=rad1=rad2=axis_ID=NULL



            ###############################
            # create the aRchi class file #
            ###############################
            aRchi = new("aRchi")

            ##########################
            # import the point cloud #
            ##########################
            if(!missing(point_cloud)){
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
            }

            ##################
            # import the QSM #
            ##################

            if(!missing(QSM)){
              # check the QSM
              if(ncol(QSM$QSM) < 7){
                stop("the provided QSM does not contains enought fields: must contain at least
           stratX | startY | startZ | endX | endY | endZ | radius")
              }
              if(QSM$model=="treeQSM"){
                aRchi@QSM = data.table::data.table(QSM$QSM)
              }else{

              # QSM is a data.table
              QSM = data.table::data.table(QSM$QSM)

              # does the axis_ID and branchind order must be computed by the finction ?
              compute_axis = TRUE
              compute_BO = TRUE

              # store original fields if needed
              if(keep_original){
                # if there are no additionnal columns, the topology is computed internally
                if(ncol(QSM) == 7){
                  warning("no original branching order and axis ID are provides, they are computed internally.")
                }

                # if there is 8 columns find which field must be computed
                if(ncol(QSM) == 8){
                  if(!colnames(QSM)[8] == "branching_order" | colnames(QSM)[8] == "axis_ID"){
                    warning("the original fields does not match the required field names. Topology is computed internally")
                  }
                  if(colnames(QSM)[8] == "branching_order"){
                    branching_orders = QSM$branching_order
                    QSM[,branching_order := NULL] # remove from table
                    cimpute_BO = FALSE
                  }
                  if(colnames(QSM)[8] == "axis_ID"){
                    axis_IDs = QSM$axis_ID # remove from table
                    QSM[,axis_ID := NULL]
                    compute_axis = FALSE
                  }
                }

                # if there are 9 columns compute both axis_ID and bronching order
                if(ncol(QSM) == 9){
                  if(!colnames(QSM)[8] == "branching_order" & colnames(QSM)[9] == "axis_ID"){
                    warning("the original fields does not match the required field names. Topology is computed internally")
                  }else{
                    branching_orders = QSM$branching_order
                    axis_IDs = QSM$axis_ID
                    QSM[,':='(axis_ID = NULL, branching_order = NULL)] # remove from table
                    compute_axis = FALSE
                    compute_BO = FALSE
                  }
                }
              }else{
                if(ncol(QSM) > 7) QSM = QSM[,1:7]
              }


              # check columns names for all required fields
              if(!identical(colnames(QSM[,1:7]),c("startX","startY","startZ","endX","endY","endZ","radius"))){
                data.table::setnames(QSM,c(c("startX","startY","startZ","endX","endY","endZ","radius")))
              }

              # store radius separately, remaining columns are cylinders coordinates
              radius = QSM$radius
              QSM[,radius := NULL]

              ############- compute topology

              # create cylinders ID
              QSM[,cyl_ID := 1:nrow(QSM)]

              # assign parent ID as the nearest cylinder end
              QSM[,parent_ID := 0]
              QSM[,parent_ID := FNN::knnx.index(data = QSM[,4:6],query = QSM[,1:3], algorithm = "kd_tree", k = 1)]
              QSM[cyl_ID == parent_ID,parent_ID := 0] # the one that is its own parent is the root of the QSMeton

              # assign successor ID as the nearest cylinder end
              QSM[,extension_ID := 0]
              QSM[,extension_ID := FNN::knnx.index(data = QSM[,1:3],query = QSM[,4:6], algorithm = "kd_tree", k = 1)]
              QSM[cyl_ID == extension_ID,extension_ID := 0] # the ones that is its own child is a tip
              QSM[! cyl_ID %in% parent_ID, extension_ID := 0]

              # compute segments length
              QSM[,length := sqrt( (startX - endX)^2 + (startY - endY)^2 + (startZ - endZ)^2)]

              # add radius
              QSM[,radius := radius]

              if(compute_axis){
                ## compute total length beared by each segment
                QSM[, bear_length := 0]

                pb <- progress::progress_bar$new(total = nrow(QSM)*2,width = 60,
                                                 format = " Segmenting axes [:bar] :percent",clear=FALSE)
                for(s in rev(sort(QSM$cyl_ID))){
                  pb$tick()
                  childs = QSM[cyl_ID == s | parent_ID == s]
                  QSM[cyl_ID == s, bear_length := length+sum(childs$bear_length)]
                }

                ## the axis follows the longest bear_length
                QSM[, axis_ID := 0]

                cur_cyl = QSM[parent_ID == 0] # start with the trunk base
                cur_ID = 1 # curent axis ID
                cur_seg = 1 # curent segment (for segment computation)

                QSM[parent_ID == 0, axis_ID := cur_ID]

                queue = c()
                while(min(QSM$axis_ID)==0){
                  pb$tick()
                  childs = QSM[cyl_ID == cur_cyl$cyl_ID,segment_ID := cur_seg]
                  childs = QSM[parent_ID == cur_cyl$cyl_ID]
                  if(nrow(childs) >= 1){
                    # if only one child -> it's in the same axis
                    if(nrow(childs) == 1){
                      QSM[cyl_ID == childs$cyl_ID, axis_ID := cur_ID]
                      cur_cyl = childs
                    }
                    # if more than one child -> the one that bare longest structure belongs
                    # to the same axis, the other ones goes to the queue
                    if(nrow(childs) > 1){
                      QSM[cyl_ID == childs$cyl_ID[which.max(childs$bear_length)], axis_ID := cur_ID]
                      cur_cyl = childs[which.max(childs$bear_length),]
                      queue = c(queue , childs$cyl_ID[-which.max(childs$bear_length)])

                      cur_seg = cur_seg + 1
                    }
                  }else{
                    # if there is no child -> pick the first one in the queue and increment cur_ID
                    cur_ID = cur_ID + 1 # increment cur_ID
                    cur_cyl = QSM[cyl_ID == queue[1]] # select curent segment
                    QSM[cyl_ID == cur_cyl$cyl_ID, axis_ID := cur_ID]

                    queue = queue[-1]

                    cur_seg = cur_seg + 1
                  }
                }
              }else{
                # axis_ID
                QSM[,axis_ID := axis_IDs]
                QSM[,segment_ID := 0]
                # segment_ID
                cur_seg = 1 # curent segment (for segment computation)
                pb <- progress::progress_bar$new(total = nrow(QSM),width = 60,
                                                 format = " Segmenting segments [:bar] :percent",clear=FALSE)
                for(i in sort(unique(QSM$axis_ID))){
                  # for each segment of each axis
                  axis = QSM[axis_ID == i]
                  for(j in 1:nrow(axis)){
                    pb$tick()
                    # the segment ID of the cylinder is the curent segment ID
                    QSM[cyl_ID == axis$cyl_ID[j],segment_ID := cur_seg]

                    # increment curent segment ID if it's a branching point
                    if(nrow(QSM[parent_ID == axis$cyl_ID[j]])>1) cur_seg = cur_seg+1
                  }

                  # increment segment ID when it's a different axis
                  cur_seg = cur_seg+1
                }
              }
              # remove the bear_length field
              QSM[,bear_length := NULL]

              if(compute_BO){
                cur_BO = QSM[QSM$axis_ID == 1] # axes of branching order 1
                QSM[axis_ID == 1,branching_order := 1]
                BO = 2 # first branching order to detect
                while(nrow(cur_BO)>0){
                  # find all child axes of the axes of the curent BO
                  child_axes=QSM[parent_ID %in% cur_BO$cyl_ID &
                                   !(axis_ID %in% unique(cur_BO$axis_ID)),c(axis_ID)]
                  # add the new BO to the child axes
                  QSM[axis_ID %in% child_axes,branching_order := BO]
                  # select the child axes for the next round
                  cur_BO = QSM[QSM$axis_ID %in% child_axes]
                  BO = BO+1
                }
              }else{
                QSM[,branching_order := branching_orders]
              }

              # right segment_ID
              QSM[,segment_ID := max(cyl_ID),by=segment_ID]

              # compute the node ID
              QSM[,node_ID := QSM$segment_ID[which(
                QSM$cyl_ID %in% parent_ID & QSM$segment_ID != segment_ID
              )],by = segment_ID]
              # Attribute 0 instead of NA for the first segment
              QSM[is.na(node_ID)]$node_ID=0
              # Bastien: sorry for that below... I did it fast ^^
              names(QSM)[which(names(QSM)=="radius")]="radius_cyl"
              QSM$volume=(pi*(QSM$radius_cyl)^2)*QSM$length
              QSM=QSM[,c("startX","startY","startZ", "endX","endY","endZ","cyl_ID","parent_ID","extension_ID","radius_cyl","length","volume","axis_ID","segment_ID","node_ID","branching_order")]

              aRchi@QSM = QSM
              }
            }

            return(aRchi)
          }


