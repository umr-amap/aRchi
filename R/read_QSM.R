#' Read a QSM
#'
#' @description Read a QSM file generated with treeQSM, simpletree, simpleforest or pypetree.
#' @param file The directory to the QSM file path
#' @param model `treeQSM`, `simpletree`, `simpleforest` or `pypetree` depending on the algorithm used to generate the QSM
#' @return  a data.table that can be used to build an aRchi object (see function [build_aRchi())])
#' @seealso [aRchi()] the aRchi class;[build_aRchi())] to build an object of class `aRchi`; [skeletonize())]; to make a QSM from a point cloud in an aRchi object.
#' @include aRchiClass.R
#' @examples
#' file=system.file("extdata","Tree_1_TreeQSM.txt",package = "aRchi")
#' QSM=read_QSM(file,model="treeQSM")

read_QSM=function(file,model){
  parent_ID=extension_ID=cyl_ID=branch_ID=PositionInBranch_ID=NULL

  data=data.table::fread(file)

  if(model=="treeQSM"){
  data.table::setnames(data,c("radius","length","startX","startY","startZ","axisX","axisY","axisZ","parent_ID","extension_ID","added","UnmodRadius","branch_ID","BranchOrder","PositionInBranch_ID"))

  ####################
  #######Step 1#######
  ####################
  # Adapt the format to have a cylinder based data instead of a circle based data. Similar to AMAPstudio-scan Format.

  data$endX=data$startX+(data$length*data$axisX) #estimates the coordinates of the end of the cylinders. X
  data$endY=data$startY+(data$length*data$axisY) #estimates the coordinates of the end of the cylinders. Y
  data$endZ=data$startZ+(data$length*data$axisZ) #estimates the coordinates of the end of the cylinders. Z
  data$cyl_ID=1:nrow(data)
  data[-1,c("startX","startY","startZ")]=data[data$parent_ID,c("endX","endY","endZ")] # Replace coordinates of the start of the children by  the end of the parents. AMAPstudio format.
  data$radius_cyl=0 # Radius of the cylinder
  data[-1,"radius_cyl"]=data[data$parent_ID,"radius"] # The radius of the cylinder is the radius of the start circle i.e the parent. Except when there is a ramification, see below
  data[data$PositionInBranch_ID==1,"radius_cyl"]=data[data$PositionInBranch_ID==1,"radius"] # When a ramification: The radius of the first cylinder of a branch is the radius of the daughter. AMAPstudio format
  names(data)[13:14]=c("axis_ID","branching_order")
  out=data[,c("startX","startY","startZ", "endX","endY","endZ","radius_cyl","axis_ID","branching_order")]
  }

  if(model == "pypetree"){
    # build initial table with the terminal nodes
    out = data.table::data.table(
      startX = 0,
      startY = 0,
      startZ = 0,
      endX = data$x,
      endY = data$y,
      endZ = data$z,
      ID = data$node_id,
      parentID = data$parent_node_id,
      rad1 = data$radius
    )

    # attribute the starting based on parent ID
    out[!is.na(parentID),':='(
      startX = out$endX[which(data$node_id == parentID)],
      startY = out$endY[which(data$node_id == parentID)],
      startZ = out$endZ[which(data$node_id == parentID)],
      rad2 = out$rad1[which(data$node_id == parentID)]
    ),by = parentID]

    # remove the first node (the starting point without parent)
    out = out[!is.na(parentID)]

    # compute the terminal radius -> the mean of the two extremities
    out[,radius := (rad1 + rad2) /2]

    # remove useless columns
    out[,':='(
      ID = NULL,
      parentID = NULL,
      rad1 = NULL,
      rad2 = NULL
    )]
  }

  if(model == "simpletree"){

    # extract the columns of interest
    out = data.table::data.table(
      startX = data$startX,
      startY = data$startY,
      startZ = data$startZ,
      endX = data$endX,
      endY = data$endY,
      endZ = data$endZ,
      radius = data$radius,
      axis_ID = data$branch_ID,
      branching_order = data$branch_Order
    )

    # adjust the axis_ID so it starts to 1
    out[axis_ID == -1, axis_ID := 0]
    out[, axis_ID := axis_ID + 1]

    warning("Using the SimpleTree axis_ID might lead to problems when using the
          functions of the aRchi package")

  }

  if(model == "simpleforest"){

    # extract the columns of interest
    out = data.table::data.table(
      startX = data$startX,
      startY = data$startY,
      startZ = data$startZ,
      endX = data$endX,
      endY = data$endY,
      endZ = data$endZ,
      radius = data$radius,
      axis_ID = data$BranchID,
      branching_order = data$BranchOrder
    )

    # adjust the axis_ID so it starts to 1
    out[axis_ID == -1, axis_ID := 0]
    out[, axis_ID := axis_ID + 1]

    warning("Using the SimpleTree axis_ID might lead to problems when using the
          functions of the aRchi package")

  }
  return(out)
}


