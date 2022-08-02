#' Read a QSM
#'
#' @description Read a QSM file generated with treeQSM, simpletree, simpleforest or pypetree.
#' @param file The directory to the QSM file path.
#' @param model \code{treeQSM}, \code{simpletree}, \code{simpleforest} or \code{pypetree} depending on the algorithm used to generate the QSM
#' @return  a list containing a data.table with the QSM and a character with the model name. This list can be used to build an aRchi object (see function \code{\link{build_aRchi}})
#' @seealso \code{\link{aRchi}} the aRchi class;\code{\link{build_aRchi}} to build an object of class \code{aRchi}
#' @details
#' For \code{treeQSM} model, .mat from treeQSM are allowed. the old format (V2.3) as well as the new format (V2.4) are allowed.
#' @include aRchiClass.R
#' @examples
#' file=system.file("extdata","Tree_1_TreeQSM.txt",package = "aRchi")
#' QSM=read_QSM(file,model="treeQSM")

read_QSM=function(file,model){
  parent_ID=extension_ID=cyl_ID=branch_ID=PositionInBranch_ID=parentID=radius=rad1=rad2=axis_ID=NULL



  if(model=="treeQSM"){
    if(stringr::str_detect(file,".mat")){
      mat=R.matlab::readMat(file)

      if(is.null(mat$OptQSM)){
        if(is.null(mat$QSM)){
          warning("No QSM found")}else{


            warning(paste0(ncol(mat$QSM[,,]), " QSMs were found. The first one only is used. Please use select_optimum function of treeQSM matlab algorithm if you want to get the optimal QSM"))
            cylinder=mat$QSM[,,1]$cylinder
            data=data.table::data.table(cbind(cylinder[,,]$radius,cylinder[,,]$length,cylinder[,,]$start,cylinder[,,]$axis,cylinder[,,]$parent,cylinder[,,]$extension,cylinder[,,]$added,cylinder[,,]$UnmodRadius,cylinder[,,]$branch,cylinder[,,]$BranchOrder,cylinder[,,]$PositionInBranch))
            data.table::setnames(data,c("radius","length","startX","startY","startZ","axisX","axisY","axisZ","parent_ID","extension_ID","added","UnmodRadius","branch_ID","BranchOrder","PositionInBranch_ID"))

          }
      }else{
        cylinder=mat$OptQSM[,,]$cylinder
        data=data.table(cbind(cylinder[,,]$radius,cylinder[,,]$length,cylinder[,,]$start,cylinder[,,]$axis,cylinder[,,]$parent,cylinder[,,]$extension,cylinder[,,]$added,cylinder[,,]$UnmodRadius,cylinder[,,]$branch,cylinder[,,]$BranchOrder,cylinder[,,]$PositionInBranch))
        data.table::setnames(data,c("radius","length","startX","startY","startZ","axisX","axisY","axisZ","parent_ID","extension_ID","added","UnmodRadius","branch_ID","BranchOrder","PositionInBranch_ID"))
      }
    }else{
      data=data.table::fread(file)
      if(ncol(data)==17){
        data=data[,-c(14,15)]
      }
      data.table::setnames(data,c("radius","length","startX","startY","startZ","axisX","axisY","axisZ","parent_ID","extension_ID","added","UnmodRadius","branch_ID","BranchOrder","PositionInBranch_ID"))

    }

  ####################
  #######Step 1#######
  ####################
  # Adapt the format to have a cylinder based data instead of a circle based data. Similar to AMAPstudio-scan Format.

  data$endX=data$startX+(data$length*data$axisX) #estimates the coordinates of the end of the cylinders. X
  data$endY=data$startY+(data$length*data$axisY) #estimates the coordinates of the end of the cylinders. Y
  data$endZ=data$startZ+(data$length*data$axisZ) #estimates the coordinates of the end of the cylinders. Z
  data$cyl_ID=1:nrow(data)

  if(nrow(data[parent_ID==0])>1){
    lost_segment=nrow(data[parent_ID==0])-1

    branch_2_remove=data[parent_ID==0]$branch_ID[-1]
    branch_2_remove=unique(c(branch_2_remove,unique(data[parent_ID%in%data[branch_ID%in%branch_2_remove]$cyl_ID]$branch_ID)))
    branch_2_remove=unique(c(branch_2_remove,unique(data[parent_ID%in%data[branch_ID%in%branch_2_remove]$cyl_ID]$branch_ID)))
    data=data[!branch_ID%in%branch_2_remove]
    branch_2_remove=unique(c(branch_2_remove,unique(data[parent_ID%in%data[branch_ID%in%branch_2_remove]$cyl_ID]$branch_ID)))
    data=data[!branch_ID%in%branch_2_remove]
    warning(paste0(length(branch_2_remove)," lost branches (not link to the rest of the tree) were removed"))
  }

  # data[-1,c("startX","startY","startZ")]=data[data$parent_ID,c("endX","endY","endZ")] # Replace coordinates of the start of the children by  the end of the parents. AMAPstudio format.
  data_parent=plyr::adply(data,1,function(x){data[cyl_ID==x$parent_ID]})
  data[-1,c("startX","startY","startZ")]=data_parent[,c("endX","endY","endZ")] # Replace coordinates of the start of the children by  the end of the parents. AMAPstudio format.
  data$radius_cyl=0 # Radius of the cylinder
  data[-1,"radius_cyl"]=data[data$parent_ID,"radius"] # The radius of the cylinder is the radius of the start circle i.e the parent. Except when there is a ramification, see below
  data[data$PositionInBranch_ID==1,"radius_cyl"]=data[data$PositionInBranch_ID==1,"radius"] # When a ramification: The radius of the first cylinder of a branch is the radius of the daughter. AMAPstudio format
  # data=data[!branch_ID%in% unique(data[is.na(data$startX)]$branch_ID)]
  ####################
  #######Step 2#######
  ####################

  #######################################################################################################
  #### a. Correct some topological errors/incoherence in terms of topology coming from treeQSM output####
  #######################################################################################################


  # Correcting a treeQSM output error: some cylinder have no extension ID (i.e 0) but they exist as parents ID => incoherence
  sub_table_cyl_error=data[parent_ID%in%data[extension_ID==0][cyl_ID%in%data$parent_ID]$cyl_ID] # a table with the last cylinder of branches that are not connected to their parents (i.e the parent extension is 0)
  if(length(unique(data$branch_ID))==1){
    segment_ID=rep(max(data$cyl_ID),nrow(data))
    Node_ID=rep(0,nrow(data))
    data=cbind(data,segment_ID,Node_ID)
    message("This data is only a trunk: no ramification")
    data$volume=(pi*(data$radius_cyl)^2)*data$length
    data$node_ID=0
    data=data[,c("startX","startY","startZ", "endX","endY","endZ","cyl_ID","parent_ID","extension_ID","radius_cyl","length","volume","branch_ID","segment_ID","node_ID","BranchOrder")]
    names(data)[c(13,16)]=c("axis_ID","branching_order")
    out=list(QSM=data,model=model)
    return(out)
  } # In case of only one branch = a trunk.
  if(nrow(sub_table_cyl_error)!=0){

    sub_table_cyl_error_no_dup=sub_table_cyl_error[duplicated(sub_table_cyl_error$parent_ID)==FALSE] # same table but only one child branch is selected when two or more children (the first one)
    data[extension_ID==0][cyl_ID%in%data$parent_ID]$extension_ID=sub_table_cyl_error_no_dup$cyl_ID # correction of the problem by assigning the good extension ID to the parent

    branch_ID_to_replace=sub_table_cyl_error_no_dup$branch_ID # This vector is the branch ID that should not exist as they follow a branch without any ramification. This second problem comes directly from the same treeQSM incoherence deteted and corrected above
    branch_ID_who_replace=data[cyl_ID%in%sub_table_cyl_error_no_dup$parent_ID]$branch_ID # this is the vector with the good branch ID that should replace the previous one. the loop below make this replacement.
    for (i in 1:length(branch_ID_to_replace)){


      if(all(is.na(match(unique(data$branch_ID),branch_ID_who_replace[i])))){
        new_branch_ID=data[cyl_ID==data[branch_ID==branch_ID_to_replace[i]]$parent_ID[1]]$branch_ID # This is the branch_ID of the parent. if the branch_ID which is supposed to replace does not exist in the branch ID of the whole QSM, it means that it has already been replaced and thus the branch_ID of the parent ID can be used.
        data[branch_ID==branch_ID_to_replace[i]]$branch_ID=rep(new_branch_ID,nrow(data[branch_ID==branch_ID_to_replace[i]]))
        data[branch_ID==new_branch_ID]$PositionInBranch_ID=seq(1,nrow(data[branch_ID==new_branch_ID]))
        next()
      }
      data[branch_ID==branch_ID_to_replace[i]]$branch_ID=rep(branch_ID_who_replace[i],nrow(data[branch_ID==branch_ID_to_replace[i]]))
      data[branch_ID==branch_ID_who_replace[i]]$PositionInBranch_ID=seq(1,nrow(data[branch_ID==branch_ID_who_replace[i]]))

    }
  }

  #################################################
  #### b. Assign a segment_ID to each cylinder ###
  ################################################

  # The nested loop below assign a segment ID to each cylinder. A segment is a succession of cylinder between two ramification points

  data$segment_ID=0

  for(b in unique(data$branch_ID)){

    tab_segm_lim_b=data[branch_ID==b&cyl_ID%in%data[PositionInBranch_ID==1,]$parent_ID,]

    if(nrow(tab_segm_lim_b)==0){
      vec_segment_ID_i=rep(max(data[branch_ID==b,]$cyl_ID),nrow(data[branch_ID==b,]))
      next
    } # if a branch does not have any ramification
    vec_segment_ID_branche_b=NULL
    for (i in 1:(nrow(tab_segm_lim_b)+1)){
      if(i==1){
        cyl_position_min=min(data[branch_ID==b,]$PositionInBranch_ID)-1
        cyl_position_max=tab_segm_lim_b[i,]$PositionInBranch_ID
        vec_segment_ID_i=rep(tab_segm_lim_b[i,]$cyl_ID,cyl_position_max-cyl_position_min)
      }
      if(i>1&i<nrow(tab_segm_lim_b)+1){
        cyl_position_min=tab_segm_lim_b[i-1,]$PositionInBranch_ID
        cyl_position_max=tab_segm_lim_b[i,]$PositionInBranch_ID
        vec_segment_ID_i=rep(tab_segm_lim_b[i,]$cyl_ID,cyl_position_max-cyl_position_min)
      }
      if(i==nrow(tab_segm_lim_b)+1){
        cyl_position_min=tab_segm_lim_b[i-1,]$PositionInBranch_ID
        cyl_position_max=max(data[branch_ID==b,]$PositionInBranch_ID)
        cyl_segm_max=max(data[branch_ID==b,]$cyl_ID)
        vec_segment_ID_i=rep(cyl_segm_max,cyl_position_max-cyl_position_min)
      }
      if(length(vec_segment_ID_i)==0){break}
      vec_segment_ID_branche_b=c(vec_segment_ID_branche_b,vec_segment_ID_i)
    }

    data[branch_ID==b,]$segment_ID=vec_segment_ID_branche_b

  }


  segment_ID_twigs=plyr::ddply(data[segment_ID==0,],~branch_ID,function(x){cbind(x,segment_ID_2=rep(max(x$cyl_ID),nrow(x)))}) # The nested loop above does not manage the twigs (i.e terminal segments). It is done here and in the two following lines
  segment_ID_twigs=segment_ID_twigs[match(data[segment_ID==0,]$cyl_ID,segment_ID_twigs$cyl_ID),]
  data[segment_ID==0,]$segment_ID=segment_ID_twigs$segment_ID_2

  #################################################
  #### b. Assign a Node_ID to each cylinder #######
  ################################################
  # A node is a group of segment (i.e a mother and her daugthers). The node ID is assign to the daugthers of a same node , the Node_ID number correspond to the ID of the mother, so the couple mother/daugthers is easy to retrace.

  QSM_ID_node=plyr::ddply(data[segment_ID!=min(data$segment_ID)],~segment_ID,function(x){
    ID_node=data[cyl_ID==x$parent_ID[1]]$segment_ID
    x=cbind(x,node_ID=rep(ID_node,nrow(x)))
  })
  QSM_first_node=cbind(data[segment_ID==min(data$segment_ID)],node_ID=rep(0,nrow(data[segment_ID==min(data$segment_ID)])))

  if(nrow(QSM_ID_node)==0){
    data=QSM_first_node
    new("aRchi",data=data)
    return(aRchi)
  } # If the tree has only one node
  data=rbind(QSM_first_node,QSM_ID_node)
  data$volume=(pi*(data$radius_cyl)^2)*data$length
  data=data[,c("startX","startY","startZ", "endX","endY","endZ","cyl_ID","parent_ID","extension_ID","radius_cyl","length","volume","branch_ID","segment_ID","node_ID","BranchOrder")]
  names(data)[c(13,16)]=c("axis_ID","branching_order")
  out=data
  out=list(QSM=out,model=model)

  }else{
    data=data.table::fread(file)
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
    out=list(QSM=out,model=model)
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
    out=list(QSM=out,model=model)

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
      axis_ID = data$branchID,
      branching_order = data$branchOrder
    )

    # adjust the axis_ID so it starts to 1
    out[axis_ID == -1, axis_ID := 0]
    out[, axis_ID := axis_ID + 1]

    warning("Using the SimpleTree axis_ID might lead to problems when using the
          functions of the aRchi package")

    out=list(QSM=out,model=model)
  }
  return(out)
}


