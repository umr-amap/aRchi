
#' Annual shoot segmentation in tree skeleton
#'
#' @description Segment the annual shoots in a tree skeleton based on the detection
#'              of the branching patterns created by acrotony.
#'
#' @param aRchi an object of class aRchi containing at least a QSM.
#' @param tree_age numeric, optional. The tree age. Helps to achieve more robust segmentation.
#' @param segment_reiterations list of numeric values. The parameters to segment traumatic
#'                             reiterations based on their age difference with the bearer
#'                             annual shoot and their elevation. The list must have the
#'                             following form: list(age_difference, elevation_angle) where
#'                             age_difference and elevation_angle are numeric vectors.
#'                             NOTE the elevation angle is defined relative to the zenith.
#'
#' @return The input aRchi file with an additional field in the QSM slot being the segmented
#'         annual shoots. NOTE that annual shoot = 1 correspond to the last growing season.
#'         If traumatic reiteration segmentation was achieved, an additional field
#'         labeling cylinders that belong to a reiteration is added.
#'
#' @references Lecigne, B., Delagrange, S., & Taugourdeau, O. (2021). Annual Shoot Segmentation and
#'             Physiological Age Classification from TLS Data in Trees with Acrotonic Growth. Forests, 12(4), 391.
#'             https://doi.org/10.3390/f12040391
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
#' plot(aRchi,color="annual_shoots",bg = "white")
#' }

setGeneric("segment_annual_shoots",
           function(aRchi,tree_age,segment_reiterations){standardGeneric("segment_annual_shoots")}
)
#' @rdname segment_annual_shoots
#' @export

setMethod("segment_annual_shoots",
          signature = "aRchi",
          function(aRchi,tree_age,segment_reiterations){

            # to pass CRAN check
            AS=AS_segmentation_flag=L2=N_child=V1=axis_ID=child_length=cluster=
              cumL=cutree=cyl_ID=dist=dist_tip=div=endX=endY=endZ=hclust=is_first=
              is_new=is_reiteration=newX=newY=newZ=new_parent=orig_ID=parent_ID=
              sorted=split_p=start=startX=startY=startZ=to_div=update_length=NULL

            if(missing(segment_reiterations)){
              reit = FALSE
              reit_param = NA
            }else{
              reit = TRUE
              reit_param = segment_reiterations
            }

            # the skeleton
            data = aRchi@QSM

            # compute the cylinder distance to axis tip for each axis
            data[,dist_tip := rev(cumsum(length)),by=axis_ID]

            # number of child axes
            data[,N_child := length(which(data$parent_ID == cyl_ID))-1,by=cyl_ID]
            data[N_child == -1, N_child := 0]

            # add the a field for reiteration segmentation if needed
            if(reit) data[, is_reiteration := 1]

            ##### compute L2 metric
            # length of child axes of each axis
            data[,child_length := sum(
              data$dist_tip[                                                  # distance from tip (i.e. total length of the axis)
                which(data$parent_ID %in% data$cyl_ID[data$axis_ID == axis_ID]# of the first child cylinders
                      & data$axis_ID != axis_ID)                              # that are not from the parent axis
              ]), by = axis_ID]

            ax_len = data[,unique(child_length),by=axis_ID] # summary child length by axis

            # L2 metric is the sum of the length of child axis and of its child axes
            data[,L2:=0]
            data[N_child >= 1,L2 := (sum(
              data$length[which( data$axis_ID %in% data$axis_ID[which(data$parent_ID == cyl_ID & data$axis_ID != axis_ID)] ) ] # child axes
            ) + sum(
              ax_len$V1[which( ax_len$axis_ID %in% data$axis_ID[which(data$parent_ID == cyl_ID & data$axis_ID != axis_ID)] ) ] # child of childs
            )) / N_child
            , by = cyl_ID ]


            #################################################
            #####- first segmentation of annual shoots -#####
            #################################################

            axes = unique(data$axis_ID)

            for(i in axes){
              axe = data[axis_ID == i]
              # number of unique length, i.e. the number of potential clusters
              ndiff = length(unique(axe$L2))

              # if more than 1 length -> clustering approcah
              if(ndiff > 1){
                ### cylinder clustering based on L2
                # distance matrix for clustering
                dist_mat = stats::dist((axe$L2), method = 'euclidean')
                # clustering
                cl = stats::hclust(dist_mat,"ward.D")


                # if there are 3 or more potential clusters -> estimate the optimal number of cluster
                #(i.e. all non branched segments are in the same group)
                if(ndiff >= 3){
                  for(k in 3:ndiff){
                    # cut the tree with k clusters
                    clusters = stats::cutree(cl, k = k)
                    # check if all non branches cylinders are in the same cluster
                    dat_cl = data.table::data.table(N_child = axe$N_child,clusters = clusters)
                    test = dat_cl[,sum(N_child),by = clusters]

                    # if all non branched are in the same cluster -> this is the optimal number of clusters
                    if(min(test$V1) == 0){break}
                  }

                  # attribute clusters so that they are sorted by L2
                  dat_cl[,L2 := axe$L2]

                  # table with sorted clusters
                  mean_cl = dat_cl[,mean(L2),by=clusters]
                  mean_cl = mean_cl[order(V1)]
                  mean_cl[,sorted := 1:max(clusters)]

                  # match new clusters to old clusters
                  dat_cl[,clusters := mean_cl$sorted[mean_cl$clusters == clusters], by = clusters]
                  axe[,cluster := dat_cl$clusters]

                  # if two groups greater than 1 are close from each other -> group them to avoid surnumemary segmentation
                  for(i in 2:(nrow(axe)-1)){
                    if(axe$cluster[i-1] == axe$cluster[i+1] & axe$cluster[i-1] != 1) axe$cluster[c(i,i+1)] = axe$cluster[i-1]
                  }
                }else{
                  # if only two different length
                  axe[,cluster:=1]
                  axe[N_child > 0,cluster:=2]
                }
              }else{
                # only one cluster
                axe[,cluster := 1]
              }

              if(ndiff > 1){
                ### if there are multiple lengths use the evolutive threshold to segment AS

                # the threshold
                th = 2

                # current annual shoot
                cur_AS = 1

                # the cylinder at the tip of the axis is of year 1
                axe[nrow(axe), AS := 1]
                for(j in (nrow(axe)-1):1){ # from the tip to the base

                  # th increases each time the cluster of the current cylinder in greater than th
                  if(axe$cluster[j] > th)th = th + 1

                  # criteria to segment a new annual shoot -> cluster is greater than th
                  # and than cluster of the previous cylinder
                  if(axe$cluster[j] > axe$cluster[j+1] & axe$cluster[j] >= th) cur_AS = cur_AS+1

                  # add the curent AS to the axis
                  axe[j,AS := cur_AS]
                }
                # add axis ASs to the data
                data[axis_ID == unique(axe$axis_ID), AS := axe$AS]
              }else{
                ### if there is only one length, there is no AS to segment
                data[axis_ID == unique(axe$axis_ID), AS := 1]
              }
            }

            # remove useless tables
            rm(axe)
            rm(dat_cl)
            rm(mean_cl)

            #################################################
            #####- annual shoot segmentation correction -#####
            #################################################
            A1 = TRUE # is it the trunk ? To avoid logical test in loop

            # done by increasing branching order
            for(i in sort(unique(data$branching_order))){
              #done for each axis
              for(j in sort(unique(data$axis_ID[data$branching_order == i]))){

                # select one axis
                axe = data[axis_ID == j]
                axe[,to_div := 0] # will further store the number of divisions required for each cylynder

                # parent annual shoot age
                if(A1){
                  # if this is the trunk -> parent AS age is the tree age
                  if(missing(tree_age)){
                    # if tree age not defined -> no corrections
                    parent_age = max(axe$AS)
                    ta = NA
                  }else{
                    parent_age = tree_age + 1
                    ta = tree_age
                  }
                  A1 = FALSE
                }else{
                  # if not the trunk -> find the age of the bearer annual shoot
                  parent_age = data$AS[which(data$cyl_ID %in% axe$parent_ID & data$axis_ID != j)]
                }

                # difference between actual age and expected age:
                # if = or > 1: the axis is too young
                # if = 0: the axis is well segmented
                # if = or < -1: the axis is too old
                age_diff = parent_age - max(axe$AS) - 1

                ############- reiteration segmentation -> if this is a reiteration -> no correction
                if(reit & j > 1){
                  # if the age difference and the elevation match the criteria or that the parent is a reiteration
                  # -> this is a reiteration
                  elevation = VoxR::axis_angle(data.frame(
                    axe$endX[nrow(axe)] - axe$startX[1],
                    axe$endY[nrow(axe)] - axe$startY[1],
                    axe$endZ[nrow(axe)] - axe$startZ[1]
                  ),axis = "Z")

                  if(
                    (age_diff >= reit_param[[1]] & elevation <= reit_param[[2]])
                    | data$is_reiteration[which(data$cyl_ID %in% axe$parent_ID & data$axis_ID != j)] == 2
                  ){
                    data[axis_ID == j, is_reiteration := 2]
                    next() # go to the next axis
                  }
                }

                # if too young -> segment longest annual shoot
                if(age_diff >= 1){

                  if(length(unique(axe$AS))>1){
                    ## if the axis was already segmented

                    # length of each AS
                    AS_length = axe[,sum(length),by=AS]
                    data.table::setnames(AS_length,old = "V1", new = "length")
                    AS_length[,':='(update_length = length,div = 1)]

                    # find how many time each AS must be devided
                    for(k in 1:age_diff){
                      # divide the longest
                      AS_length[which.max(update_length),div := div+1]

                      # update length
                      AS_length[,update_length := length/div]
                    }

                    AS_length = AS_length[div > 1]
                    # identify which segments to divide

                    first = TRUE
                    for(k in AS_length$AS){
                      # keep the target AS
                      in_AS = axe[AS == k]

                      # compute the distance of cylinder tips to the axis bae
                      in_AS[,cumL := cumsum(length)]

                      # identify the distances to add divisions
                      D_div = AS_length$update_length[AS_length$AS == k] * c(1:(AS_length$div[AS_length$AS == k]-1))

                      # add a column to store the number of division for each cylinder
                      in_AS[,to_div := 0]

                      # find which cylinders to split
                      for(l in D_div){
                        if(l > min(in_AS$cumL)){
                          in_AS[which( abs(in_AS$cumL[in_AS$cumL > l] - l) == min(abs(in_AS$cumL[in_AS$cumL > l] - l)))+max(which(in_AS$cumL < l)),to_div := to_div + 1]
                        }else{
                          in_AS[1,to_div := to_div + 1]
                        }
                      }


                      ## build the splitting table
                      # keep only the cylinders to split
                      in_AS = in_AS[to_div > 0]

                      # duplicate the cylinders Nsplit times
                      in_AS = in_AS[rep(c(1:nrow(in_AS)), to_div)]

                      # where is the split point located in % of the cylinder length
                      in_AS[, split_p := (D_div-(cumL-length))/length]

                      # produce the split table containing each cylinder to split and the indication about where to split
                      if(first){
                        split_tab = in_AS
                        first = FALSE
                      }else{
                        split_tab = data.table::rbindlist(list(split_tab,in_AS))
                      }
                    }
                  }else{
                    ## if the axis was not previously segmented
                    # compute the distance of cylinder tips to the axis bae
                    axe[,cumL := cumsum(length)]

                    # the split distances are equally distributed along the axis
                    D_div = c(1:age_diff) * (sum(axe$length)/(age_diff+1))

                    # find which cylinders to split
                    for(l in D_div){
                      if(l > min(axe$cumL)){
                        axe[which( abs(axe$cumL[axe$cumL > l] - l) == min(abs(axe$cumL[axe$cumL > l] - l))) + max(which(axe$cumL < l)),to_div := to_div + 1]
                      }else{
                        axe[1,to_div := to_div + 1]
                      }
                    }

                    ## build the spliting tsable
                    # keep only the cylinders to split
                    axe = axe[to_div > 0]

                    # duplicate the cylinders Nsplit times
                    axe = axe[rep(c(1:nrow(axe)), to_div)]

                    # where is the split point located in % of the cylinder length
                    split_tab = axe[, split_p := (D_div-(cumL-length))/length]
                  }

                  ## split cylinders : interpolate coordinates based on the split table infos

                  # add a cylinder to the tip of each splitting section
                  # duplicate all
                  split_tab = split_tab[rep(1:nrow(split_tab),each=2)]
                  # identify duplicates
                  split_tab[,is_new := rep(c(1,2),nrow(split_tab)/2)]
                  # add the ID of the next cylinder
                  split_tab[,test := data.table::shift(cyl_ID,type = "lead")]
                  # if a cylinder is new but does not have the same ID -> it's a tip
                  split_tab[is.na(test), is_new := 1] # the last one is always a tip
                  split_tab = split_tab[!(test == cyl_ID & is_new == 2)]

                  # the cylinder at the tip keeps the same coordinates
                  split_tab[test != cyl_ID | is.na(test), split_p := 1]

                  # interpolate coordinates
                  split_tab[,':='(
                    endX = endX * split_p + startX * (1-split_p),
                    endY = endY * split_p + startY * (1-split_p),
                    endZ = endZ * split_p + startZ * (1-split_p)
                  )]

                  # for the cylinders that are not the first -> the start coordinates are the end of the parent
                  split_tab[,':='(
                    newX = data.table::shift(endX),
                    newY = data.table::shift(endY),
                    newZ = data.table::shift(endZ),
                    start = 1:nrow(split_tab)
                  )]
                  # identify the first cylinder of each split section to not change its coordinates
                  split_tab[,is_first := start == min(start),by=cyl_ID]
                  # change coordinates
                  split_tab[is_first ==  FALSE,':='(
                    startX = newX,
                    startY = newY,
                    startZ = newZ
                  )]

                  # keep the original cylinder ID for matching
                  split_tab[,orig_ID := cyl_ID]

                  # update the table before adding to the original table
                  # the last cylinder of the split section keeps the same ID to avoid problems in the hierarchy
                  split_tab[cyl_ID == test, cyl_ID := 1:nrow(split_tab[cyl_ID == test])+max(data$cyl_ID)]
                  # the first cylinder of the split section keeps the same parent ID
                  split_tab[,new_parent := data.table::shift(cyl_ID)]
                  split_tab[is_first == FALSE, parent_ID := new_parent]


                  # add the new segments to the original table
                  # remove useless columns
                  split_tab[,':='(
                    to_div = NULL,
                    cumL = NULL,
                    split_p = NULL,
                    is_new = NULL,
                    test = NULL,
                    newX = NULL,
                    newY = NULL,
                    newZ = NULL,
                    start = NULL,
                    is_first = NULL,
                    new_parent = NULL
                  )]

                  for(k in sort(unique(split_tab$orig_ID),decreasing = TRUE)){
                    before = data[cyl_ID < k] # all before the target cylinder
                    after = data[cyl_ID > k] # all after target cylinder
                    to_add = split_tab[orig_ID == k] # traget cylinder

                    # update tables
                    # those in the same axis but before are one year older than in original segmentation
                    before[axis_ID == j & cyl_ID < k, AS := AS + (nrow(to_add)-1)]
                    # the first cylinder is one year older than the second
                    to_add[(1:nrow(to_add))-1, AS := AS+((nrow(to_add)-1):1)]
                    # also update the splitting table
                    split_tab[, AS := AS+(nrow(to_add)-1)]

                    # remove before merging
                    to_add[, orig_ID := NULL]

                    # insert the new segments in the original data
                    data = data.table::rbindlist(list(before,to_add,after))
                  }
                }
                # if too old -> merge shortest annual shoot
                if(age_diff <= -1){

                  # remove a useless column
                  axe[,to_div := NULL]

                  for(k in 1:abs(age_diff)){
                    # compute AS length in the axis
                    AS_length = axe[,sum(length),by=AS]

                    # merge the shortest to the previous except if it's the first
                    if(which.min(AS_length$V1) == 1){
                      axe[AS == AS_length$AS[1], AS := AS - 1]
                    }else{
                      axe[AS >= AS_length$AS[which.min(AS_length$V1)-1], AS := AS - 1]
                    }
                  }
                  # replace the curent axis
                  data[axis_ID == j, AS := axe$AS]
                }
              }
              data[,AS_segmentation_flag := 0]
              data[AS < 1, ':='(AS = 1,AS_segmentation_flag = 1)]
            }

            data[,length := sqrt( (startX-endX)^2 + (startY-endY)^2 + (startZ-endZ)^2 )]

            #remove useless variables
            data[,':='(
              dist_tip = NULL,
              N_child = NULL,
              child_length = NULL,
              L2 = NULL
            )]

            aRchi@QSM = data

            aRchi@operations$segment_annual_shoots = list(tree_age = ta, segment_reiterations = reit_param)

            return(aRchi)
          }
)
