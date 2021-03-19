#' Make Node
#'
#' Compute several node metrics (i.e Leonardo's ratio and WBE parameters) from an object aRchi at different distance from the node.
#'
#' @export
#' @docType methods
#' @rdname Make_Node
#' @param aRchi a file of class aRchi containing at least a QSM
#' @param all_combination logical (default = \code{FALSE}). Should the node metrics be computed at each combination of distance from the node ? (see details).
#' @return The aRchi file with the Nodes slot filled.
#' @details The Nodes slot contains a list of two data.table (Absolute_positions and Relative_positions). Each data.tables contains for each node several values of Leonardo's rule ratio (R_ratio), radius scaling exponent (alpha), length scaling exponent (beta) and the estimated metabolic rate (theta) which are parameters of the West Brown and Enquist metabolic theory (see Bentley et al. 2013, Lau et al 2019 and Martin-Ducup et al. 2020).
#'
#' For a given node, each R_ratio and alpha correspond to its value estimated with the branch radius (for alpha) or the cross section area (for R_ratio) at a given position from the node for the parent (pos_parent) and for the daughters (pos_daugthers). The positions (i.e pos_parent and pos_daughters) are the distances in meters (for Absolute_position data.table) or in percentage (for Relative_position data.table) from the node (i.e the ramification point) to the point of radius or cross section area estimation. These positions are given at a 10 cm step and cannot be higher than the length of the shorter segment among the couple daughters/parent for Absolute_position and at a 10\% step for the Relative_position.
#'
#' Kriging models are used to estimate radius (and thus cross section area) along the segment positions. For example pos_parent = 0.2 and pos_daugthers = 0.2 for Absolute_position data.table means that the parameters (i.e R_ratio and alpha) have been estimated at 20cm from the node position for the parent and 20 cm from the node position for the daughters. For the Relative_position data.table, this position is given in proportion to the total length of the segment. For example, pos_parent = 0.5 and pos_daughter = 0.5 means that the parameters are estimated at mid_length of the segments for both the parent and the daughters).
#'
#' beta values are repeated for each node as it depend on segment length only and not on radius position.
#'
#' If \code{all_combination = TRUE} all the possible combination for pos_parent and pos_daughter are computed (e.g pos_parent = 0.3 and pos_daughter = 0.5) but note that this processing might take several minutes and is \code{FALSE} by default.
#'
#' @references
#'
#' 	Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' 	Lau, A. et al. Estimating architecture-based metabolic scaling exponents of tropical trees using terrestrial LiDAR and 3D modelling. Forest Ecology and Management 439, 132–145 (2019).
#'
#' 	Bentley, L. P. et al. An empirical assessment of tree branching networks and implications for plant allometric scaling models. Ecology Letters 16, 1069–1078 (2013).
#'
#' @seealso \code{\link{WBEparameters}} to estimates WBE parameters at different level; \code{\link{LeonardoRatio}} to estimates Leonardo Da Vinci's ratio  at different level.
#' @include aRchiClass.R
#' @examples
#' \dontrun{
#' # Read a QSM file
#' file=system.file("extdata","Tree_1_TreeQSM.txt",package = "aRchi")
#' QSM=read_QSM(file,model="treeQSM")
#' # Build an object of class aRchi
#' Tree1_aRchi<-build_aRchi(QSM=QSM)
#' Tree1_aRchi
#' # Make the node table
#' Tree1_aRchi<-Make_Node(Tree1_aRchi)
#' Tree1_aRchi
#' # WBE parameters at the tree level
#' LeonardoRatio(Tree1_aRchi)
#' }
#'
setGeneric("Make_Node",
           function(aRchi,all_combination=F){standardGeneric("Make_Node")}
)

#' @rdname Make_Node
#' @export


setMethod("Make_Node",
          signature = "aRchi",

          function(aRchi,all_combination){
            segment_ID=cyl_ID=node_ID=position=pos_node=alpha=extension_ID=.=NULL

            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            QSM=aRchi@QSM
            if(all_combination){
              message("Computing all the possible combination of position (i.e all_combination=TRUE) might take several minutes.")
            }


            node_table=NULL
            node_table_rel=NULL
          pb <- progress::progress_bar$new(total =length( unique(QSM$node_ID)),width = 60,
                                           format = "Kriging model per segment [:bar] :percent",clear=FALSE)
            for (i in unique(QSM$node_ID)) {
              # skip the first node_ID as it is the trunk and it does not have a mother
              pb$tick()
              if(i==0){next}

              # Parent of node i
              sub_parent=QSM[segment_ID==QSM[cyl_ID==i]$segment_ID]

              # Daughters of node i
              daugthers=QSM[node_ID ==i]

              # Start the computation of node metrics at different position from the node (every 10cm from the insertion for the absolute table and every 10\% fr the relative table)
              Rs_levels_i=NULL
              Krig_daughters=NULL
              Krig_daughters_rel=NULL
              for (j in unique(daugthers$segment_ID)) {
                sub_br=daugthers[segment_ID==j]

                # Construct a table for the kriging model. A kriging model is used to estimate the diameter at any position in a segment.
                if(unique(sub_br$axis_ID)==unique(sub_parent$axis_ID)){
                  tab_Krig=data.table::data.table(cbind(x=cumsum(c(0,sub_br$length)),y=c(sub_br$radius_cyl,sub_br$radius_cyl[nrow(sub_br)]-0.01)))
                }else{
                  tab_Krig=data.table::data.table(cbind(x=cumsum(c(0,sub_br$length)),y=c(sub_br$radius_cyl[1]+0.01,sub_br$radius_cyl)))
                }

                if(nrow(tab_Krig)==2){

                  tab_Krig=rbind(tab_Krig[1],cbind(x=mean(tab_Krig$x),y=mean(tab_Krig$y)),tab_Krig[2])

                }


                model=NULL
                it=0
                while (class(model)[1]!="km") {
                  it=it+1
                  if(it>10000) stop("kriging model for diameter estimation at absolute position did not converge.")
                  lower=0.20
                  if(all(tab_Krig$x<0.20)){lower=0}
                   sink(paste0(tempdir(),"\\sink-examp.txt"))
                  try(model <- pkgcond::suppress_messages(DiceKriging::km(formula=~x,design=data.frame(x=tab_Krig$x), response=data.frame(y=tab_Krig$y),covtype="matern5_2",lower=lower)),silent=T)
                   sink()
                }

                # Absolute distance from the node
                ListX=seq(min(tab_Krig$x),max(tab_Krig$x),0.1) # Absolute => seq 0.10 m
                ListY <- stats::predict(model, newdata=data.frame(x=ListX), type="UK")$mean #
                Krig_sub_br=data.table::data.table(cbind(x=ListX,y=as.vector(ListY)))
                Krig_br=data.table::data.table(cbind(segment_ID=j,position=Krig_sub_br$x,radius=Krig_sub_br$y,length_tot=sum(sub_br$length)))
                Krig_daughters=data.table::data.table(rbind(Krig_daughters,Krig_br))

                # Relative distance from the node
                ListX_rel=seq(0,1,0.1)*max(tab_Krig$x) # Absolute => seq 10%
                ListY_rel <- stats::predict(model, newdata=data.frame(x=ListX_rel), type="UK")$mean #
                Krig_sub_br_rel=data.table::data.table(cbind(x=ListX_rel,y=as.vector(ListY_rel)))
                Krig_br_rel=data.table::data.table(cbind(segment_ID=j,position=seq(0,1,0.1),radius=Krig_sub_br_rel$y,length_tot=sum(sub_br$length)))
                Krig_daughters_rel=data.table::data.table(rbind(Krig_daughters_rel,Krig_br_rel))
              }



              # Estimate the daugthers area at the different absolute position
              Krig_daughters$radius=abs(Krig_daughters$radius)
              Krig_daughters$area=pi*((Krig_daughters$radius)^2)

              # Estimate the daugthers area at the different relative position
              Krig_daughters_rel$radius=abs(Krig_daughters_rel$radius)
              Krig_daughters_rel$area=pi*((Krig_daughters_rel$radius)^2)

              # Branch number
              Krig_daughters$N_branch=length(unique(Krig_daughters$segment_ID))
              Krig_daughters_rel$N_branch=length(unique(Krig_daughters_rel$segment_ID))

              # Start the computation of area at different position (every 10cm from the insertion) for the parent (same approach as above for daughters).
              tab_Krig_par=data.table::data.table(cbind(x=cumsum(c(0,sub_parent$length)),y=c(sub_parent$radius_cyl[1]+0.01,sub_parent$radius_cyl)))

              if(nrow(tab_Krig_par)==2){
                tab_Krig_par=rbind(tab_Krig_par[1],cbind(x=mean(tab_Krig_par$x),y=mean(tab_Krig_par$y)+0.1),tab_Krig_par[2])
              }


              model=NULL
              it=1
              while (class(model)[1]!="km") {
                it=it+1
                if(it>10000) stop("kriging model for diameter estimation at regular position did not converge.")
                lower=0.20
                if(all(tab_Krig_par$x<0.20)){lower=0}
                 sink(paste0(tempdir(),"\\sink-examp.txt"))
                 try(model<-pkgcond::suppress_messages(DiceKriging::km(formula=~x,design=data.frame(x=tab_Krig_par$x), response=data.frame(y=tab_Krig_par$y),covtype="matern5_2",lower=lower)),silent=T)
                 sink()
              }


              ListX=seq(min(tab_Krig_par$x),max(tab_Krig_par$x),0.1)
              ListY <- stats::predict(model, newdata=data.frame(x=ListX), type="UK")$mean
              Krig_parent_tab=data.table::data.table(cbind(x=ListX,y=as.vector(ListY)))
              Krig_parent=data.table::data.table(cbind(Segment_ID=i,position=Krig_parent_tab$x,radius=Krig_parent_tab$y,length_tot=sum(sub_parent$length)))


              # Relative distance from the node
              ListX_rel=seq(0,1,0.1)*max(tab_Krig_par$x) # Absolute => seq 10%
              ListY_rel <- stats::predict(model, newdata=data.frame(x=ListX_rel), type="UK")$mean #
              Krig_parent_tab_rel=data.table::data.table(cbind(x=ListX_rel,y=as.vector(ListY_rel)))
              Krig_parent_rel=data.table::data.table(cbind(Segment_ID=i,position=seq(0,1,0.1),radius=Krig_parent_tab_rel$y,length_tot=sum(sub_parent$length)))



              # Estimate parent area at the different positions
              Krig_parent$radius=abs(Krig_parent$radius)
              Krig_parent$area=pi*((Krig_parent$radius)^2)

              Krig_parent_rel$radius=abs(Krig_parent_rel$radius)
              Krig_parent_rel$area=pi*((Krig_parent_rel$radius)^2)

              # Reverse the position => distance from the insertion (the node)
              Krig_parent$pos_node=round(abs(Krig_parent$position-max(Krig_parent$position)),2)

              Krig_parent_rel$pos_node=round(abs(Krig_parent_rel$position-max(Krig_parent_rel$position)),2)

              # Maximum distance to the node to compute the ratio
              pos_lim=min(c(Krig_daughters[,max(position),.(segment_ID)]$V1,max(Krig_parent$pos_node)))
              # Apply the maximum distance to daughters
              Krig_daughters=Krig_daughters[position<=pos_lim]

              Krig_daughters$position=round(Krig_daughters$position,2)
              Krig_daughters_rel$position=round(Krig_daughters_rel$position,2)
              # Sum of daughters area at every position (until the maximum position)
              Krig_daughters_R_ratio=plyr::ddply(Krig_daughters,("position"),function(x){sum_area=sum(x$area)})
              names(Krig_daughters_R_ratio)=c("position","sum_area")
              Krig_daughters_R_ratio=data.table::data.table(Krig_daughters_R_ratio)

              Krig_daughters_R_ratio_rel=plyr::ddply(Krig_daughters_rel,("position"),function(x){sum_area=sum(x$area)})
              names(Krig_daughters_R_ratio_rel)=c("position","sum_area")
              Krig_daughters_R_ratio_rel=data.table::data.table(Krig_daughters_R_ratio_rel)

              # Apply the maximum distance to the parent
              Krig_parent=Krig_parent[pos_node<=pos_lim]




              # Combination of all possible position (i.e distance to the node) for daugthers and parent
              if(all_combination){
              combinaison_position=data.table::data.table(expand.grid(Krig_daughters_R_ratio$position,Krig_daughters_R_ratio$position))
              names(combinaison_position)=c("pos_parent","pos_daugthers")

              combinaison_position_rel=round(data.table::data.table(expand.grid(Krig_daughters_R_ratio_rel$position,Krig_daughters_R_ratio_rel$position)),2)
              names(combinaison_position_rel)=c("pos_parent","pos_daugthers")
              }

              if(all_combination==FALSE){
              combinaison_position_rel=data.table::data.table(cbind(pos_parent=seq(0,1,0.1),pos_daugthers=seq(0,1,0.1)))
              combinaison_position=data.table::data.table(cbind(pos_parent=Krig_daughters_R_ratio$position,pos_daugthers=Krig_daughters_R_ratio$position))
              }

              # Compute the Leonardo's ratio for all these positions
              #absolute
              Rs_levels_i=plyr::adply(combinaison_position,1,function(x){cbind(Krig_parent[pos_node==x$pos_parent,3],R_ratio=Krig_daughters_R_ratio[position==x$pos_daugthers]$sum_area/Krig_parent[pos_node==x$pos_parent]$area)})
              names(Rs_levels_i)[3]="radius_parent"

              #Relative
              Rs_levels_i_rel=plyr::adply(combinaison_position_rel,1,function(x){cbind(Krig_parent_rel[pos_node==x$pos_parent,3],R_ratio=Krig_daughters_R_ratio_rel[position==x$pos_daugthers]$sum_area/Krig_parent_rel[pos_node==x$pos_parent]$area)})
              names(Rs_levels_i_rel)[3]="radius_parent"

              # Compute the WBE parameters for all these positions
              #absolute
              WBE_param_i=plyr::adply(combinaison_position,1,function(x){cbind(Krig_parent[pos_node==x$pos_parent,3],Krig_daughters[position==x$pos_daugthers,1],
                                                                               alpha=-(log(Krig_daughters[position==x$pos_daugthers]$radius/Krig_parent[pos_node==x$pos_parent]$radius)/log(Krig_daughters[position==x$pos_daugthers]$N_branch)),
                                                                               beta=-(log(Krig_daughters[position==x$pos_daugthers]$length_tot/Krig_parent[pos_node==x$pos_parent]$length_tot)/log(Krig_daughters[position==x$pos_daugthers]$N_branch)))})

              WBE_param_i=WBE_param_i[,.(alpha=median(alpha),beta=median(beta)),by=c("pos_parent","pos_daugthers")]
              WBE_param_i$theta=1/((2*WBE_param_i$alpha)+WBE_param_i$beta)

              #Relative

              WBE_param_i_rel=plyr::adply(combinaison_position_rel,1,function(x){cbind(Krig_parent_rel[pos_node==x$pos_parent,3],Krig_daughters_rel[position==x$pos_daugthers,1],
                                                                                       alpha=-(log(Krig_daughters_rel[position==x$pos_daugthers]$radius/Krig_parent_rel[pos_node==x$pos_parent]$radius)/log(Krig_daughters_rel[position==x$pos_daugthers]$N_branch)),
                                                                                       beta=-(log(Krig_daughters_rel[position==x$pos_daugthers]$length_tot/Krig_parent_rel[pos_node==x$pos_parent]$length_tot)/log(Krig_daughters_rel[position==x$pos_daugthers]$N_branch)))})

              WBE_param_i_rel=WBE_param_i_rel[,.(alpha=median(alpha),beta=median(beta)),by=c("pos_parent","pos_daugthers")]
              WBE_param_i_rel$theta=1/((2*WBE_param_i_rel$alpha)+WBE_param_i_rel$beta)


              node_table_i=cbind(node_ID=i,Rs_levels_i,WBE_param_i[,3:5],Daughters_number=length(unique(daugthers$segment_ID)))
              node_table=data.table::data.table(rbind(node_table,node_table_i))


              node_table_i_rel=cbind(node_ID=i,Rs_levels_i_rel,WBE_param_i_rel[,3:5],Daughters_number=length(unique(daugthers$segment_ID)))
              node_table_rel=data.table::data.table(rbind(node_table_rel,node_table_i_rel))
            }

            Nodes=list(node_table,node_table_rel)
            names(Nodes)=c("Absolute_positions","Relative_positions")
            aRchi@Nodes=Nodes
            aRchi@operations$Make_Node=c(all_combination=all_combination)
            return(aRchi)
            close(pb)
          }
)




