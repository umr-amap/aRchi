#' Generate foliage to a QSM or skeleton with segmented annual shoots.
#'
#' @param aRchi a file of class aRchi.
#' @param aArea numeric. Allometric coefficient (a) for leaf area per annual shoot.
#' @param bArea numeric. Allometric coefficient (b) for leaf area per annual shoot.
#' @param aNl numeric. Allometric coefficient (a) for the number of leaves per annual shoot.
#' @param bNl numeric. Allometric coefficient (b) for the number of leaves per annual shoot.
#' @param aNin numeric. Allometric coefficient (a) for the number of internodes per annual shoot.
#' @param bNin numeric. Allometric coefficient (b) for the number of internodes per annual shoot.
#' @param elev numeric. A single value or a vector of length 2 giving the leaves elevation angles.
#' @param elev_error numeric. A random error added to leaves elevation angle.
#' @param az_err numeric. A random error added to le leaves azimuth.
#' @param phyllo character. The phyllotaxy used to insert the leaves along the annual shoot.
#'               Accepted values: "alternate", "spiral", "opposite", "opposite_decussate".
#' @param simple logical. Bypass the allometry for the number of leaves per annual shoot
#'               to insert the leaves on each internode.
#' @param petiole numeric. The petiole length given as a multiplier of the leaf length.
#'
#' @details Allometries for the leaf area, number of leaves and number of internodes are computed
#'          at the annual shoot level and are of the form X = a*AS_length^b. The leaves elevation angle
#'          is constant if a single value is provided for \code{elev} or is linearly interpolated based
#'          of the leaf relative height within the tree crown if two values are provided.
#'
#' @return The aRchi file now including the reconstructed foliage.
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
#' # add leaves
#' aRchi = aRchi::add_leaves(aRchi)
#'
#' plot(aRchi,leaves=TRUE,bg="white",color="chocolate4")
#' }

setGeneric("add_leaves",
           function(aRchi,aArea = 73.027,bArea = 0.826,aNl = 2.563, bNl = 0.4402,
                    aNin = 2.1337, bNin = 0.3818, elev = c(17,45), elev_error = 5,
                    az_err = 20, phyllo = "opposite_decussate", simple = T, petiole = 0.5){standardGeneric("add_leaves")}
)

#' @rdname add_leaves
#' @export

setMethod("add_leaves",
          signature = "aRchi",
          function(aRchi,aArea = 73.027,bArea = 0.826,aNl = 2.563, bNl = 0.4402,
                   aNin = 2.1337, bNin = 0.3818, elev = c(17,45), elev_error = 5,
                   az_err = 20, phyllo = "opposite_decussate", simple = T, petiole = 0.5){

            # to pass CRAN checks
            AS=axis_ID=Nin=LA=Nl=X=Y=Z=leaf_ID=rot=D_tip=cumL=to_div=split_p=endX=startX=endY=startY=endZ=startZ=
              .=intX=intY=intZ=AS_ID=Z_temp=H=elevation=radius=NULL

            # number of leaves per internode according to phyllotaxy
            if(phyllo == "alternate"){
              N_leaves = 1
            }
            if(phyllo == "spiral"){
              N_leaves = 1
            }
            if(phyllo == "opposite"){
              N_leaves = 2
            }
            if(phyllo == "opposite_decussate"){
              N_leaves = 2
            }

            # keep curent annual shoots
            AS_tab = aRchi@QSM[AS==1]

            # compute annual shoot length
            AS_car = AS_tab[,sum(length)*100,by=axis_ID] # in cm for allometries
            data.table::setnames(AS_car,c("axis_ID","length"))

            # compute the number of internodes, leaves and leaf area per annual shoot
            # based on allometric coefficients (from cm)
            if(simple){
              AS_car[, ':='(
                Nin = round(aNin*length^bNin)
              )]
              AS_car[, ':='(
                Nl = Nin*N_leaves,
                LA = aArea*length^bArea
              )]
            }else{
              AS_car[, ':='(
                Nin = round(aNin*length^bNin),
                Nl = round(aNl*length^bNl),
                LA = aArea*length^bArea
              )]
            }

            # internode length
            AS_car[,':='(
              In_length = (length/Nin)/100, # back in meter for compatibility
              Ind_LA = LA / Nl
            )]

            tot_leaves = sum(AS_car$Nl)

            # scalable leaf model
            leaf_model <- data.table::data.table(X=c(0, 2.828427, 4.242641, 2.828427), Y=c(0, 1.414214, 0, -1.414214), Z=c(0, 0, 0, 0))


            # replicate the leaf model N leaves * N leaflets times
            leaves_tab = leaf_model[rep(c(1:4),tot_leaves),]

            # add leaves ID and AS_ID
            leaves_tab[,':='(
              leaf_ID = rep(c(1:tot_leaves),each = 4),
              AS_ID = rep(AS_car$axis_ID,AS_car$Nl*4),
              LA = rep(AS_car$Ind_LA,AS_car$Nl*4)
            )]

            # scales the leaves
            leaves_tab[,':='(
              X = X*sqrt(2/3*LA)/2,
              Y = Y*sqrt(2/3*LA)/2,
              Z = Z*sqrt(2/3*LA)/2,
              LA = NULL
            )]

            # back to meters
            leaves_tab[,':='(
              X = X/100,
              Y = Y/100,
              Z = Z/100
            )]

            # add petiole length
            leaves_tab[,X := X+max(X)*petiole, by = leaf_ID]

            # rotate leaves according to phyllotaxy
            if(phyllo == "alternate"){
              leaves_tab[, rot := rep(c(0,0,0,0,pi,pi,pi,pi),(tot_leaves))[1:nrow(leaves_tab)]]
            }
            if(phyllo == "spiral"){
              leaves_tab[, rot := leaf_ID * pracma::deg2rad(90)]
            }
            if(phyllo == "opposite"){
              leaves_tab[, rot := rep(c(0,0,0,0,pi,pi,pi,pi),(tot_leaves))[1:nrow(leaves_tab)] ]
            }
            if(phyllo == "opposite_decussate"){
              leaves_tab[, rot := rep(c(0,0,0,0,pi,pi,pi,pi,-pi/2,-pi/2,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/2),(tot_leaves))[1:nrow(leaves_tab)] ]
            }

            #############- compute internode position for each AS

            # identify the internodes distances from tip (for internodes with leaves only)
            N_occ = ceiling(AS_car$Nl / N_leaves) # Number of internodes with leaves
            D_div_tab = data.table::data.table(AS = rep(AS_car$axis_ID,N_occ),D_tip = rep(AS_car$In_length,N_occ))
            D_div_tab[,D_tip := cumsum(D_tip),by=AS]

            # does the number of leaves in each internode differs from N_leaves ?
            Nl_test = !(AS_car$Nl / N_leaves - trunc(AS_car$Nl / N_leaves)) == 0

            # matrix to store internodes length
            AS_out = matrix(ncol=4,nrow=sum(N_occ))

            start = 1 # current empty line in the matrix
            for(i in 1:nrow(AS_car)){
              # keep the target AS
              in_AS = AS_tab[axis_ID == AS_car$axis_ID[i]]

              # compute the distance of cylinder tips to the axis base
              in_AS[,cumL := rev(cumsum(rev(length)))]

              # distances for division
              D_div = D_div_tab[AS == AS_car$axis_ID[i],D_tip]
              D_div = D_div - min(D_div)

              # add a column to store the number of division for each cylinder
              in_AS[,to_div := 0]

              # set index for splitting identification
              in_AS[,index := 1:nrow(in_AS)]

              # find which cylinders to split
              for(l in D_div){
                in_AS[which.max( in_AS$index[in_AS$cumL > l] ),to_div := to_div + 1]
              }

              ## build the splitting table
              # keep only the cylinders to split
              in_AS = in_AS[to_div > 0]

              # duplicate the cylinders Nsplit times
              in_AS = in_AS[rep(c(1:nrow(in_AS)), to_div)]

              # where is the split point located in % of the cylinder length
              in_AS[, split_p := (rev(D_div)-(cumL-length))/length]

              # interpolate coordinates
              in_AS[,':='(
                intX = endX * (1-split_p) + startX * split_p,
                intY = endY * (1-split_p) + startY * split_p,
                intZ = endZ * (1-split_p) + startZ * split_p
              )]

              # add the number of leaves for each internode
              in_AS[,N_leaves := N_leaves]

              if(Nl_test[i]){
                in_AS[N_occ[i],N_leaves := 1]
              }

              AS_out[start:(start+nrow(in_AS)-1),1:4] = as.matrix(in_AS[,.(intX,intY,intZ,N_leaves)])
              start = start+nrow(in_AS)
            }
            # replivate the internodes coordinates
            AS_out = AS_out[rep(1:nrow(AS_out),AS_out[,4]*4),1:3]

            ################- Compute leaves geometry

            # add random azimuth by AS
            leaves_tab[,rot := rot + sample(seq(-pi,pi,0.0001),1),by = AS_ID]

            # add random azimuth by leaf
            leaves_tab[,rot := rot + pracma::deg2rad(sample(seq(-az_err/2,az_err/2,0.0001),1)),by = leaf_ID]

            ## compute elevation angle
            # relative height
            leaves_tab[,Z_temp := AS_out[,3]]
            leaves_tab[,H := (Z_temp-min(Z_temp))/(max(Z_temp)-min(Z_temp))]

            # elevation angle
            leaves_tab[,elevation := pracma::deg2rad(elev[1]*(1-H)+elev[2]*H)]

            # add random error to elevation angle
            leaves_tab[,elevation := elevation + pracma::deg2rad(sample(seq(-elev_error/2,elev_error/2,0.001),1)),by=leaf_ID]

            # perform rotation for each leaf
            temp_leaves = as.matrix(leaves_tab[,.(X,Y,Z)])
            for(i in 1:nrow(temp_leaves)){
              temp_leaves[i,] = rgl::rotate3d(temp_leaves[i,], leaves_tab$elevation[i],0,1,0)
              temp_leaves[i,] = rgl::rotate3d(temp_leaves[i,], leaves_tab$rot[i],0,0,1)
            }

            leaves_tab[,':='(
              X = temp_leaves[,1],
              Y = temp_leaves[,2],
              Z = temp_leaves[,3],
              rot = NULL,
              Z_temp = NULL,
              H = NULL,
              elevation = NULL
            )]


            # translate the leaves to the right position
            leaves_tab[,':='(
              X = X + AS_out[,1],
              Y = Y + AS_out[,2],
              Z = Z + AS_out[,3]
            )]


            ###### Build triangular mesh
            leaves_tab = leaves_tab[order(leaf_ID)]
            leaves=leaves_tab[,.(X,Y,Z)]

            # devide the quads into triangular faces
            nl=nrow(leaves)/4
            list = rep(c(1,2,3,1,3,4),nl)
            list = list + rep(seq(0,(nl-1)*4,4),each=6)
            leaves = leaves[list,]

            # build the mesh
            mesh = rgl::as.mesh3d(x=leaves$X,y=leaves$Y,z=leaves$Z,type = "triangles",smooth = T,merge=T)

            # add infos to the mesh
            mesh$meshColor = "faces"
            mesh$leaves_ID = leaves_tab$leaf_ID[which(c(1:nrow(leaves_tab))%%2==0)]
            mesh$AS_ID = leaves_tab$AS_ID[which(c(1:nrow(leaves_tab))%%2==0)]

            aRchi@leaves = list(mesh = mesh, AS_properties = AS_car[order(axis_ID)])

            aRchi@operations$add_leaves = list(aArea = aArea,bArea = bArea,aNl = aNl, bNl = bNl,
                                               aNin = aNin, bNin = bNin, elev = elev, elev_error = elev_error,
                                               az_err = az_err, phyllo = phyllo, simple = simple, petiole = petiole)

            return(aRchi)
          }
)


