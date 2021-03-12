#' Compute Moment of force
#'
#' Compute the moment of gravity force Mf and the moment of gravity force relative to cylinder radius Mf_r from an object of class aRchi.
#'
#' @export
#' @docType methods
#' @rdname Compute_Mf
#' @param aRchi an object of class aRchi with at least the QSM and the paths table.
#' @param WoodDensity a numeric or a data.table. A single wood density value for the whole tree or one value per cylinder in kg/m3. If wood density is given for each cylinder a data.table with two column (i.e cyl_ID and WoodDensity) must be given.
#' @return The aRchi file with the QSM slot having three new columns: the biomass upstream the cylinder \code{sub_tree_biomass}, the moment of gravity force \code{Mf} and the moment of gravity force relative to cylinder radius \code{Mf_r}.
#' @details
#'
#' The moment of gravity force (i.e \code{Mf}) is calculated at each cylinder position. \code{Mf} can be seen as a proxy of the mechanical loading history due to gravity at a given position of a tree.
#' This quantity is defined by the following the equation: Mf=R*F where R is the lever arm, which is the norm of the horizontal vector between the position where \code{Mf} is measured (i.e a cylinder) and the position where the force is applied (i.e., the center of mass, G, of the whole structure upstream a cylinder: a subtree).
#'
#' The mass of the cylinders are needed to calculate the center of mass and are estimated using their volume and the wood density provided in argument \code{WoodDensity}. Finally, F is the weight of the subtree: F=g*M	with g the standard acceleration due to gravity (9.81 m.s-²).
#'
#' The moment of gravity force relative (i.e \code{Mf_r}) to cylinder radius r is also computed following the formula: Mf_r = Mf/r^3
#'
#'
#' @include aRchiClass.R
#'
#' @examples
#' # Read an aRchi file with a QSM and paths tables.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' # Compute the moment of force for each cylinder
#' Tree1_aRchi=Compute_Mf(Tree1_aRchi,WoodDensity=550)
#'
#' # show the QSM data.table with the three new columns sub_tree_biomass, MF and Mf_r)
#' get_QSM(Tree1_aRchi)
#'
setGeneric("Compute_Mf",
           function(aRchi,WoodDensity){standardGeneric("Compute_Mf")}
)

#' @rdname Compute_Mf
#' @export

setMethod("Compute_Mf",
          signature = "aRchi",

          function(aRchi,WoodDensity=NULL){
            cyl_ID=ID_Path=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(WoodDensity)) stop("Incorrect WooDensity argument. Please fill the WoodDensity argument with a single wood density value or with a data.table with two column (i.e cyl_ID and WoodDensity).")
            if(is.numeric(WoodDensity)==F&is.data.frame(WoodDensity)==F) stop("Incorrect WoodDensity argument. Please fill the WoodDensity argument with a single wood density value or with a data.table with two column (i.e cyl_ID and WoodDensity).")
            if(is.null(aRchi@Paths)) stop("The archi file does not contains Paths")

            QSM=aRchi@QSM
            Paths=aRchi@Paths

            if(is.numeric(WoodDensity)){
              # Wood density kg.m^3
              QSM$WD=WoodDensity
              QSM$Biomass=QSM$volume*QSM$WD
            }
            if(is.data.frame(WoodDensity)){
              QSM=merge(QSM,WoodDensity,by="cyl_ID")
              # Wood density  to kg.m^3
              QSM$WD=QSM$WoodDensity
              QSM$Biomass=QSM$volume*QSM$WD
            }


            # Compute sub_tree_biomass, Mf and sigma for each cylinder
            MC_tab= plyr::ddply(QSM,("cyl_ID"),function(x){
              cylinder_of_interest=QSM[cyl_ID==x$cyl_ID]
              coord_cylinder_of_interest=cylinder_of_interest[,c("startX","startY","startZ")]
              sub_tree=QSM[cyl_ID%in%unique(Paths[ID_Path%in%Paths[cyl_ID==x$cyl_ID]$ID_Path&cyl_ID>=x$cyl_ID]$cyl_ID)]
              sub_tree$mean_x=(sub_tree$startX+sub_tree$endX)/2
              sub_tree$mean_Y=(sub_tree$startY+sub_tree$endY)/2
              sub_tree$mean_Z=(sub_tree$startZ+sub_tree$endZ)/2

              coord_CM=c(
                x=sum(sub_tree$mean_x*sub_tree$Biomass)/sum(sub_tree$Biomass),
                y=sum(sub_tree$mean_Y*sub_tree$Biomass)/sum(sub_tree$Biomass),
                z=sum(sub_tree$mean_Z*sub_tree$Biomass)/sum(sub_tree$Biomass))

              # Lever arm
              H=sqrt((coord_cylinder_of_interest$startX-coord_CM[1])^2+(coord_cylinder_of_interest$startY-coord_CM[2])^2)

              sub_tree_biomass=sum(sub_tree$Biomass)
              # Weight of the sub-tree
              P=9.81*sub_tree_biomass

              # Moment of force of the sub-tree
              Mf=H*P

              # Sigma
              Mf_r = Mf / cylinder_of_interest$radius_cyl^3

              MC_tab=data.table::data.table(cbind(cyl_ID=x$cyl_ID,sub_tree_biomass,Mf,Mf_r))
              return(MC_tab)
            },.progress = "text")

            # Merge QSM and MC table
            QSM=merge(QSM,MC_tab,by="cyl_ID")
            QSM=QSM[,-("WD")]
            aRchi@QSM=QSM
            aRchi@operations$"Compute_Mf" = c(WoodDensity=WoodDensity)
            return(aRchi)
          }
)

