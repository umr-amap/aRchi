#' Tree biomass estimation of the woody part
#'
#' Compute the tree wood biomass at different level of organization
#'
#' @export
#' @docType methods
#' @rdname TreeBiomass
#' @param aRchi an object of class aRchi with at least a QSM
#' @param WoodDensity a numeric or a data.table. A single wood density value for the whole tree or one value per cylinder in kg/m3. If wood density is given for each cylinder a data.table with two column (i.e cyl_ID and WoodDensity) must be given.
#' @param level character. The level at which the wood biomass is computed. \code{Tree}, \code{branching_order} or \code{Axis}.
#' @return a numeric or data.table. The wood biomass in Kg at the requested level
#' @include aRchiClass.R
#' @examples
#' # Read an aRchi file with at least a QSM
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Compute the whole tree wood biomass.
#' TreeBiomass(Tree1_aRchi,WoodDensity=550)

setGeneric("TreeBiomass",
           function(aRchi,WoodDensity,level="Tree"){standardGeneric("TreeBiomass")}
)

#' @rdname TreeBiomass
#' @export
#'
setMethod("TreeBiomass",
          signature = "aRchi",

          function(aRchi,WoodDensity=NULL,level){
            Biomass=axis_ID=.=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(is.null(WoodDensity)) stop("Incorrect WooDensity argument. Please fill the WoodDensity argument with a single wood density value or with a data.table with two column (i.e cyl_ID and WoodDensity).")
            if(is.numeric(WoodDensity)==FALSE&is.data.frame(WoodDensity)==FALSE) stop("Incorrect WoodDensity argument. Please provide a WoodDensity argument with a single wood density value or with a data.table with two column (i.e cyl_ID and WoodDensity).")
            if(level %in% c("Tree","branching_order","Axis")==FALSE) stop("Incorrect level argument")

            QSM=aRchi@QSM
            if(is.numeric(WoodDensity)){
              QSM$WD=WoodDensity
              QSM$Biomass=QSM$volume*QSM$WD
            }
            if(is.data.frame(WoodDensity)){
              QSM=merge(QSM,WoodDensity,by="cyl_ID")
              QSM$WD=QSM$WoodDensity
              QSM$Biomass=QSM$volume*QSM$WD
            }
            if(level=="Tree"){return(Biomass=sum(QSM$Biomass))}
            if(level=="branching_order"){return(QSM[,.(Biomass=sum(Biomass)),by=c("branching_order")])}
            if(level=="Axis"){Biomass_branch=data.table::setorder(QSM[,.(Biomass=sum(Biomass)),by=c("axis_ID")],axis_ID)
            return(Biomass_branch)
            }
          }
)


