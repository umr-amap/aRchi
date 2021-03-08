#' Estimate the unrolled wood surface from a QSM
#'
#' Estimate the unrolled wood surface from a QSM by summing the area of all cylinders at several level of organization.
#'
#' @export
#' @docType methods
#' @rdname WoodSurface
#' @param aRchi a file of class aRchi
#' @param level text. The level at which the wood surface is computed. "Tree", "branching_order" or "Axis".
#' @return a numeric or data.table. The wood surface in m2
#' @include aRchiClass.R
#' @examples
#' # Read an aRchi file with at least a QSM.
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Compute the wood surface at branching order level.
#' WoodSurface(Tree1_aRchi,level="branching_order")
#'
setGeneric("WoodSurface",
           function(aRchi,level="Tree"){standardGeneric("WoodSurface")}
)

setMethod("WoodSurface",
          signature = "aRchi",

          function(aRchi,level){
            axis_ID=NULL


            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(level %in% c("Tree","branching_order","Branch")==F) stop("Incorrect level argument")
            QSM=aRchi@QSM

            # Compute the Wood Surface
            QSM$WoodSurface=2*pi*(QSM$radius_cyl/2)*QSM$length

            if(level=="Tree"){return(WoodSurface=sum(QSM$WoodSurface))}
            if(level=="branching_order"){return(QSM[,.(WoodSurface=sum(WoodSurface)),by=c("branching_order")])}
            if(level=="Axis"){Wood_surface_branch=data.table::setorder(QSM[,.(WoodSurface=sum(WoodSurface)),by=c("axis_ID")],axis_ID)
            return(Wood_surface_branch)
            }
          }
)



