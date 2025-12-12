#' Tree volume estimation of the woody part
#'
#' Compute the tree wood volume at different level of organization
#'
#' @export
#' @docType methods
#' @rdname Treevolume
#' @param aRchi an object of class aRchi with at least a QSM
#' @param level character. The level at which the wood volume is computed. \code{Tree}, \code{branching_order} or \code{Axis}.
#' @return a numeric or data.table. The wood volume in m3 at the requested level
#' @include aRchiClass.R
#' @examples
#' # Read an aRchi file with at least a QSM
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#' # Compute the whole tree wood biomass.
#' Treevolume(Tree1_aRchi)

setGeneric("Treevolume",
           function(aRchi,level="Tree"){standardGeneric("Treevolume")}
)

#' @rdname Treevolume
#' @export
#'
setMethod("Treevolume",
          signature = "aRchi",

          function(aRchi,level){
            volume=axis_ID=.=NULL


            if(inherits(aRchi,"aRchi")==F) stop("The provided data is not of class aRchi")
            if(is.null(aRchi@QSM)) stop("The archi file does not contains a QSM")
            if(level %in% c("Tree","branching_order","Axis")==FALSE) stop("Incorrect level argument")

            QSM=aRchi@QSM



            if(level=="Tree"){return(volume=sum(QSM$volume))}
            if(level=="branching_order"){return(QSM[,.(volume=sum(volume)),by=c("branching_order")])}
            if(level=="Axis"){volume_branch=data.table::setorder(QSM[,.(volume=sum(volume)),by=c("axis_ID")],axis_ID)
            return(volume_branch)
            }
          }
)


