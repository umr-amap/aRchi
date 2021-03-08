#' Compute the path fraction
#'
#' Compute from an object of class aRchi the path fraction
#' @export
#' @docType methods
#' @rdname PathFraction
#' @param aRchi an object of class aRchi with at least the QSM and the Paths table
#' @return The path fraction
#' @details The path fraction is the ratio between the mean path length and the maximum path length.
#' @references
#'
#' Martin-Ducup, O. et al. Terrestrial laser scanning reveals convergence of tree architecture with increasingly dominant crown canopy position. Functional Ecology (2020).
#'
#' Smith, D. D. et al. Deviation from symmetrically self-similar branching in trees predicts altered hydraulics, mechanics, light interception and metabolic scaling. New Phytologist 201, 217–229 (2014).
#'
#' @seealso [Make_Path())] to compute the paths table.
#' @include aRchiClass.R
#' @examples
#' # Read an aRchi file with at least the QSM and the paths table
#' file=system.file("extdata","Tree_1_aRchi.aRchi",package = "aRchi")
#' Tree1_aRchi=read_aRchi(file)
#'
#' PathFraction(Tree1_aRchi)
#'
setGeneric("PathFraction",
           function(aRchi){standardGeneric("PathFraction")}
)

setMethod("PathFraction",
          signature = "aRchi",
          function(aRchi){
            ID_Path=NULL
            Paths=aRchi@Paths

            if(class(aRchi) != "aRchi") stop("The provided data is not of class aRchi")
            if(is.null(aRchi@Paths)) stop("The aRchi file does not contains Paths")

            PathLength=Paths[,sum(length),by=ID_Path]
            Pf=mean(PathLength$V1)/max(PathLength$V1)
            return(Pf)
          }
)

