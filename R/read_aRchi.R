#' Read an aRchi file
#'
#' @docType methods
#' @rdname read_aRchi
#' @description Read an aRchi file
#' @param file The directory to the .aRchi file.
#' @include aRchiClass.R
#' @seealso \code{\link{write_aRchi}}

setGeneric("read_aRchi",
           function(file){standardGeneric("read_aRchi")}
)

#' @rdname read_aRchi
#' @export

setMethod("read_aRchi",
          signature = "character",
          function(file){
            readRDS(file)
          }
)

