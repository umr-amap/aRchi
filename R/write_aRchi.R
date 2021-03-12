#' Write an aRchi file
#'
#' @docType methods
#' @rdname write_aRchi
#' @description write an aRchi file
#' @param aRchi the object of class aRchi to write
#' @param file The directory where to write the .aRchi file.
#' @include aRchiClass.R
#' @seealso \code{\link{read_aRchi}}

setGeneric("write_aRchi",
           function(aRchi,file){standardGeneric("write_aRchi")}
)

#' @rdname write_aRchi
#' @export

setMethod("write_aRchi",
          signature(aRchi="aRchi",file="character"),
          function(aRchi,file){
            saveRDS(aRchi,file=file)
          }
)

