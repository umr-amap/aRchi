#' Write an aRchi file
#'
#' @docType methods
#' @rdname write_aRchi
#' @description write an aRchi file
#' @param aRchi the object of class aRchi to write
#' @param file The directory where to write the `.aRchi` file.
#' @include aRchiClass.R
#' @seealso [read_aRchi()]
#' @examples
#' # Create an empty object of class aRchi
#'  empty_aRchi=aRchi()
#'  # Write the aRchi object
#' write_aRchi(empty_aRchi,file="my_empty_aRchi_file.aRchi")
#' read_aRchi(file = "my_empty_aRchi_file.aRchi")

setGeneric("write_aRchi",
           function(x,file){standardGeneric("write_aRchi")}
)
setMethod("write_aRchi",
          signature(x="aRchi",file="character"),
          function(x,file){
            saveRDS(x,file=file)
          }
)

