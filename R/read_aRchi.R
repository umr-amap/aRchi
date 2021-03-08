#' Read an aRchi file
#'
#' @docType methods
#' @rdname read_aRchi
#' @description Read an aRchi file
#' @param file The directory to the `.aRchi` file.
#' @include aRchiClass.R
#' @seealso [write_aRchi()]
#' @examples
#' # Create an empty object of class aRchi
#'  empty_aRchi=aRchi()
#'  # Write the aRchi object
#' write_aRchi(empty_aRchi,file="my_empty_aRchi_file.aRchi")
#' read_aRchi(file = "my_empty_aRchi_file.aRchi")
#'
setGeneric("read_aRchi",
           function(file){standardGeneric("read_aRchi")}
)

setMethod("read_aRchi",
          signature = "character",
          function(file){
            readRDS(file)
          }
)

