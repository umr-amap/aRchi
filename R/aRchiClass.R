#' @name aRchi
#' @title aRchi
#' @description Class containing files to compute and display in three dimensions tree architectural metrics at different level of organization
#' @docType class
#' @slot QSM a data.table containing QSM information according to read_QSM function format
#' @slot pointcloud a data.table containing the point cloud used to generated the QSM
#' @slot Paths a data.table of Paths according to Make_Path function (see [Make_Path()])
#' @slot Nodes Metrics computed at the node scale (see [Make_node()])
#' @slot operations Record all the operations realized on the object.

setClassUnion("nullOrDatatable", c("NULL", "data.table"))
setClassUnion("nullOrLASOrDatatable", c("NULL", "LAS", "data.table"))
setClassUnion("nullOrlist", c("NULL", "list"))
setClassUnion("nullOrnumeric", c("NULL", "numeric"))
aRchi=setClass("aRchi",slots=c(QSM = "nullOrDatatable",
                               Paths="nullOrDatatable",
                               pointcloud = "nullOrLASOrDatatable",
                               foliage = "nullOrDatatable",
                               Nodes = "nullOrlist",
                               operations = "nullOrlist"
                               ))


setMethod("show",
          "aRchi",
          function(object) {
            if(!is.null(object@QSM)){
              cat("QSM:",nrow(object@QSM),"cylinders,",length(unique(object@QSM$segment_ID)), "segments and", length(unique(object@QSM$node_ID)-1), "node\n")
            }else{
              cat("QSM is empty\n")
            }
            if(!is.null(object@Paths)){
              cat("Paths:", length(unique(object@Paths$ID_Path)),"Paths \n")
            }else{
              cat("Paths is empty\n")
            }
            if(!is.null(object@pointcloud)){
              if(is.data.frame(object@pointcloud)){cat("Point cloud:", nrow(object@pointcloud),"points \n")}
              if(is.data.frame(object@pointcloud)==F){cat("Point cloud:", nrow(object@pointcloud@data),"points \n")}
            }else{
              cat("Point cloud is empty\n")
            }
            if(!is.null(object@Nodes)){
              cat("List with two elements: Relative and Absolute nodes tables\n")
            }else{
              cat("Node table is empty\n")
            }
            if(!is.null(object@operations)){
              cat("Operations: ", paste(names(object@operations),collapse="; "))
            }else{
              cat("No operation realized\n")
            }
          }
)



