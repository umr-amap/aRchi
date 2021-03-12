#' @title nullOrDatatable
#' @description Class union
#' @export nullOrDatatable

setClassUnion("nullOrDatatable", c("NULL", "data.table"))
