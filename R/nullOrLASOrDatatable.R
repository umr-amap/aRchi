#' @title nullOrLASOrDatatable
#' @description Class union
#' @export nullOrLASOrDatatable

setClassUnion("nullOrLASOrDatatable", c("NULL", "LAS", "data.table"))
