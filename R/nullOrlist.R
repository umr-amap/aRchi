#' @title nullOrlist
#' @description Class union
#' @export nullOrlist

setClassUnion("nullOrlist", c("NULL", "list"))
