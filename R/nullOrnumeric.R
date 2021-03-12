#' @title nullOrnumeric
#' @description Class union
#' @export nullOrnumeric

setClassUnion("nullOrnumeric", c("NULL", "numeric"))
