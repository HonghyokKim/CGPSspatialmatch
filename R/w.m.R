#' w.m. function
#'
#' This internal function calculates weighted average.
#' @param x is a variable to be calculated; w is a weight.
#' w.m()
w.m <- function(x, w) sum(x*w)/sum(w)
