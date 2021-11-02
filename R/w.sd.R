#' w.sd. function
#'
#' This internal function calculates weighted standard deviation
#' @param x is a variable to be calculated; w is a weight.
#' w.sd()
w.sd <- function(x, w) sqrt(sum(((x - w.m(x,w))^2)*w)/(sum(w)*(length(x)-1)/length(x)))
