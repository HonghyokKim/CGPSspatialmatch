#' %!in% function
#'
#' This internal function is a reverse %in%
#' @export
#' @examples 
#' %!in%()

'%!in%' <- function(ppp,qqq)!('%in%'(ppp,qqq))