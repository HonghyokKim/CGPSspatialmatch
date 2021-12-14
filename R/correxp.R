#' correxp function
#'
#' This function calcualtes standardized mean difference in variable(s) in a matched dataset.
#' @param data a dataset object.
#' @param bexp a character string indicating the name of the binary exposure. Use apostrophe like "VariableName"
#' @param cexp a character string indicating the name of the continuous exposure. Use apostrophe like "VariableName"
#' @param exp.status a numeric vector indicating the value indicating exposed units. Defalut=1
#' @param varinames a vector of variable names for which you wish to calculate standardized mean difference. List variable names as a vector like c("VariableA","VariableB")
#' @param weightname the name of the weight variable. This weight is calculated by cgps.match. See cgps.match function. If NULL, weight is not given (i.e., all observations are equally weighted)
#' @param method a character string indicating which correlation coefficient is to be computed. These include "Pearson" (default), "Spearman", "Polychoric", or "Polyserial". For tetrachoric use "Polychoric" and for biserial use "Polyserial". This relies on wCorr::weightedCorr
#' @param aftermatch a character string indicating whether this function is used for a matched sample. Default="no"
#' @export
#' @examples correxp(test.dat,"exposure","exp_cons",exp.status=1,varinames=c("X1","X2","X3","U"),weightname="weight",method="Pearson")
#' correxp()
correxp<-function(data,bexp,cexp,exp.status=1,varinames=NULL,weightname=NULL,method="Pearson",aftermatch="no") {
  
  if(aftermatch=="no"){
  exposed <- data[data[,bexp] == exp.status,]
  }
  if(aftermatch=="yes"){
  exposed <- data
  }
  
  if( nrow(exposed)>0) {
  if(is.null(weightname)==TRUE) {
    weightname<-"W"
    exposed$W<-1
  }
  result<-mapply(FUN=wCorr::weightedCorr,x=data.frame(exposed[,cexp]),y=data.frame(exposed[,varinames]),weights=data.frame(exposed[,weightname]),method=method)
  }
  else {
    result<-rep(NA,length(varinames))
  }
  names(result)<-varinames
  return(result)
}