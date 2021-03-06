#' smd function
#'
#' This function calcualtes standardized mean difference in variable(s) in a matched dataset.
#' @param data a dataset object.
#' @param bexp a character string indicating the name of binary exposure. Use apostrophe like "VariableName"
#' @param exp.status a numeric vector indicating the value indicating exposed units. Defalut=1
#' @param varinames a character vector indicating variable names for which you wish to compute standardized mean difference. List variable names as a vector like c("VariableA","VariableB")
#' @param weightname the name of the weight variable. This weight is calculated by cgps.match. See cgps.match function. If NULL, weight is not given (i.e., all observations are equally weighted)
#' @export
#' @examples smd(test.dat,"exposure",exp.status=1,varinames=c("X1","X2","X3","U"))
#' smd()
smd<-function(data,bexp,exp.status=1,varinames=NULL, weightname=NULL) {
  exposed <- data[data[,bexp] == exp.status,]
  unexposed <- data[data[,bexp] != exp.status,]
  
  if( nrow(exposed) * nrow(unexposed)>0) {
    if(is.null(weightname)==TRUE) {
  exposed$W<-1
  unexposed$W<-1
  exposed.w.m<-sapply(data.frame(exposed[,varinames]),FUN=w.m,w=exposed$W)
  unexposed.w.m<-sapply(data.frame(unexposed[,varinames]),FUN=w.m,w=unexposed$W)
  exposed.w.sd<-sapply(data.frame(exposed[,varinames]),FUN=w.sd,w=exposed$W)
  unexposed.w.sd<-sapply(data.frame(unexposed[,varinames]),FUN=w.sd,w=unexposed$W)
  result<-(exposed.w.m-unexposed.w.m)/sqrt((exposed.w.sd^2+unexposed.w.sd^2)/2)
    }
    
    if(is.null(weightname)==FALSE) {
      exposed$W<-exposed[,weightname]
      unexposed$W<-unexposed[,weightname]
      exposed.w.m<-sapply(data.frame(exposed[,varinames]),FUN=w.m,w=exposed$W)
      unexposed.w.m<-sapply(data.frame(unexposed[,varinames]),FUN=w.m,w=unexposed$W)
      exposed.w.sd<-sapply(data.frame(exposed[,varinames]),FUN=w.sd,w=exposed$W)
      unexposed.w.sd<-sapply(data.frame(unexposed[,varinames]),FUN=w.sd,w=unexposed$W)
      result<-(exposed.w.m-unexposed.w.m)/sqrt((exposed.w.sd^2+unexposed.w.sd^2)/2)
    }
    
  }
  else {
    result<-rep(NA,length(varinames))
  }
  names(result)<-varinames
  return(result)
}