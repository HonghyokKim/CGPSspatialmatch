#' cgpsmatch function
#'
#' This function computes generalized propensity score for exposed and unexposed units using the product of the propensity score and the conditional generalized propensity score and/or conducts matching by this generalized propensity score.
#' @param data a dataset object.
#' @param bexp a character string indicating the name of the binary exposure variable. Use apostrophe like "VariableName"
#' @param cexp a character string indicating the name of the continuous exposure variable. Use apostrophe like "VariableName"
#' @param ps a character string indicating the name of propensity score. 
#' @param model.exponly a regression model object that is a generalized propensity score regression model fitted using only exposed units.
#' @param exp.status the value indicating exposed units. Default is 1
#' @param method a character string indicating the matching method used to conduct matching. Default is NULL. If NULL, only generalized propensity score is computed and matching is not done. Options include "nearest" (nearest neighbor matching) and "caliper" (caliper matching)
#' @param caliper_bw a numeric vector indicating caliper bandwidth. Default is 0.1. If method is "nearest", this parameter is ignored.
#' @param replace an indicator of whether matching is done with replacement. Default is TRUE. If FALSE, matching is done without replacement. If FALSE, note that the output of this function may differ by the order of observation units in the original dataset.
#' @param weight.cutoff a numeric vector indicating the cutoff value of the weight. Default is 10.
#' @export
#' @examples 
#' cgpsmatch()

cgpsmatch<-function(data,bexp,cexp,ps,model.exponly,expstatus=1,method=NULL,caliper_bw=0.1,replace=TRUE,weight.cutoff=10) {
  exponly <- data[data[,bexp]==expstatus,]
  unexponly <- data[data[,bexp] != expstatus,]
  
  if(sum(grepl("xgboost",model.exponly$call))==1) {
    exponly.pred <- data.matrix(exponly[,c(model.exponly$feature_names)])
    unexponly.pred <- data.matrix(unexponly[,c(model.exponly$feature_names)])
  }
  if(sum(grepl("xgboost",model.exponly$call))==0) {
    exponly.pred <- exponly
    unexponly.pred <- unexponly
  }
  
  GPS_exponly <- dnorm(exponly[,cexp],mean=predict(model.exponly,newdata=exponly.pred),sd=sqrt(mean( (predict(model.exponly,newdata=exponly.pred)-exponly[,cexp])^2 )))
  GPS_exponly_Cstab <- dnorm(exponly[,cexp],mean=mean(exponly[,cexp],na.rm=T),sd=sd(exponly[,cexp],na.rm=T))
  
 #cf_unexposed <- rnorm(nrow(unexponly.pred),mean=predict(model.exponly,newdata=unexponly.pred),sd=sqrt(mean( (predict(model.exponly,newdata=exponly.pred)-exponly[,cexp])^2)) )
  cf_unexposed <- predict(model.exponly,newdata=unexponly.pred)+rnorm(nrow(unexponly.pred),mean=0,sd=sqrt(mean( (predict(model.exponly,newdata=exponly.pred)-exponly[,cexp])^2 )))
  GPS_unexponly <- dnorm(cf_unexposed,mean=predict(model.exponly,newdata=unexponly.pred),sd=sqrt(mean( (predict(model.exponly,newdata=exponly.pred)-exponly[,cexp])^2 )))
  GPS_unexponly_Cstab <- dnorm(cf_unexposed,mean=mean(cf_unexposed,na.rm=T),sd=sd(cf_unexposed,na.rm=T))
  
  exponly[,paste0(cexp,"_cf")]<- exponly[,cexp]
  unexponly[,paste0(cexp,"_cf")]<- cf_unexposed
  
  gpsname<-paste0(ps,"_GPS")
  exponly[,gpsname]<-exponly[,ps]*GPS_exponly
  unexponly[,gpsname]<-unexponly[,ps]*GPS_unexponly
  
  GPS_exponly_Bstab <- sum(exponly[,ps])/nrow(exponly)
  GPS_unexponly_Bstab <- sum(unexponly[,ps])/nrow(unexponly)
  
  exponly$weight<-  (1/exponly[,gpsname])*GPS_exponly_Cstab*GPS_exponly_Bstab
  unexponly$weight<-  (1/unexponly[,gpsname])*GPS_unexponly_Cstab*GPS_unexponly_Bstab
  
  exponly$weight <- ifelse(exponly$weight>weight.cutoff,weight.cutoff,exponly$weight)
  unexponly$weight <- ifelse(unexponly$weight>weight.cutoff,weight.cutoff,unexponly$weight)
  
  result<-rbind(exponly,unexponly)
  if(is.null(method)==FALSE) {
    if(method=="nearest") {
      result<-match.nearest(result,bexp,gpsname,replace2=replace)
    }
    if(method=="caliper") {
      result<-match.caliper(result,bexp,gpsname,caliper_bw)
    }
    if(method=="optimal") {
      result<-match.optimal(result,bexp,gpsname,caliper_bw)
    }
  }
  if(replace==FALSE) {
    print("if replacement is FALSE, NA may come out because unexposed units ran out")
  }
  if(is.null(result)==TRUE) {
    print("Matching is failed maybe because unexposed units satisfying the matching criteria do not exist")
    
  }
  if(is.null(result)==FALSE) {
    print(paste0("Matching is sucessfully done: / matched exposed units: ", nrow(result[result[,"index"]==1,])," out of ",nrow(data[data[,"index"]==1,])))
    
  }
  return(result)
}