#' match.caliper function
#'
#' This internal function is
#' @export
#' @examples 
#' match.caliper()
match.caliper <-function(data,bexp,gpsname,caliper_bw,exp.status=1,index.status=1,replace2=TRUE) {
  data[,"GPS"]<-data[,gpsname]
  num_exp<-length(unique(data[,"strata_matchdist"]))
  strata<-vector("list",num_exp)
  PS_caliper <- caliper_bw*sd(data[,"GPS"],na.rm=T)
  
  exposed<-data[data[,"index"]==index.status,c("strata_matchdist","GPS")]
  
  exposed[,"GPS_range1"] <- exposed[,"GPS"]-PS_caliper
  exposed[,"GPS_range2"] <- exposed[,"GPS"]+PS_caliper
  exposed<-exposed[ , -which(names(exposed) %in% c("GPS"))]
  data<-merge(data,exposed,by="strata_matchdist",all.x=T)
  
  data[,"selected"] <- ifelse(data[,"GPS"]>=data[,c("GPS_range1")],1,0) & ifelse(data[,"GPS"]<=data[,c("GPS_range2")],1,0)
  data<-data[data[,"selected"]==TRUE,]
  
  check<-table(data[,"strata_matchdist"],data[,"index"])
  if(ncol(check) !=2) {
    result<-NULL
  }
  else {
    matched.units.exist<-as.numeric(check[,as.numeric(colnames(check)) != 1])
    matched.units.exist.order<-ifelse(matched.units.exist>0,1,0)
    matched.units.exist<-which(matched.units.exist.order==1)
    result<-data[data[,"strata_matchdist"] %in% matched.units.exist,]
    result<-result[order(result[,"index"],decreasing=TRUE),]
    result<-result[order(result[,"strata_matchdist"],rev(result[,"index"])),]
  }
}