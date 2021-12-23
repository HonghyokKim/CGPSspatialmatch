#' match.nearestcaliper function
#'
#' This internal function matches units by the generalized propensity score computed by cgpsmatch function. This function is used only within cgpsmatch function.
#' @export
#' @examples 
#' match.caliper()
match.nearestcaliper <-function(data,bexp,gpsname,caliper_bw,exp.status=1,index.status=1,replace2=TRUE) {
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
  
  data<-data[order(data[,"index"],decreasing=TRUE),]
  data<-data[order(data[,"strata_matchdist"], rev(data[,"index"])),]
  data<-data %>%group_by(strata_matchdist)%>%mutate(GPSdiff=abs(first(GPS,default=first(GPS))-GPS))
  data<-data.frame(data)
  num_exp<-length(unique(data[,"strata_matchdist"]))
  index.units<-data[data[,"index"] == index.status,]
  matched.units<-data[data[,"index"] != index.status,]
  selected<-matched.units %>% 
    group_by(strata_matchdist) %>% 
    slice(which.min(GPSdiff))
  if(replace2==FALSE) {
    selected.nonreplace<-selected[1,]
    selected.FID<-selected[1,"FID"]
    matched.units<-matched.units[matched.units[,"FID"] %!in% selected.FID,]
    matched.units<-matched.units[matched.units[,"strata_matchdist"] != 1,]
    for (ii in seq(2,num_exp)) {
      ss<-matched.units %>%
        group_by(strata_matchdist) %>%
        slice(which.min(GPSdiff))
      selected.nonreplace<-rbind(selected.nonreplace,ss[1,])
      selected.FID<-ss[1,"FID"]
      matched.units<-matched.units[matched.units[,"FID"] %!in% selected.FID,]
      matched.units<-matched.units[matched.units[,"strata_matchdist"] != ii,]
    }
    selected<-selected.nonreplace
  }
  
  data<-rbind(index.units,selected)
  
  
  check<-table(data[,c("strata_matchdist","index")])
  if(ncol(check) !=2) {
    result<-NULL
  }
  else {
    
    
    data<-data[data[,"strata_matchdist"] %!in% as.numeric(names(check[,1][check[,1]==0])),]
    matched.units.exist<-as.numeric(check[,as.numeric(colnames(check)) != 1])
    matched.units.exist.order<-ifelse(matched.units.exist>0,1,0)
    matched.units.exist<-which(matched.units.exist.order==1)
    result<-data[data[,"strata_matchdist"] %in% matched.units.exist,]
    result<-result[order(result[,"index"],decreasing=TRUE),]
    result<-result[order(result[,"strata_matchdist"],rev(result[,"index"])),]
    
    
    ##IF Caliper...
    #### INDEX 0's weight-> 1/(Number of Index=0)
    
  }
}