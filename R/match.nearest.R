#' match.nearest function
#'
#' This internal function matches units by the generalized propensity score computed by cgpsmatch function. This function is used only within cgpsmatch function.
#' @export
#' @examples 
#' match.nearest()
match.nearest <-function(data,bexp,gpsname,index.status=1,replace2=TRUE) {
  data[,"GPS"]<-data[,gpsname]
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
  result<-rbind(index.units,selected)
  result<-result[order(result[,"strata_matchdist"],rev(result[,"index"]),result[,gpsname]),]
  result[,paste0(gpsname,"_diff")] <-result[,"GPSdiff"]
  result<-result[ , -which(names(result) %in% c("GPS","GPSdiff"))]
  return(result)
  
}