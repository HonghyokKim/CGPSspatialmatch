#' matchdist function
#'
#' This function matches exposed units and unexposed units by a pre-specified buffer distance (Euclidean distance)
#' @param data a dataset object.
#' @param bexp a character string indicating the name of the binary exposure. Use apostrophe like "VariableName"
#' @param long a character string indicating the name of the longitude variable of observation units
#' @param lat a character string indicating the name of the latitude variable of observation units
#' @param exp.status a numeric vector indicating the value indicating exposed units. Defalut=1
#' @param distbuf a numeric vector indicating the buffer distance by which exposed units and unexposed units are matched
#' @param exp.included an indicator of whether exposed units are matched with not only unexposed units but also other exposed units. Defalut is TRUE. If FALSE, exposed units are matched with only unexposed units. See details
#' @param replace an indicator of whether matching is done with replacement. Default is TRUE. If FALSE, matching is done without replacement. If FALSE, note that the output of this function may differ by the order of observation units in the original dataset.
#' @export
#' @examples 
#' matchdist()

matchdist<-function(data,bexp,long,lat,exp.status=1,distbuf=0.1,exp.included=TRUE,replace=TRUE) {
  start.time<-Sys.time()
  NUM_EXP <- sum(data[,bexp]==exp.status)
  data[,paste0(bexp,"2")] <- data[,bexp]
  data[,paste0(bexp,"2")][data[,paste0(bexp,"2")]==exp.status] <-seq(NUM_EXP)
  strata<-vector("list",NUM_EXP)
  data[,"FID"] <-seq(nrow(data))
  not.matched.exp<-NULL
  matched.exp<-NULL
  for (ii in seq(NUM_EXP)) {
    
    if(exp.included==FALSE) {
      match.sample<-rbind(data[data[,paste0(bexp,"2")]==ii,],
                          data[data[,bexp] != exp.status,])
    }
    if(exp.included==TRUE) {
      match.sample<-rbind(data[data[,paste0(bexp,"2")]==ii,],
                          data[data[,paste0(bexp,"2")] != ii,])
    }
    match.sample$nearestdist<-0
    match.sample$matchdist<-0
    match.sample$multiple_selected<-0
    # dist.mat<-as.matrix(dist(cbind(match.sample[,long], match.sample[,lat])))
    # diag(dist.mat)<-NA
    # dist.mat<-dist.mat[1,]
    dist.mat <- match.sample %>% mutate(diff=sqrt( (first(eval(as.name(long)),default=first(eval(as.name(long))))-eval(as.name(long)))^2 
                                                   +(first(eval(as.name(lat)),default=first(eval(as.name(lat))))-eval(as.name(lat)))^2)
    )
    dist.mat<-dist.mat$diff
    dist.mat[1]<-NA
    match_num<-which(dist.mat<distbuf)
    matching.dist<-dist.mat[which(dist.mat<distbuf)]
    matched.nonexp.row<-NULL
    
    if(length(match_num)>0L) {
      nearest_match<-as.integer(rev(names(dist.mat[dist.mat==min(dist.mat,na.rm=T)]))[1])
      match.sample[match.sample[,bexp]==1,"nearestdist"] <-1
      match.sample[nearest_match,"nearestdist"] <-1
      match.sample[match_num,"matchdist"] <- matching.dist
      
      index<-match.sample[match.sample[,paste0(bexp,"2")]==ii,]
      ref<-match.sample[match_num,]
      index[,"index"]<-1
      ref[,"index"]<-0
      strata[[ii]]<-rbind(index,
                          ref)
      strata[[ii]]$strata_matchdist<-ii
      #print(match.sample[match_num,"FID"])
      matched.units.row<-match.sample[match_num,"FID"]
      matched.exp<-rbind(matched.exp,match.sample[match.sample[,paste0(bexp,"2")]==ii,])
    }
    else {
      not.matched.exp<-rbind(not.matched.exp,match.sample[match.sample[,paste0(bexp,"2")]==ii,])
    }
    if(replace==FALSE) {
      if(is.null(matched.units.row) != TRUE) {
        data<-subset(data, FID %!in% matched.units.row)
        if(exp.included==TRUE) {
          selected<-match.sample[match_num,]
          exp.selected <- selected[selected[,bexp]==exp.status,]
          matched.exp <- rbind(matched.exp,exp.selected)
        }
      }
      else {
        data<-data
      }
    }
  }
  if(is.null(not.matched.exp) ==FALSE) {
    not.matched.exp<-not.matched.exp[order(not.matched.exp[,paste0(bexp,"2")]),]
  }
  if(is.null(matched.exp) ==FALSE) {
    matched.exp<-matched.exp[order(matched.exp[,paste0(bexp,"2")]),]
  }
  
  matched<-do.call(rbind,strata)
  end.time<-Sys.time()
  if(is.null(matched)) {
    print(paste0("FAILED: units satisfying the matching criteria could not be found (perhaps too short buffer distance) / Elapsed time: ",round(end.time-start.time,3),
                 " ",attr(end.time-start.time,"units")
    ))
  }
  else {
    matched$ID_matchdist<-seq(nrow(matched))
    matched_exp_count<-length(unique(matched[,"strata_matchdist"]))
    
    if(exp.included==TRUE & replace==FALSE) {
      matched[matched[,bexp]==exp.status,"index"]<-1
      matched_exp_count<-sum(matched[matched[,bexp]==exp.status,"index"])
      cat("Recommend not to use this for GPS matching
            If you use this for GPS matching, GPS of onlythe first row of exposed units in strata where more than one exposed units exist 
          is used for matching criteria
          ")
    }
    
    if(length(which(duplicated(matched$FID)==TRUE))>0L) {
      matched[matched$FID %in% which(duplicated(matched$FID)==TRUE),"multiple_selected"]<-1
    }
    print(paste0("Matching by distance DONE / Total count of exposure: ",NUM_EXP," / Matched exposure count: ", matched_exp_count,
                 "/ Elapsed time: ",round(end.time-start.time,3),
                 " ",attr(end.time-start.time,"units")
    ))
    result<-list(distance=distbuf,replace=replace,total.number.exp=NUM_EXP,total.number.matched.exp=matched_exp_count,
                 matched.exp=matched.exp,not.matched.exp=not.matched.exp,matched.dataset=matched)
    return(result)
  }
}