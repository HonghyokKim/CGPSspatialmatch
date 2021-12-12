#' estimateATT function
#'
#' This function matches exposed units and unexposed units by a pre-specified buffer distance (Euclidean distance)
#' @param dataset a dataset object.
#' @param bexp a character string indicating the name of the binary exposure. Use apostrophe like "VariableName"
#' @param exp.status a numeric vector indicating the value indicating exposed units. Defalut=1
#' @param cexp a character string indicating the name of the continuous exposure. Use apostrophe like "VariableName"
#' @param fmethod.replace an indicator of whether matching by distance is done with replacement. Default=TRUE. If FALSE, matching is done without replacement. If FALSE, note that the output of this function may differ by the order of observation units in the original dataset.
#' @param distbuf a numeric vector indicating the buffer distance by which exposed units and unexposed units are matched
#' @param exp.included an indicator of whether exposed units are matched with not only unexposed units but also other exposed units. Defalut is TRUE. If FALSE, exposed units are matched with only unexposed units. See details
#' @param long a character string indicating the name of the longitude variable of observation units
#' @param lat a character string indicating the name of the latitude variable of observation units
#' @param PS.method a character string or a vector of variable names, indicating the method of propensity score estimation. Options include "mgcv.GAM" (Generalized additive model in mgcv package), "xgboost" (extreme gradient boosting in xgboost package), and "xgboost.cv" (xgboost with cross-validation). Default="mgcv.GAM"
#' @param PS.formula a character string indicating the formula of propensity score estimation. For PSmethod="mgcv.GAM", this string must be like "the binary exposure variable ~ variableA+variableB+variableC+s(long,lat)". For PSmethod="xgboost", this must be a vector like c("variableA", "variableB", "variableC", "long", "lat")
#' @param PS.max_depth (xgboost only) a numeric vector indicating maximum depth of a tree. Default=5
#' @param PS.eta (xgboost only) a numeric vector indicating the learning rate. Default=0.1
#' @param PS.nthread (xgboost only) a numeric vector indicating the number of thread. Default=1. See xgb.train
#' @param PS.eval_metric (xgboost only) a character string indicating evaluation metrics for validation data. Default="auc". See xgb.train
#' @param PS.objective (xgboost only) a character string indicating the objective function. Default="binary:logistic". See xgb.train
#' @param PS.nrounds (xgboost only) a numeric vector indicating the number of rounds. Default=100. See xgb.train.
#' @param PS.cv.nround (xgboost.cv only) a numeric vector indicating the maximum number of rounds.
#' @param PS.cv.nfold (xgboost.cv only) a numeric vector indicating N-fold cross-validation.
#' @param PS.cv.objective (xgboost.cv only) a character string indicating the objective function. Default="binary:logistic" See xgb.train
#' @param PS.cv.max_depth (xgboost.cv only) a numeric vector indicating maximum depth of a tree.
#' @param PS.cv.eta (xgboost.cv only) a numeric vector indicating the learning rate. 
#' @param PS.cv.nthread (xgboost.cv only) a numeric vector indicating the number of thread. See xgb.train
#' @param PS.cv.subsample (xgboost.cv only) a numeric vector indicating subsample ratio of the training instance. Default=0.5. See xgb.train
#' @param PS.cv.eval_metric (xgboost.cv only) a character string indicating evaluation metrics for validation data. Default="auc". See xgb.train
#' @param PS.cv.colsample_bytree (xgboost.cv only) a numeric vector indicating subsample ratio of columns when constructing each tree. Default=1. See. xgb.train
#' @param PS.cv.min_child_weight (xgboost.cv only) Default=1. See. xgb.train
#' @param PS.cv.min.burnin (xgboost.cv only) a numeric vector indicating the minimum number of rounds. Default=99, which indicates the minimum number of rounds is 100. This must be smaller than PS.early_stopping_rounds.
#' @param PS.early_stopping_rounds (xgboost.cv only) a numeric vector indicating when xgboost stops. If the evaluation metric did not decrease until when (code)PS.early_stopping_rounds, xgboost stops. This saves time. Default=100.
#' @param PS.cv.local.N (xgboost.cv only) a numeric vector indicating how many times xgboost will search (local) mini-ma of the evaluation metric function. Default=10
#' @param CGPS.method a character string or a vector of variable names, indicating the method of conditional propensity score estimation. Options include "mgcv.GAM", "xgboost", and "xgboost.cv". Default="mgcv.GAM".
#' @param CGPS.formula a character string indicating the formula of conditional generalized propensity score estimation
#' @param CGPS.max_depth (xgboost only) a numeric vector indicating maximum depth of a tree. Default=5
#' @param CGPS.eta (xgboost only) a numeric vector indicating the learning rate. Default=0.1
#' @param CGPS.nthread (xgboost only) a numeric vector indicating the number of thread. Default=1. See xgb.train
#' @param CGPS.eval_metric (xgboost only) a character string indicating evaluation metrics for validation data. Default="rmse". See xgb.train
#' @param CGPS.objective (xgboost only) a character string indicating the objective function. Default="reg:squarederror". See xgb.train
#' @param CGPS.nrounds (xgboost only) a numeric vector indicating the number of rounds. Default=100. See xgb.train.
#' @param CGPS.cv.nround (xgboost.cv only) a numeric vector indicating the maximum number of rounds.
#' @param CGPS.cv.nfold (xgboost.cv only) a numeric vector indicating N-fold cross-validation.
#' @param CGPS.cv.objective (xgboost.cv only) a character string indicating the objective function. Default="binary:logistic" See xgb.train
#' @param CGPS.cv.max_depth (xgboost.cv only) a numeric vector indicating maximum depth of a tree.
#' @param CGPS.cv.eta (xgboost.cv only) a numeric vector indicating the learning rate.
#' @param CGPS.cv.nthread (xgboost.cv only) a numeric vector indicating the number of thread. Default=1. See xgb.train
#' @param CGPS.cv.subsample (xgboost.cv only) a numeric vector indicating subsample ratio of the training instance. See xgb.train
#' @param CGPS.cv.eval_metric (xgboost.cv only) a character string indicating evaluation metrics for validation data. Default="auc". See xgb.train
#' @param CGPS.cv.colsample_bytree (xgboost.cv only) a numeric vector indicating subsample ratio of columns when constructing each tree. Default=1. See. xgb.train
#' @param CGPS.cv.min_child_weight (xgboost.cv only) Default=1. See. xgb.train
#' @param CGPS.cv.min.burnin (xgboost.cv only) a numeric vector indicating the minimum number of rounds. Default=99, which indicates the minimum number of rounds is 100. This must be smaller than PS.early_stopping_rounds.
#' @param CGPS.early_stopping_rounds (xgboost.cv only) a numeric vector indicating when xgboost stoCGPS. If the evaluation metric did not decrease until when (code)CGPS.early_stopping_rounds, xgboost stops This saves time. Default=100.
#' @param CGPS.cv.local.N (xgboost.cv only) a numeric vector indicating how many times xgboost will search (local) mini-ma of the evaluation metric function. Default=10
#' @param smethod method a character string indicating the matching method used to conduct matching by GCGPS. Default="caliper". Options include "nearest" (nearest neighbor matching) and "caliper" (caliper matching)
#' @param caliper_bw a numeric vector indicating caliper bandwidth. Default=0.1. If method is "nearest", this parameter is ignored.
#' @param smethod.replace an indicator of whether matching by GPS is done with replacement. Default=TRUE. If FALSE, matching is done without replacement. If FALSE, note that the output of this function may differ by the order of observation units in the original dataset.
#' @param weight.cutoff a numeric vector indicating the cutoff value of the weight. Default=10.
#' @param formulaDisease a character string indicating the formula of the disease model.
#' @param family a character string indicating the error distribution and link function to be used in the disease model.
#' @param bs.N a numeric vector indicating the number of bootstrapping samples. If bs.N=1, then bootstrapping is not used and bs.replace is ignored.
#' @param bs.replace a character string indicating whether bootstrapping is done with replacement. Default=TRUE.
#' @param varilist a character vector indicating variable names for which you wish to compute standardized mean difference. List variable names as a vector like c("VariableA","VariableB")
#' @param corrmethod a character string indicating which correlation coefficient is to be computed. These include "Pearson" (default), "Spearman", "Polychoric", or "Polyserial". For tetrachoric use "Polychoric" and for biserial use "Polyserial". This relies on wCorr::weightedCorr
#' @export
#' @examples 
#' estimateATT()

estimateATT<-function(dataset,bexp,exp.status=1,cexp,fmethod.replace=TRUE,distbuf=0.1,exp.included=FALSE,long,lat,
                     PS.method="mgcv.GAM",
                     PS.method.data="original",
                     PS.formula,
                     PS.max_depth=5, PS.eta=0.1, PS.nthread=1, PS.eval_metric="auc", PS.objective="binary:logistic", PS.nrounds=100,
                     PS.cv.nround=10000,PS.cv.nfold=5,PS.cv.min.burnin=99,PS.early_stopping_rounds=100,
                     PS.cv.objective="binary:logistic",PS.cv.max_depth=c(4,5),PS.cv.eta=c(0.01,0.05,0.1),PS.cv.nthread=1,PS.cv.subsample=1,PS.cv.eval_metric="auc",PS.cv.colsample_bytree=1,PS.cv.min_child_weight=1,PS.cv.local.N=10,
                     CGPS.method="mgcv.GAM",CGPS.formula,
                     CGPS.cv.nround=10000,CGPS.cv.nfold=5,CGPS.cv.min.burnin=99,CGPS.early_stopping_rounds=100,
                     CGPS.max_depth=5, CGPS.eta=0.1, CGPS.nthread=1, CGPS.eval_metric="rmse", CGPS.objective="reg:squarederror", CGPS.nrounds=100,
                     CGPS.cv.objective="reg:squarederror",CGPS.cv.max_depth=c(4,5),CGPS.cv.eta=c(0.01,0.05,0.1),CGPS.cv.nthread=1,CGPS.cv.subsample=1,CGPS.cv.eval_metric="rmse",CGPS.cv.colsample_bytree=1,CGPS.cv.min_child_weight=1,CGPS.cv.local.N=10,
                     smethod="nearest",caliper_bw=NULL,smethod.replace=FALSE,weight.cutoff=10,
                     formulaDisease,family,
                     bs.N,bs.replace=TRUE,
                     varilist,corrmethod="Pearson") {
  start.time<-Sys.time()
  PSerror<-0
  CGPSerror<-0
  
  distmatched<-CGPSspatialmatch::matchdist(dataset,bexp,long=long,lat=lat,exp.status=exp.status,distbuf=distbuf,exp.included=TRUE,replace=FALSE)$matched.dataset
  
  if(bs.N>1) {
  bootsp<-replicate(bs.N,dplyr::sample_n(dataset,nrow(dataset),replace=bs.replace),simplify=FALSE)
  }
  if(bs.N==1) {
    bootsp<-list(dataset)
  }
  message(">>>>>>>>STEP 1: Matching by distance initiated")
  bootsp.m <- lapply(bootsp, 
                     function(data){
                       CGPSspatialmatch::matchdist(data,bexp,long=long,lat=lat,exp.status=exp.status,distbuf=distbuf,exp.included=exp.included,replace=fmethod.replace)$matched.dataset
                     })
  message(">>>>>>>>STEP 1: Matching by distance sucessfully done")
  
  
  message(">>>>>>>>STEP 2: PS estimation initiated")
  
  if(PS.method=="mgcv.GAM") {
      tryCatch(expr={
        if(PS.method.data=="original") {
          f1 <- as.formula(
            paste(PS.formula))
          PSmodel<-eval(bquote(mgcv::gam(.(f1), data=distmatched, family="binomial")))
          PS.m<-lapply(bootsp.m,function(data){
            data$PS<-predict(PSmodel,newdata=data,type="response")
            data
          })
        }
        
        if(PS.method.data=="aftermatch") {
          f1 <- as.formula(
            paste0(PS.formula,"+as.factor(strata_matchdist)"))
          PSmodel<-lapply(bootsp.m,function(data) {
            eval(bquote(mgcv::gam(.(f1), data=data, family="binomial")))
          })
          
          PS.m<-mapply(data=bootsp.m,model=PSmodel,function(data,model){
            data$PS<-predict(model,newdata=data,type="response")
            data
          },SIMPLIFY=FALSE)
        }
        
        message(">>>>>>>>STEP 2: PS estimation (GAM) sucessfully done")
        message(">>>>>>>>STEP 3: CGPS estimation initiated")
      }
      ,
      error=function(e) {e;message("PS estimation (GAM) failed: check the dataset and/or parameterization"); PSerror<<-1;CGPSerror<<-1},
      warning=function(w) {w;message("PS estimation (GAM) may have failed: check the dataset and/or parameterization")}
      )
    }
    

  if(PS.method=="xgboost") {
  tryCatch(expr={
    boost.fitdat<-data.matrix(distmatched[,PS.formula])
    boost.dat<-xgboost::xgb.DMatrix(boost.fitdat, label = distmatched[,bexp])
    param <- list(max_depth = PS.max_depth, eta = PS.eta, nthread = PS.nthread,
                  objective = PS.objective, eval_metric = PS.eval_metric)
    PSmodel <- xgboost::xgb.train(param=param,boost.dat,nrounds=PS.nrounds)
    
    pred.dat<-lapply(bootsp.m, function(data) {
      data.matrix(data[,PS.formula])
    })
    
    PS.m<-mapply(data=bootsp.m,pred.dat=pred.dat,function(data,pred.dat){
      data$PS<-predict(PSmodel,newdata=pred.dat,type="response")
      data
    },SIMPLIFY = FALSE)
    
    message(">>>>>>>>STEP 2: PS estimation (xgboost) sucessfully done")
    message(">>>>>>>>STEP 3: CGPS estimation initiated")
  }
  ,
  error=function(e) {e;message("PS estimation (xgboost) failed: check the dataset and/or parameterization"); PSerror<<-1;CGPSerror<<-1},
  warning=function(w) {w;message("PS estimation (xgboost) may have failed: check the dataset and/or parameterization")}
  )
  }
  
  
  if(PS.method=="xgboost.cv") {
    tryCatch(expr={
      boost.fitdat<-data.matrix(distmatched[,PS.formula])
      boost.dat<-xgboost::xgb.DMatrix(boost.fitdat, label = distmatched[,bexp])
      
      PSmodel <- CGPSspatialmatch::xgb.model.fit.cv(data=boost.dat,
                                                    cv.nround=PS.cv.nround,
                                                    cv.nfold=PS.cv.nfold,
                                                    cv.objective=PS.cv.objective,
                                                    cv.max_depth=PS.cv.max_depth,
                                                    cv.eta=PS.cv.eta,
                                                    cv.nthread=PS.cv.nthread,
                                                    cv.subsample=PS.cv.subsample,
                                                    cv.eval_metric=PS.cv.eval_metric,
                                                    cv.colsample_bytree=PS.cv.colsample_bytree,
                                                    cv.min_child_weight=PS.cv.min_child_weight,
                                                    cv.min.burnin=PS.cv.min.burnin,
                                                    early_stopping_rounds=PS.early_stopping_rounds,
                                                    cv.local.N=PS.cv.local.N)

      pred.dat<-lapply(bootsp.m, function(data) {
        data.matrix(data[,PS.formula])
      })
      
      PS.m<-mapply(data=bootsp.m,pred.dat=pred.dat,function(data,pred.dat){
        data$PS<-predict(PSmodel,newdata=pred.dat,type="response")
        data
      },SIMPLIFY = FALSE)
      
      message(">>>>>>>>STEP 2: PS estimation (xgboost.cv) sucessfully done")
      message(">>>>>>>>STEP 3: CGPS estimation initiated")
    }
    ,
    error=function(e) {e;message("PS estimation (xgboost.cv) failed: check the dataset and/or parameterization"); PSerror<<-1;CGPSerror<<-1},
    warning=function(w) {w;message("PS estimation (xgboost.cv) may have failed: check the dataset and/or parameterization")}
    )
  }
  
  
  
  
  if(PSerror !=1) {
    
    if(CGPS.method=="mgcv.GAM") {
    f2 <- as.formula(
      paste(CGPS.formula))
    tryCatch(expr={
      CGPS.model <- 
        eval(bquote(mgcv::gam(.(f2), data=dataset[dataset[,bexp]==1,])))

      message(">>>>>>>>STEP 3: CGPS estimation (GAM) sucessfully done")
      message(">>>>>>>>STEP 4: Matching by GPS initiated")
    }
    ,
    error=function(e) {e;message("CGPS estimation (GAM) failed: check the dataset and/or parameterization"); CGPSerror<<-1},
    warning=function(w) {w;message("CGPS estimation (GAM) may have failed: check the dataset and/or parameterization")}
    )
    }
    
    if(CGPS.method=="xgboost.cv") {

      tryCatch(expr={
        boost.fitdat<-data.matrix(dataset[dataset[,bexp]==1,CGPS.formula])
        boost.dat<-xgboost::xgb.DMatrix(boost.fitdat, label = dataset[dataset[,bexp]==1,cexp])
        CGPS.model <- xgb.model.fit.cv(data=boost.dat,
                                       cv.nround=CGPS.cv.nround,
                                       cv.nfold=CGPS.cv.nfold,
                                       cv.objective=CGPS.cv.objective,
                                       cv.max_depth=CGPS.cv.max_depth,
                                       cv.eta=CGPS.cv.eta,
                                       cv.nthread=CGPS.cv.nthread,
                                       cv.subsample=CGPS.cv.subsample,
                                       cv.eval_metric=CGPS.cv.eval_metric,
                                       cv.colsample_bytree=CGPS.cv.colsample_bytree,
                                       cv.min_child_weight=CGPS.cv.min_child_weight,
                                       cv.min.burnin=CGPS.cv.min.burnin,
                                       early_stopping_rounds=CGPS.early_stopping_rounds,
                                       cv.local.N=CGPS.cv.local.N)
        message(">>>>>>>>STEP 3: CGPS estimation (xgboost) sucessfully done")
        message(">>>>>>>>STEP 4: Matching by GPS initiated")
      }
      ,
      error=function(e) {e;message("CGPS estimation (xgboost) failed: check the dataset and/or parameterization"); CGPSerror<<-1},
      warning=function(w) {w;message("CGPS estimation (xgboost) may have failed: check the dataset and/or parameterization")}
      )
    }
    
    if(CGPS.method=="xgboost") {
      
      tryCatch(expr={
        boost.fitdat<-data.matrix(dataset[dataset[,bexp]==1,CGPS.formula])
        boost.dat<-xgboost::xgb.DMatrix(boost.fitdat, label = dataset[dataset[,bexp]==1,cexp])
        param <- list(max_depth = CGPS.max_depth, eta = CGPS.eta, nthread = CGPS.nthread,
                      objective = CGPS.objective, eval_metric = CGPS.eval_metric)
        CGPS.model <- xgboost::xgb.train(param=param,boost.dat,nrounds=CGPS.nrounds)
        
        message(">>>>>>>>STEP 3: CGPS estimation (xgboost.cv) sucessfully done")
        message(">>>>>>>>STEP 4: Matching by GPS initiated")
      }
      ,
      error=function(e) {e;message("CGPS estimation (xgboost.cv) failed: check the dataset and/or parameterization"); CGPSerror<<-1},
      warning=function(w) {w;message("CGPS estimation (xgboost.cv) may have failed: check the dataset and/or parameterization")}
      )
    }
  }
  
  
  if(CGPSerror !=1) {
    f3 <- as.formula(
      paste(formulaDisease))
    
    findat<-lapply(PS.m,
                   function(data) {
                     CGPSspatialmatch::cgpsmatch(data,bexp,cexp,"PS",CGPS.model,expstatus=exp.status,method=smethod,caliper_bw=caliper_bw,replace=smethod.replace,weight.cutoff=weight.cutoff)
                   })
    message(">>>>>>>>STEP 4: Matching by GPS successfully done")
    message(">>>>>>>>STEP 5: Disease model estimation initiated")
    
    tryCatch(expr={
      modelfit<-lapply(findat,function(data){
        eval(bquote(gnm::gnm(.(f3), data=data,family=family,eliminate=as.factor(strata_matchdist))))
      })
      message(">>>>>>>>STEP 5: Disease model estimation sucessfully done")
      message(">>ATT has been successfully estimated by CGPS spatial matching procedure with bootstrapping. Check the distribution of bootstrapped estimates")
    }
    ,
    error=function(e) {e;message("Disease model estimation failed: check the dataset and/or parameterization")},
    warning=function(w) {w;message("Disease model estimation may have failed: check the dataset and/or parameterization")}
    )
    
    coefest<-sapply(modelfit,function(object){
      summary(object)$coef[1]
    })
    
    seest<-sapply(modelfit,function(object){
      summary(object)$coef[2]
    })
    result<-c(mean(coefest),quantile(coefest,0.025),quantile(coefest,0.975))
    names(result)<-c("Mean","95% Lower eCI","95% Upper eCI")
    result2<-c(mean(coefest),mean(coefest)-1.96*mean(seest),mean(coefest)+1.96*mean(seest))
    names(result2)<-c("Mean","95% Lower CI","95% Upper CI")
    
    smd.org<-CGPSspatialmatch::smd(dataset,bexp,exp.status=exp.status,varinames=varilist)
    correxp.org<-CGPSspatialmatch::correxp(dataset,bexp,cexp,exp.status=exp.status,varinames=varilist,method=corrmethod)
    smd.m<-lapply(findat,function(data) {CGPSspatialmatch::smd(data,bexp,exp.status=exp.status,varinames=varilist)})
    correxp.m<-lapply(findat,function(data) {CGPSspatialmatch::correxp(data,bexp,paste0(cexp,"_cf"),exp.status=exp.status,varinames=varilist,method=corrmethod,weightname="weight")})
    smd.m<-do.call("rbind",smd.m)
    correxp.m<-do.call("rbind",correxp.m)
    
    total.exp.num<-lapply(bootsp,function(data){
      nrow(data[data[,bexp]==exp.status,])
    })
    matched.exp.num<-lapply(findat,function(data){
      length(unique(data[,"strata_matchdist"]))
    })
    total.exp.num<-do.call(rbind,total.exp.num)
    matched.exp.num<-do.call(rbind,matched.exp.num)
    
    match.info<-data.frame(total.exp.N=total.exp.num,
               matched.exp.N=matched.exp.num,
               matched.percentage=matched.exp.num/total.exp.num*100)

    
    
    rm(CGPSerror)
    rm(PSerror)
    end.time<-Sys.time()
    print(paste0("Elapsed time: ",round(end.time-start.time,3),attr(end.time-start.time,"units")))
    return(list(match.info=match.info,summary.empirical=result,summary=result2,distribution=coefest,modelfit=modelfit,
                smd.org=smd.org,smd.matched=colMeans(smd.m,na.rm=T),
                correxp.org=correxp.org,correxp.matched=colMeans(correxp.m,na.rm=T),
                smd.matched.bs=smd.m,correxp.matched.bs=correxp.m))
  }
}