#' estimateATT function
#'
#' This function matches exposed units and unexposed units by a pre-specified buffer distance (Euclidean distance)
#' @param dataset a dataset object.
#' @param PSdataset a dataset for PS estimation.
#' @param bexp a character string indicating the name of the binary exposure. Use apostrophe like "VariableName"
#' @param exp.status a numeric vector indicating the value indicating exposed units. Defalut=1
#' @param cexp a character string indicating the name of the continuous exposure. Use apostrophe like "VariableName"
#' @param fmethod.replace an indicator of whether one-to-n distance-matching will be performed with replacement or without replacement. Default=TRUE. If FALSE, matching is done without replacement, which of the performance has not been tested. If FALSE, note that the output of this function may differ by the order of observation units in the original dataset.
#' @param distbuf a numeric vector indicating the buffer distance by which exposed units and unexposed units are matched
#' @param exp.included an indicator of whether exposed units are matched with not only unexposed units but also other exposed units. Defalut is TRUE. If FALSE, exposed units are matched with only unexposed units. See details
#' @param long a character string indicating the name of the longitude variable of observation units
#' @param lat a character string indicating the name of the latitude variable of observation units
#' @param PS.method a character string or a vector of variable names, indicating the method of propensity score estimation. Options include "mgcv.GAM" (Generalized additive model in mgcv package), "xgboost" (extreme gradient boosting in xgboost package), and "xgboost.cv" (xgboost with cross-validation). Default="mgcv.GAM"
#' @param PS.formula a character string indicating the formula of propensity score estimation. For PSmethod="mgcv.GAM", this string must be like "the binary exposure variable ~ variableA+variableB+variableC+s(long,lat)". For PSmethod="xgboost", this must be a vector like c("variableA", "variableB", "variableC", "long", "lat")
#' @param PS.max_depth (xgboost only) a numeric vector indicating maximum depth of a tree. Default=5. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.eta (xgboost only) a numeric vector indicating the learning rate. Default=0.1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.nthread (xgboost only) a numeric vector indicating the number of thread. Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.eval_metric (xgboost only) a character string indicating evaluation metrics for validation data. Default="logloss". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.objective (xgboost only) a character string indicating the objective function. Default="binary:logistic". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.nrounds (xgboost only) a numeric vector indicating the number of rounds. Default=100. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.nround (xgboost.cv only) a numeric vector indicating the maximum number of rounds. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.nfold (xgboost.cv only) a numeric vector indicating N-fold cross-validation. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.objective (xgboost.cv only) a character string indicating the objective function. Default="binary:logistic" This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.max_depth (xgboost.cv only) a numeric vector indicating maximum depth of a tree. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.eta (xgboost.cv only) a numeric vector indicating the learning rate. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.nthread (xgboost.cv only) a numeric vector indicating the number of thread. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.subsample (xgboost.cv only) a numeric vector indicating subsample ratio of the training instance. Default=0.5. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.eval_metric (xgboost.cv only) a character string indicating evaluation metrics for validation data. Default="logloss". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.colsample_bytree (xgboost.cv only) a numeric vector indicating subsample ratio of columns when constructing each tree. Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.min_child_weight (xgboost.cv only) Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.lambda (xgboost.cv only) L2 Regularization term on weights. Default=0. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.lambda_bias (xgboost.cv only) L2 Regularization term on bias. Default=0. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.alpha (xgboost.cv only) L1 Regularization term on weights. Default=0. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.scale_pos_weight (xgboost.cv only) Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param PS.cv.min.burnin (xgboost.cv only) a numeric vector indicating the minimum number of rounds. Default=99, which indicates the minimum number of rounds is 100. This must be smaller than PS.early_stopping_rounds. 
#' @param PS.early_stopping_rounds (xgboost.cv only) a numeric vector indicating when xgboost stops. If the evaluation metric did not decrease until when (code)PS.early_stopping_rounds, xgboost stops. This saves time. Default=100.
#' @param PS.cv.local.N (xgboost.cv only) a numeric vector indicating how many times xgboost will search (local) mini-ma of the evaluation metric function. Default=10. 
#' @param CGPS.method a character string or a vector of variable names, indicating the method of conditional propensity score estimation. Options include "mgcv.GAM", "xgboost", and "xgboost.cv". Default="mgcv.GAM".
#' @param CGPS.formula a character string indicating the formula of conditional generalized propensity score estimation.
#' @param CGPS.max_depth (xgboost only) a numeric vector indicating maximum depth of a tree. Default=5. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.eta (xgboost only) a numeric vector indicating the learning rate. Default=0.1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.nthread (xgboost only) a numeric vector indicating the number of thread. Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.eval_metric (xgboost only) a character string indicating evaluation metrics for validation data. Default="rmse". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.objective (xgboost only) a character string indicating the objective function. Default="reg:squarederror". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.nrounds (xgboost only) a numeric vector indicating the number of rounds. Default=100. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.nround (xgboost.cv only) a numeric vector indicating the maximum number of rounds. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.nfold (xgboost.cv only) a numeric vector indicating N-fold cross-validation. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.objective (xgboost.cv only) a character string indicating the objective function. Default="binary:logistic". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.max_depth (xgboost.cv only) a numeric vector indicating maximum depth of a tree. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.eta (xgboost.cv only) a numeric vector indicating the learning rate. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.nthread (xgboost.cv only) a numeric vector indicating the number of thread. Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.subsample (xgboost.cv only) a numeric vector indicating subsample ratio of the training instance. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.eval_metric (xgboost.cv only) a character string indicating evaluation metrics for validation data. Default="rmse". This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.colsample_bytree (xgboost.cv only) a numeric vector indicating subsample ratio of columns when constructing each tree. Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.min_child_weight (xgboost.cv only) Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.lambda (xgboost.cv only) L2 Regularization term on weights. Default=0. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.lambda_bias (xgboost.cv only) L2 Regularization term on bias. Default=0. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.alpha (xgboost.cv only) L1 Regularization term on weights. Default=0. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.scale_pos_weight (xgboost.cv only) Default=1. This is a wrapper for xgboost::xgb.train. See https://xgboost.readthedocs.io/en/latest/parameter.html
#' @param CGPS.cv.min.burnin (xgboost.cv only) a numeric vector indicating the minimum number of rounds. Default=99, which indicates the minimum number of rounds is 100. This must be smaller than PS.early_stopping_rounds.
#' @param CGPS.early_stopping_rounds (xgboost.cv only) a numeric vector indicating when xgboost stoCGPS. If the evaluation metric did not decrease until when (code)CGPS.early_stopping_rounds, xgboost stops This saves time. Default=100.
#' @param CGPS.cv.local.N (xgboost.cv only) a numeric vector indicating how many times xgboost will search (local) mini-ma of the evaluation metric function. Default=10
#' @param smethod method a character string indicating the matching method used to conduct matching by GPS Default="nearest". Options include "nearest" (nearest neighbor matching), "nearestcaliper" (nearest neighbor caliper matching), "caliper" (caliper matching)
#' @param caliper_bw a numeric vector indicating caliper bandwidth. Default=0.1. If method is "nearest", this parameter is ignored.
#' @param smethod.replace an indicator of whether matching by GPS is done with replacement. Default=FALSE. If FALSE, matching is done without replacement. If FALSE, note that the output of this function may differ by the order of observation units in the original dataset. Note: Matching with replacement has not been tested for its performance.
#' @param formulaDisease a character string indicating the formula of the disease model.
#' @param family a character string indicating the error distribution and link function to be used in the disease model.
#' @param bs.N a numeric vector indicating the number of bootstrapping samples. If bs.N=1, then bootstrapping is not used and bs.replace is ignored.
#' @param bs.replace a character string indicating whether bootstrapping is done with replacement. Default=TRUE.
#' @param varilist a character vector indicating variable names for which you wish to compute standardized mean difference. List variable names as a vector like c("VariableA","VariableB")
#' @param corrmethod a character string indicating which correlation coefficient is to be computed. These include "Pearson" (default), "Spearman", "Polychoric", or "Polyserial". For tetrachoric use "Polychoric" and for biserial use "Polyserial". This is a wrapper for wCorr::weightedCorr
#' @export
#' @examples 
#' estimateATT()

estimateATT<-function(dataset,PSdataset,bexp,exp.status=1,cexp,fmethod.replace=TRUE,distbuf=0.1,exp.included=FALSE,long,lat,
                     PS.method="mgcv.GAM",
                     PS.method.data="original",
                     PS.formula,
                     PS.max_depth=5, PS.eta=0.1, PS.nthread=1, PS.eval_metric="logloss", PS.objective="binary:logistic", PS.nrounds=50,
                     PS.cv.nround=100,PS.cv.nfold=5,PS.cv.min.burnin=49,PS.early_stopping_rounds=50,
                     PS.cv.objective="binary:logistic",PS.cv.max_depth=c(4,5),PS.cv.eta=c(0.01,0.05,0.1),PS.cv.gamma=c(0),
                     PS.cv.nthread=1,PS.cv.subsample=1,PS.cv.eval_metric="logloss",
                     PS.cv.colsample_bytree=1,PS.cv.min_child_weight=1,PS.cv.lambda=c(0),PS.cv.lambda_bias=c(0),PS.cv.alpha=c(0),PS.cv.scale_pos_weight=1,
                     PS.cv.local.N=10,
                     CGPS.method="mgcv.GAM",CGPS.formula,
                     CGPS.max_depth=5, CGPS.eta=0.1, CGPS.nthread=1, CGPS.eval_metric="rmse", CGPS.objective="reg:squarederror", CGPS.nrounds=50,
                     CGPS.cv.nround=100,CGPS.cv.nfold=5,CGPS.cv.min.burnin=49,CGPS.early_stopping_rounds=50,
                     CGPS.cv.objective="reg:squarederror",CGPS.cv.max_depth=c(4,5),CGPS.cv.eta=c(0.01,0.05,0.1),CGPS.cv.gamma=c(0),
                     CGPS.cv.nthread=1,CGPS.cv.subsample=1,CGPS.cv.eval_metric="rmse",
                     CGPS.cv.colsample_bytree=1,CGPS.cv.min_child_weight=1,CGPS.cv.lambda=c(0),CGPS.cv.lambda_bias=c(0),CGPS.cv.alpha=c(0),CGPS.cv.scale_pos_weight=1,
                     CGPS.cv.local.N=10,
                     smethod="nearest",caliper_bw=NULL,smethod.replace=FALSE,
                     formulaDisease,family,
                     bs.N,bs.replace=TRUE,corrmethod="Pearson",
                     varilist,modelinfo=FALSE) {
  start.time<-Sys.time()
  PSerror<-0
  CGPSerror<-0
  

  bootsp<-list(dataset)
  
  message(">>>>>>>>STEP 1: Matching by distance initiated")
  bootsp.m <- lapply(bootsp, 
                     function(data){
                       matchdist(data,bexp,long=long,lat=lat,exp.status=exp.status,distbuf=distbuf,exp.included=exp.included,replace=fmethod.replace)$matched.dataset
                     })
  message(">>>>>>>>STEP 1: Matching by distance sucessfully done")
  
  
  message(">>>>>>>>STEP 2: PS estimation initiated")
  
  if(PS.method=="mgcv.GAM") {
      tryCatch(expr={
        if(PS.method.data=="original") {
          f1 <- as.formula(
            paste(PS.formula))
          PSmodel<-eval(bquote(mgcv::gam(.(f1), data=PSdataset, family="binomial")))
          PS.m<-lapply(bootsp.m,function(data){
            data$PS<-predict(PSmodel,newdata=data,type="response")
            data
          })
        }
        
        if(PS.method.data=="aftermatch") {
          f1 <- as.formula(
            paste0(PS.formula,"+as.factor(strata_matchdist)"))
          PSmodel<-lapply(bootsp.m,function(data) {
            eval(bquote(mgcv::gam(.(f1), data=PSdataset, family="binomial")))
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
    boost.fitdat<-data.matrix(PSdataset[,PS.formula])
    boost.dat<-xgboost::xgb.DMatrix(boost.fitdat, label = PSdataset[,bexp])
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
      boost.fitdat<-data.matrix(PSdataset[,PS.formula])
      boost.dat<-xgboost::xgb.DMatrix(boost.fitdat, label = PSdataset[,bexp])
      
      PSmodel <- xgb.model.fit.cv(data=boost.dat,
                                                    cv.nround=PS.cv.nround,
                                                    cv.nfold=PS.cv.nfold,
                                                    cv.objective=PS.cv.objective,
                                                    cv.max_depth=PS.cv.max_depth,
                                                    cv.eta=PS.cv.eta,
                                                    cv.gamma=PS.cv.gamma,
                                                    cv.lambda=PS.cv.lambda,
                                                    cv.lambda_bias=PS.cv.lambda_bias,
                                                    cv.alpha=PS.cv.alpha,
                                                    cv.nthread=PS.cv.nthread,
                                                    cv.subsample=PS.cv.subsample,
                                                    cv.eval_metric=PS.cv.eval_metric,
                                                    cv.colsample_bytree=PS.cv.colsample_bytree,
                                                    cv.min_child_weight=PS.cv.min_child_weight,
                                                    cv.scale_pos_weight=PS.cv.scale_pos_weight,
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
                                       cv.gamma=CGPS.cv.gamma,
                                       cv.lambda=CGPS.cv.lambda,
                                       cv.lambda_bias=CGPS.cv.lambda_bias,
                                       cv.alpha=CGPS.cv.alpha,
                                       cv.nthread=CGPS.cv.nthread,
                                       cv.subsample=CGPS.cv.subsample,
                                       cv.eval_metric=CGPS.cv.eval_metric,
                                       cv.colsample_bytree=CGPS.cv.colsample_bytree,
                                       cv.min_child_weight=CGPS.cv.min_child_weight,
                                       cv.scale_pos_weight=CGPS.cv.scale_pos_weight,
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
                     cgpsmatch(data,bexp,cexp,"PS",CGPS.model,expstatus=exp.status,method=smethod,caliper_bw=caliper_bw,replace=smethod.replace,weight.cutoff=100)
                   })
    message(">>>>>>>>STEP 4: Matching by GPS successfully done")
    message(">>>>>>>>STEP 5: Disease model estimation initiated")
    
    if(bs.N>1) {
      findat.bootsp<-replicate(bs.N,dplyr::sample_n(findat[[1]],nrow(findat[[1]]),replace=bs.replace),simplify=FALSE)
    }
    if(bs.N==1) {
      findat.bootsp<-list(findat)
    }
    
    tryCatch(expr={
      modelfit<-lapply(findat.bootsp,function(data){
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
    
    smd.org<-smd(dataset,bexp,exp.status=exp.status,varinames=varilist)
    correxp.org<-correxp(dataset,bexp,cexp,exp.status=exp.status,varinames=varilist,method=corrmethod)
    smd.m<-lapply(findat,function(data) {smd(data,bexp,exp.status=exp.status,varinames=varilist)})
    smd.m<-do.call("rbind",smd.m)
    
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
    
    if(modelinfo==TRUE) {
    return(list(match.info=match.info,summary=result2,coef_dist=coefest,se_dist=seest,modelfit=modelfit,
                smd.org=smd.org,smd.matched=colMeans(smd.m,na.rm=T),
                correxp.org=correxp.org,
                smd.matched.bs=smd.m))
    }
    if(modelinfo==FALSE) {
      return(list(match.info=match.info,summary=result2,coef_dist=coefest,se_dist=seest,
                  smd.org=smd.org,smd.matched=colMeans(smd.m,na.rm=T),
                  correxp.org=correxp.org,
                  smd.matched.bs=smd.m))
    }    
  }
}