#' xgb.model.fit.cv function
#'
#' This internal function explores optimal parameters in xgb.train using cross-validation
#' xgb.model.fit.cv()

xgb.model.fit.cv<-function(data,cv.objective,cv.nround=10000,cv.nfold=5,cv.max_depth=c(4,5),cv.eta=c(0.01,0.05,0.1),cv.nthread=1,cv.subsample=1,cv.eval_metric,cv.colsample_bytree=1,cv.min_child_weight=1,early_stopping_rounds=50,cv.local.N=10) {
  CV.start.time<-Sys.time()
  hyperparam<-expand.grid(cv.max_depth,cv.eta,cv.subsample,cv.local.N)
  colnames(hyperparam)<-c("max_depth","subsample","eta","cv.local.N")
  CVrun<-mapply(CGPSspatialmatch::xgb.model.cv,hyperparam[,4],hyperparam[,1],hyperparam[,2],hyperparam[,3],SIMPLIFY = FALSE)
  CV.result<-do.call(rbind,CVrun)
  CV.result<-CV.result[CV.result$test_auc_mean==min(CV.result$test_auc_mean),]         
  
  nround = CV.result$iter
  param <- list(objective = cv.objective,
                max_depth = CV.result$max_depth,
                eta = CV.result$eta,
                nthread=cv.nthread,
                subsample=CV.result$subsample,
                eval_metric=cv.eval_metric,
                colsample_bytree = cv.colsample_bytree, 
                min_child_weight = cv.min_child_weight
  )
  set.seed(CV.result$seedn)
  return(xgboost::xgb.train(data=data, params=param, nrounds=nround))
  CV.end.time<-Sys.time()
  print(paste0(round(CV.end.time-CV.start.time,3),attr(CV.end.time-CV.start.time,"units")))
}
