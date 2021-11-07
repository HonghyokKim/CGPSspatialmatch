#' xgb.model.fit.cv function
#'
#' This internal function explores optimal parameters in xgb.train using cross-validation
#' xgb.model.fit.cv()

xgb.model.fit.cv<-function(data,cv.objective,cv.nround=10000,cv.nfold=5,cv.max_depth=c(4,5),cv.eta=c(0.01,0.05,0.1),cv.nthread=1,cv.subsample=1,cv.eval_metric,cv.colsample_bytree=1,cv.min_child_weight=1,early_stopping_rounds=50,cv.local.N=10) {
  CV.start.time<-Sys.time()
  hyperparam<-expand.grid(cv.max_depth,cv.eta,cv.subsample,cv.local.N)
  colnames(hyperparam)<-c("max_depth","subsample","eta","cv.local.N")
  CVrun<-mapply(hyperparam[,4],hyperparam[,1],hyperparam[,2],hyperparam[,3],
                FUN=function(cv.local.N,max_depth,eta,subsample) {
                  param <- list(objective = cv.objective,
                                max_depth = max_depth,
                                eta = eta,
                                nthread=cv.nthread,
                                subsample=subsample,
                                eval_metric=cv.eval_metric,
                                colsample_bytree = cv.colsample_bytree, 
                                min_child_weight = cv.min_child_weight
                  )
                  seedn = sample.int(100000, 1)
                  set.seed(seedn)
                  cvfit <- xgboost::xgb.cv(data=data, params = param, 
                                           nfold=cv.nfold, nrounds=cv.nround,
                                           verbose = F,early_stopping_rounds=early_stopping_rounds,maximize=FALSE)
                  
                  min.evalmetric<-cvfit$evaluation_log[cvfit$best_iteration,4]
                  min.evalmetric.index<-cvfit$evaluation_log[cvfit$best_iteration,1]
                  return(data.frame(max_depth=max_depth,eta=eta,subsample=subsample,
                                    seedn=seedn,iter=min.evalmetric.index,evalmetric=min.evalmetric))
                },
                
                SIMPLIFY = FALSE)
  CV.result<-do.call(rbind,CVrun)
  CV.result<-CV.result[CV.result[,6]==min(CV.result[,6]),]         
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
