#' xgb.model.cv function
#'
#' This internal function explores optimal parameters in xgb.train using cross-validation
#' xgb.model.cv()

xgb.model.cv<-function(cv.local.N,max_depth,eta,subsample) {
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
  cvfit <- xgboost::xgb.cv(data=boost.dat, params = param, 
                           nfold=cv.nfold, nrounds=cv.nround,
                           verbose = F,early_stopping_rounds=early_stopping_rounds,maximize=FALSE)
  
  min.evalmetric<-cvfit$evaluation_log[cvfit$best_iteration,4]
  min.evalmetric.index<-cvfit$evaluation_log[cvfit$best_iteration,1]
  return(data.frame(max_depth=max_depth,eta=eta,subsample=subsample,
                    seedn=seedn,iter=min.evalmetric.index,evalmetric=min.evalmetric))
}
