#' xgb.model.cv function
#'
#' This internal function explores optimal parameters in xgb.train using cross-validation
#' xgb.model.cv()

xgb.model.cv<-function(data,cv.objective,cv.nround=10000,cv.nfold=5,cv.max_depth=5,cv.eta=0.1,cv.nthread=1,cv.subsample=0.5,cv.gamma=0.1,cv.eval_metric,cv.colsample_bytree=1,cv.min_child_weight=1,early_stopping_rounds=50,cv.local.N=100) {

best.seedn = NA
best.param = list()
best.evalmetric = Inf
best.evalmetric.index = 0
for (aaaa in 1:cv.local.N) {
  cv.nround = cv.nround
  cv.nfold = cv.nfold
  param <- list(objective = cv.objective,
                max_depth = cv.max_depth,
                eta = cv.eta,
                nthread=cv.nthread,
                subsample=cv.subsample,
                eval_metric=cv.eval_metric,
                colsample_bytree = cv.colsample_bytree, 
                min_child_weight = cv.min_child_weight
  )
  seedn = sample.int(100000, 1)
  set.seed(seedn)
  cvfit <- xgboost::xgb.cv(data=data, params = param, 
                 nfold=cv.nfold, nrounds=cv.nround,
                 verbose = F,early_stopping_rounds=early_stopping_rounds,maximize=FALSE)
  
  min.evalmetric = min(cvfit$evaluation_log[,4])
  min.evalmetric.index = which(cvfit$evaluation_log[,4]==min.evalmetric)
  if (min.evalmetric < best.evalmetric) {
    best.evalmetric = min.evalmetric
    best.evalmetric.index = min.evalmetric.index
    best.seedn = seedn
    best.param = param
  }
}
nround = best.evalmetric.index
set.seed(best.seedn)
return(xgboost::xgb.train(data=data, params=best.param, nrounds=nround))
}
