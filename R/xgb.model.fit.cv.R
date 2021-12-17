#' xgb.model.fit.cv function
#'
#' This internal function explores optimal parameters in xgb.train using cross-validation
#' xgb.model.fit.cv()

xgb.model.fit.cv<-function(data,cv.objective,cv.nround,cv.nfold,cv.min.burnin,cv.max_depth,cv.eta,cv.gamma,cv.nthread,cv.subsample,cv.eval_metric,cv.colsample_bytree,cv.min_child_weight,cv.lambda,cv.lambda_bias,cv.alpha,cv.scale_pos_weight,early_stopping_rounds,cv.local.N) {
  CV.start.time<-Sys.time()
  hyperparam<-expand.grid(cv.max_depth,cv.eta,cv.gamma,cv.subsample,cv.colsample_bytree,cv.min_child_weight,cv.lambda,cv.lambda_bias,cv.alpha,cv.scale_pos_weight,cv.local.N)
  colnames(hyperparam)<-c("max_depth","eta","gamma","subsample","colsample_bytree","min_child_weight","lambda","lambda_bias","alpha","cv.local.N")
  CVrun<-mapply(hyperparam[,1],hyperparam[,2],hyperparam[,3],hyperparam[,4],hyperparam[,5],hyperparam[,6],hyperparam[,7],hyperparam[,8],hyperparam[,9],hyperparam[,10],hyperparam[,11],
                FUN=function(max_depth,eta,gamma,subsample,colsample_bytree,min_child_weight,lambda,lambda_bias,alpha,scale_pos_weight,cv.local.N) {
                  param <- list(objective = cv.objective,
                                max_depth = max_depth,
                                eta = eta,
                                gamma = gamma,
                                nthread=cv.nthread,
                                eval_metric=cv.eval_metric,
                                subsample=subsample,
                                colsample_bytree = colsample_bytree, 
                                min_child_weight = min_child_weight,
                                lambda=lambda,
                                lambda_bias=lambda_bias,
                                alpha=alpha,
                                scale_pos_weight=scale_pos_weight
                                )
                  seedn = sample.int(100000, 1)
                  set.seed(seedn)
                  cvfit <- xgboost::xgb.cv(data=data, params = param, 
                                           nfold=cv.nfold, nrounds=cv.nround,
                                           verbose = F,early_stopping_rounds=early_stopping_rounds,maximize=FALSE)
                  cvlog<-cvfit$evaluation_log
                  cvlog<-subset(cvlog,iter>cv.min.burnin)
                  cvlog<-data.frame(cvlog)
                  if(cv.eval_metric %in% c("rmse","logloss")) {
                    cvlog<-cvlog[cvlog[,4]==min(cvlog[,4]),]
                  }
                  
                  if(cv.eval_metric %in% c("auc","map")) {
                    cvlog<-cvlog[cvlog[,4]==max(cvlog[,4]),]
                  }
                  
                  return(data.frame(max_depth=max_depth,eta=eta,gamma=gamma,subsample=subsample,
                                    colsample_bytree=colsample_bytree,min_child_weight=min_child_weight,lambda=lambda,lambda_bias=lambda_bias,alpha=alpha,scale_pos_weight=scale_pos_weight,
                                    seedn=seedn,iter=cvlog[,1],evalmetric=cvlog[,4]))
                },
                
                SIMPLIFY = FALSE)
  CV.result<-do.call(rbind,CVrun)
  if(cv.eval_metric %in% c("rmse","logloss")) {
  CV.result<-CV.result[CV.result[,13]==min(CV.result[,13]),]
  }
  if(cv.eval_metric %in% c("auc","map")) {
    CV.result<-CV.result[CV.result[,13]==max(CV.result[,13]),]
  }
  
  if(nrow(CV.result)>1) {
    CV.result<-CV.result[1,]
  }
  
  print(CV.result)
  nround = CV.result$iter
  param <- list(objective = cv.objective,
                max_depth = CV.result$max_depth,
                eta = CV.result$eta,
                gamma = CV.result$gamma,
                nthread=cv.nthread,
                subsample=CV.result$subsample,
                eval_metric=cv.eval_metric,
                colsample_bytree = CV.result$colsample_bytree, 
                min_child_weight = CV.result$min_child_weight,
                lambda = CV.result$lambda,
                lambda_bias = CV.result$lambda_bias,
                alpha = CV.result$alpha,
                scale_pos_weight=CV.result$scale_pos_weight
  )
  set.seed(CV.result$seedn)
  return(xgboost::xgb.train(data=data, params=param, nrounds=nround))
  CV.end.time<-Sys.time()
  print(paste0(round(CV.end.time-CV.start.time,3),attr(CV.end.time-CV.start.time,"units")))
}
