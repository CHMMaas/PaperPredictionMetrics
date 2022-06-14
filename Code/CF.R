#######
####### CAUSAL FOREST
#######
train.CF <- function(Y.train=NULL, X.train=NULL, W.train=NULL, results=FALSE, tune=FALSE){
  if (tune){tune.par <- "all"}    # tune when training model
  else{tune.par <- "none"}        # don't tune in case of simulating treatment effect

  ####### ESTIMATE E[W|X]
  forest.W <- grf::regression_forest(X.train, W.train, tune.parameters = tune.par)
  W.hat <- predict(forest.W)$predictions

  ####### ESTIMATE E[Y|X]
  forest.Y <- grf::regression_forest(X.train, Y.train, tune.parameters = tune.par)
  Y.hat <- predict(forest.Y)$predictions

  ####### VARIABLE SELECTION
  # forest.Y.varimp <- grf::variable_importance(forest.Y)
  # selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.05)
  # cat('selected variables:', selected.vars, '\n')

  # split data set on treatment
  Y.train <- 1-Y.train  # take reverse of Y, because we want to compute E[Y(0) - Y(1) | X=x]

  ####### CAUSAL FOREST
  # TODO: for small sample: increase honesty.fraction=0.8, honesty.prune.leaves=FALSE, and increase num.trees in training
  # tune.parameters: "sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty"
  num.trees <- 2000
  tau.forest <- grf::causal_forest(X.train, Y.train, W.train,   # use selected variables, Y and W
                              W.hat = W.hat, Y.hat = Y.hat,  # pre-estimates of Y and W
                              seed = 1,                      # reproducability
                              num.trees = num.trees,         # increase (double) number of trees for small samples
                              tune.parameters = tune.par,    # tune all or no parameters
                              tune.num.trees = 1000)         # increase number of trees for parameter tuning for small samples

  # print(tau.forest$tuning.output)                          # print tuning output

  tau.hat <- predict(tau.forest)$predictions
  mu.x <- Y.hat - (W.train - W.hat) * tau.hat                      # E[Y | X, W] = Y.hat + (W - W.hat) * (- tau.hat)

  if (results == TRUE){
    # check propensity scores bounded away from 0 and 1
    hist(W.hat, xlab="Check propensity score bounded away from 0 and 1", main="")

    # print variable importance
    print(as.data.frame(cbind(labels(X.train)[[2]], round(variable_importance(tau.forest), 5))))

    # Estimate treatment effects for the training data using out-of-bag prediction.
    tau.hat.oob <- predict(tau.forest)
    hist(tau.hat.oob$predictions, xlab="Predictions probability of diabetes on OOB", main="")

    # Estimate the conditional average treatment effect on the full sample (CATE).
    CATE <- grf::average_treatment_effect(tau.forest, target.sample = "all")
    cat('average treatment effect on full training sample:    ', CATE[1], 'sd:', CATE[2], '\n')

    # Estimate the conditional average treatment effect on the treated sample (CATT).
    CATT <- grf::average_treatment_effect(tau.forest, target.sample = "treated")
    cat('average treatment effect on treated training sample: ', CATT[1], 'sd:', CATT[2], '\n')
  }

  return(list(tau.forest=tau.forest, num.trees=num.trees, forest.Y=forest.Y, forest.W=forest.W, mu.x=mu.x))
}

#######
####### OBTAIN TREATMENT EFFECT PREDICTIONS FOR CAUSAL FOREST
#######
predictions.CF <- function(X.test=NULL, final.model=NULL){
  # selected.vars <- final.model$selected.vars
  tau.hat <- predict(final.model$tau.forest, newdata=X.test, estimate.variance = FALSE)$predictions

  Y.hat <- predict(final.model$forest.Y, newdata=X.test)$predictions
  W.hat <- predict(final.model$forest.W, newdata=X.test)$predictions

  p.0 <- as.numeric(Y.hat + W.hat * tau.hat)       # E[Y | X, W = 0] = Y.hat - W.hat * (- tau.hat)
  p.1 <- as.numeric(Y.hat - (1 - W.hat) * tau.hat) # E[Y | X, W = 1] = Y.hat + (1 - W.hat) * (- tau.hat)

  return(list(tau.hat=tau.hat, p.0=p.0, p.1=p.1))
}
