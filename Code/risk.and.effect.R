#######
####### TRAIN RISK MODEL
#######
train.risk <- function(treatment.arm=NULL, Y.train=NULL, X.train=NULL,
                       W.train=NULL, alpha.reg=0.5, folds=5, spline=FALSE){
  # fit model with only prognostic factors
  lp.model <- stats::glm(Y.train ~ X.train, family="binomial")
  
  # obtain linear predictors
  lp.train <- as.vector(X.train%*%lp.model$coefficients[-1]) # omit intercept
  lp.train.center <- lp.train - mean(lp.train)  # center LP
  data.train <- data.frame(Y=Y.train, W=W.train, lp=lp.train.center)

  # train model
  if (spline) { # smooth linear predictors
    final.model <- stats::glm(Y ~ W*rms::rcs(lp, 3), family="binomial", data=data.train) # 3 knots, 2 degrees of freedom
  } else if (!spline){
    final.model <- stats::glm(Y ~ W + I(1-W):lp + W:lp, family="binomial", data=data.train)
  }

  # print(coef(final.model))

  return(list(lp.model=lp.model, final.model=final.model))
}

#######
####### OBTAIN PREDICTIONS RISK MODEL
#######
predictions.risk <- function(X.test=NULL, lp.model=NULL, final.model=NULL){
  # dataframes for treatment effect predictions
  lp.test <- as.vector(X.test%*%lp.model$coefficients[-1])
  lp.test.center <- lp.test - mean(lp.test)
  X.0.test <- data.frame(W=rep(0, length(lp.test)), lp=lp.test.center)
  X.1.test <- data.frame(W=rep(1, length(lp.test)), lp=lp.test.center)

  # obtain predictions on test set
  p.0 <- predict(final.model, newdata=X.0.test, type="response")
  p.1 <- predict(final.model, newdata=X.1.test, type="response")
  tau.hat <- p.0 - p.1

  return(list(tau.hat=tau.hat, p.0=p.0, p.1=p.1, X.0.test=X.0.test,
              X.1.test=X.1.test, lp.test=lp.test))
}

#######
####### CREATE DATAFRAME FOR EFFECT MODEL
#######
create.X <- function(treatment.arm=NULL, X.PF=NULL, W=NULL){
  # merge independent variables
  X <- cbind(W, X.PF)

  # add cross terms
  crossterm.labels <- c()
  labels.PF <- colnames(as.data.frame(X.PF))
  for (i in 1:(length(labels.PF))){
    X <- cbind(X,  W * X.PF[, i])
    crossterm.labels <- c(crossterm.labels, paste(treatment.arm, ".", labels.PF[i], sep=""))
  }
  colnames(X) <- c(treatment.arm, labels.PF, crossterm.labels)

  return(X)
}

#######
####### DO CROSS-VALIDATION N TIMES
#######
train.effect <- function(treatment.arm=NULL, Y.train=NULL, X.train=NULL, W.train=NULL,
                         folds=5, alpha.reg=alpha.reg, penalize=FALSE, print=FALSE){
  # create dataframe with cross terms
  X.full <- create.X(treatment.arm=treatment.arm, X.PF=X.train, W=W.train)

  if (penalize){
    # perform regularization with cross-validation once
    penalty.vec <- c(0, rep(1, ncol(X.full)-1)) # penalty vector indicating which coefficients shrinkage is applied to
    # no shrinkage on treatment arm indicator
    result <- glmnet::cv.glmnet(X.full, Y.train, alpha=alpha.reg, family="binomial", nfolds=folds,
                        lambda=seq(0, 0.01, 0.0001), penalty.factor=penalty.vec)

    # return minimum cross-validated lambda
    best.lambda <- result$lambda.min

    # obtain results with best lambda, thus not cross-validation again
    final.model <- glmnet::glmnet(X.full, Y.train, alpha=alpha.reg, family="binomial",
                             lambda=best.lambda, penalty.factor=penalty.vec)

    # show coefficients with best lambda
    if (print){
      cat('best lambda: ', final.model$lambda, '\n')
      result.beta <- data.frame(Names=row.names(coef(final.model)), Coefficients=coef(final.model)[,1])
      print(result.beta[order(result.beta[,"Coefficients"]),])
    }
  } else{
    final.model <- stats::glm(Y.train~X.full, family="binomial")
  }

  return(list(final.model=final.model, X.full=X.full))
}

#######
####### OBTAIN TREATMENT EFFECT PREDICTIONS FOR EFFECT MODEL
#######
predictions.effect <- function(treatment.arm=NULL, X=NULL, penalized=TRUE, final.model=NULL){
  # obtain treatment effect predictions
  # type="response" gives the fitted probabilities for family="binomial"
  X.test.0 <- create.X(treatment.arm=treatment.arm, X.PF=X, W=rep(0, nrow(X)))
  X.test.1 <- create.X(treatment.arm=treatment.arm, X.PF=X, W=rep(1, nrow(X)))
  if (penalized){
    p.0 <- as.numeric(predict(final.model, newx=X.test.0, type="response"))
    p.1 <- as.numeric(predict(final.model, newx=X.test.1, type="response"))
  } else{
    p.0 <- as.numeric(predict(final.model, newdata=as.data.frame(X.test.0), type="response"))
    p.1 <- as.numeric(predict(final.model, newdata=as.data.frame(X.test.1), type="response"))
  }
  tau.hat <- as.numeric(p.0 - p.1)
  return(list(X.test.0=X.test.0, X.test.1=X.test.1, p.0=p.0, p.1=p.1, tau.hat=tau.hat))
}

#######
####### SUBOPTIMAL MODEL
#######
create.suboptimal.model <- function(true.model=NULL, X.test=NULL,
                                   coef.W=1, coef.LP=1, coef.W.LP=1,
                                   constant=0){
  optimal.model <- true.model$final.model
  suboptimal.model <- optimal.model
  suboptimal.model$coefficients["(Intercept)"] <- coef(optimal.model)["(Intercept)"]
  suboptimal.model$coefficients["W"] <- coef(optimal.model)["W"]*coef.W+constant  # increase treatment effect; add constant to equate ATT
  suboptimal.model$coefficients["I(1 - W):lp"] <- coef(optimal.model)["I(1 - W):lp"]*coef.LP   # increase discriminative ability
  suboptimal.model$coefficients["W:lp"] <- coef(optimal.model)["W:lp"]*coef.W.LP         # discriminative ability of treated patients is overestimated

  p.0.suboptimal <- predict(suboptimal.model, newdata=X.test$X.0.test, type="response")
  p.1.suboptimal <- predict(suboptimal.model, newdata=X.test$X.1.test, type="response")
  tau.hat.suboptimal <- p.0.suboptimal - p.1.suboptimal

  # cat('ATT: ', round(mean(tau.hat.suboptimal)*100, 1), '\n')

  return(list(tau.hat=tau.hat.suboptimal, p.0=p.0.suboptimal, p.1=p.1.suboptimal))
}
