#######
####### COMMON FUNCTIONS
#######
library(HTEPredictionMetrics) # new metrics
library(dplyr)                # bootstrap matched pairs slice()
library(gridExtra)            # table in figure

source('./Code/data.loading.R')          # Load data
source('./Code/risk.and.effect.R')       # risk and effect model
source('./Code/CF.R')                    # CF

run.analysis <- function(CF=FALSE, effect=FALSE, R=1, name.data.set=NULL,
                         treatment.arm=NULL, plot.cal=FALSE,
                          B=1, folds=3, alpha.reg=0.5,
                         simulation.results=NULL, application.results=NULL){
  # load data
  data <- load.data()
  dat.cc.select <- select.on.treatment(treatment.arm=treatment.arm, data=data$dat.cc, scale=FALSE) # TODO: scale data?

  # dataframe for metrics
  metrics.df <- as.data.frame(c("Calibration-in-the-large", "Eavg-for-benefit", "E50-for-benefit", "E90-for-benefit",
                                "Log-loss-for-benefit", "Brier-for-benefit",
                                "C-for-benefit"))
  names(metrics.df) <- "Metric names"

  # SIMULATION STUDY
  if (R > 0){
    results <- simulation.study(original.data=dat.cc.select, alpha.reg=alpha.reg,
                                folds=folds, R=R, plot.cal=plot.cal,
                                treatment.arm=treatment.arm,
                                metrics.df=metrics.df,
                                saved.results=simulation.results)
  }

  # CASE STUDY
  if (B > 0){
    results <- application(original.data=dat.cc.select,
                           treatment.arm=treatment.arm,
                           folds=folds, B=B,
                           metrics.df=metrics.df,
                           effect=effect, CF=CF,
                           saved.results=application.results)
  }

  return(results)
}

#######
####### PERFORM SIMULATION STUDY
#######
simulation.study <- function(original.data=NULL, alpha.reg=0.5, folds=5, R=1,
                             plot.cal=FALSE, treatment.arm=NULL,
                             metrics.df=NULL, saved.results=NULL){
  # names of models in simualtion study
  list.model.names <- list('optimal', 'suboptimal.1', 'suboptimal.2', 'suboptimal.3')

  # OBTAIN MATCHED PAIRS
  if (is.null(saved.results$dup.matched.patients)){
    true.model <- train.risk(treatment.arm=treatment.arm, Y.train=original.data$Y,
                               X.train=original.data$X, W.train=original.data$W,
                               alpha.reg=alpha.reg, folds=folds, spline=FALSE)

    # STEP 2: obtain P[Y=1]
    original.probabilities <- as.numeric(predict(true.model$final.model, type="response"))

    # STEP 3: obtain tau.hat predictions using true model
    pred.true <- predictions.risk(X.test=original.data$X,
                                    lp.model=true.model$lp.model,
                                    final.model=true.model$final.model)

    # STEP 4: perturbate models
    X.test <- c()
    X.test$X.0.test <- pred.true$X.0.test
    X.test$X.1.test <- pred.true$X.1.test
    X.test$lp.test <- pred.true$lp.test
    pred.optimal <- create.suboptimal.model(true.model=true.model,
                                            X.test=X.test,
                                            coef.W=1, coef.LP=1, coef.W.LP=1,
                                            constant=0)
    pred.suboptimal.1 <- create.suboptimal.model(true.model=true.model,
                                                 X.test=X.test,
                                                 coef.W=2, coef.LP=1, coef.W.LP=1,
                                                 constant=0)
    # correction of model 2 and 3
    if (treatment.arm=='life'){
      optimal.2 <- -0.14
      optimal.3 <- 0.53
    }
    else if (treatment.arm=='met'){
      optimal.2 <- -0.02
      optimal.3 <- 0.375
    }
    pred.suboptimal.2 <- create.suboptimal.model(true.model=true.model,
                                                 X.test=X.test,
                                                 coef.W=1, coef.LP=2, coef.W.LP=2,
                                                 constant=optimal.2)
    pred.suboptimal.3 <- create.suboptimal.model(true.model=true.model,
                                                 X.test=X.test,
                                                 coef.W=1, coef.LP=2, coef.W.LP=0.5,
                                                 constant=optimal.3)

    # save ATT of each model
    ATT.df <- data.frame(c(mean(pred.optimal$tau.hat),
                           mean(pred.suboptimal.1$tau.hat),
                           mean(pred.suboptimal.2$tau.hat),
                           mean(pred.suboptimal.3$tau.hat))*100)
    rownames(ATT.df) <- c("optimal model", "suboptimal model 1", "suboptimal model 2", "suboptimal model 3")

    # save metrics
    utils::write.table(ATT.df, file=paste0('./Results/', treatment.arm, '/ATT.', treatment.arm, '.suboptimal.optimal.txt'),
                col.names=FALSE, sep=",")
    print(ATT.df)

    # plot suboptimal models
    cat('Plot log odds \n')
    for (model.name in list.model.names){
      pred <- eval(parse(text=paste0('pred.', model.name)))
      plot.log.odss(treatment.arm=treatment.arm, model=model.name,
                    lp.test=X.test$lp.test, pred=pred)
    }

    # STEP 5: match patient pairs by distance of covariates
    cat('Match patient pairs \n')
    out.matching <- match.patients(Y=original.data$Y, W=original.data$W,
                                    X=original.data$X,
                                    p.0=pred.optimal$p.0,
                                    p.1=pred.optimal$p.1,
                                    tau.hat=pred.optimal$tau.hat)
    matched.pairs <- out.matching$df.matched.patients
    for (name.optimal in c("p.0", "p.1", "tau.hat")){
      colnames(matched.pairs)[which(colnames(matched.pairs)==name.optimal)] <- paste0(name.optimal, ".optimal")
      colnames(matched.pairs)[which(colnames(matched.pairs)==paste0("matched.", name.optimal))] <- paste0("matched.", name.optimal, ".optimal")
    }

    # add predictions of suboptimal model to dataframe
    discarded <- out.matching$discarded

    # add probabilities
    matched.pairs <- matched.pairs[order(matched.pairs$match.id), ]
    matched.pairs$p.y.1 <- original.probabilities[-discarded]

    # add probabilities for perturbed models
    for (model.name in list.model.names[-1]){
      # sort on match.id
      matched.pairs <- matched.pairs[order(matched.pairs$match.id), ]
      matched.pairs[, paste0('p.0.', model.name)] <- eval(parse(text=paste0('pred.', model.name, '$p.0[-discarded]')))
      matched.pairs[, paste0('p.1.', model.name)] <- eval(parse(text=paste0('pred.', model.name, '$p.1[-discarded]')))
      matched.pairs[, paste0('tau.hat', model.name)] <- eval(parse(text=paste0('pred.', model.name, '$tau.hat[-discarded]')))

      # sort on subclass and W
      matched.pairs <- matched.pairs[with(matched.pairs, order(subclass, 1-W)), ]

      # matched p.0 = P[Y = 1| W = 0] so the probability of an outcome given no treatment of the untreated patient
      matched.p.0 <- (1-matched.pairs$W)*matched.pairs[, paste0('p.0.', model.name)]
      matched.pairs[, paste0('matched.p.0.', model.name)] <- rep(matched.p.0[matched.p.0!=0], each=2)

      # matched p.1 = P[Y = 1| W = 1] so the probability of an outcome given no treatment of the treated patient
      matched.p.1 <- matched.pairs$W*matched.pairs[, paste0('p.1.', model.name)]
      matched.pairs[, paste0('matched.p.1.', model.name)] <- rep(matched.p.1[matched.p.1!=0], each=2)

      # matched treatment effect
      matched.pairs[, paste0('matched.tau.hat.', model.name)] <- matched.pairs[, paste0('matched.p.0.', model.name)] - matched.pairs[, paste0('matched.p.1.', model.name)]
    }

    # STEP 6: duplicate matched pairs
    dup.indices <- rep(1:nrow(matched.pairs), each=R)
    dup.matched.patients <- as.data.frame(matched.pairs[dup.indices,])

    # create unique ID's for duplicated matched pairs
    dup.matched.patients$subclass <- dup.matched.patients$subclass + (0:(R-1))*max(as.numeric(dup.matched.patients$subclass))
    dup.matched.patients$match.id <- dup.matched.patients$match.id + (0:(R-1))*max(as.numeric(dup.matched.patients$match.id))

    # STEP 7: simulate outcomes
    cat('Simulate Y \n')
    dup.matched.patients$Y.simulated <- stats::rbinom(length(dup.matched.patients$p.y.1), 1, dup.matched.patients$p.y.1)

    # STEP 8: obtain new observed treatment effect for (sub)optimal models
    cat('Update observed treatment effect of (sub)optimal model \n')
    # sort on subclass and W
    dup.matched.patients <- dup.matched.patients[with(dup.matched.patients, order(subclass, 1-W)), ]
    # update matched observed treatment effect
    updated.observed.TE <- stats::aggregate(dup.matched.patients, list(dup.matched.patients$subclass), diff)$Y.simulated
    dup.matched.patients$matched.updated.tau.obs <- rep(updated.observed.TE, each=2)
    saved.results$dup.matched.patients <- dup.matched.patients
  }
  else{
    dup.matched.patients <- saved.results$dup.matched.patients
  }

  # calibration plot limits
  limits.benefit <- list(xmax=1, ymax=1)
  if (treatment.arm=="life"){
    limits.benefit$xmin <- -0.6
    limits.benefit$ymin <- limits.benefit$xmin
  } else if (treatment.arm=="met"){
    limits.benefit$xmin <- -0.8
    limits.benefit$ymin <- limits.benefit$xmin
  }
  
  # put table in right bottom corner
  limits.table <- list(xmax=1, ymax=0, ymin=limits.benefit$ymin)
  if (treatment.arm=="life"){
    limits.table$xmin <- 0.2
  } else if (treatment.arm=="met"){
    limits.table$xmin <- 0.1
    limits.table$ymax <- -0.1
  }
  # OBTAIN MODEL METRICS AND PLOT CALIBRATION
  panel.nr.df <- data.frame(name=c("optimal", "suboptimal.1", "suboptimal.2", "suboptimal.3"), panel=c("A", "B", "C", "D"))
  cat('Obtaining model metrics... \n')
  for (model.name in list.model.names){
    cat('For', model.name, 'model \n')

    # matched df
    matched.df <- data.frame(subclass=dup.matched.patients$subclass,
                             matched.tau.hat=eval(parse(text=paste0('dup.matched.patients$matched.tau.hat.', model.name))),
                             matched.tau.obs=eval(parse(text=paste0('dup.matched.patients$matched.updated.tau.obs'))),
                             matched.p.0=eval(parse(text=paste0('dup.matched.patients$matched.p.0.', model.name))),
                             matched.p.1=eval(parse(text=paste0('dup.matched.patients$matched.p.1.', model.name))))

    # new metrics
    overall.cal.measure <- mean(matched.df$matched.tau.obs) - mean(matched.df$matched.tau.hat)
    out.E <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE, replace=FALSE)
    out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE, replace=FALSE)
    out.OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.df,
                             CI=FALSE, message=FALSE, replace=FALSE)
    metrics <- c(overall.cal.measure, out.E$Eavg.for.benefit,
                 out.E$E50.for.benefit, out.E$E90.for.benefit,
                 out.OP$Log.Loss.for.Benefit, out.OP$Brier.for.Benefit,
                 out.C$c.for.benefit)

    # save metrics
    assign(paste0(model.name, '.model.metrics'), metrics)

    # PLOT CALIBRATION PLOT
    if (plot.cal){
      cal.plot <- calibration.plot(matched.patients=out.E$matched.patients, g=5,
                                   plot.CI=FALSE, show=FALSE)
      cat("Quantiles: ", round(cal.plot$quantiles, 3), '\n')
      metric.table <- cbind(metrics.df[1:7, 1],
                            sprintf("%.3f", as.numeric(metrics)))
      plot <- cal.plot$build.plot+ggplot2::annotation_custom(gridExtra::tableGrob(metric.table,
                                                                   theme=ttheme_default(core=list(fg_params=list(hjust=1, x=1, fontsize=14),
                                                                                bg_params=list(fill=c("lightgrey", 'white'))))),
                                                                   xmin=limits.table$xmin, xmax=limits.table$xmax, ymin=limits.table$ymin, ymax=limits.table$ymax)+
        ggplot2::scale_y_continuous(labels=round(seq(from=min(limits.benefit$ymin, -1), to=max(limits.benefit$ymax, 1), length.out=5), 1),
                                    breaks=round(seq(from=min(limits.benefit$ymin, -1), to=max(limits.benefit$ymax, 1), length.out=5), 1),
                                    limits=c(limits.benefit$ymin, limits.benefit$ymax))+
        ggplot2::scale_x_continuous(labels=round(seq(from=min(limits.benefit$xmin, 0), to=max(limits.benefit$xmax, 1), length.out=5), 1),
                                    breaks=round(seq(from=min(limits.benefit$xmin, 0), to=max(limits.benefit$xmax, 1), length.out=5), 1),
                                    limits=c(limits.benefit$xmin, limits.benefit$xmax))+
        ggplot2::theme(plot.title=element_text(hjust=0.5))+
        ggplot2::annotate(geom="label", x=limits.benefit$xmin, y=limits.benefit$ymax, size=15, fontface=2, fill="white", label.size=NA,
                                                          label=panel.nr.df[panel.nr.df$name==model.name, "panel"])+
        ggplot2::theme_light(base_size=25)+                   # increase font size
        ggplot2::theme(axis.title.x=element_blank(), axis.title.y=element_blank())
      save(plot, file=paste0('./Results/', treatment.arm, '/', model.name, '.simulation.calibration.plot.Rdata'))
    }
  }

  metric.values <- cbind(c(metrics.df[1:7,1]),
                           c(optimal.model.metrics),
                           c(suboptimal.1.model.metrics),
                           c(suboptimal.2.model.metrics),
                           c(suboptimal.3.model.metrics))
  utils::write.table(metric.values,
              file=paste0('./Results/', treatment.arm, '/metrics.', treatment.arm, '.suboptimal.optimal.txt'),
              row.names=FALSE, sep=",")

  return(saved.results)
}

#######
####### PLOT LOG ODDS
#######
plot.log.odss <- function(treatment.arm=NULL, model=NULL, lp.test=NULL, pred=NULL){
  plot <- ggplot(data=data.frame(lp.test=lp.test, p.0=qlogis(pred$p.0),
                                 p.1=qlogis(pred$p.1)), aes(x=lp.test), show.legend=TRUE)
  x.lab.plot <- ""
  y.lab.plot <- ""
  plot <- plot+ggplot2::labs(x=x.lab.plot, y=y.lab.plot, color=" ")   # axis names
  plot <- plot+ggplot2::geom_line(aes(y=p.0), color="blue", size=1)   # draw p.0
  plot <- plot+ggplot2::geom_line(aes(y=p.1), color="red", size=1)    # draw p.1
  plot <- plot+ggplot2::theme_light(base_size=25)                     # increase font size
  plot <- plot+ggplot2::theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  panel.nr.df <- data.frame(name=c("optimal", "suboptimal.1", "suboptimal.2", "suboptimal.3"), panel=c("A", "B", "C", "D"))
  if (treatment.arm=="life"){
    ymin <- -10
    ymax <- 10
    xmin <- 47
    xmax <- 57
  }
  else if (treatment.arm=="met"){
    ymin <- -10
    ymax <- 10
    xmin <- 45
    xmax <- 53
  }
  plot <- plot+ggplot2::scale_y_continuous(labels=round(seq(ymin, ymax, 2), 0), 
                                           breaks=round(seq(ymin, ymax, 2), 0), 
                                           limits=c(ymin, ymax)) # axis limits
  plot <- plot+ggplot2::scale_x_continuous(labels=round(seq(xmin, xmax, 2), 0),
                                           breaks=round(seq(xmin, xmax, 2), 0), 
                                           limits=c(xmin, xmax)) # axis limits
  plot <- plot+ggplot2::annotate(geom="label", x=xmin, y=ymax, size=15, fontface=2, fill="white", label.size=NA,
                        label=panel.nr.df[panel.nr.df$name==model, "panel"])
  save(plot, file=paste0('./Results/', treatment.arm, '/', 'log.odds.', model, '.Rdata'))
}

######
###### OBTAIN METRIC VALUES FOR CASE STUDY
######
metric.values.three.models <- function(boot=0, pred=NULL,
                                       treatment.arm=NULL, metrics.df=NULL,
                                       Y=NULL, X=NULL, W=NULL, matched.test=NULL){
  # limits calibration plot
  limits.benefit <- list(xmin=-0.5,
                              xmax=1,
                              ymin=-1.5,
                              ymax=1.5)
  xmin.fact <- 1.57
  ymax.fact <- 1.3
  limits.table <- list(xmin=limits.benefit$xmin-limits.benefit$xmin*xmin.fact,
                       xmax=limits.benefit$xmax,
                       ymin=limits.benefit$ymin,
                       ymax=limits.benefit$ymax-limits.benefit$ymax*ymax.fact)

  ##### OBTAIN METRIC VALUES
  for (method in list('risk', 'effect', 'CF')){
    pred.method <- data.frame(tau.hat=eval(parse(text=paste0('pred$', method, '.tau.hat'))),
                              p.0=eval(parse(text=paste0('pred$', method, '.p.0'))),
                              p.1=eval(parse(text=paste0('pred$', method, '.p.1'))))

    if (boot > 0){
      # bootstrap matched pairs
      matched.patients <- eval(parse(text=paste0('matched.test$matched.', method)))

      # duplicate subclass.id
      subclass.IDs <- unique(matched.patients$subclass)
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      dup.subclass.IDs <- c()
      for (i in sample.subclass){
        dup.subclass.IDs <- c(dup.subclass.IDs, matched.patients[matched.patients$subclass==i, 'match.id'])
      }
      matched.df <- dplyr::slice(matched.patients, dup.subclass.IDs)
    }
    else{
      # match patient pairs by distance of covariates
      matched.df <- match.patients(Y=Y, W=W,
                                     X=X,
                                     p.0=pred.method$p.0,
                                     p.1=pred.method$p.1,
                                     tau.hat=pred.method$tau.hat)$df.matched.patients
    }
    assign(paste0('matched.', method), matched.df)

    overall.cal.measure <- mean(matched.df$matched.tau.obs) - mean(matched.df$matched.tau.hat)
    out.E <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE, replace=FALSE)
    out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE, replace=FALSE)
    out.OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.df,
                             CI=FALSE, message=FALSE, replace=FALSE)
    metric.values <- c(overall.cal.measure, out.E$Eavg.for.benefit,
                       out.E$E50.for.benefit, out.E$E90.for.benefit,
                       out.OP$Log.Loss.for.Benefit, out.OP$Brier.for.Benefit,
                       out.C$c.for.benefit)
    assign(paste0('metric.values.', method), metric.values)

    # save metric values
    if (boot == -1){
      data.name <- "train"
      results.metrics <- metrics.df
      results.metrics[, ncol(results.metrics)+1] <- metric.values
      colnames(results.metrics)[ncol(results.metrics)] <- 'original train'
    }
    else if (boot == 0){
      data.name <- "test"
      results.metrics <- metrics.df
      results.metrics[, ncol(results.metrics)+1] <- metric.values
      colnames(results.metrics)[ncol(results.metrics)] <- 'original test'
    }
    else if (boot > 0){
      results.metrics <- as.data.frame(read.table(file=paste0('./Results/', treatment.arm, '/metrics.', treatment.arm, '.', method, '.txt'),
                                                  header = TRUE, sep = ",", dec = "."))
      results.metrics[, ncol(results.metrics)+1] <- metric.values
      colnames(results.metrics)[ncol(results.metrics)] <- paste0('CB.', boot)
    }
    utils::write.table(results.metrics,
                file=paste0('./Results/', treatment.arm, '/metrics.', treatment.arm, '.', method, '.txt'),
                row.names=FALSE, sep=",")

    # plot calibration for training and test set
    if (boot == -1 | boot == 0){
      cal.plot <- calibration.plot(matched.patients=out.E$matched.patients, g=5,
                                   plot.CI=TRUE, show=FALSE)
      assign(paste0('cal.plot.', method), cal.plot$build.plot)
    } else{
      cal.plot.risk <- NA
      cal.plot.effect <- NA
      cal.plot.CF <- NA
    }
  }
  metric.values <- data.frame(risk=metric.values.risk, effect=metric.values.effect, CF=metric.values.CF)

  return(list(matched.risk=matched.risk, matched.effect=matched.effect,
              matched.CF=matched.CF, metric.values=metric.values,
              cal.plot.risk=cal.plot.risk, cal.plot.effect=cal.plot.effect,
              cal.plot.CF=cal.plot.CF))
}

#####
##### APPLICATION
#####
application <- function(original.data=NULL, treatment.arm=NULL, folds=5, B=B,
                        metrics.df=NULL, effect=FALSE, CF=FALSE, saved.results=NULL){
  # split training and test set
  n <- length(original.data$Y)

  # split train and test set
  split.ratio <- 0.7
  ind.train <- 1:(floor(n*split.ratio))
  ind.test <- (floor(n*split.ratio)+1):n

  # obtain training set
  training.data <- original.data$data[ind.train,]
  select.train <- select.on.treatment(treatment.arm=treatment.arm, data=training.data, scale=FALSE)
  Y.train <- select.train$Y
  X.train <- select.train$X
  W.train <- select.train$W

  #####
  ##### TRAIN MODELS on training set
  #####
  cat('ORIGINAL SAMPLE \n')
  trained <- train.three.models(fold=fold, treatment.arm=treatment.arm,
                                Y=Y.train, X=X.train, W=W.train,
                                alpha.reg=alpha.reg, folds=folds, effect=effect, CF=CF)

  # test models on training fold
  test.on.train <- test.three.models(X=X.train, trained.risk=trained$risk,
                                     treatment.arm=treatment.arm,
                                     trained.effect=trained$effect,
                                     trained.CF=trained$CF,
                                     effect=effect, CF=CF)

  matched.train <- metric.values.three.models(boot=-1, pred=test.on.train,
                                                        treatment.arm=treatment.arm,
                                                        metrics.df=metrics.df,
                                                        Y=Y.train, X=X.train,
                                                        W=W.train,
                                                        matched.test=NULL)

  #####
  ##### TEST MODELS
  #####
  # obtain training set
  test.data <- original.data$data[ind.test,]
  select.test <- select.on.treatment(treatment.arm=treatment.arm, data=test.data, scale=FALSE)
  Y.test <- select.test$Y
  X.test <- select.test$X
  W.test <- select.test$W
  cat("Number of observations for test set: ", length(Y.test), "\n")

  # test models on test fold
  test <- test.three.models(X=X.test, trained.risk=trained$risk,
                            treatment.arm=treatment.arm,
                            trained.effect=trained$effect,
                            trained.CF=trained$CF,
                            effect=effect, CF=CF)
  # obtain metric values
  matched.test <- metric.values.three.models(boot=0, pred=test,
                                             treatment.arm=treatment.arm,
                                             metrics.df=metrics.df,
                                             Y=Y.test, X=X.test,
                                             W=W.test,
                                             matched.test=NULL)

  #####
  ##### BOOTSTRAP TEST SET
  #####
  cat('BOOTSTRAP SAMPLES \n')
  for (boot in 1:B){
    # use matched test set (don't match again)
    bootstrap <- metric.values.three.models(boot=boot, pred=test,
                               treatment.arm=treatment.arm, metrics.df=metrics.df,
                               Y=Y.test, X=X.test, W=W.test,
                               matched.test=matched.test)
  }

  # obtain CI
  results.matrix <- as.data.frame(metrics.df[1:7, ])
  for (method in c('risk', 'effect', 'CF')){
    results.metrics <- as.matrix(read.table(file=paste0('./Results/', treatment.arm, '/metrics.', treatment.arm, '.', method, '.txt'),
                                            header = TRUE, row.names=1, sep = ",", dec = "."))
    original.metrics <- as.vector(results.metrics[, grepl("original.test", colnames(results.metrics))])
    CB.metrics <- as.matrix(results.metrics[, grepl("CB.", colnames(results.metrics))])
    assign(paste0('original.', method), original.metrics)
    assign(paste0('CB.', method), CB.metrics)
    if (ncol(CB.metrics) > 0){
      # construct confidence interval
      CI.lower <- c()
      CI.upper <- c()
      for (row in 1:length(original.metrics)){
        CI.lower.new <- as.numeric(quantile(CB.metrics[row,], 0.025))
        CI.upper.new <- as.numeric(quantile(CB.metrics[row,], 0.975))
        CI.lower <- c(CI.lower, CI.lower.new)
        CI.upper <- c(CI.upper, CI.upper.new)
      }

      # add to results table
      results.matrix <- cbind(results.matrix, as.vector(original.metrics),
                              paste0('[', sprintf('%.3f', CI.lower),
                                    '; ', sprintf('%.3f', CI.upper), ']'))
    }
  }
  colnames(results.matrix) <- c("Metric names", "risk", "risk.CI", "effect", "effect.CI", "CF", "CF.CI")

  print(results.matrix)

  # plot calibration
  for (method in c('risk', 'effect', 'CF')){
    for (data.name in c('train', 'test')){
      # limits calibration plot
      limits.benefit <- list(xmin=-0.5,
                             xmax=1,
                             ymin=-1.5,
                             ymax=1.5)
      cal.plot <- eval(parse(text=paste0('matched.', data.name, '$cal.plot.', method)))
      cal.plot <- cal.plot+ggplot2::theme_light(base_size=25)+                   # increase font size
        ggplot2::theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        ggplot2::scale_y_continuous(labels=seq(limits.benefit$ymin, limits.benefit$ymax, length.out=7),
                                    breaks=seq(limits.benefit$ymin, limits.benefit$ymax, length.out=7),
                                    limits=c(limits.benefit$ymin, limits.benefit$ymax+0.1))+
        ggplot2::scale_x_continuous(labels=seq(limits.benefit$xmin, limits.benefit$xmax, length.out=7),
                                    breaks=seq(limits.benefit$xmin, limits.benefit$xmax, length.out=7),
                                    limits=c(limits.benefit$xmin, limits.benefit$xmax+0.1))+
        ggplot2::theme(plot.title=element_text(hjust=0.5))
        
      if (data.name== 'test'){
        metric.table <- cbind(metrics.df[1:7, 1],
                              sprintf("%.3f", as.numeric(results.matrix[, method])),
                              results.matrix[, paste0(method, ".CI")])
        cal.plot <- cal.plot+ggplot2::annotation_custom(gridExtra::tableGrob(metric.table,
                                               theme=ttheme_default(core=list(fg_params=list(hjust=1, x=1, fontsize=14),
                                                                              bg_params=list(fill=c("lightgrey", 'white'))))),
                                               xmin=0.05, xmax=limits.benefit$xmax, 
                                               ymin=limits.benefit$ymin, ymax=-0.55)
      }

      if (data.name=='train'){
        panel.nr.df <- data.frame(name=c("risk", "effect", "CF"), panel=c("A", "C", "E"))
        cal.plot.appendix <- cal.plot+ggplot2::annotate(geom="label", x=limits.benefit$xmin, y=limits.benefit$ymax, size=15, fontface=2, fill="white", label.size=NA,
                                      label=panel.nr.df[panel.nr.df$name==method, "panel"])
      } else if (data.name=='test'){
        panel.nr.df <- data.frame(name=c("risk", "effect", "CF"), panel=c("B", "D", "F"))
        cal.plot.appendix <- cal.plot+ggplot2::annotate(geom="label", x=limits.benefit$xmin, y=limits.benefit$ymax, size=15, fontface=2, fill="white", label.size=NA,
                                      label=panel.nr.df[panel.nr.df$name==method, "panel"])
      }
      plot <- cal.plot.appendix
      save(plot, file=paste0('./Results/', treatment.arm, '/', data.name, '.', method, '.case.study.appendix.calibration.plot.Rdata'))

      if (data.name=='test'){
        panel.nr.df <- data.frame(name=c("risk", "effect", "CF"), panel=c("A", "B", "C"))
        cal.plot <- cal.plot+ggplot2::annotate(geom="label", x=limits.benefit$xmin, y=limits.benefit$ymax, size=15, fontface=2, fill="white", label.size=NA,
                                      label=panel.nr.df[panel.nr.df$name==method, "panel"])
      }
      plot <- cal.plot
      save(plot, file=paste0('./Results/', treatment.arm, '/', data.name, '.', method, '.case.study.calibration.plot.Rdata'))
    }
  }

  return(saved.results)
}

#####
##### TRAIN THREE MODELS
#####
train.three.models <- function(fold=0, treatment.arm=NULL, Y=NULL, X=NULL, W=NULL,
                               alpha.reg=NULL, folds=NULL, effect=FALSE, CF=FALSE){
  ##### TRAIN RISK MODEL
  # cat('train risk model    on fold number', fold, '\n')
  trained.risk <- train.risk(treatment.arm=treatment.arm, Y.train=Y,
                             X.train=X, W.train=W,
                             alpha.reg=alpha.reg, folds=folds, spline=TRUE)

  ##### TRAIN EFFECT MODEL
  if (effect){
    # cat('train effect model  on fold number', fold, '\n')
    trained.effect <- train.effect(treatment.arm=treatment.arm, Y.train=Y, X.train=X, W.train=W,
                                   folds=folds, alpha.reg=alpha.reg, penalize=TRUE, print=FALSE)
  }
  else{
    trained.effect <- trained.risk
  }

  ##### TRAIN CAUSAL FOREST
  if (CF){
    # cat('train causal forest on fold number', fold, '\n')
    trained.CF <- train.CF(Y.train=Y, X.train=X, W.train=W, results=FALSE, tune=TRUE)
  }
  else{
    trained.CF <- trained.risk
  }

  return(list(risk=trained.risk, effect=trained.effect, CF=trained.CF))
}

#####
##### TEST THREE MODELS
#####
test.three.models <- function(X=X, trained.risk=NULL, treatment.arm=NULL,
                              trained.effect=NULL, trained.CF=NULL,
                              effect=FALSE, CF=FALSE){
  ##### OBTAIN PREDICTIONS FOR RISK MODEL
  pred.risk <- predictions.risk(X.test=X, lp.model=trained.risk$lp.model,
                                final.model=trained.risk$final.model)

  if (effect){
    ##### OBTAIN PREDICTIONS FOR EFFECT MODEL
    pred.effect <- predictions.effect(treatment.arm=treatment.arm, X=X,
                                      penalized=TRUE,
                                      final.model=trained.effect$final.model)
  }
  else{
    pred.effect <- pred.risk
  }

  if (CF){
    ##### OBTAIN PREDICTIONS FOR CAUSAL FOREST
    pred.CF <- predictions.CF(X.test=X, final.model=trained.CF)
  }
  else{
    pred.CF <- pred.risk
  }

  return(list(risk.tau.hat=pred.risk$tau.hat, risk.p.0=pred.risk$p.0, risk.p.1=pred.risk$p.1,
              effect.tau.hat=pred.effect$tau.hat, effect.p.0=pred.effect$p.0, effect.p.1=pred.effect$p.1,
              CF.tau.hat=pred.CF$tau.hat, CF.p.0=pred.CF$p.0, CF.p.1=pred.CF$p.1))
}
