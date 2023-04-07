#######
####### COMMON FUNCTIONS
#######
source('./Code/data.loading.R')          # Load data
source('./Code/risk.and.effect.R')       # risk and effect model
source('./Code/CF.R')                    # CF
library(dplyr)
options(dplyr.summarise.inform=FALSE)

run.analysis <- function(CF=FALSE, effect=FALSE, R=1, name.data.set=NULL,
                         treatment.arm=NULL, plot.cal=FALSE, 
                         random.matching=FALSE, match.on.benefit=FALSE, 
                         match.on.covariates=TRUE,
                         B=1, folds=3, alpha.reg=0.5,
                         simulation.results=NULL, application.results=NULL){
  # load data
  data <- load.data()
  
  # SIMULATION STUDY
  if (R > 0){
    # dataframe for metrics
    metrics.df <- data.frame(names=c("RMSE", "CitL", "Eavg-B", 
                                  "E50-B", "E90-B",
                                  "CE-B", "Brier-B",
                                  "C-B"))
    
    results <- simulation.study(data=data$dat.cc, alpha.reg=alpha.reg,
                                folds=folds, R=R, plot.cal=plot.cal,
                                random.matching=random.matching,
                                match.on.benefit=match.on.benefit,
                                match.on.covariates=match.on.covariates,
                                treatment.arm=treatment.arm,
                                metrics.df=metrics.df,
                                saved.results=simulation.results)
  }

  # CASE STUDY
  if (B > 0){
    # dataframe for metrics
    metrics.df <- data.frame(names=c("CitL", "Eavg-B", "E50-B", "E90-B",
                                  "CE-B", "Brier-B",
                                  "C-B"))
    
    results <- application(data=data$dat.cc,
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
simulation.study <- function(data=NULL, alpha.reg=0.5, folds=5, R=1,
                             plot.cal=FALSE, random.matching=FALSE,
                             match.on.benefit=FALSE, match.on.covariates=TRUE,
                             treatment.arm=NULL,
                             metrics.df=NULL, saved.results=NULL){
  # load data
  original.data <- select.on.treatment(treatment.arm=treatment.arm, data=data, scale=TRUE)
  
  # names of models in simualtion study
  list.model.names <- list('optimal', 'suboptimal.1', 'suboptimal.2', 'suboptimal.3', 'suboptimal.4')

  # OBTAIN MATCHED PAIRS
  if (is.null(saved.results$dup.matched.patients)){
    # STEP 1: train effect model penalized LASSO
    # print(summary(stats::glm(Y ~ W+X, data=original.data, family="binomial")))
    true.model <- train.effect(treatment.arm=treatment.arm, 
                               Y.train=original.data$Y,
                               X.train=original.data$X, 
                               W.train=original.data$W,
                               alpha.reg=0, folds=folds, 
                               print=FALSE, simulation=TRUE)
    
    # STEP 2: obtain P[Y=1]
    original.probabilities <- as.numeric(predict(true.model$final.model, 
                                                 newx=true.model$X.full, 
                                                 type="response"))
    
    # STEP 3: obtain tau.hat predictions using true model
    pred.true <- predictions.effect(treatment.arm=treatment.arm,
                                    X=original.data$X,
                                    penalized=TRUE,
                                    final.model=true.model$final.model)
    # train.risk: final.X.test=original.data$X, lp.model=true.model$lp.model
    
    # TODO: how to do perturbe models for LASSO?
    # STEP 4: perturbate models
    pred.optimal <- create.suboptimal.model(true.model=true.model,
                                            X=original.data$X,
                                            coef.W=1, coef.X=1, coef.X.W=1,
                                            c=0, c.W=0)
    pred.suboptimal.1 <- create.suboptimal.model(true.model=true.model,
                                                 X=original.data$X,
                                                 coef.W=2, coef.X=1, coef.X.W=1,
                                                 c=0, c.W=0)
    pred.suboptimal.2 <- create.suboptimal.model(true.model=true.model,
                                                 X=original.data$X,
                                                 coef.W=0.5, coef.X=1, coef.X.W=1,
                                                 c=0, c.W=0)
    # correction of model 3 and 4
    if (treatment.arm=='life'){
      c.3 <- 0
      c.W.3 <- -0.195
      c.W.4 <- -0.19
    }
    else if (treatment.arm=='met'){
      c.3 <- 0
      c.W.3 <- -0.085
      c.W.4 <- -0.16
    }
    pred.suboptimal.3 <- create.suboptimal.model(true.model=true.model,
                                                 X=original.data$X,
                                                 coef.W=1, coef.X=2, coef.X.W=1,
                                                 c=c.3, c.W=c.W.3)
    pred.suboptimal.4 <- create.suboptimal.model(true.model=true.model,
                                                 X=original.data$X,
                                                 coef.W=1, coef.X=1, coef.X.W=3,
                                                 c=0, c.W=c.W.4)

    # save ATT of each model
    ATT.df <- data.frame(c(mean(pred.optimal$tau.hat),
                           mean(pred.suboptimal.1$tau.hat),
                           mean(pred.suboptimal.2$tau.hat),
                           mean(pred.suboptimal.3$tau.hat),
                           mean(pred.suboptimal.4$tau.hat))*100)
    rownames(ATT.df) <- c("optimal model", "suboptimal model 1", 
                          "suboptimal model 2", "suboptimal model 3",
                          "suboptimal model 4")
    
    print(data.frame(names=c("optimal model", "suboptimal model 1", 
                             "suboptimal model 2", "suboptimal model 3",
                             "suboptimal model 4"),
                     values= c(mean(pred.optimal$p.0),
                                mean(pred.suboptimal.1$p.0),
                                mean(pred.suboptimal.2$p.0),
                                mean(pred.suboptimal.3$p.0),
                                mean(pred.suboptimal.4$p.0))*100))

    # save ATT
    utils::write.table(ATT.df, file=paste0('./Results/', treatment.arm, '/ATT.', treatment.arm, '.suboptimal.optimal.txt'),
                col.names=FALSE, sep=",")
    print(ATT.df)
    
    # plot suboptimal models (TODO: plot log odds for LASSO?)
    cat('Plot log odds \n')
    if (match.on.covariates){
      for (model.name in list.model.names){
        pred <- eval(parse(text=paste0('pred.', model.name)))
        
        # plot on probability scale
        plot.p <- plot.perturbations(treatment.arm=treatment.arm,
                                    model=model.name, pred=pred, 
                                    log.scale=FALSE)
        assign(paste0("p.", model.name), plot.p)
        
        # plot on log odds scale
        plot.log.odds <- plot.perturbations(treatment.arm=treatment.arm,
                                    model=model.name, pred=pred, 
                                    log.scale=TRUE)
        assign(paste0("log.odds.", model.name), plot.log.odds)
      }
      ggplot2::ggsave(filename=paste0('./Results/Simulation/p.', treatment.arm, '.png'),
                      plot=ggpubr::annotate_figure(ggpubr::ggarrange(p.optimal, p.suboptimal.1,
                                                                     p.suboptimal.3, p.suboptimal.4,
                                                                     nrow=2, ncol=2, align="h", labels=c("A", "B", "C", "D", "E")),
                                                   left=ggpubr::text_grob("Probability of potential outcome under treatment", rot=90, size=20),
                                                   bottom=ggpubr::text_grob("Probability of potential outcome under control treatment", size=20)),
                      width=10, height=10, dpi=300)
      ggplot2::ggsave(filename=paste0('./Results/Simulation/log.odds.', treatment.arm, '.png'),
                       plot=ggpubr::annotate_figure(ggpubr::ggarrange(log.odds.optimal, log.odds.suboptimal.1,
                                                                      log.odds.suboptimal.3, log.odds.suboptimal.4,
                                                                      nrow=2, ncol=2, align="h", labels=c("A", "B", "C", "D", "E")),
                                                    left=ggpubr::text_grob("Log odds of potential outcome under treatment", rot=90, size=20),
                                                    bottom=ggpubr::text_grob("Log odds of potential outcome under control treatment", size=20)),
                       width=10, height=10, dpi=300)
    }
    
    # STEP 5: match patient pairs randomly, by distance of treatment effect or by distance of covariates
    cat('Match patient pairs \n')
    if (random.matching){
      # randomly match patients
      out.matching <- HTEPredictionMetrics::match.patients(Y=original.data$Y, W=original.data$W,
                                                           X=rep(1, length(original.data$Y)),
                                                           p.0=pred.optimal$p.0,
                                                           p.1=pred.optimal$p.1,
                                                           tau.hat=pred.optimal$tau.hat)
    } else if (match.on.benefit){
      # match predicted treatment effect
      out.matching <- HTEPredictionMetrics::match.patients(Y=original.data$Y, W=original.data$W,
                                                           X=pred.optimal$tau.hat,
                                                           p.0=pred.optimal$p.0,
                                                           p.1=pred.optimal$p.1,
                                                           tau.hat=pred.optimal$tau.hat)
    } else if (match.on.covariates){
      # match using covariates
      parameter.set <- list("nearest", "mahalanobis")
      out.matching <- HTEPredictionMetrics::match.patients(Y=original.data$Y, W=original.data$W,
                                                           X=original.data$X,
                                                           p.0=pred.optimal$p.0,
                                                           p.1=pred.optimal$p.1,
                                                           tau.hat=pred.optimal$tau.hat,
                                                           measure=parameter.set[[1]],
                                                           distance=parameter.set[[2]])
    }
    matched.pairs <- out.matching$df.matched.patients
    
    for (name.optimal in c("p.0", "p.1", "tau.hat")){
      colnames(matched.pairs)[which(colnames(matched.pairs)==name.optimal)] <- paste0(name.optimal, ".optimal")
      colnames(matched.pairs)[which(colnames(matched.pairs)==paste0("matched.", name.optimal))] <- paste0("matched.", name.optimal, ".optimal")
    }
    
    # add predictions of suboptimal model to dataframe
    discarded <- out.matching$discarded
    
    # omit patients who were discarded
    matched.pairs <- matched.pairs[order(matched.pairs$match.id), ]
    
    # add probabilities
    matched.pairs$p.y.1 <- original.probabilities[-discarded]
    
    # add probabilities for perturbed models
    for (model.name in list.model.names[-1]){
      # sort on match.id
      matched.pairs <- matched.pairs[order(matched.pairs$match.id), ]
      matched.pairs[, paste0('p.0.', model.name)] <- eval(parse(text=paste0('pred.', model.name, '$p.0[-discarded]')))
      matched.pairs[, paste0('p.1.', model.name)] <- eval(parse(text=paste0('pred.', model.name, '$p.1[-discarded]')))
      matched.pairs[, paste0('tau.hat.', model.name)] <- eval(parse(text=paste0('pred.', model.name, '$tau.hat[-discarded]')))

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
  if (treatment.arm=="life"){
    limits.benefit <- list(xmin=-0.6, xmax=1, 
                           ymin=-0.6, ymax=1)
  } else if (treatment.arm=="met"){
    limits.benefit <- list(xmin=-0.8, xmax=1, 
                           ymin=-0.8, ymax=1)
  }
  
  # put table in corner
  if (treatment.arm=="life"){
    limits.table <- list(xmin=-0.8, xmax=0,
                         ymin=0.4, ymax=1)
  } else if (treatment.arm=="met"){
    limits.table <- list(xmin=-1.2, xmax=0, 
                         ymin=0.3, ymax=1)
  }
  
  # OBTAIN MODEL METRICS AND PLOT CALIBRATION
  cat('Obtaining model metrics... \n')
  quantiles.calibration.plot <- c()
  for (model.name in list.model.names){
    cat('For', model.name, 'model \n')

    # matched df
    matched.df <- data.frame(subclass=dup.matched.patients$subclass,
                             Y=eval(parse(text=paste0('dup.matched.patients$Y'))),
                             W=eval(parse(text=paste0('dup.matched.patients$W'))),
                             matched.tau.hat=eval(parse(text=paste0('dup.matched.patients$matched.tau.hat.', model.name))),
                             matched.tau.obs=eval(parse(text=paste0('dup.matched.patients$matched.updated.tau.obs'))),
                             matched.p.0=eval(parse(text=paste0('dup.matched.patients$matched.p.0.', model.name))),
                             matched.p.1=eval(parse(text=paste0('dup.matched.patients$matched.p.1.', model.name))))
    matched.df$P <- (1-matched.df$W)*matched.df$matched.p.0+matched.df$W*matched.df$matched.p.1
    
    # conventional metrics
    # loess.calibrate <- loess(Y ~ P, data=matched.df)
    # P.calibrate <- predict(loess.calibrate, newdata=matched.df$P)
    # ICI <- mean(abs(P.calibrate - matched.df$P))
    # E50 <- median(abs(P.calibrate - matched.df$P))
    # E90 <- quantile(abs(P.calibrate - matched.df$P), probs=0.9)
    # CE <- -1/nrow(matched.df)*sum(matched.df$Y*log(matched.df$P)+
    #                                 (1-matched.df$Y)*log(1-matched.df$P))
    # Brier <- sum((matched.df$Y-matched.df$P)^2)/nrow(matched.df)
    # Cindex <- Hmisc::rcorr.cens(matched.df$P, matched.df$Y)["C Index"]
    
    # new metrics
    RMSE <- sqrt(mean((pred.true$tau.hat - eval(parse(text=paste0('pred.', model.name, '$tau.hat'))))^2))
    overall.cal.measure <- mean(matched.df$matched.tau.obs) - mean(matched.df$matched.tau.hat)
    
    out.E <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE)
    out.OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.df,
                             CI=FALSE, message=FALSE)
    out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                                                 CI=FALSE, message=FALSE)
    metrics <- c(RMSE, overall.cal.measure, 
                 out.E$Eavg.for.benefit, 
                 out.E$E50.for.benefit, 
                 out.E$E90.for.benefit,
                 out.OP$Cross.entropy.for.benefit, 
                 out.OP$Brier.for.benefit,
                 out.C$C.for.benefit)
    
    # save metrics
    assign(paste0(model.name, '.model.metrics'), metrics)

    # PLOT CALIBRATION PLOT
    if (plot.cal){
      cal.plot <- HTEPredictionMetrics::calibration.plot(matched.patients=out.E$matched.patients,
                                                         g=5, plot.CI=FALSE, 
                                                         show=FALSE)
      cat("Quantiles: ", round(cal.plot$quantiles, 3), '\n')
      quantiles.calibration.plot <- rbind(quantiles.calibration.plot, 
                                          c(model.name, round(cal.plot$quantiles, 3)))
      metric.table <- cbind(metrics.df[, "names"], sprintf("%.3f", as.numeric(metrics)))
      
      # make histogram table
      matched.df$colors <- ifelse(matched.df$matched.tau.obs==-1, "harm", 
                                  ifelse(matched.df$matched.tau.obs==0, "noeffect", "benefit"))
      matched.df <- matched.df[order(matched.df$matched.tau.hat), ]
      matched.df$tau.obs.cut <- cut(matched.df$matched.tau.hat, 
                                    breaks=seq(min(matched.df$matched.tau.hat), 
                                               max(matched.df$matched.tau.hat), 
                                               by=0.01))
      
      # count number of observations in each group
      data.cut <- matched.df %>%
        dplyr::group_by(colors, tau.obs.cut) %>%
        dplyr::summarise(n=n()) %>%
        dplyr::mutate(tau.obs.cut=ifelse(is.na(tau.obs.cut), 0, tau.obs.cut))

      # set new x-axis range for histogram
      nr.empty.neg <- (limits.benefit$xmin-min(matched.df$matched.tau.hat))/0.01
      nr.empty.pos <- max(data.cut$tau.obs.cut)+(limits.benefit$xmax-max(matched.df$matched.tau.hat))/0.01
      
      # make histogram figure
      dist <- ggplot2::ggplot(data=data.cut, 
                           mapping=ggplot2::aes(x=.data$tau.obs.cut,
                                                y=.data$n,
                                                fill=.data$colors))+
             ggplot2::geom_bar(stat="identity", alpha=1, show.legend=FALSE, 
                               position="dodge2")+
             ggplot2::scale_x_continuous(breaks=c(nr.empty.neg, nr.empty.pos),
                                         labels=c(nr.empty.neg, nr.empty.pos),
                                         limits=c(nr.empty.neg, nr.empty.pos))+
             ggplot2::theme_void()+
             ggplot2::scale_fill_manual("legend", values=c("harm"="red", "noeffect"="darkgrey", "benefit"= "blue"))
      assign(paste0("simulation.dist.", model.name), dist)
      
      # add metrics to calibration plot
      plot <- cal.plot$build.plot+ggplot2::annotation_custom(gridExtra::tableGrob(metric.table,
                                                                   theme=gridExtra::ttheme_default(core=list(padding=ggplot2::unit(c(1.5,1.5), "mm"),
                                                                               fg_params=list(hjust=1, x=1, fontsize=12),
                                                                               bg_params=list(fill="white")))),
                                                             xmin=limits.table$xmin, 
                                                             xmax=limits.table$xmax, 
                                                             ymin=limits.table$ymin, 
                                                             ymax=limits.table$ymax)+
        ggplot2::scale_y_continuous(labels=round(seq(from=limits.benefit$ymin, to=limits.benefit$ymax, length.out=5), 1),
                                    breaks=round(seq(from=limits.benefit$ymin, to=limits.benefit$ymax, length.out=5), 1),
                                    limits=c(limits.benefit$ymin, limits.benefit$ymax))+
        ggplot2::scale_x_continuous(labels=round(seq(from=limits.benefit$xmin, to=limits.benefit$xmax, length.out=5), 1),
                                    breaks=round(seq(from=limits.benefit$xmin, to=limits.benefit$xmax, length.out=5), 1),
                                    limits=c(limits.benefit$xmin, limits.benefit$xmax))+
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
        ggplot2::theme_light(base_size=25)                   # increase font size
      assign(paste0("simulation.plot.", model.name), plot)
    }
  }

  # combine calibration plots
  ggplot2::ggsave(file=paste0('./Results/Simulation/simulation.', treatment.arm,
                              ifelse(random.matching, ".random", ""),
                              ifelse(match.on.benefit, ".benefit", ""),
                              ifelse(match.on.covariates, ".covariates", ""), '.png'),
                  plot=ggpubr::annotate_figure(ggpubr::ggarrange(simulation.plot.optimal+
                                                                   ggpubr::rremove("xlab")+
                                                                   ggpubr::rremove("ylab"),
                                                                 simulation.plot.suboptimal.1+
                                                                   ggpubr::rremove("xlab")+
                                                                   ggpubr::rremove("ylab"),
                                                                 simulation.dist.optimal,
                                                                 simulation.dist.suboptimal.1,
                                                                 simulation.plot.suboptimal.3+
                                                                   ggpubr::rremove("xlab")+
                                                                   ggpubr::rremove("ylab"),
                                                                 simulation.plot.suboptimal.4+
                                                                   ggpubr::rremove("xlab")+
                                                                   ggpubr::rremove("ylab"),
                                                                 simulation.dist.suboptimal.3,
                                                                 simulation.dist.suboptimal.4,
                                                                 nrow=4, ncol=2, align="hv",
                                                                 labels=c("A", "B", "", "", "C", "D", "", ""),
                                                                 heights=c(1, 0.4, 1, 0.4)),
                                               left=ggpubr::text_grob("Observed pairwise treatment effect", rot=90, size=20),
                                               bottom=ggpubr::text_grob("Predicted pairwise treatment effect", size=20)),
                  width=10, height=12, dpi=300)

  # save plot of calibration for model underestimating ATE
  ggplot2::ggsave(file=paste0('./Results/Simulation/appendix.simulation.underestimateATE.', treatment.arm,
                              ifelse(random.matching, ".random", ""),
                              ifelse(match.on.benefit, ".benefit", ""),
                              ifelse(match.on.covariates, ".covariates", ""), '.png'),
                  plot=ggpubr::annotate_figure(ggpubr::ggarrange(simulation.plot.suboptimal.2+
                                           ggplot2::theme_light(base_size=14)+
                                           ggpubr::rremove("xlab"),
                                         simulation.dist.suboptimal.2,
                                         nrow=2, ncol=1, align="hv",
                                         heights=c(1, 0.4)),
                              bottom=ggpubr::text_grob("Predicted pairwise treatment effect", size=14)),
                  width=5, height=6, dpi=300)
  
  # save quantiles
  utils::write.table(quantiles.calibration.plot, 
                     file=paste0('./Results/', treatment.arm, 
                                 '/quantiles.calibration.plot.', 
                                 treatment.arm, '.txt'),
                     col.names=FALSE, sep=",")
  
  # save metric values
  metric.values <- cbind(c(metrics.df[, "names"]),
                           c(optimal.model.metrics),
                           c(suboptimal.1.model.metrics),
                           c(suboptimal.2.model.metrics),
                           c(suboptimal.3.model.metrics),
                           c(suboptimal.4.model.metrics))
  utils::write.table(metric.values,
             file=paste0('./Results/', treatment.arm, '/metrics.', treatment.arm, 
                         '.suboptimal.optimal',
                         ifelse(random.matching, ".random", ""),
                         ifelse(match.on.benefit, ".benefit", ""),
                         ifelse(match.on.covariates, ".covariates", ""), '.txt'),
             row.names=FALSE, sep=",")
  
  return(saved.results)
}

#######
####### PLOT PERTURBATIONS
#######
plot.perturbations <- function(treatment.arm=NULL, model=NULL, pred=NULL, log.scale=FALSE){
  if (log.scale){
    data.plot <- data.frame(p.0=qlogis(pred$p.0),
                            p.1=qlogis(pred$p.1))
  } else{
    data.plot <- data.frame(p.0=pred$p.0,
                            p.1=pred$p.1)
  }
  
  plot <- ggplot2::ggplot(data=data.plot, 
                          ggplot2::aes(x=p.0, y=p.1), show.legend=TRUE)
  x.lab.plot <- ""
  y.lab.plot <- ""
  plot <- plot+ggplot2::labs(x=x.lab.plot, y=y.lab.plot, color=" ")# axis names
  plot <- plot+ggplot2::geom_abline(intercept=0, slope=1, col="red")               # diagonal line
  plot <- plot+ggplot2::geom_point(ggplot2::aes(y=p.1), size=1)    # draw p.1
  plot <- plot+ggplot2::theme_light(base_size=25)                     # increase font size
  plot <- plot+ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())
  
  if (log.scale){
    limit <- 10
    plot <- plot+ggplot2::scale_y_continuous(limits=c(-limit, limit)) # axis limits
    plot <- plot+ggplot2::scale_x_continuous(limits=c(-limit, limit)) # axis limits
  } else{
    plot <- plot+ggplot2::scale_y_continuous(limits=c(0, 1)) # axis limits
    plot <- plot+ggplot2::scale_x_continuous(limits=c(0, 1)) # axis limits
  }
  
  return(plot)
}

######
###### OBTAIN METRIC VALUES FOR CASE STUDY
######
metric.values.three.models <- function(boot=0, pred=NULL,
                                       treatment.arm=NULL, metrics.df=NULL,
                                       Y=NULL, X=NULL, W=NULL, matched.test=NULL){
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
      matched.df <- HTEPredictionMetrics::match.patients(Y=Y, W=W,
                                     X=X,
                                     p.0=pred.method$p.0,
                                     p.1=pred.method$p.1,
                                     tau.hat=pred.method$tau.hat)$df.matched.patients
      # cat("ATE:", mean(pred.method$tau.hat), "\n")
      # cat("ATE matched:", mean(matched.df$matched.tau.hat), "\n")
    }
    assign(paste0('matched.', method), matched.df)
    
    overall.cal.measure <- mean(matched.df$matched.tau.obs) - mean(matched.df$matched.tau.hat)
    out.E <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE)
    out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                           CI=FALSE, message=FALSE)
    out.OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.df,
                             CI=FALSE, message=FALSE)
    metric.values <- c(overall.cal.measure, out.E$Eavg.for.benefit,
                       out.E$E50.for.benefit, out.E$E90.for.benefit,
                       out.OP$Cross.entropy.for.benefit, out.OP$Brier.for.benefit,
                       out.C$C.for.benefit)
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
      cal.plot <- HTEPredictionMetrics::calibration.plot(matched.patients=out.E$matched.patients, g=5,
                                   plot.CI=TRUE, show=FALSE)
      assign(paste0('cal.plot.', method), cal.plot$build.plot)
      
      # temp
      if (boot==-1){
        ggplot2::ggsave(file=paste0('./Results/', treatment.arm, '/case.study.', method, ".", data.name, '.png'),
                        plot=cal.plot$build.plot, 
                        width=10, height=10, dpi=300)   
      }
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
##### SENSITIVITY ANALYSIS FOR MATCHING
#####
SA.matching <- function(pred=NULL, treatment.arm=NULL, metrics.df=NULL,
                        Y=NULL, X=NULL, W=NULL){
  for (parameter.set in list(list("nearest", "mahalanobis", FALSE), # currently used
                             list("nearest", "robust_mahalanobis", FALSE),
                             list("optimal", "mahalanobis", FALSE))){
    cal.plots <- list()
    i <- 1
    # differ measure, distance, caliper, replace
    for (method in c("risk", "effect", "CF")){
      pred.method <- data.frame(tau.hat=eval(parse(text=paste0('pred$', method, '.tau.hat'))),
                                p.0=eval(parse(text=paste0('pred$', method, '.p.0'))),
                                p.1=eval(parse(text=paste0('pred$', method, '.p.1'))))
      matched.df <- HTEPredictionMetrics::match.patients(Y=Y, W=W,
                                                         X=X,
                                                         p.0=pred.method$p.0,
                                                         p.1=pred.method$p.1,
                                                         tau.hat=pred.method$tau.hat,
                                                         measure=parameter.set[[1]],
                                                         distance=parameter.set[[2]])$df.matched.patients
      
      ATE <- mean(pred.method$p.0*(1-W)+pred.method$p.1*W)
      paired.ATE <- mean(matched.df$matched.tau.hat)
      overall.cal.measure <- mean(matched.df$matched.tau.obs) - mean(matched.df$matched.tau.hat)
      
      # proposed metrics
      out.E <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.df,
                                                   CI=FALSE, message=FALSE)
      out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                                                   CI=FALSE, message=FALSE)
      out.OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.df,
                                                     CI=FALSE, message=FALSE)
      
      metrics.values <- c(ATE, paired.ATE, 
                          overall.cal.measure,
                          out.E$Eavg.for.benefit,
                          out.E$E50.for.benefit,
                          out.E$E90.for.benefit,
                          out.OP$Cross.entropy.for.benefit,
                          out.OP$Brier.for.benefit,
                          out.C$C.for.benefit)
      metrics.df[, paste0(method, ".",
                          parameter.set[[1]], ".",
                          parameter.set[[2]], ".",
                          parameter.set[[3]])] <- metrics.values
      
      # calibration plot
      cal.plot <- HTEPredictionMetrics::calibration.plot(matched.patients=out.E$matched.patients, 
                                                         g=5, plot.CI=TRUE, 
                                                         show=FALSE)
      cal.plots[[i]] <- cal.plot$build.plot+
        ggplot2::scale_x_continuous(limits=c(-0.5, 1))+
        ggplot2::scale_y_continuous(limits=c(-2, 2))+
        ggplot2::xlab("")+ggplot2::ylab("")+
        ggplot2::annotation_custom(gridExtra::tableGrob(d=data.frame(names=metrics.df[, "names"],
                                                                   values=sprintf("%.3f", metrics.values)), 
                                                       rows=NULL, cols=NULL,
                                                       theme=gridExtra::ttheme_default(core=list(padding=ggplot2::unit(c(1.5,1.5), "mm"),
                                                                                                 fg_params=list(hjust=1, x=1, fontsize=12),
                                                                                                 bg_params=list(fill=c("white"))))),
                                  ymin=0.5, ymax=2, xmin=-0.5, xmax=0)
      i <- i + 1
    }
    ggplot2::ggsave(file=paste0('./Results/Application/SA/case.study.',
                                treatment.arm, ".", parameter.set[[1]], ".",
                                parameter.set[[2]], ".", parameter.set[[3]],
                                '.png'),
                    plot=ggpubr::annotate_figure(ggpubr::ggarrange(cal.plots[[1]],
                                                                   cal.plots[[2]],
                                                                   cal.plots[[3]],
                                                                   nrow=2, ncol=2, align="h",
                                                                   labels=c("A", "B", "C", "")),
                                                 left=ggpubr::text_grob("Observed pairwise treatment effect", rot=90, size=20),
                                                 bottom=ggpubr::text_grob("Predicted pairwise treatment effect", size=20)),
                    width=10, height=10, dpi=300)
  }
  openxlsx::write.xlsx(as.data.frame(metrics.df), file=paste0(paste0('./Results/Application/SA/SA.matching.', treatment.arm, '.xlsx')))
  return(metrics.df)
} 

#####
##### APPLICATION
#####
application <- function(data=NULL, treatment.arm=NULL, folds=5, B=B,
                        metrics.df=NULL, effect=FALSE, CF=FALSE, saved.results=NULL){
  # select data
  original.data <- select.on.treatment(treatment.arm=treatment.arm, data=data, scale=FALSE)
  
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

  # cat("ATE:", mean(test.on.train$risk.tau.hat), "\n")
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
  Y.test <- original.data$Y[ind.test]
  X.test <- original.data$X[ind.test,]
  W.test <- original.data$W[ind.test]
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
  ##### SENSITIVITY ANALYSIS FOR MATCHING
  #####
  metrics.SA <- SA.matching(pred=test, treatment.arm=treatment.arm,
                              metrics.df=rbind(data.frame(names=c("ATE", "Pairwise ATE")), metrics.df),
                              Y=Y.test, X=X.test, W=W.test)
  
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
    results.metrics <- as.matrix(read.table(file=paste0('./Results/', treatment.arm, 
                                                        '/metrics.', treatment.arm, '.', method, '.txt'),
                                            header = TRUE, row.names=1, 
                                            sep = ",", dec = "."))
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
      
      # limits calibration plot
      limits.table <- list(xmin=0.25,
                           xmax=1,
                           ymin=-1.5,
                           ymax=-0.3)
      
      # make histogram table
      matched.df <- eval(parse(text=paste0('matched.', data.name, '$matched.', method)))
      matched.df$colors <- ifelse(matched.df$matched.tau.obs==-1, "harm", 
                                  ifelse(matched.df$matched.tau.obs==0, "noeffect", "benefit"))
      matched.df <- matched.df[order(matched.df$matched.tau.hat), ]
      matched.df$tau.obs.cut <- cut(matched.df$matched.tau.hat, 
                                    breaks=seq(min(matched.df$matched.tau.hat), 
                                               max(matched.df$matched.tau.hat), 
                                               by=0.01))
      # count number of observations in each group
      data.cut <- matched.df %>%
        dplyr::group_by(colors, tau.obs.cut) %>%
        dplyr::summarise(n=n()) %>%
        dplyr::mutate(tau.obs.cut=ifelse(is.na(tau.obs.cut), 0, tau.obs.cut))
      
      # set new x-axis range for histogram
      nr.empty.neg <- (limits.benefit$xmin-min(matched.df$matched.tau.hat))/0.01
      nr.empty.pos <- max(data.cut$tau.obs.cut)+(limits.benefit$xmax-max(matched.df$matched.tau.hat))/0.01
      
      # make histogram figure
      dist <- ggplot2::ggplot(data=data.cut, 
                              mapping=ggplot2::aes(x=.data$tau.obs.cut, 
                                                   y=.data$n, 
                                                   fill=.data$colors))+
        ggplot2::geom_bar(stat="identity", alpha=1, show.legend=FALSE, 
                          position="dodge2")+
        ggplot2::scale_x_continuous(breaks=c(nr.empty.neg, nr.empty.pos),
                                    labels=c(nr.empty.neg, nr.empty.pos),
                                    limits=c(nr.empty.neg, nr.empty.pos))+
        ggplot2::theme_void()+
        ggplot2::scale_fill_manual("legend", values=c("harm"="red", "noeffect"="darkgrey", "benefit"= "blue"))
      assign(paste0("case.study.dist.", method, ".", data.name), dist)
      
      # create calibration plots
      cal.plot <- eval(parse(text=paste0('matched.', data.name, '$cal.plot.', method)))
      cal.plot <- cal.plot+ggplot2::theme_light(base_size=20)+                   # increase font size
        ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())+
        ggplot2::scale_y_continuous(labels=round(seq(limits.benefit$ymin, limits.benefit$ymax, length.out=7), 1),
                                    breaks=seq(limits.benefit$ymin, limits.benefit$ymax, length.out=7),
                                    limits=c(limits.benefit$ymin, limits.benefit$ymax+0.1))+
        ggplot2::scale_x_continuous(labels=round(seq(limits.benefit$xmin, limits.benefit$xmax, length.out=7), 1),
                                    breaks=seq(limits.benefit$xmin, limits.benefit$xmax, length.out=7),
                                    limits=c(limits.benefit$xmin, limits.benefit$xmax))+
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
      
      if (data.name== 'test'){
        metric.table <- cbind(metrics.df[1:7, "names"],
                              sprintf("%.3f", as.numeric(results.matrix[, method])),
                              results.matrix[, paste0(method, ".CI")])
        cal.plot <- cal.plot+ggplot2::annotation_custom(gridExtra::tableGrob(metric.table,
                                               theme=gridExtra::ttheme_default(core=list(padding=ggplot2::unit(c(1.5,1.5), "mm"),
                                                                              fg_params=list(hjust=1, x=1, fontsize=12),
                                                                              bg_params=list(fill=c("white"))))),
                                               xmin=limits.table$xmin, 
                                               xmax=limits.table$xmax, 
                                               ymin=limits.table$ymin, 
                                               ymax=limits.table$ymax)
      }
      assign(paste0("case.study.plot.", method, ".", data.name), cal.plot)
    }
  }

  # manuscript figure
  ggplot2::ggsave(file=paste0('./Results/Application/case.study.', treatment.arm, '.png'),
                  plot=ggpubr::annotate_figure(ggpubr::ggarrange(case.study.plot.risk.test,
                                                                 case.study.plot.effect.test,
                                                                 case.study.dist.risk.test,
                                                                 case.study.dist.effect.test,
                                                                 case.study.plot.CF.test,
                                                                 ggplot2::ggplot(data.frame())+ggplot2::theme_void(),
                                                                 case.study.dist.CF.test,
                                                                 nrow=4, ncol=2, align="hv", 
                                                                 heights=c(1, 0.4, 1, 0.4),
                                                                 labels=c("A", "B", "", "", "C", "", "", "")), 
                                               left=ggpubr::text_grob("Observed pairwise treatment effect", rot=90, size=20),
                                               bottom=ggpubr::text_grob("Predicted pairwise treatment effect", size=20)), 
                  width=10, height=10, dpi=300)
  
  # appendix figure
  ggplot2::ggsave(file=paste0('./Results/Application/appendix.case.study.', treatment.arm, '.png'),
                  plot=ggpubr::annotate_figure(ggpubr::ggarrange(case.study.plot.risk.train,
                                                                 case.study.plot.risk.test,
                                                                 case.study.dist.risk.train,
                                                                 case.study.dist.risk.test,
                                                                 case.study.plot.effect.train,
                                                                 case.study.plot.effect.test,
                                                                 case.study.dist.effect.train,
                                                                 case.study.dist.effect.test,
                                                                 case.study.plot.CF.train,
                                                                 case.study.plot.CF.test,
                                                                 case.study.dist.CF.train,
                                                                 case.study.dist.CF.test,
                                                                 nrow=6, ncol=2, align="hv", 
                                                                 heights=c(1, 0.4, 1, 0.4, 1, 0.4),
                                                                 labels=c("A", "B", "", "", "C", "D", "", "", "E", "F",  "", "")), 
                                               left=ggpubr::text_grob("Observed pairwise treatment effect", rot=90, size=20),
                                               bottom=ggpubr::text_grob("Predicted pairwise treatment effect", size=20)), 
                  width=10, height=15, dpi=300)
  
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
    trained.effect <- train.effect(treatment.arm=treatment.arm, Y.train=Y, 
                                   X.train=X, W.train=W,
                                   folds=folds, alpha.reg=alpha.reg, 
                                   penalize=TRUE, print=FALSE, simulation=FALSE)
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
