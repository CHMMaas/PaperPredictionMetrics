#######
####### SIMULATION EXAMPLE
#######
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) grDevices::dev.off()
cat('\014')

# library(rms)
# library(gridExtra)            # table in figure
# library(png)                  # read png to combine plots
# remotes::install_github("CHMMaas/HTEPredictionMetrics")
# library(HTEPredictionMetrics) # new metrics
# remotes::install_github("CHMMaas/PredictionTools")
library(PredictionTools)    # calibration plot

set.seed(1)

######
###### SIMULATE DATA
######
n.population <- 100000 # 100000
m.X <- 12

OR.control <- c(rep(1,3),rep(1.2,3),rep(1.5,3),rep(2,3))
OR.treated <- c(rep(1,3),c(1.4,1.2,1),c(2,1.5,1),c(2.5,2,1.5))
# 6 interactions: sum(OR.control!=OR.treated)

intercept.control <- -2
OR.treatment <- 0.8

beta.control <- c(intercept.control, 0, log(OR.control))
beta.treated <- c(intercept.control, log(OR.treatment), log(OR.treated))

X <- matrix(rbinom(n = n.population * m.X, 1, 0.2), n.population, m.X)
X <- cbind(1,c(rep(0,n.population/2),rep(1,n.population/2)),X)

# c = intercept
# z = treatment
colnames(X) <- c("c", "z", paste("x", 1:12, sep=""))
head(X)

# create linear predictor of covariates and coefficients
logit.y <- rep(NA,n.population)
logit.y[X[,"z"]==0] <- X[X[,"z"]==0,] %*% beta.control
logit.y[X[,"z"]==1] <-  X[X[,"z"]==1,] %*% beta.treated

# draw outcomes based on the probability plogis(logit.y)
y <- rbinom(n.population,1,plogis(logit.y))

mean(y[X[,"z"]==0])
mean(y[X[,"z"]==1])

f.population <- rms::lrm(y~z*(x1+x2+x3+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12),data = data.frame(X))

######
###### DRAW SAMPLE
######
n.sample <- 3600
i.sample <- sample(n.population,n.sample,replace=FALSE)

f.sample <- rms::lrm(y~z*(x1+x2+x3+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12),data = data.frame(X),sub=i.sample)

pred <- plogis(predict(f.sample,newdata=data.frame(X)))

######
###### CALIBRATION PLOT
######
lp.val <- predict(f.sample,newdata=data.frame(X), type="lp")
y.val <- y
PredictionTools::val.prob.mi(lp.mi=lp.val,y=y.val,g=5,main="Calibration of outcome predictions")
outcome.calibration <- PredictionTools::val.prob.mi(lp.mi=lp.val,y=y.val,g=5,main="Calibration of outcome predictions")

p.val <- predict(f.sample,newdata=data.frame(X), type="fitted")
outcome.cross.entropy <- -(sum(y.val*log(p.val))+sum((1-y.val)*log(1-p.val)))/length(y.val)
outcome.brier <- sum((y.val-p.val)^2)/length(y.val)

y.smoothed <- predict(loess(y~pred,degree=2))
df.outcome <- data.frame(x=pred, y=y.smoothed)
segments.outcome <- data.frame(x=outcome.calibration$p.mi,
                               y=outcome.calibration$obs.mi,
                               y.lower=outcome.calibration$obs.mi.lower,
                               y.upper=outcome.calibration$obs.mi.upper)
outcome.plot <- ggplot2::ggplot(data=df.outcome, ggplot2::aes(x=x, y=y),
                                show.legend=TRUE)+                     # set data
  ggplot2::theme_light(base_size=22)+               # increase font size
  ggplot2::geom_line(data=df.outcome, ggplot2::aes(y=y), col="blue", size=1)+   # blue LOESS line
  ggplot2::geom_segment(data=segments.outcome, 
                        mapping=ggplot2::aes(x=x, y=y.lower, xend=x, yend=y.upper), size=1)+
  ggplot2::geom_point(data=segments.outcome, mapping=ggplot2::aes(x=x, y=y), size=2)+
  ggplot2::geom_abline(intercept=0, linetype="dashed")+ # 45-degree line
  ggplot2::labs(x="Predicted outcome", y="Observed outcome", color=" ")+   # axis names
  ggplot2::ylim(0, 0.5)+
  ggplot2::xlim(0, 0.5)+
  ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0.5, 0, 0), "cm"))
metric.table.outcome <- cbind(c("Eavg", "E90", "Cross-entropy", "Brier score", "C-index"),
                            sprintf("%.3f", c(outcome.calibration$e.avg, 
                                              outcome.calibration$e.90, 
                                              outcome.cross.entropy,
                                              outcome.brier,
                                              outcome.calibration$cindex)))
outcome.plot <- outcome.plot+ggplot2::annotation_custom(gridExtra::tableGrob(metric.table.outcome,
                                                                             theme=gridExtra::ttheme_default(core=list(fg_params=list(hjust=1, x=1, fontsize=14),
                                                                             bg_params=list(fill=c("white"))))),
                                                                             xmin=0.35, xmax=0.5, ymin=0, ymax=0.1)+
  ggplot2::ggtitle("A. Calibration of outcome")+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))

######
###### BENEFIT PREDICTION
######
X.0 <- X.1 <- X
X.0[,"z"] <- 0
X.1[,"z"] <- 1

pred.0 <- stats::plogis(predict(f.sample,newdata=data.frame(X.0)))
pred.1 <- stats::plogis(predict(f.sample,newdata=data.frame(X.1)))

pred.benefit <- pred.0 - pred.1

#####
##### MATCH PATIENTS
#####
out.matching <- HTEPredictionMetrics::match.patients(Y=y, W=X[, "z"],
                               X=X[, -c(1, 2)],
                               p.0=pred.0,
                               p.1=pred.1,
                               tau.hat=pred.benefit)
save(out.matching, file="./Results/Illustration/illustration.matched.RData")
load(file="./Results/Illustration/illustration.matched.RData")

#####
##### C-INDEX
#####
matched.df <- out.matching$df.matched.patients
out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                       CI=FALSE, message=FALSE, replace=FALSE)

#####
##### GROUPED CALIBRATION PLOT
#####
grouped.calibration <- function(Y=NULL, W=NULL, tau.hat=NULL, g=5, limits=NULL){
  # quantile cuts based on unmatched treatment benefit
  unmatched.df <- data.frame(Y=Y, W=W, tau.hat=tau.hat)
  quantiles <- c(min(unmatched.df$tau.hat)-0.01,
                 as.numeric(quantile(unmatched.df$tau.hat, (1:(g-1)/g))),
                 max(unmatched.df$tau.hat))

  # order data frame
  ordered.df <- unmatched.df[order(unmatched.df$tau.hat), ]

  # define quantiles
  ordered.df$quantile.nr <- cut(ordered.df$tau.hat, breaks=quantiles, labels=FALSE)

  # x-axis value is mean of unmatched treatment benefit in the quantiles
  x.of.quant <- stats::aggregate(ordered.df, list(ordered.df$quantile.nr), mean)$tau.hat

  # untreated df
  df.0 <- ordered.df[ordered.df$W==0,]
  mean.0 <- stats::aggregate(df.0, list(df.0$quantile.nr), mean)$Y

  # treated df
  df.1 <- ordered.df[ordered.df$W==1,]
  mean.1 <- stats::aggregate(df.1, list(df.1$quantile.nr), mean)$Y

  # y-axis value is the difference of the mean prevalence between control and treated patients in a quantile
  y.of.quant <- mean.0 - mean.1

  # confidence interval
  sd.0 <- stats::aggregate(df.0, list(df.0$quantile.nr), sd)$Y
  sd.1 <- stats::aggregate(df.1, list(df.1$quantile.nr), sd)$Y
  n.0 <- rep(0, g)
  n.1 <- rep(0, g)
  for (q.nr in 1:g){
    n.0[q.nr] <- nrow(df.0[df.0$quantile.nr==q.nr,])
    n.1[q.nr] <- nrow(df.1[df.1$quantile.nr==q.nr,])
  }
  sd <- sqrt(sd.0^2/n.0+sd.1^2/n.1)
  y.lower <- y.of.quant-1.96*sd
  y.upper <- y.of.quant+1.96*sd
  df.grouped <- data.frame(x=x.of.quant, y=y.of.quant, y.lower=y.lower, y.upper=y.upper)

  # create plot
  build.plot <- ggplot2::ggplot(data=df.grouped, ggplot2::aes(x=x, y=y),
                                show.legend=TRUE)+                     # set data
    ggplot2::theme_light(base_size=22)+               # increase font size
    ggplot2::geom_segment(data=df.grouped,
                          mapping=ggplot2::aes(x=x, y=y.lower,
                          xend=x, yend=y.upper), size=1)+
    ggplot2::geom_point(data=df.grouped, mapping=ggplot2::aes(x=x, y=y), size=2)+
    ggplot2::geom_abline(intercept=0, linetype="dashed")+# 45-degree line
    ggplot2::labs(x="Predicted treatment effect",
                  y="Observed treatment effect", color=" ")+   # axis names
    ggplot2::ylim(limits$ymin, limits$ymax)+
    ggplot2::xlim(limits$xmin, limits$xmax)+
    ggplot2::ggtitle("B. Grouped calibration of benefit")+
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))

  return(build.plot)
}
g <- 5
limits <- list(xmin=-0.2, xmax=0.3, ymin=-0.2, ymax=0.3)
grouped.plot <- grouped.calibration(Y=y, W=X[, "z"], tau.hat=pred.benefit, g=g, limits=limits)
Cindex <- c("C-for-benefit", sprintf("%.3f", out.C$C.for.benefit))
grouped.plot <- grouped.plot+ggplot2::annotation_custom(gridExtra::tableGrob(Cindex,
                                                         theme=gridExtra::ttheme_default(core=list(fg_params=list(hjust=1, x=1, fontsize=14),
                                                                                        bg_params=list(fill=c("white"))))),
                                               xmin=0.22, xmax=0.3, ymin=-0.2, ymax=-0.19)
grouped.plot <- grouped.plot+ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0.5, 0, 0), "cm"))

#####
##### CALCULATE NEW METRICS
#####
overall.cal.measure <- mean(matched.df$matched.tau.obs) - mean(matched.df$matched.tau.hat)
out.E <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.df,
                       CI=FALSE, message=FALSE, replace=FALSE)
out.OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.df,
                         CI=FALSE, message=FALSE, replace=FALSE)
out.C <- HTEPredictionMetrics::C.for.Benefit(matched.patients=matched.df,
                                             CI=FALSE, message=FALSE, replace=FALSE)
metrics <- c(overall.cal.measure, out.E$Eavg.for.benefit,
             out.E$E50.for.benefit, out.E$E90.for.benefit,
             out.OP$Cross.entropy.for.benefit, out.OP$Brier.for.benefit,
             out.C$C.for.benefit)

# plot calibration
cal.plot <- HTEPredictionMetrics::calibration.plot(matched.patients=out.E$matched.patients, g=g,
                             limits=limits, plot.CI=FALSE, show=FALSE)
metric.table <- cbind(c("Calibration-in-the-large", "Eavg-for-benefit", "E50-for-benefit", "E90-for-benefit",
                        "Cross-entropy-for-benefit", "Brier-for-benefit",
                        "C-for-benefit"),
                      c(sprintf("%.3f", as.numeric(metrics[1:7]))))
limits.table <- list(xmin=limits$xmin-limits$xmin*1.23,
                     xmax=limits$xmax,
                     ymin=limits$ymin,
                     ymax=limits$ymax-limits$ymax*1.15)
cal.plot$build.plot <- cal.plot$build.plot+ggplot2::annotation_custom(gridExtra::tableGrob(metric.table,
                                                                       theme=gridExtra::ttheme_default(core=list(fg_params=list(hjust=1, x=1, fontsize=14),
                                                                                                      bg_params=list(fill=c("white"))))),
                                                             xmin=limits.table$xmin, xmax=limits.table$xmax, ymin=limits.table$ymin, ymax=limits.table$ymax)+
  ggplot2::scale_y_continuous(labels=round(seq(from=limits$ymin, to=limits$ymax, length.out=6), 1),
                             breaks=round(seq(from=limits$ymin, to=limits$ymax, length.out=6), 1),
                             limits=c(limits$ymin, limits$ymax))+
  ggplot2::scale_x_continuous(labels=round(seq(from=limits$xmin, to=limits$xmax, length.out=6), 1),
                              breaks=round(seq(from=limits$xmin, to=limits$xmax, length.out=6), 1),
                              limits=c(limits$xmin, limits$xmax))+
  ggplot2::ggtitle("C. Calibration of benefit")+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
  ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0.5, 0, 0), "cm"))

# conventional metrics
conv <- rms::val.prob(pred,y,g=4)
# val.prob.ci.2 ?
conv.metrics <- as.numeric(c(conv["Intercept"], conv["Eavg"], NA, conv["E90"], NA, conv["Brier"], conv["C (ROC)"]))
conv.new <- cbind(metric.table, round(conv.metrics, 3))
colnames(conv.new) <- c("Metric names", "New", "Conventional")
conv.new

# combine calibration plots into one
ggplot2::ggsave(filename="./Results/Illustration/illustration.png", 
       plot=ggpubr::ggarrange(outcome.plot, grouped.plot, cal.plot$build.plot, ncol=3),
       width=20, height=7, dpi=300)
