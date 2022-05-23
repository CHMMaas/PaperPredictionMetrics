library(rms)
library(Hmisc)

# Rubin's rules
Rubin.combine<-function(est,se){
  m<-length(est)
  est.mi<-mean(est)
  var.w<-mean(se^2)
  var.b<-0
  if (m>1) var.b<-sum((est - est.mi)^2)/(m - 1)
  se.mi<-sqrt(var.w+(1+1/m)*var.b)
  return(list(est=est.mi,se=se.mi))

}

# Model-based concordance
mb.c <- function(p.hat){
  n<-length(p.hat)
  ord<-order(p.hat)
  p.hat<-p.hat[ord]
  q.hat<-1-p.hat
  V1<-(p.hat*(cumsum(q.hat)-q.hat)+q.hat*(sum(p.hat)-cumsum(p.hat)))/(n-1)
  V2<-(p.hat*(sum(q.hat)-q.hat)+q.hat*(sum(p.hat)-p.hat))/(n-1)
  mb.c<-sum(V1)/sum(V2)
  return(mb.c)
}


#To test code
#lp.mi=lp.val
#y=y.val
#g=5
#main="All patients"
#save_plots = FALSE


# val.prob function, including multiple (imputation) linear predictors
val.prob.mi<-function(lp.mi,y,g,main="",dist=FALSE,save_plots=FALSE)
{

  n<-length(y)
  m.imp.val<-ncol(lp.mi)
  if (is.null(m.imp.val)){
    m.imp.val <- 1
  }

  cindex<-rep(0,m.imp.val)
  cindex.se<-rep(0,m.imp.val)
  slope<-rep(0,m.imp.val)
  slope.se<-rep(0,m.imp.val)
  int<-rep(0,m.imp.val)
  int.se<-rep(0,m.imp.val)
  cindex<-rep(0,m.imp.val)
  cindex.se<-rep(0,m.imp.val)
  E<-rep(0,m.imp.val)
  E.90<-rep(0,m.imp.val)
  mbc<-rep(0,m.imp.val)
  sm.y<-NULL

  p.groups<-array(rep(0,g*m.imp.val),dim=c(m.imp.val,g),dimnames=list(1:m.imp.val,1:g))
  y.groups<-array(rep(0,2*g*m.imp.val),dim=c(m.imp.val,g,2),dimnames=list(1:m.imp.val,1:g,c("obs","se")))

  for (i in 1:m.imp.val)
  {
    if (m.imp.val!=1){
      lp.val<-lp.mi[,i]
    }

    f.val<-stats::glm(y~lp.val,family='binomial')
    f.val.offset<-stats::glm(y~offset(lp.val),family='binomial')

    #cindex[i]<-f.val$stats["C"]

    rc<-Hmisc::rcorr.cens(lp.val,y)
    cindex[i]<-rc["C Index"]
    cindex.se[i]<-rc["S.D."]/2

    slope[i]<-f.val$coefficients[[2]]
    slope.se[i]<-sqrt(vcov(f.val)[[2,2]])

    int[i]<-f.val.offset$coefficients[[1]]
    int.se[i]<-sqrt(vcov(f.val.offset)[[1,1]])

    p.val<-stats::plogis(lp.val)
    quants<-stats::quantile(p.val,(1:(g-1))/g)
    cuts<-cut(p.val,breaks=c(0,quants,1))
    p.groups[i,]<-tapply(p.val,cuts,mean)
    for (j in 1:g)
    {
      sub<-(cuts==levels(cuts)[j])
      if (sum(y[sub])>0){
        f.0<-stats::glm(y~1,family = binomial,subset=sub)
        y.groups[i,j,1]<-f.0$coef
        y.groups[i,j,2]<-sqrt(vcov(f.0)[1,1])} else {y.groups[i,j,]<- -Inf}
    }

    Sm <- stats::loess(y~p.val,degree=2)
    sm.y<-cbind(sm.y,predict(Sm,newdata=data.frame(p.val=(0:100)/100)))
    E[i]<-mean(abs(Sm$x-Sm$fitted))
    E.90[i]<-stats::quantile(abs(Sm$x-Sm$fitted),.9)

    mbc[i]<-mb.c(p.val)

  }

  p.mi<-colMeans(p.groups)
  obs.mi<-rep(0,g)
  obs.mi.lower<-rep(0,g)
  obs.mi.upper<-rep(0,g)
  for (j in 1:g)
  {
    RC<-Rubin.combine(y.groups[,j,1],y.groups[,j,2])
    obs.mi[j]<-plogis(RC$est)
    obs.mi.lower[j]<-stats::plogis(RC$est+qnorm(.025)*RC$se)
    obs.mi.upper[j]<-stats::plogis(RC$est+qnorm(.975)*RC$se)
  }

  if (save_plots) grDevices::dev.new(grDevices::windows.options(width = 5, height = 5))
  par(mar = c(5,5,2,1))

  lim<-c(0,0.5)
  par(mar = c(5,5,2,1))
  plot(lim,lim,type='l',xlab="Predicted probability",ylab="Observed frequency",main=main,lwd=1,bty='n')

  lines(lim,lim)
  abline(v=quants,col="darkgrey",lwd=1,lty=2)

  sm.y.mi<-rowMeans(sm.y)
  lines((0:100)/100,sm.y.mi,col="dark gray",lwd=2,lty=2)

  segments(p.mi,obs.mi.lower,p.mi,obs.mi.upper)
  points(p.mi,obs.mi,pch=20)

  int.mi<-Rubin.combine(int,int.se)
  slope.mi<-Rubin.combine(slope,slope.se)
  cindex.mi<-Rubin.combine(cindex,cindex.se)
  mbc.mi<-mean(mbc) #mb.c(plogis(as.vector(lp.mi)))
  E.mi<-mean(E) ## Standard errors unclear
  E.90.mi<-mean(E.90) ## Standard errors unclear

  legend(lim[1], lim[2], c(paste("n =",format(n,big.mark=",")),
                           paste("a =",format(round(int.mi$est,2),nsmall=2)),
                           paste("b =",format(round(slope.mi$est,2),nsmall=2)),
                           paste("c =",format(round(cindex.mi$est,2),nsmall=2)),
                           paste("mb.c =",format(round(mbc.mi,2),nsmall=2)),
                           paste("e =",format(round(E.mi,3),nsmall=3)),
                           paste("e.90 =",format(round(E.90.mi,3),nsmall=3))),
         box.col="white",  bg = "white",cex=1)

  # Histogram of risk distribution
  dist=FALSE
  line.bins <- 0.05
  length.seg <- 1
  dist.label <- 0.04
  dist.label2 <- 0.03
  d0lab <- 0
  d1lab <- 1
  cex.d01 <- 0.7

  if (dist){
    x <- rowMeans(plogis(lp.mi))
    bins <- seq(0, min(1,max(lim[2])), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0	<-table(cut(x[y==0],bins))
    f1	<-table(cut(x[y==1],bins))
    j0	<-f0 > 0
    j1	<-f1 > 0
    bins0 <-(bins[-101])[j0]
    bins1 <-(bins[-101])[j1]
    f0	<-f0[j0]
    f1	<-f1[j1]
    maxf <-max(f0,f1)
    f0	<-(0.1*f0)/maxf
    f1	<-(0.1*f1)/maxf

    segments(bins1,line.bins,bins1,length.seg*f1+line.bins)
    segments(bins0,line.bins,bins0,length.seg*-f0+line.bins)
    lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(line.bins,line.bins))
    text(max(bins0,bins1)+dist.label,line.bins+dist.label2,d1lab,cex=cex.d01)
    text(max(bins0,bins1)+dist.label,line.bins-dist.label2,d0lab,cex=cex.d01)
    }

  if (save_plots) {savePlot(main,type="emf"); grDevices::dev.off()}
  return(list(main=main,
              n=n,quants=quants,
              p.mi=p.mi,
              obs.mi=obs.mi,
              obs.mi.lower=obs.mi.lower,
              obs.mi.upper=obs.mi.upper,
              int=int.mi$est,
              int.lower=int.mi$est+qnorm(.025)*int.mi$se,
              int.upper=int.mi$est+qnorm(.975)*int.mi$se,
              slope=slope.mi$est,
              slope.lower=slope.mi$est+qnorm(.025)*slope.mi$se,
              slope.upper=slope.mi$est+qnorm(.975)*slope.mi$se,
              cindex=cindex.mi$est,
              cindex.lower=cindex.mi$est+qnorm(.025)*cindex.mi$se,
              cindex.upper=cindex.mi$est+qnorm(.975)*cindex.mi$se,
              mb.c=mbc.mi,
              e=E.mi,
              e.90=E.90.mi
              ))

}




# setwd("F:\\COVID\\analysis\\PCORI")
# save_plots<-FALSE
#
# # Read COPE validation  dataset, choose between 1st wave (apparent) and 2nd wave (temporal),
# #  and between death, need for icu, or death conditional on icu
# data.val<-read.csv("Wave1_death_28days.csv")
# data.val<-read.csv("Wave2_death_28days.csv")
#
# data.val<-read.csv("Wave1_icu_28days.csv")
# data.val<-read.csv("Wave2_icu_28days.csv")
#
# data.val<-read.csv("Apparent_death_icu_28days.csv")
# data.val<-read.csv("Temporal_death_icu_28days.csv")
#
# data.val<-read.csv("copeV2GeographicValidation.csv")
#
# head(data.val)
#
#
# # Choose between NOCOS and COPE predictions
#
# ### COPE
# lp.val<-qlogis(as.matrix(data.val[,2:6]))
# y.val<-data.val$y.val
#
# ### NOCOS
# lp.val<-qlogis(as.matrix(data.val[,7:11]))
# y.val<-data.val$y.val
#
# ### NOCOS  in Northwell data
# #lp.val<-qlogis(as.matrix(data.val[,c(4:6,8)]))  # Leave out imputation 4: set g=4
# lp.val<-qlogis(as.matrix(data.val[,c(4:8)]))
# y.val<-1-data.val$y28
# lp.val[lp.val==-Inf]<- -10  # Probability 0 ???
#
# sum(is.na(lp.val))
# colMeans(data.val[,4:9])
#
# # Overall validation
# val.prob.mi(lp.mi=lp.val,y=y.val,g=5,main="All patients",save_plots = save_plots)
#
# # Stratified by hospital 1-4
# hospital<-data.val[,1]
# hospitals<-unique(hospital)
# for (h in hospitals)
# {
#   sel.val<-hospital==h
#   print(val.prob.mi(lp.mi=lp.val[sel.val,],y=y.val[sel.val],g=5,main=paste("Hospital",h),save_plots = save_plots))
# }

