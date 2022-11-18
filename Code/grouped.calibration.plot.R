###
### GROUPED CALIBRATION PLOT
### Author: C.H.M. Maas
### 

# create a dataframe
df <- data.frame(tau.hat=tau.hat, 
                 tau.obs=tau.obs)

# plot quantiles
g <- 10

# define the edges of the quantiles
quantiles <- c(min(df$tau.hat)-0.01,
               as.numeric(quantile(df$tau.hat, (1:(g-1)/g))),
               max(df$tau.hat))

# sort dataframe on tau.hat
ordered.df <- df[order(df$tau.hat), ]

# split the ordered dataframe in g groups by creating a variable quantile.nr that indicates the quantile number
ordered.df$quantile.nr <- cut(ordered.df$tau.hat, breaks=quantiles, labels=FALSE)

# x-axis value is mean of predicted treatment effect in the quantiles
x.of.quant <- aggregate(ordered.df, list(ordered.df$quantile.nr), mean)$tau.hat

# y-axis value is mean of observed treatment effect in the quantiles
y.of.quant <- aggregate(ordered.df, list(ordered.df$quantile.nr), mean)$tau.obs

# determine standard deviation in group
sd <- aggregate(ordered.df, list(ordered.df$quantile.nr), sd)$tau.obs

# calculate the number of observations in each quantile
n <- rep(NA, g)
for (q.nr in 1:g){
  # n is the number of pairs in each group
  n[q.nr] <- nrow(unique(ordered.df[ordered.df$quantile.nr==q.nr,]))
}

# calculate the standard deviation in each group
sd <- sqrt(sd^2/n)

# determine confidence interval
y.lower <- y.of.quant-1.96*sd
y.upper <- y.of.quant+1.96*sd

# plot
data.plot <- data.frame(x=x.of.quant, y.lower=y.lower, y.upper=y.upper)
plot <- ggplot2::ggplot(data=data.plot, 
                        ggplot2::aes(x= .data$tau.hat, y=.data$tau.obs), 
                        show.legend=TRUE)+          # set data
  ggplot2::theme_light(base_size=22)+                     # increase font size
  ggplot2::geom_abline(intercept=0, linetype="dashed")+   # 45-degree line
  ggplot2::labs(x="Predicted grouped treatment effect",
                y="Observed grouped treatment effect", color=" ")+
  ggplot2::scale_x_continuous(limits=c(-1, 1))+
  ggplot2::scale_y_continuous(limits=c(-1, 1))+
  ggplot2::geom_segment(data=data.plot,
                        mapping=ggplot2::aes(x = x.of.quant, 
                                             y = y.lower, 
                                             xend = x.of.quant, 
                                             yend = y.upper), size=1)+
  ggplot2::geom_point(data=data.plot,
                      mapping=ggplot2::aes(x=x.of.quant, y=y.of.quant), size=2)
show(plot)