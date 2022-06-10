#######
####### RUN SCRIPTS
#######
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) dev.off()
cat('\014')

####### Install packages
## install txBenefit:
# install.packages("devtools")
# devtools::install_github("msadatsafavi/txBenefit")
# install.packages('mice')
# install.packages('glmnet')
# install.packages('grf')
# install.packages('MatchIt')
# install.packages('optmatch')
# install.packages('fANCOVA')
# install.packages('remotes')
# remotes::install_github("CHMMaas/HTEPredictionMetrics")

source('./Code/common.functions.R')

# reproducability
set.seed(1)

# set parameters
alpha.reg <- 0                # Ridge (alpha.reg=0), Lasso (alpha.reg=1), Elastic Net (alpha.reg=0.5) with no shrinkage on treatment arm indicator
folds <- 5                    # number of folds to do cross validation on training/test and training/validation set (folds=5)
treatment.arm <- 'life'       # life or met
compute.new.results <- TRUE   # compute new results or not

#####
##### Simulations
#####
R <- 0                         # repetitions in simulation (R = 100)
plot.cal <- TRUE              # if TRUE calibration plot of simulation is made

#####
##### Application
#####
B <- 100                      # number of bootstrap samples (additional to original sample computations) (B = 100)
effect.indicator <- TRUE      # use penalized treatment effect model to estimate treatment effect (TRUE)
CF.indicator <- TRUE          # use causal forest to estimate treatment effect (TRUE)

#####
##### Run analysis
#####
if (R > 0 & B > 0){
  cat("Set either R or B to zero \n")
  break
}

if (compute.new.results){
  simulation.results <- NULL
  application.results <- NULL
} else{
  if (R > 0 & B == 0){
    simulation.results <- get(load(paste("./Results/Simulation/simulation.results.", treatment.arm, ".RData", sep="")))
    application.results <- NULL
  } else if (R == 0 & B > 0){
    simulation.results <- NULL
    application.results <- get(load(paste("./Results/Application/application.results.", treatment.arm, ".RData", sep="")))
  }
}

results <- run.analysis(CF=CF.indicator, effect=effect.indicator, R=R,
                name.data.set='DPP', treatment.arm=treatment.arm,
                plot.cal=plot.cal, B=B, folds=folds, alpha.reg=alpha.reg,
                simulation.results=simulation.results, application.results=application.results)

if (compute.new.results){
  cat("Overwriting old files..\n")
  if (R > 0 & B == 0){
    save(results, file=paste("./Results/Simulation/simulation.results.", treatment.arm, ".RData", sep=""))
  } else if (R == 0 & B > 0){
    save(results, file=paste("./Results/Application/application.results.", treatment.arm, ".RData", sep=""))
  }
}
