#######
####### LOAD DATA
#######
library(mice)         # imputation of missing values

load.data <- function(){
  # read original data
  dat.orig <- read.csv('./Data/DPP_DPPOS.csv')

  # omit release_id, DIABV, diabf_OS, diabv_OS, months_OS
  dat.cc <- dat.orig[,-c(1, 3, 6, 7, 11)]

  # omit na values
  # dat.cc <- na.omit(dat.cc)
  # cat("number of outliers removed: ", nrow(dat.orig) - nrow(dat.cc), '\n')

  # impute NA values
  m <- 1 # single imputation
  imputed.data <- mice::mice(dat.cc, m, print=FALSE, seed=500)
  dat.cc <- mice::complete(imputed.data, m)

  # set new id
  id <- 1:nrow(dat.cc)
  dat.cc <- cbind(id, dat.cc)

  # data transformations
  dat.cc$fpg <- log(dat.cc$fpg)
  dat.cc$MAQMETHR <- log(1+dat.cc$MAQMETHR) # avoid log(0)

  return(list(dat.orig=dat.orig, dat.cc=dat.cc))
}

#######
####### SELECT PROGNOSTIC FACTORS
#######
select.X <- function(data=NULL){
  # omit dat$waist, dat$whr, and dat$MAQMETHR
  X <- as.matrix(cbind(data$TRIG, data$HBA1, data$female, data$age, data$black, data$hispa, data$bmi,
                          data$hxhbg, data$FAMHXDB, data$smk, data$htn, data$fpg, data$GDM))
  colnames(X) <- c('TRIG', 'HBA1', 'female', 'age', 'black', 'hispa', 'bmi',
                      'hxhbg', 'FAMHXDB', 'smk', 'htn', 'fpg', 'GDM')

  return(X)
}

#######
####### SELECT ON TREATMENT
#######
select.on.treatment <- function(treatment.arm=NULL, data=NULL, scale=FALSE){
  if (treatment.arm == 'life'){
    # select data for which met == 0 to compare life with placebo
    ind.tot <- which(data$met==0)  # indices

    # only select met == 0 if met == 1 is present in data
    if (length(ind.tot) == 0){
      dat <- data
      ind.tot <- 1:nrow(data)
    }
    else{
      dat <- data[data$met==0,]
    }
    W <- dat$life                   # select treatment assignment indicator
  }
  else if (treatment.arm == 'met'){
    # select data for which life == 0 to compare life with placebo
    ind.tot <- which(data$life==0)

    # only select life == 0 if life == 1 is present in data
    if (length(ind.tot) == 0){
      dat <- data
      ind.tot <- 1:nrow(data)
    }
    else{
      dat <- data[data$life==0,]
    }
    W <- dat$met                    # select treatment assignment indicator
  }

  Y <- dat$DIABF                    # event indicator
  time <- dat$months                # time indicator
  id <- dat$id                      # id of patients

  # select prognostic factors
  X <- select.X(data=dat)
  if (scale){
    X <- scale(X) # standardize prognostic factors
  }

  # create data frame with selected treatment arm
  data <- as.data.frame(cbind(id, Y, time, W, X))
  colnames(data) <- c(c('id', 'DIABF', 'months', treatment.arm), colnames(X))

  return(list(Y=Y, X=X, W=W, time=time, ind.tot=ind.tot, id=id, data=data))
}
