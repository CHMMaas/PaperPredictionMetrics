library(openxlsx)
dat.orig <- as.data.frame(read.csv('./Data/DPP_DPPOS.csv'))

results.table <- c()
life <- dat.orig[dat.orig$life==1,-c(8, 9)]
met <- dat.orig[dat.orig$met==1,-c(8, 9)]
control <- dat.orig[dat.orig$life!=1&dat.orig$met!=1, -c(8, 9)]
dat.orig <- dat.orig[, -c(8, 9)]

n <- nrow(dat.orig)
results.table <- rbind(results.table, c(n, 0, 0, nrow(life), 0, nrow(met), 0, nrow(control), 0))
for (column in c(2, 10, 12, 13, 15, 19, 20, 21, 23)){
  results.table <- rbind(results.table, c(sum(dat.orig[!is.na(dat.orig[, column]), column]),
                                          round(sum(dat.orig[!is.na(dat.orig[, column]), column])/n*100, 1), sum(is.na(dat.orig[, column])),
                                          sum(life[!is.na(life[, column]), column]),
                                          round(sum(life[!is.na(life[, column]), column])/n*100, 1),
                                          sum(met[!is.na(met[, column]), column]),
                                          round(sum(met[!is.na(met[, column]), column])/n*100, 1),
                                          sum(control[!is.na(control[, column]), column]),
                                          round(sum(control[!is.na(control[, column]), column])/n*100, 1)))
}
for (column in c(11, 14, 4, 5, 22)){
  results.table <- rbind(results.table, c(stats::median(dat.orig[!is.na(dat.orig[, column]), column]),
                                          paste0('[', as.numeric(quantile(dat.orig[!is.na(dat.orig[, column]), column], 0.25)), '; ',
                                                 as.numeric(quantile(dat.orig[!is.na(dat.orig[, column]), column], 0.75)), ']'),
                                          sum(is.na(dat.orig[, column])),
                                          stats::median(life[!is.na(life[, column]), column]),
                                          paste0('[', as.numeric(quantile(life[!is.na(life[, column]), column], 0.25)), '; ',
                                                 as.numeric(quantile(life[!is.na(life[, column]), column], 0.75)), ']'),
                                          stats::median(met[!is.na(met[, column]), column]),
                                          paste0('[', as.numeric(quantile(met[!is.na(met[, column]), column], 0.25)), '; ',
                                                 as.numeric(quantile(met[!is.na(met[, column]), column], 0.75)), ']'),
                                          stats::median(control[!is.na(control[, column]), column]),
                                          paste0('[', as.numeric(quantile(control[!is.na(control[, column]), column], 0.25)), '; ',
                                                 as.numeric(quantile(control[!is.na(control[, column]), column], 0.75)), ']')))
}
colnames(results.table) <- c("Total", "%", "Missing", "life", "%", "met", "%", "control", "%")
results.table <- rbind(results.table[1:3,], rep(0, 9), results.table[4:nrow(results.table),])
rownames(results.table) <- c("Sample size", "Diabetes", "Females",
                             "Ethnicity", "Black", "Hispanic", "Hxhbg", "FAMHXDB", "Smoking",
                             "Hyptertension", "GDM", "Age", "BMI", "Triglycerides", "Hemoglobin", "Fpg")
results.table[results.table==0] <- ""
openxlsx::write.xlsx(x=as.data.frame(results.table), file="./Results/descriptives.xlsx", asTable=FALSE, overwrite=TRUE, rowNames=TRUE)
results.table
