#######
####### MAKE FOREST PLOT
#######
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) dev.off()
cat('\014')

library(png)
library(ggplot2)
library(grid)
library(gridExtra)

# simulation study log odds figures
for (treatment.arm in c("life", "met")){
  for (model in c("optimal", "suboptimal.1", "suboptimal.2", "suboptimal.3")){
    load(paste0('./Results/', treatment.arm, "/log.odds.", model, ".Rdata"))
    assign(paste0("log.odds.", model), plot)
  }
  png(file=paste0('./Results/Simulation/log.odds.', treatment.arm, '.png'), width=960, height=960, units="px")
  bottom.title <- grid::textGrob("Prognostic index", gp=gpar(fontsize=30))
  bottom.title <- gridExtra::arrangeGrob(bottom.title, ggplot2::zeroGrob(), widths=unit(1, 'npc'),
                              heights=unit(c(0.5, 1), c('cm', 'npc')),
                              as.table=FALSE)
  gridExtra::grid.arrange(log.odds.optimal, log.odds.suboptimal.1,
               log.odds.suboptimal.2, log.odds.suboptimal.3,
               bottom=bottom.title,
               left=textGrob("Log odds", gp=gpar(fontsize=30), rot=90))
  grDevices::dev.off()
}

# simulation study figures
for (treatment.arm in c("life", "met")){
  for (model in c("optimal", "suboptimal.1", "suboptimal.2", "suboptimal.3")){
    load(paste0('./Results/', treatment.arm, "/", model, ".simulation.calibration.plot.Rdata"))
    assign(paste0("simulation.", model), plot)
  }
  png(file=paste0('./Results/Simulation/simulation.', treatment.arm, '.png'), width=960, height=960, units="px")
  bottom.title <- grid::textGrob("Predicted treatment effect", gp=gpar(fontsize=30))
  bottom.title <- gridExtra::arrangeGrob(bottom.title, ggplot2::zeroGrob(), widths=unit(1, 'npc'),
                              heights=unit(c(0.5, 1), c('cm', 'npc')),
                              as.table=FALSE)
  gridExtra::grid.arrange(simulation.optimal, simulation.suboptimal.1,
               simulation.suboptimal.2, simulation.suboptimal.3,
               bottom=bottom.title,
               left=textGrob("Observed treatment effect", gp=gpar(fontsize=30), rot=90))
  grDevices::dev.off()
}

# case study figures
for (treatment.arm in c("life", "met")){
  for (data.name in c("test", "train")){
    for (model in c("risk", "effect", "CF")){
      load(paste0('./Results/', treatment.arm, "/", data.name, ".", model, ".case.study.calibration.plot.Rdata"))
      assign(paste0("case.study.", model, ".", data.name), plot)
    }
  }
  png(file=paste0('./Results/Application/case.study.', treatment.arm, '.png'), width=1000, height=1000, units="px")
  bottom.title <- grid::textGrob("Predicted treatment effect", gp=gpar(fontsize=30))
  bottom.title <- gridExtra::arrangeGrob(bottom.title, ggplot2::zeroGrob(), widths=unit(1, 'npc'),
                              heights=unit(c(0.5, 1), c('cm', 'npc')),
                              as.table=FALSE)
  gridExtra::grid.arrange(case.study.risk.test, case.study.effect.test,
               case.study.CF.test,
               bottom=bottom.title, ncol=2,
               left=textGrob("Observed treatment effect", gp=gpar(fontsize=30), rot=90))
  grDevices::dev.off()
}

# case study figures
for (treatment.arm in c("life", "met")){
  for (data.name in c("test", "train")){
    for (model in c("risk", "effect", "CF")){
      load(paste0('./Results/', treatment.arm, "/", data.name, ".", model, ".case.study.appendix.calibration.plot.Rdata"))
      assign(paste0("appendix.case.study.", model, ".", data.name), plot)
    }
  }
  png(file=paste0('./Results/Application/appendix.case.study.', treatment.arm, '.png'), width=1000, height=1500, units="px")
  bottom.title <- grid::textGrob("Predicted treatment effect", gp=gpar(fontsize=30))
  bottom.title <- gridExtra::arrangeGrob(bottom.title, ggplot2::zeroGrob(), widths=unit(1, 'npc'),
                              heights=unit(c(0.5, 1), c('cm', 'npc')),
                              as.table=FALSE)
  gridExtra::grid.arrange(appendix.case.study.risk.train,
               appendix.case.study.risk.test,
               appendix.case.study.effect.train,
               appendix.case.study.effect.test,
               appendix.case.study.CF.train,
               appendix.case.study.CF.test,
               bottom=bottom.title, ncol=2,
               left=textGrob("Observed treatment effect", gp=gpar(fontsize=30), rot=90))
  grDevices::dev.off()
}
