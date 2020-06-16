require(asreml)
require(asremlPlus)
library(grid)

############## Results
facs <- fac.gen(list(Design = 0:5, simul = 1000))

correl <- list()
correl[["1"]] <- cbind(facs, 
                       data.frame(corr = c(AnalyseDes0_04[[2]],AnalyseDes1_04[[2]],
                                           AnalyseDes2_04[[2]],AnalyseDes3_04[[2]],
                                           AnalyseDes4_04[[2]],AnalyseDes5_04[[2]])))

correl[["2"]] <-  cbind(facs, 
                        data.frame(corr = c(AnalyseDes0_04[[2]],AnalyseDes1_04a[[2]],
                                            AnalyseDes2_04a[[2]],AnalyseDes3_04a[[2]],
                                            AnalyseDes4_04a[[2]],AnalyseDes5_04a[[2]])))


correl[["4"]] <- cbind(facs, 
                       data.frame(corr = c(AnalyseDes0_25[[2]],AnalyseDes1_25[[2]],
                                           AnalyseDes2_25[[2]],AnalyseDes3_25[[2]],
                                           AnalyseDes4_25[[2]],AnalyseDes5_25[[2]])))

correl[["3"]] <- cbind(facs, 
                       data.frame(corr = c(AnalyseDes0_25[[2]],AnalyseDes1_25a[[2]],
                                           AnalyseDes2_25a[[2]],AnalyseDes3_25a[[2]],
                                           AnalyseDes4_25a[[2]],AnalyseDes5_25a[[2]])))


RGG <- list()
RGG[["1"]] <- cbind(facs, 
                       data.frame(RGG = c(AnalyseDes0_04[[4]],AnalyseDes1_04[[4]],
                                           AnalyseDes2_04[[4]],AnalyseDes3_04[[4]],
                                           AnalyseDes4_04[[4]],AnalyseDes5_04[[4]])))

RGG[["2"]] <-  cbind(facs, 
                        data.frame(RGG = c(AnalyseDes0_04[[4]],AnalyseDes1_04a[[4]],
                                            AnalyseDes2_04a[[4]],AnalyseDes3_04a[[4]],
                                            AnalyseDes4_04a[[4]],AnalyseDes5_04a[[4]])))


RGG[["4"]] <- cbind(facs, 
                       data.frame(RGG = c(AnalyseDes0_25[[4]],AnalyseDes1_25[[4]],
                                           AnalyseDes2_25[[4]],AnalyseDes3_25[[4]],
                                           AnalyseDes4_25[[4]],AnalyseDes5_25[[4]])))

RGG[["3"]] <- cbind(facs, 
                       data.frame(RGG = c(AnalyseDes0_25[[4]],AnalyseDes1_25a[[4]],
                                           AnalyseDes2_25a[[4]],AnalyseDes3_25a[[4]],
                                           AnalyseDes4_25a[[4]],AnalyseDes5_25a[[4]])))


APEV <- list()
APEV[["1"]] <- cbind(facs, 
                    data.frame(PEV = c(AnalyseDes0_04[[3]],AnalyseDes1_04[[3]],
                                        AnalyseDes2_04[[3]],AnalyseDes3_04[[3]],
                                        AnalyseDes4_04[[3]],AnalyseDes5_04[[3]])))

APEV[["2"]] <-  cbind(facs, 
                     data.frame(PEV = c(AnalyseDes0_04[[3]],AnalyseDes1_04a[[3]],
                                         AnalyseDes2_04a[[3]],AnalyseDes3_04a[[3]],
                                         AnalyseDes4_04a[[3]],AnalyseDes5_04a[[3]])))


APEV[["4"]] <- cbind(facs, 
                    data.frame(PEV = c(AnalyseDes0_25[[3]],AnalyseDes1_25[[3]],
                                        AnalyseDes2_25[[3]],AnalyseDes3_25[[3]],
                                        AnalyseDes4_25[[3]],AnalyseDes5_25[[3]])))

APEV[["3"]] <- cbind(facs, 
                    data.frame(PEV = c(AnalyseDes0_25[[3]],AnalyseDes1_25a[[3]],
                                        AnalyseDes2_25a[[3]],AnalyseDes3_25a[[3]],
                                        AnalyseDes4_25a[[3]],AnalyseDes5_25a[[3]])))


simil <- list()
simil[["1"]] <- cbind(facs, 
                         data.frame(simil = c(AnalyseDes0_04[[5]]/3*10,AnalyseDes1_04[[5]]/3*10,
                                            AnalyseDes2_04[[5]]/3*10,AnalyseDes3_04[[5]]/3*10,
                                            AnalyseDes4_04[[5]]/3*10,AnalyseDes5_04[[5]]/3*10)))

simil[["2"]] <-  cbind(facs, 
                      data.frame(simil = c(AnalyseDes0_04[[5]]/3*10,AnalyseDes1_04a[[5]]/3*10,
                                         AnalyseDes2_04a[[5]]/3*10,AnalyseDes3_04a[[5]]/3*10,
                                         AnalyseDes4_04a[[5]]/3*10,AnalyseDes5_04a[[5]]/3*10)))


simil[["4"]] <- cbind(facs, 
                     data.frame(simil = c(AnalyseDes0_25[[5]]/3*10,AnalyseDes1_25[[5]]/3*10,
                                        AnalyseDes2_25[[5]]/3*10,AnalyseDes3_25[[5]]/3*10,
                                        AnalyseDes4_25[[5]]/3*10,AnalyseDes5_25[[5]]/3*10)))

simil[["3"]] <- cbind(facs, 
                     data.frame(simil = c(AnalyseDes0_25[[5]]/3*10,AnalyseDes1_25a[[5]]/3*10,
                                        AnalyseDes2_25a[[5]]/3*10,AnalyseDes3_25a[[5]]/3*10,
                                        AnalyseDes4_25a[[5]]/3*10,AnalyseDes5_25a[[5]]/3*10)))


simul.anal <- function(measure, results)
{
  cat("\n\n#### Results of mixed model analysis\n\n")
  fix.mod <- as.formula(paste(measure, "~ Design"))
  simul.asr <- do.call(asreml, 
                       list(fix.mod, random = ~ simul, data = results))
  print(summary(simul.asr)$varcomp)
  wald.tab <- wald(simul.asr, denDF = "numeric")$Wald
  print(wald.tab)
  cat("\n\n#### Diagnostice checking")
  fit <- fitted(simul.asr)
  res <- residuals(simul.asr)
  plot(x= fit, y = res)
  QTLRel::qqPlot(y = res, sd = sqrt(simul.asr$sigma2))
  cat("#### Predictions for the designs\n\n")
  diffs <- predictPlus(simul.asr, classify = "Design", wald.tab = wald.tab, 
                       tables = "predictions", which.predictions = c("title", "table"),
                       error.intervals = "halfLeast")
  ggplot(diffs$predictions, aes(y = predicted.value, x = Design)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower.halfLeastSignificant.limit,
                      ymax = upper.halfLeastSignificant.limit))
  plotPvalues(diffs, show.sig = TRUE)
  cat("#### Design medians\n\n")
  print(aggregate(results[measure], by = results["Design"], FUN = median, na.rm = TRUE))
  invisible(diffs)
}

### Analyse
correl.diffs <- lapply(correl, simul.anal, measure  = "corr")
RGG.diffs <- lapply(RGG, simul.anal, measure  = "RGG")
PEV.diffs <- lapply(APEV, simul.anal, measure  = "PEV")
simil.diffs <- lapply(simil, simul.anal, measure  = "simil")


predPlots <- function(diffs, y.lab, row, col)
{
  preds <- do.call(rbind, lapply(diffs, 
                                 function(diffs)
                                 {
                                   diffs$predictions
                                 }))
  preds <- cbind(fac.gen(list(DesignRho = c("+","-"), ModelRho = c("+", "-")), each = 6),
                 preds)
  plt <- ggplot(preds, aes(y = predicted.value, x = Design)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower.halfLeastSignificant.limit,
                      ymax = upper.halfLeastSignificant.limit)) + 
    facet_grid(rows = vars(ModelRho), cols = vars(DesignRho), 
               labeller = label_bquote(rows = Model~rho[r]:~.(levels(ModelRho)[ModelRho]),
                                       cols = Design~rho[r]:~.(levels(DesignRho)[DesignRho]))) + 
    ylab(y.lab) + theme_bw() + theme(strip.text = element_text(size = 12))
  print(plt, vp=viewport(layout.pos.row=row, layout.pos.col=col))
  
  return(preds)  
}

#+ SimulResults
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))

preds.correl <- predPlots(correl.diffs, y.lab =  "r_a", 1, 1)
preds.RGG <- predPlots(RGG.diffs, y.lab = "RGG", 1, 2)
preds.PEV <- predPlots(PEV.diffs, y.lab = "APEV", 2, 1)
preds.simil <- predPlots(simil.diffs, y.lab = "Similarity (%)", 2, 2)
