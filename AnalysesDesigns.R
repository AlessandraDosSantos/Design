require(asreml)
require(dae)

##############################################################################
####################   Analysis ##############################################
##############################################################################
banded <- function(order, kappa) {
  H <- mat.banded(c(0,1),order,order)
  V <- kappa^H - (mat.J(order)-mat.banded(c(1,1),order,order))
  ## derivative
  dV <- H*(kappa^(H-1))
  return(list(V, dV))
}

### for design 0, 1, 3 and 4

fitted.simulation <- function(opDesign,environmental){
  convergeModel <- 0
  correlation <- rep(NA,1000)
  APEV <- rep(NA,1000)
  relativeGain <- rep(NA,1000)
  similarity <- rep(0,1000)
  parametros <- matrix(NA,ncol=6,nrow = 1000)
  criteria<- list(convergeModel,correlation, APEV, relativeGain,
                  similarity,parametros)
  Z <- model.matrix(~ -1+opDesign$New)[,-c(415)]
  for(i in 1:1000)
  {
    opDesign$Gen <- Z%*%genetic[,i]
    opDesign$Y <- opDesign$Gen+environmental[,i]
    model_op <- asreml(Y ~ Control/Check, random =~ New + Column + Row,
                       residual = ~ ar1(Column):own(Row,"banded",0.1,"R"),
                       aom=T,trace=F,data = opDesign,maxit=30,gammaPar=TRUE)
    if(model_op$converge == TRUE)
    {
      criteria[[1]][1] <- criteria[[1]][1] +1
      BESTX <- tapply(opDesign$Gen,opDesign$New,mean)[-415]
      selec <- names(BESTX[order(BESTX,decreasing=T)][1:30])
      EBLUP1 <- coef(model_op)$random[46:459,1]
      # correlation
      criteria[[2]][i] <- cor(BESTX, EBLUP1) 
      # APEV
      criteria[[3]][i] <- sum((BESTX - EBLUP1)^2)/414
      # relative gain
      criteria[[4]][i] <- mean(EBLUP1[order(EBLUP1,decreasing=T)][1:30])
      # parametrs
      criteria[[6]][i,] <- summary(model_op)$varcomp[,1]
      pred<-predict(model_op,classify = "New",
                    present = c("Control","Check","New"),data=opDesign)$pvals[1:414,]
      selec1 <- pred[order(pred[,2],decreasing=T),][1:30,1]
      # similarity
      {
        for(a in 1:30)
        {
          for(b in 1:30)
          {
            if(selec[a] == selec1[b])
              criteria[[5]][i] <- criteria[[5]][i]+1
            else
              criteria[[5]][i] <- criteria[[5]][i]
          }
        }
      }
    } 
  }
  return(criteria)
}


### Design 0
summary(syst)
AnalyseDes0_04 <- fitted.simulation(opDesign=syst,environmental=uni1)
AnalyseDes0_25 <- fitted.simulation(opDesign=syst,environmental=uni2)

### Design 1
summary(des1_04)
AnalyseDes1_04 <- fitted.simulation(opDesign=des1_04,environmental=uni1)
AnalyseDes1_04a <- fitted.simulation(opDesign=des1_04,environmental=uni2)
summary(des1_25)
AnalyseDes1_25 <- fitted.simulation(opDesign=des1_25,environmental=uni2)
AnalyseDes1_25a <- fitted.simulation(opDesign=des1_25,environmental=uni1)

### Design 3
summary(des3_04)
des3_04$Control <- fac.recode(des3_04$Check,c(1,1,1,2),labels=c("Test","New"))
AnalyseDes3_04 <- fitted.simulation(opDesign=des3_04,environmental=uni1)
AnalyseDes3_04a <- fitted.simulation(opDesign=des3_04,environmental=uni2)
summary(des3_25)
des3_25$Control <- fac.recode(des3_25$Check,c(1,1,1,2),labels=c("Test","New"))
AnalyseDes3_25 <- fitted.simulation(opDesign=des3_25,environmental=uni2)
AnalyseDes3_25a <- fitted.simulation(opDesign=des3_25,environmental=uni1)

### Design 4
summary(des4_04)
des4_04$Control <- fac.recode(des4_04$Check,c(1,1,1,2),labels=c("Test","New"))
AnalyseDes4_04 <- fitted.simulation(opDesign=des4_04,environmental=uni1)
AnalyseDes4_04a <- fitted.simulation(opDesign=des4_04,environmental=uni2)
summary(des4_25)
des4_25$Control <- fac.recode(des4_25$Check,c(1,1,1,2),labels=c("Test","New"))
AnalyseDes4_25 <- fitted.simulation(opDesign=des4_25,environmental=uni2)
AnalyseDes4_25a <- fitted.simulation(opDesign=des4_25,environmental=uni1)
