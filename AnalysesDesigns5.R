require(asreml)
require(dae)
require(parallel)
require(doParallel)
require(foreach)
require(asremlPlus)

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

### for design 5

fitted5simul <- function(opDesign,environmental){
  convergeModel <- 0
  correlation <- rep(NA,1000)
  APEV <- rep(NA,1000)
  relativeGain <- rep(NA,1000)
  similarity <- rep(0,1000)
  parametros <- matrix(NA,ncol=6,nrow = 1000)
  criteria<- list(convergeModel,correlation, APEV, relativeGain,
                  similarity,parametros)
  Z <- model.matrix(~ -1+opDesign$Treat5)
  for(i in 1:1000)
  {
    opDesign$Gen <- Z%*%genetic[,i]
    opDesign$Y <- opDesign$Gen+environmental[,i]
    model_op <- asreml(Y ~ 1, random =~ Treat5 + Column + Row,
                       residual = ~ ar1(Column):own(Row,"banded",0.1,"R"),
                       aom=T,trace=F,data = opDesign,maxit=30,gammaPar=TRUE)
    if(model_op$converge == TRUE)
    {
      criteria[[1]][1] <- criteria[[1]][1] +1
      BESTX <- tapply(opDesign$Gen,opDesign$Treat5,mean)
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
      pred<-predict(model_op,classify = "Treat5",data=opDesign)$pvals
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

summary(des5_04)
AnalyseDes5_04 <- fitted5simul(opDesign=des5_04,environmental=uni1)
AnalyseDes5_04a <- fitted2simul(opDesign=des5_04,environmental=uni2)
summary(des5_25)
AnalyseDes5_25 <- fitted5simul(opDesign=des5_25,environmental=uni2)
AnalyseDes5_25a <- fitted5simul(opDesign=des5_25,environmental=uni1)