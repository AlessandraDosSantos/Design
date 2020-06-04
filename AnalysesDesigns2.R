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

### for design 0, 1, 3 and 4

fitted2simul <- function(opDesign,environmental){
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
    model_op <- asreml(Y ~ Check, random =~ New + Column + Row,
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
                    present = c("Check","New"),data=opDesign)$pvals[1:414,]
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

summary(des2_04)
AnalyseDes2_04 <- fitted2simul(opDesign=des2_04,environmental=uni1)
AnalyseDes2_04a <- fitted2simul(opDesign=des2_04,environmental=uni2)
summary(des2_25)
AnalyseDes2_25 <- fitted2simul(opDesign=des2_25,environmental=uni2)
AnalyseDes2_25a <- fitted2simul(opDesign=des2_25,environmental=uni1)



################# another function
cl <- makeCluster(6)
#registerDoParallel(cl)
registerDoSEQ()
anal.out.S1<- NULL
bootsamp = 1:1000
anal.out.S1 <- foreach(i = bootsamp, .inorder=FALSE, 
                       .packages = c("asreml"))  %dopar%
  {
    cat("#### Simulation ",i,"\n")
    Z <- model.matrix(~ -1+des2_04$New)[,-c(415)]
    des2_04$Gen <- Z%*%genetic[,i]
    des2_04$Y <- des2_04$Gen+uni1[,i]
    asreml(Y ~ Check, random =~ New + Column + Row,
           residual = ~ ar1(Column):own(Row,"banded",0.1,"R"),
           aom=T,trace=F,data = des2_04,maxit=30,gammaPar=TRUE)
  }
stopCluster(cl)


CV.S3 <- rbind(summary(anal.out.S1[[1]])$varcomp$component,
               summary(anal.out.S1[[2]])$varcomp$component,
               summary(anal.out.S1[[3]])$varcomp$component,
               summary(anal.out.S1[[4]])$varcomp$component,
               summary(anal.out.S1[[5]])$varcomp$component)
colnames(CV.S3) <- c("sima2r", "sigma2c", "sigma2g", "sigma2", "phic", "band1r", "band2r")
CV.S3 <- as.data.frame(CV.S3)

CV.S3.df <- data.frame(
  Amostra = rep(1:5, 7),
  Parametro = rep(colnames(CV.S3), each = 5),
  Estimativa = c(CV.S3[,1],
                 CV.S3[,2],
                 CV.S3[,3],
                 CV.S3[,4],
                 CV.S3[,5],
                 CV.S3[,6],
                 CV.S3[,7]),
  Verdadeiro = c(rep(c(0.3,0.3,0.8,1,0.5,-0.4,0.3),
                     each=5))
)

CV.S3.df
