require(dae)
require(MASS)
########################################################################
#################### Simulation Effects ################################
########################################################################
## for environmental
unif.trial<- function(gmu, gamac,gamar, rho.c,rho.r,sim)
{
  Row <- factor(c(rep(1:20, times=25)))
  Column <- factor(c(rep(1:25, each=20)))
  variancia <- fac.vcmat(Column,gamac)+ fac.vcmat(Row,gamar)+
    kronecker(mat.ar1(rho.c,25),mat.banded(c(1,rho.r),20,20))
  Y <- matrix(0,ncol=sim,nrow=500)
  for(i in 1:sim)
  {
    Y[,i] <- mvrnorm(1,matrix(gmu,nrow=500),variancia)
  }
  return(Y)
}

uni1 <- unif.trial(gmu=90, gamac=0.3,gamar=0.3, rho.c=0.5,rho.r=0.3,sim=1000)
uni2 <- unif.trial(gmu=90, gamac=0.3,gamar=0.3, rho.c=0.5,rho.r=-0.3,sim=1000)


##### for genetic effects
gen.trial <- function(gamag,sim)
{
  Ug <- matrix(0, ncol=sim, nrow=414)
  for(i in 1:sim)
  {
    Ug[,i] <- rnorm(414,0,gamag)
  }
  return(Ug)
}

genetic <- gen.trial(gamag=1,sim=1000)
