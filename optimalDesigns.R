library(knitr)
library(dae)
packageVersion("dae")
library(od)
packageVersion("od")
require(MASS)
require(asreml)
packageVersion("asreml")
require(asremlPlus)
packageVersion("asremlPlus")
require(grid)
require(ggplot2)


############# Set up for all designs
Column<- factor(rep(1:25,each=20))
Row<- factor(rep(1:20,25))
###########################################################################
####################### Construct the six designs  #######################
###########################################################################

#############################  Design 0  #################################
## Setup the data.frame syst for the systematic design factors: Treat, New, Control and Check
treat <- c(1:500)
a <- c(2,76,103,156,203,278,324,398,423,477)
b <- c(3,77,105,159,205,276,325,399,422,479)
d <- c(4,78,104,158,204,277,323,397,425,478)
treat[a] <- 2
treat[b] <- 3 
treat[d] <- 4
k <- c(0:19)
m <- c(1+21*k,10+21*(k[1:11]),19,40,181+21*k[1:16],361+21*k[1:7])
treat[m] <- 1
Treat <- as.factor(treat)
New <- factor(ifelse(Treat =="1" | Treat=="2" | Treat =="3" |
                       Treat =="4","Check",Treat))
Control <- factor(ifelse(Treat =="1" | Treat=="2" | Treat =="3" | 
                           Treat =="4","Check","New"))
Check <- factor(ifelse(Treat =="1" | Treat=="2" | Treat =="3" | 
                         Treat =="4",Treat,"New"))
syst <- data.frame(Column,Row,Control,Check,New,Treat)
exp.ini <- od(fixed = ~ Check,
              random=~New + Column + Row,
              residual=~ar1(Column):ar1(Row), 
              start.values = TRUE, data=syst)
vp.tab <- exp.ini$vparameters.table

#########################  optimal design function ##########################
od.options(P = 0.01, localSearch = 10000)
design<- function(corR,design.used){
    Ginit <- c(1, 0.3, 0.3)
    Rinit <- c(0.5, corR)
    vp.rand <- vp.tab
    vp.rand$Value = c(Ginit,1,Rinit)
    exp.rand <- od(fixed = ~ Check,
                   random=~New + Column + Row,
                   residual=~ar1(Column):ar1(Row), 
                   permute=~ New | Check, maxit=1000,
                   search= "tabu+rw",
                   G.param = vp.rand, R.param = vp.rand, 
                   data=design.used)
  return(exp.rand$design)
}

######## for design 1 #################################
des1_04<- design(corR=0.4,design.used=syst)
des1_25<- design(corR=-0.25,design.used=syst)

########  for design 2  #################################
New <- sample(as.factor(c(1:414,rep("T",86))))
Check <- factor(ifelse(New =="T","Check","New"))
Des2 <- data.frame(Column,Row,Check, New)

des2_04<- design(corR=0.4,design.used=Des2)
des2_25<- design(corR=-0.25,design.used=Des2)

########  for design 3  #################################
Treat3 <- c(sample(c(rep(c("A","B","C"),c(29,29,28)),1:414)))
Check <- ifelse(Treat3=="A"| Treat3 == "B"| Treat3 == "C", Treat3,"New")
New <- ifelse(Treat3=="A"| Treat3 == "B"| Treat3 == "C", "T",Treat3)
Des3 <- data.frame(Column,Row,Check, New)

des3_04<- design(corR=0.4,design.used=Des3)
des3_25<- design(corR=-0.25,design.used=Des3)

########  for design 4  #################################
Treat4 <- c(sample(c(rep(c("A","B","C"),5),1:235)),
            sample(c(rep(c("A","B","C"),5),1:56,236:414)))
Check <- ifelse(Treat4=="A"| Treat4 == "B"| Treat4 == "C", Treat4,"New")
New <- ifelse(Treat4=="A"| Treat4 == "B"| Treat4 == "C", "T",Treat4)
Des4 <- data.frame(Column,Row,Check, New)

des4_04<- design(corR=0.4,design.used=Des4)
des4_25<- design(corR=-0.25,design.used=Des4)

########  for design 5  #################################
Treat5<-factor(c(sample(1:250),sample(c(1:86,251:414))))
Des5<- data.frame(Column,Row,Treat5)

design5<- function(corR){
    Ginit <- c(1, 0.3, 0.3)
    Rinit <- c(0.5, corR)
    vp.rand5 <- vp.tab
    vp.rand5$Value = c(Ginit,1,Rinit)
    exp.rand5 <- od(fixed = ~ 1, random=~Treat5 + Column + Row,
                    residual=~ar1(Column):ar1(Row), 
                    permute=~ Treat5, maxit=1000,
                    search= "tabu+rw",G.param = vp.rand5,
                    R.param = vp.rand5,data=Des5)
  return(exp.rand5$design)
}

des5_04<- design5(corR=0.4)
des5_25<- design5(corR=-0.25)

######## A function to calculate APEVs and AVPDs for one or more values of rho_Row for a design
calcAvals <- function(rho_R, used.design)
{
  Criteria <- c(rep(NA,4))
  Vp <- mat.Vpredicts(target = ~ New - 1, Gt = 1,
                      fixed = ~ Check - 1,
                      random = ~ Column + Row -1 , 
                      G = list(Columns = 0.3, Row = 0.3),
                      R = kronecker(mat.ar1(0.5,25),mat.ar1(rho_R,20)),
                      design = used.design)
  Criteria[1] <- mean(diag(Vp))
  Criteria[2] <- mean(diag(Vp[,-415]))
  Criteria[3] <- designAmeasures(Vp)
  Criteria[4] <- designAmeasures(Vp, groupsizes = 414)
  return(Criteria)
}

######## Calculate the APEVs and AVPDs for all design
Measures <- matrix(NA,12,4)
colnames(Measures) <- c("APEV_New+Check","APEV_New",
                        "AVPD_New+Check","AVPD_New")
rownames(Measures) <- c(rep(c("rho_r=0.4","rho_r=-0.25"),each=6))
Measures[1,] <- calcAvals(rho_R = 0.4, used.design = syst)
Measures[2,] <- calcAvals(rho_R = 0.4, used.design = des1_04)
Measures[3,] <- calcAvals(rho_R = 0.4, used.design = des2_04)
Measures[4,] <- calcAvals(rho_R = 0.4, used.design = des3_04)
Measures[5,] <- calcAvals(rho_R = 0.4, used.design = des4_04)
Vp5 <- mat.Vpredicts(target = ~ Treat5 - 1, Gt = 1,
                     random = ~ Column + Row -1 , 
                     G = list(Columns = 0.3, Row = 0.3),
                     R = kronecker(mat.ar1(0.5,25),mat.ar1(0.4,20)),
                     design = des5_04)
Measures[6,] <- c(NA,mean(diag(Vp5)), NA, designAmeasures(Vp5))
# for -0.25
Measures[7,] <- calcAvals(rho_R = -.25, used.design = syst)
Measures[8,] <- calcAvals(rho_R = -.25, used.design = des1_25)
Measures[9,] <- calcAvals(rho_R = -.25, used.design = des2_25)
Measures[10,] <- calcAvals(rho_R =-.25, used.design = des3_25)
Measures[11,] <- calcAvals(rho_R =-.25, used.design = des4_25)
Vp5a <- mat.Vpredicts(target = ~ Treat5 - 1, Gt = 1,
                      random = ~ Column + Row -1 , 
                      G = list(Columns = 0.3, Row = 0.3),
                      R = kronecker(mat.ar1(0.5,25),mat.ar1(-.25,20)),
                      design = des5_25)
Measures[12,] <- c(NA,mean(diag(Vp5a)), NA, designAmeasures(Vp5a))
Measures

###### Comparing to Design 0
Effic0 <- matrix(NA,12,4)
for(i in 1:6)
{
  Effic0[i,] <- Measures[1,]/Measures[i,]
  Effic0[i+6,] <- Measures[7,]/Measures[i+6,]
}
Effic0

####### Layout for design 1
des1_04$Cplot <- fac.recode(des1_04$Check,c(1:5),labels=c("T","A","B","C",""))
des1_25$Cplot <- fac.recode(des1_25$Check,c(1:5),labels=c("T","A","B","C",""))

designGGPlot(design = des1_25, labels = "Cplot",
             row.factors = "Row", column.factors = "Column",
             colour.values = c("lightblue","lightcoral","yellow",
                               "lightgreen","white"),celllinesize=0.3,
             title="", axis.text.size = 10,title.size=10,
             cellfillcolour.column = "Check", cellalpha = 0.5)

####### Layout for design 2
des2_04$Cplot <- fac.recode(des2_04$Check,c(1:2),labels=c("T",""))
des2_25$Cplot <- fac.recode(des2_25$Check,c(1:2),labels=c("T",""))

designGGPlot(design = des2_04, labels = "Cplot",
             row.factors = "Row", column.factors = "Column",
             colour.values = c("lightblue","white"),celllinesize=0.3,
             title="", axis.text.size = 10,title.size=10,
             cellfillcolour.column = "Check", cellalpha = 0.5)

####### Layout for design 3
des3_04$Cplot <- fac.recode(des3_04$Check,c(1:4),labels=c("A","B","C",""))
des3_25$Cplot <- fac.recode(des3_25$Check,c(1:4),labels=c("A","B","C",""))

designGGPlot(design = des3_04, labels = "Cplot",
             row.factors = "Row", column.factors = "Column",
             colour.values = c("lightcoral","yellow","lightgreen","white"),
             celllinesize=0.3,title="", axis.text.size = 10,title.size=10,
             cellfillcolour.column = "Check", cellalpha = 0.5)

####### Layout for design 4
clones <- droplevels(subset(des4_04,des4_04$New != "T"))
sorted_labels <- paste(sort(as.integer(levels(droplevels(clones)$New))))
clones$label <- factor(clones$New, levels = c(sorted_labels))
clones$label <- fac.recode(clones$label,c(1:56,rep(57,358)),labels=c(1:56,""))
Test <- subset(des4_04,des4_04$New == "T")
Test$label <- droplevels(Test)$Check
dados <- rbind(clones,Test)
dados <- dados[order(dados$Column,dados$Row),]

designGGPlot(design = dados, labels = "label",
             row.factors = "Row", column.factors = "Column",
             colour.values = c(rep("orange",56),"white","lightblue","lightcoral","yellow","lightgreen"),
             celllinesize=0.3,title="", axis.text.size = 10,title.size=10,
             cellfillcolour.column = "label", cellalpha = 0.5)

####### Layout for design 5
des5_25$Nplot <- fac.recode(des5_25$Treat5,c(1:86,rep(87,328)),
                            labels=c(1:86,""))

designGGPlot(design = des5_25, labels = "Nplot",
             row.factors = "Row", column.factors = "Column",
             colour.values = c(rep("orange",86),
                               rep("white",328)),celllinesize=0.3,
             title="", axis.text.size = 10,title.size=10,
             cellfillcolour.column = "Nplot", cellalpha = 0.5)
