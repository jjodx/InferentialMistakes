library(dplyr)
library(ggpubr)
library(ggplot2)
library(Hmisc)
library(MASS)
library(ppcor)
library(boot)

CorFunc <- function(d,indices){
  R <- rcorr(d$X[indices],d$Y[indices])
  return(R$r[1,2])
}

# unused example of difference in significance but no difference in effect.
set.seed(1235)  # for reproducibility
A <- rnorm(20,0.5,1)# simulate random data
t.test(A, mu=0)
B <- rnorm(20,0.5,1.5)# simulate random data with larger variance
t.test(B, mu=0)
AB=as.matrix(c(A,B))
group = factor(rep(c('A','B'),each=20));
dt = data.frame(AB,group)

t.test(A,B,var.equal = TRUE)

sumup = data.frame("Gr" = c('A','B'))

sumup<-cbind(sumup,Mean=with(dt, tapply(AB, group, mean)))
sumup<-cbind(sumup,SD=with(dt, tapply(AB, group, sd)))
sumup<-cbind(sumup,N=with(dt, tapply(AB, group, length)))
sumup<-cbind(sumup,SE =sumup$SD/sqrt(sumup$N))  # Calculate standard error of the mean
ciMult <- qt(.95/2 + .5, sumup$N-1)
sumup<-cbind(sumup,CI = sumup$SE * ciMult)

EffNieuw <- ggplot(sumup, aes(x=Gr, y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-CI, ymax=Mean+CI), width=NA,position = position_nudge(x = -0.1)) +
  geom_point(position = position_nudge(x = -0.1)) +
  geom_jitter(data=dt,aes(x=group,y=AB),color = "gray",width = 0.05)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(EffNieuw)

# figure on Nieuwenhuis with correlations
Rval <- 0.5
Sigma <- matrix(c(1,Rval,Rval,1),2,2)
for (i in c(10,26)){ #,11,8,18,22,24 --> other possible examples
  print(i)
  set.seed(1234+i)  # for reproducibility
  D <- mvrnorm(n = 20, rep(0, 2), Sigma, empirical = FALSE)
  X <- D[,1]
  Y <- D[,2]
  Add <- data.frame(X,Y)
  R<-rcorr(Add$X,Add$Y,type=c("pearson"))
  RValue = R$r[1,2]
  print(R$r[1,2])
  print(R$P[1,2])
  # boostrap CI --> CorFunc defined above
  Rbci <- boot.ci(boot(data = Add, statistic = CorFunc, R = 10000),type="norm")
  
  txt = toString(paste("Pearson R:", format(RValue,digits=2), "CI: [",format( Rbci$normal[2],digits=2),",",format(Rbci$normal[3],digits=2) ,"]"))
  NiewenCorr<-ggplot(Add, aes(X, Y)) +
    geom_point() + 
    labs(x = "X", y = "Y") +
    geom_smooth(method="lm",color = "black",linetype = 2, size=0.5, se=FALSE) +
    coord_cartesian(xlim = c(-3,3),ylim = c(-3,3)) +
    annotate("text", x = 0, y = 2, label = txt,size=2)+
    theme(legend.position="none")
  print(NiewenCorr)
}


###############################
# figure on underpower
###############################
P=NA;R=NA;Cond=NA;
LowP = data.frame(P,R,Cond)
HighP= data.frame(P,R,Cond)

for (i in 1:10000){
A <- rnorm(15,0.5,1)
B <- rnorm(15,0.5,1)
R<-rcorr(A,B)
LowP <- rbind(LowP,c(R$P[1,2], R$r[1,2],R$P[1,2]<0.05))

A <- rnorm(100,0.5,1)
B <- rnorm(100,0.5,1)
R<-rcorr(A,B)
HighP <- rbind(HighP,c(R$P[1,2], R$r[1,2],R$P[1,2]<0.05))
}
LowP$Cond = factor(LowP$Cond)
HighP$Cond = factor(HighP$Cond)

L2 <-filter(LowP, Cond == "1")
R2 <-filter(HighP, Cond == "1")
  ggplot(L2, aes(x=R)) + 
    geom_histogram(aes(x=R),binwidth=.025)+
    geom_histogram(data=R2,aes(x=R),binwidth=.025,fill="red")+
    xlim(-1, 1)

###############################
# figure on correlation
###############################

set.seed(1234)  # for reproducibility
PP = list()

# effect of outlier
X <- rnorm(20,0,1)
Y <- rnorm(20,0,1)
Z <- rep(1,20)
Dcor = data.frame(X,Y,Z)
for (i in c(1,3,5)){
  print(i)
  Add <- rbind(Dcor,c(i,i,2))
  R<-rcorr(Add$X,Add$Y,type=c("pearson"))
  RValue = R$r[1,2]
  print(R$P[1,2])
  # boostrap CI
  Rbci <- boot.ci(boot(data = Add, statistic = CorFunc, R = 10000),type="norm")

  txt = toString(paste("Pearson R:", format(RValue,digits=2), "CI: [",format( Rbci$normal[2],digits=2),",",format(Rbci$normal[3],digits=2) ,"]"))
  PP[[length(PP) + 1]]<-ggplot(Add, aes(X, Y),fill=factor(Z)) +
    scale_color_manual(values=c("black", "red")) +
    scale_fill_manual(values=c(NA, "red")) +
    geom_point(aes(colour = factor(Z),fill=factor(Z)),shape=21) + 
    labs(x = "X", y = "Y") +
    geom_smooth(method="lm",color = "black",linetype = 2, size=0.5) +
    coord_cartesian(xlim = c(-2.5,6.5),ylim = c(-2.5,6.5)) +
    annotate("text", x = 0, y = 6, label = txt,size=2)+
    theme(legend.position="none")
}


# effect of subgroups
set.seed(1234)  # for reproducibility
QQ = list()
Xbase <- rnorm(20,0,1)
Ybase <- rnorm(20,0,1)
Z <- c(rep(0,10),rep(1,10))

for (i in c(0,2,4)){
  print(i)
  X = Xbase+i*Z
  Y = Ybase+i*Z
  Add = data.frame(X,Y,Z)
  R<-rcorr(Add$X,Add$Y)
  RCor = R$r[1,2]
  print(R$P[1,2])
  # boostrap CI
  Rbci <- boot.ci(boot(data = Add, statistic = CorFunc, R = 10000),type="norm")
  
  txt = toString(paste("Pearson R:", format(RCor,digits=2), "CI: [",format( Rbci$normal[2],digits=2),",",format(Rbci$normal[3],digits=2) ,"]"))
  QQ[[length(QQ) + 1]]<-ggplot(Add, aes(X, Y),fill=factor(Z)) +
    scale_color_manual(values=c("black", "red")) +
    scale_fill_manual(values=c(NA, "red")) +
    geom_point(aes(colour = factor(Z),fill=factor(Z)),shape=21) + 
    labs(x = "X", y = "Y") +
    geom_smooth(method="lm",color = "black",linetype = 2, size=0.5) +
    coord_cartesian(xlim = c(-2.5,6.5),ylim = c(-2.5,6.5)) +
    annotate("text", x = 0, y = 6, label = txt,size=2)+
    theme(legend.position="none")
}

ggarrange(PP[[1]],PP[[2]],PP[[3]],QQ[[1]],QQ[[2]],QQ[[3]],ncol = 3, nrow = 2,labels = c("A","B","C","D","E","F"))
