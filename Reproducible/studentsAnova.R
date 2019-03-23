#to be run on cluster (anticipated running time: 8 hrs)
#edit Aug 2018: added blocking simulation
#ANOVA, J=3
#compute X2 and B, special case of ANOVA, every 1000 calls of this fn is about 2s
PermDistANOVA=function(Yobs,countsf,ns,shuffle){
  if(shuffle)
    Yobs=sample(Yobs)
  J=length(ns)
  temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x))) #calculate mean and variance of values of Y subsetted by levels of countsf
  Ybar=temp$Y[,1]
  Sob2=temp$Y[,2]/ns
  Q=1/Sob2
  Yobw=weighted.mean(Ybar,Q)
  X2=sum(Q*(Ybar-Yobw)^2)
  B=(J*sum(Ybar^2)-sum(Ybar)^2)/(J-1)/sum(Sob2) #2x faster to get F stat by direct calculation than thru mod=aov(Yobs~countsf); B=summary(mod)[[1]]$F[1]
  return(c(X2,B)) }

#compute X2 and F, in blocking design, fn is very specialized, still takes more than 2x as long to run as its "no blocking" counterpart
CalcStatBlock=function(Yo1,Yo2,countsf,ns,shuffle){
  if(shuffle){
    Yo1=sample(Yo1); Yo2=sample(Yo2) }
  t1=aggregate(Yo1~countsf,FUN=function(x) c(mean(x),var(x))) #calculate mean and variance of values of Y subsetted by levels of countsf
  t2=aggregate(Yo2~countsf,FUN=function(x) c(mean(x),var(x)))
  mY1=t1$Y[,1]; mY2=t2$Y[,1]
  vY1=t1$Y[,2]/ns; vY2=t2$Y[,2]/ns
  Yps=(mY1+mY2)/2
  Q=4/(vY1+vY2)
  Yobw=weighted.mean(Yps,Q)
  X2=sum(Q*(Yps-Yobw)^2)
  mmY1=mean(mY1)/2; mmY2=mean(mY2)/2
  dF1=Yo1-rep(Yps,times=ns)-mmY1+mmY2
  dF2=Yo2-rep(Yps,times=ns)-mmY2+mmY1
  Fstat=232*ns[1]*var(Yps)/(sum(dF1^2)+sum(dF2^2)) #leading const is (J-1)*(N+1-J-H) in general
  return(c(X2,Fstat)) }

#function to be replicated
TBR=function(Y,countsf,ns){
  Y=Y[sample(60),]; Y2=Y[sample(60),]
  Yo1=Y[cbind(1:60,countsf)] #factor changed to numeric
  Yo2=Y2[cbind(1:60,countsf)]
  St=PermDistANOVA(Yo1,countsf,ns,F)
  result=replicate(n=2500,PermDistANOVA(Yo1,countsf,ns,T))
  pX2=mean(St[1]<result[1,])
  pF=mean(St[2]<result[2,])
  Stb=CalcStatBlock(Yo1,Yo2,countsf,ns,F)
  result=replicate(n=2500,CalcStatBlock(Yo1,Yo2,countsf,ns,T))
  pX2b=mean(Stb[1]<result[1,])
  pFb=mean(Stb[2]<result[2,])
  return(c(St,pX2,pF,Stb,pX2b,pFb)) }

set.seed(1) #added Aug 2018 for reproducibility
temp=rnorm(60)
temp=temp-mean(temp)
Y=cbind(temp,2*temp,3*temp)
ns=rep(20,3)
counts=rep(1:3,times=ns)
countsf=factor(counts)
result=replicate(n=2000,TBR(Y,countsf,ns))
save(result,file='Anova.RData')
