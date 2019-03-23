# setwd("Downloads/Causal/Coding")
# par(family='serif')

#File to create all figures after running simulations on cluster
load('Anova.RData')
brk=c(0:10/50,3:10/10)
X2=result[1,]; Fstat=result[2,]
pX2=1-pchisq(X2,df=2); pF=1-pf(Fstat,df1=2,df2=57)
pX2f=result[3,]; pFf=result[4,]
hist(pFf,breaks=brk,col='gray',lty='blank',xlab='p values',main='ANOVA Simulation') #BW plot
hist(pX2f,breaks=brk,col=NULL,add=T)
abline(h=1,lty=2)
legend(x=0.1,y=3/2,legend=c("Box-type stat",expression(paste(X^2," stat"))),fill=c('gray','white'),bty='n') #save as BrunnerAnova.pdf, 5in x 6in, landscape 
X2b=result[5,]; Fb=result[6,]
pX2b=1-pchisq(X2,df=2); pFb=1-pf(Fstat,df1=2,df2=117)
pX2bf=result[7,]; pFbf=result[8,]
hist(pFbf,breaks=brk,col='gray',lty='blank',xlab='p values',main='Blocking Simulation')
hist(pX2bf,breaks=brk,col=NULL,add=T)
abline(h=1,lty=2)
legend(x=0.1,y=3/2,legend=c("F stat",expression(paste(X^2," stat"))),fill=c('gray','white'),bty='n') #save as Block.pdf
plot(pX2,pX2f,col=rgb(0,0,0,0.5),pch=16,xlim=c(0,0.1),ylim=c(0,0.15),main="p values",xlab=expression(paste(chi^2," approximation")),ylab='FRT')
abline(a=0,b=1,lty=2) #save as pvscatterFact.pdf

load('Fact.RData')
hist(result[2,],breaks=brk,col='gray',lty='blank',xlab='p values',main='Factorial Simulation')
hist(result[1,],breaks=brk,col=NULL,add=T)
abline(h=1,lty=2)
legend(x=0.1,y=3/2,legend=c("Box-type Stat",expression(paste(X^2," stat"))),fill=c('gray','white'),bty='n') #save as BrunnerFact.pdf

load('raw.RData')
ns=rep(10,4)
counts=rep(1:4,times=ns)
countsf=factor(counts)
Yobs=Y[cbind(1:40,counts)]
temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x)))
Ybar=temp$Y[,1]
Sob2=temp$Y[,2]
C=matrix(c(-1,-1,1,1,-1,1,-1,1,1,-1,-1,1,1,1,1,1),nrow=4)/2
tao=C%*%Ybar
x1=-20:34/32; x2=-25:29/32
xx1=rep(x1,each=55); xx2=rep(x2,times=55)
#some asymptotic quantities
k1=qchisq(0.95,df=1); k2=qchisq(0.95,df=2) #=2*log(20)
b0=sqrt(k1*sum(Sob2)/40) #(half-)width of 1D CI's
b1=sqrt(k2*(Sob2[1]+Sob2[4])/40) #"size" of c(1,1) axis
b2=sqrt(k2*(Sob2[2]+Sob2[3])/40) #size of c(1,-1) axis
# act11=(x1-tao[1])^2/sum(Sob2); act12=(x2-tao[2])^2/sum(Sob2) #observed X2, except multiplier of N
# act2=(xx1+xx2-tao[1]-tao[2])^2/(Sob2[1]+Sob2[4])+(xx1-xx2-tao[1]+tao[2])^2/(Sob2[2]+Sob2[3])
# Yimp2=mapply(FUN=imputeY,xx1,xx2,MoreArgs=list(Yobs=Yobs),SIMPLIFY=F,USE.NAMES=F)
# Yimp11=Yimp2[55*1:55]; Yimp12=Yimp2[1:55]
load('CI.RData') #actual FRT stuff sent to cluster
tx1=expression(x[1]); tx2=expression(x[2])
plot(x=x1,y=hyp11/2500,pch=16,ylim=c(0,1.5),col='gray',xlab=tx1,ylab=expression(paste('p value of test ',tau[1],'=',x[1])))
title(main=expression(paste('95% CI for ',tau[1])))
abline(v=c(0,tao[1]),col='gray',lty=2) #truth and estimate
text(x=c(0.2,-0.3),y=c(0.1,0.6),labels=c('estimate','truth'),pos=4,col='gray')
lines(x=c(-19,49)/64,y=c(1.1,1.1),lwd=2)
text(x=1/4,y=1.2,labels='CI from FRT',pos=4)
lines(x=c(tao[1]-b0,tao[1]+b0),y=c(1.3,1.3),lwd=2)
text(x=1/4,y=1.4,labels='CI from asymptotics',pos=4) #save as 1D-CI1.pdf
plot(x=x2,y=hyp12/2500,pch=16,ylim=c(0,1.5),col='gray',xlab=tx2,ylab=expression(paste('p value of test ',tau[2],'=',x[2])))
title(main=expression(paste('95% CI for ',tau[2])))
abline(v=c(0,tao[2]),col='gray',lty=2) #truth and estimate
text(x=c(1/20,-1/4),y=c(0.1,0.1),labels=c('estimate','truth'),pos=4,col='gray')
lines(x=c(-31,37)/64,y=c(1.1,1.1),lwd=2)
text(x=0.1,y=1.2,labels='CI from FRT',pos=4)
lines(x=c(tao[2]-b0,tao[2]+b0),y=c(1.3,1.3),lwd=2)
text(x=0.1,y=1.4,labels='CI from asymptotics',pos=4) #save as 1D-CI2.pdf
A2=matrix(cut(hyp2,breaks=c(-1,24,124,249,2500),labels=F),nrow=55,byrow=T) #cut changes hyp2 to 1 for sig at 0.01, 2 for sig at 0.05, 3 for sig at 0.1, 4 for not sig
Ae=matrix(c(b1,b1,b2,-b2),nrow=2)
t=pi*(0:100)/50
elipse=Ae%*%rbind(cos(t),sin(t))
image(x=c(x1,35:52/32),y=x2,z=rbind(A2,matrix(0,nrow=18,ncol=55)),xlab=tx1,ylab=tx2,col=c('white','white','gray','gray'),asp=1)
title(main='95% Confidence Region')
points(elipse[1,]+tao[1],elipse[2,]+tao[2],type='l') #asymptotic CI outline
points(x=c(0,tao[1]),y=c(0,tao[2]),pch=c(16,17)) #add estimate and truth
text(x=tao[1],y=tao[2],labels='estimate',pos=4)
text(x=0,y=0,labels='truth',pos=2)
legend(x=4/5,y=7/8,legend=c('from FRT','asymptotic'),fill=c('gray','white'),border=c('gray','black')) #save as 2D-FRTvsAsymp.pdf

#some functions required to analyze actual datasets
#compute X2 and F, special case of 1D hypotheses. For treatment-control, use C=c(1,-1)
PermDist1D=function(Yobs,countsf,ns,C,shuffle){
  if(shuffle)
    Yobs=sample(Yobs)
  temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x)))
  Ybar=temp$Y[,1]; Sob2=temp$Y[,2]
  m2=sum(C*Ybar)^2
  C2=C^2
  Fstat=m2/weighted.mean(Sob2,ns-1)/sum(C2/ns)
  Sob2=Sob2/ns
  X2=m2/sum(Sob2*C2)
  return(c(X2,Fstat)) }
#compute X2,B, and/or F, special case of 2-level factorial design, testing for no main effects
PermDistFact=function(Yobs,countsf,ns,shuffle){
  if(shuffle)
    Yobs=sample(Yobs)
  temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x))) #calculate mean and variance of values of Y subsetted by levels of countsf
  Ybar=temp$Y[,1]; Sob2=temp$Y[,2]
  m14=(Ybar[1]-Ybar[4])^2; m23=(Ybar[2]-Ybar[3])^2
  Fstat=(m14/(1/ns[1]+1/ns[4])+m23/(1/ns[2]+1/ns[3]))/weighted.mean(Sob2,ns-1)
  Sob2=Sob2/ns
  X2=m14/(Sob2[1]+Sob2[4])+m23/(Sob2[2]+Sob2[3])
  #B=(m14+m23)/sum(Sob2)
  return(c(X2,Fstat)) }
#compute X2 and B, special case of ANOVA
PermDistANOVA=function(Yobs,countsf,ns,shuffle){
  if(shuffle)
    Yobs=sample(Yobs)
  J=length(ns)
  temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x))) #calculate mean and variance of values of Y subsetted by levels of countsf
  Ybar=temp$Y[,1]; Sob2=temp$Y[,2]
  gm=weighted.mean(Ybar,ns)
  Fstat=sum(ns*(Ybar-gm)^2)/weighted.mean(Sob2,ns-1)
  Q=ns/Sob2
  Yobw=weighted.mean(Ybar,Q)
  X2=sum(Q*(Ybar-Yobw)^2)
  #B=(J*sum(Ybar^2)-sum(Ybar)^2)/(J-1)/sum(Sob2)
  return(c(X2,Fstat)) }

set.seed(1) #reproducibility, even though there should be very little difference if this line is not run
#Application to exercise data
IE=read.csv("IncentiveExercise.csv")
IEdiff=IE$After-IE$Before+rnorm(120,sd=1e-3) #jitter data
ns=rep(40,3)
counts=rep(1:3,times=ns); countsf=factor(counts)
temp=aggregate(IEdiff~countsf,FUN=function(x) c(mean(x),var(x)))
temp #print means and variances
#hypthesis 2,-1,-1
TS=PermDist1D(IEdiff,countsf,ns,c(2,-1,-1),F)
1-pchisq(TS[1],df=1); 1-pf(TS[2],df1=1,df2=117)
FRT=replicate(n=1e4,PermDist1D(IEdiff,countsf,ns,c(2,-1,-1),T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])
#hypothesis ANOVA
TS=PermDistANOVA(IEdiff,countsf,ns,F)
1-pchisq(TS[1],df=2); 1-pf(TS[2]/2,df1=2,df2=117) #p values from theoretical quantiles
FRT=replicate(n=1e4,PermDistANOVA(IEdiff,countsf,ns,T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2]) #p values from FRT
mod=lm(IEdiff~countsf); summary(mod) #sanity check to make sure F stat and p value computed correctly
#hypothesis 1 vs 2,3
ns=c(40,80)
counts=rep(1:2,times=ns); countsf=factor(counts)
TS=PermDist1D(IEdiff,countsf,ns,c(1,-1),F)
1-pchisq(TS[1],df=1); 1-pf(TS[2],df1=1,df2=118)
FRT=replicate(n=1e4,PermDist1D(IEdiff,countsf,ns,c(1,-1),T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])
#hypothesis 1 vs 2
ns=c(40,40)
counts=rep(1:2,times=ns); countsf=factor(counts)
TS=PermDist1D(IEdiff[1:80],countsf,ns,c(1,-1),F) #Fstat=X2
1-pchisq(TS[1],df=1); 1-pf(TS[2],df1=1,df2=78)
FRT=replicate(n=1e4,PermDist1D(IEdiff[1:80],countsf,ns,c(1,-1),T))
mean(FRT[1,]>TS[1]) #don't need to re-compute for F b/c of equivalence

#Application to college data
stardata=read.csv("star.csv")
stardata$treat=1
stardata[stardata$sfp==1,7]=2
stardata[stardata$ssp==1,7]=3
stardata[stardata$sfsp==1,7]=4
stardata=stardata[order(stardata$treat),] #sort by treat
gr=na.omit(stardata[,c(3,7)])
counts=gr[,2]; countsf=factor(counts)
gr=gr[,1]
ns=c(854,219,212,119) #from summary(countsf), N=1404
temp=aggregate(gr~countsf,FUN=function(x) c(mean(x),var(x)))
temp
#hypothesis no effect 1 (ssp, services)
TS=PermDist1D(gr,countsf,ns,c(-1,-1,1,1),F)
1-pchisq(TS[1],df=1); 1-pf(TS[2],df1=1,df2=1400)
FRT=replicate(n=1e4,PermDist1D(gr,countsf,ns,c(-1,-1,1,1),T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])
#hypothesis no effect 2 (sfp, incentives)
TS=PermDist1D(gr,countsf,ns,c(-1,1,-1,1),F)
1-pchisq(TS[1],df=1); 1-pf(TS[2],df1=1,df2=1400)
FRT=replicate(n=1e4,PermDist1D(gr,countsf,ns,c(-1,1,-1,1),T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])
#hypothesis no main effects
TS=PermDistFact(gr,countsf,ns,F)
1-pchisq(TS[1],df=2); 1-pf(TS[2]/2,df1=2,df2=1400)
FRT=replicate(n=1e4,PermDistFact(gr,countsf,ns,T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])
#hypothesis no interaction
TS=PermDist1D(gr,countsf,ns,c(1,-1,-1,1),F)
1-pchisq(TS[1],df=1); 1-pf(TS[2],df1=1,df2=1400)
FRT=replicate(n=1e4,PermDist1D(gr,countsf,ns,c(1,-1,-1,1),T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])
#hypothesis ANOVA
TS=PermDistANOVA(gr,countsf,ns,F)
1-pchisq(TS[1],df=3); 1-pf(TS[2]/3,df1=3,df2=1400)
FRT=replicate(n=1e4,PermDistANOVA(gr,countsf,ns,T))
mean(FRT[1,]>TS[1]); mean(FRT[2,]>TS[2])

# par(pty='s')
# plot(x=elipse[1,]+tao[1],y=elipse[2,]+tao[2],xlim=c(-10,17)/16,ylim=c(-25,29)/32,type='l',xlab=tx1,ylab=tx2)
# title(main='Asymptotic region')
# abline(v=c(tao[1]+b0,tao[1]-b0),col='gray')
# abline(h=c(tao[2]+b0,tao[2]-b0),col='gray')
# points(x=tao[1],y=tao[2],pch=15) #save as ApproxShape.pdf