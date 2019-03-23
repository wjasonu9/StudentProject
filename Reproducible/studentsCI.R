#function for imputing in K=2 factorial design
imputeY=function(x1,x2,Yobs){ #need to append x1 and x2 to the output
  Yobs=c(0,Yobs)
  A0=c(x1,x2,0,0)
  A1=matrix(c(0,x2,x1,x1+x2),nrow=10,ncol=4,byrow=T)
  A2=matrix(c(-x2,0,x1-x2,x1),nrow=10,ncol=4,byrow=T)
  A3=matrix(c(-x1,x2-x1,0,x2),nrow=10,ncol=4,byrow=T)
  A4=matrix(c(-x1-x2,-x1,-x2,0),nrow=10,ncol=4,byrow=T)
  return(cbind(Yobs,Yobs,Yobs,Yobs)+rbind(A0,A1,A2,A3,A4)) }
#compute X2 (except scaling by N), 2-level factorial design, given a permutation.
PermFact=function(Y,perm,effect){
  x1=Y[1,1]; x2=Y[1,2] #1st row of Y must contain x1,x2
  counts=rep(1:4,each=10)
  countsf=factor(counts)
  Yobs=Y[cbind(1+perm,counts)]
  temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x)))
  Ybar=temp$Y[,1]
  Sob2=temp$Y[,2]
  tao=matrix(c(-1,-1,-1,1,1,-1,1,1),nrow=2)%*%Ybar/2
  if(effect<1.1)
    return((tao[1]-x1)^2/sum(Sob2))
  else if(effect<2.1)
    return((tao[2]-x2)^2/sum(Sob2))
  else
    return((x1+x2-tao[1]-tao[2])^2/(Sob2[1]+Sob2[4])+(x1-x2-tao[1]+tao[2])^2/(Sob2[2]+Sob2[3])) }

load('raw.RData')
set.seed(1) #for reproducibility
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
act11=(x1-tao[1])^2/sum(Sob2); act12=(x2-tao[2])^2/sum(Sob2) #observed X2, except multiplier of N
act2=(xx1+xx2-tao[1]-tao[2])^2/(Sob2[1]+Sob2[4])+(xx1-xx2-tao[1]+tao[2])^2/(Sob2[2]+Sob2[3])
Yimp2=mapply(FUN=imputeY,xx1,xx2,MoreArgs=list(Yobs=Yobs),SIMPLIFY=F,USE.NAMES=F)
Yimp11=Yimp2[55*1:55]; Yimp12=Yimp2[1:55]
hyp2=rep(0,3025); hyp11=rep(0,55); hyp12=hyp11
for(i in 1:2500){
  p1=sample(40)
  s2=sapply(Yimp2,FUN=PermFact,perm=p1,effect=3)>act2
  s11=sapply(Yimp11,FUN=PermFact,perm=p1,effect=1)>act11
  s12=sapply(Yimp12,FUN=PermFact,perm=p1,effect=2)>act12
  hyp2=hyp2+s2; hyp11=hyp11+s11; hyp12=hyp12+s12 }
save(hyp2,hyp11,hyp12,file='CI.RData')
print(1)
