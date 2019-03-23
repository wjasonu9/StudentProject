#to be run on cluster
#edit Aug 2017
#factorial, K=2, test for no main effects
#compute X2 and B, special case of 2-level factorial design, testing for no main effects
PermDistFact=function(Yobs,countsf,ns,shuffle){
  if(shuffle)
    Yobs=sample(Yobs)
  temp=aggregate(Yobs~countsf,FUN=function(x) c(mean(x),var(x))) #calculate mean and variance of values of Y subsetted by levels of countsf
  Ybar=temp$Y[,1]
  Sob2=temp$Y[,2]/ns
  X2=(Ybar[1]-Ybar[4])^2/(Sob2[1]+Sob2[4])+(Ybar[2]-Ybar[3])^2/(Sob2[2]+Sob2[3])
  B=((Ybar[1]-Ybar[4])^2+(Ybar[2]-Ybar[3])^2)/sum(Sob2)
  return(c(X2,B)) }

#function to be replicated
TBR=function(Y,counts,countsf,ns){
  Y=Y[sample(80),]
  Yobs=Y[cbind(1:80,counts)]
  X2B=PermDistFact(Yobs,countsf,ns,F)
  result1=replicate(n=2500,PermDistFact(Yobs,countsf,ns,T))
  pvX2=mean(X2B[1]<result1[1,])
  pvB=mean(X2B[2]<result1[2,])
  return(c(pvX2,pvB)) }

set.seed(1) #added Aug 2018 for reproducibility
temp=rnorm(80)
temp=temp-mean(temp)
Y=cbind(3*temp,temp,temp,3*temp)
ns=rep(20,4)
counts=rep(1:4,times=ns)
countsf=factor(counts)
result=replicate(n=2000,TBR(Y,counts,countsf,ns))
save(result,file='Fact.RData')
