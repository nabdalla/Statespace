library(nimble)
library(rjags)
library(plotrix)
library(compositions)
library(mvtnorm)
library(expm)
library(DPpackage)
library(ggvis)
library(MSGARCH)
library(data.table)

library(mltools)
library(scoringRules)
library(loo)
n=100
Qprime=13.8
Gprime=351.54
V=3.8
sigma<-0.1
aomega<-2
bomega<-1
kl=0.1
set.seed(2017)
wt<-rgamma(n,aomega,bomega)
yt<-NULL
c<-NULL
h=0.01
for(i in 1:(n-1)){
  c[1]<-1
  c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V
}
#exact solution
c2<-NULL
for(i in 1:(n-1)){
  c2[1]<-1
  c2[i+1]<-exp(-h*i*(Qprime+kl*V)/V)*c[1]+1/(-(Qprime+kl*V)/V)*
    (exp(-h*i*(Qprime+kl*V)/V)-1)*Gprime/V
}
c2
plot(c2)
plot(c2,c)
abline(0,1)
Qprime/V+kl
set.seed(000)
vt<-rnorm(n,0,sigma)

logyt<-log(c2)+(vt)
plot(c2)
points(exp(logyt),col="green")
hist(exp(logyt))
#different sims
ytsim1<-vtsim1<-vector("list")
sigmasim1<-NULL
for(i in 1:100){
sigmasim1[i]<-sigma*exp(rnorm(1,0,1))
vtsim1[[i]]<-rnorm(n,0,sigmasim1[i])
ytsim1[[i]]<-log(c2)+(vtsim1[[i]])
}
plot(c2)
points(exp(ytsim1[[20]]),col="green")
#jags
# cat("model{
#     for (i in 1:N){  
#     omega[i+1]~dgamma(aomega,bomega)
#     c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+(omega[i+1])
#     logyt[i+1]~dnorm(log(c[i+1]),sigmav[z2[i+1]])
#     z1[i+1]~dcat(p1[])
#     z2[i+1]~dcat(p2[])
#     }
#     p1[1]<-r1[1]
#     p2[1]<-r2[1]
#     for(j in 2:(m-1)){p1[j]<-r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1]
#     p2[j]<-r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1]}
#     for(l in 1:(m-1)){r1[l]~dbeta(1,alpha)
#     r2[l]~dbeta(1,alpha)}
#     ps1<-sum(p1[1:(m-1)])
#     p1[m]<-1-ps1
#     ps2<-sum(p2[1:(m-1)])
#     p2[m]<-1-ps2
#     for(l in 1:m){
#     muv[l]~dnorm(basemu, basetau)
#     sigmav[l]~dgamma(2,1)
#     }
#     basemu<-0
#     basetau<-10
#     omega[1]~dgamma(aomega,bomega)
#     c[1]<-1
#     logyt[1]~dnorm(log(c[1]), prec)
#     alpha~dunif(0,20)
#     aomega~dunif(0.5,3)
#     bomega~dunif(0.5,3)
#     sigma~dgamma(2,0.01)
#     prec<-1/sigma
#     Qprime~dunif(11,17)
#     Gprime~dunif(281,482)
#     kl~dunif(0,0.8)
#     }",file="modeljagsdp.jag")
# jags.dp<-jags.model("modeljagsdp.jag",data=list("logyt"=logyt, N=99,"h"=0.01, "V"=3.8,"m"=5))
# update(jags.dp, 1000)
# mcmcjags.dp<-jags.samples(jags.dp,
#                           c('Qprime','Gprime','kl','sigma','c','aomega','bomega','alpha'),
#                           1000)
# mcmcsamplesjags.dp<-coda.samples(jags.dp,
#                                  c('Qprime','Gprime','kl','sigma','c','aomega','bomega','alpha'),
#                                  1000)
# 
# m1.mcmc.dp<-(as.mcmc(mcmcsamplesjags.dp))
# m1.mat.dp<-as.matrix(mcmcsamplesjags.dp)
# m1.dat.dp<-as.data.frame(m1.mat.dp)
# aomega.post1zone.dp<-m1.dat.dp$aomega
# bomega.post1zone.dp<-m1.dat.dp$bomega
# qprime.post1zone.dp<-m1.dat.dp$Qprime
# gprime.post1zone.dp<-m1.dat.dp$Gprime
# kl.post1zone.dp<-m1.dat.dp$kl
# quantile(qprime.post1zone.dp,c(0.025,0.5,0.975))
# quantile(gprime.post1zone.dp,c(0.025,0.5,0.975))
# quantile(kl.post1zone.dp,c(0.025,0.5,0.975))
# quantile(aomega.post1zone.dp,c(0.025,0.5,0.975))
# quantile(bomega.post1zone.dp,c(0.025,0.5,0.975))
# c.hat.dp <- apply(mcmcjags.dp$c, 1, mean)
# c.hatlower.dp<- apply(mcmcjags.dp$c, 1, function(x) quantile(x,0.025))
# c.hatupper.dp<-apply(mcmcjags.dp$c, 1, function(x) quantile(x,0.975))
# plotCI(c, c.hat.dp,ui=c.hatupper.dp, li=c.hatlower.dp)
# abline(c(0,1))
# 
# quartz()
# plot(c,col="blue")
# points(c.hat.dp, col="red")
# points(exp(logyt),col="green")
# legend("bottomrigh", legend=c("True","Predicted","Observed"),col=c("blue","red","green"),
#        pch=15,cex=1)
# #interactive plot
# d1<-data.frame(cbind(t=seq(1,100,1),c,c.hat.dp,y=exp(logyt)))
# quartz()
# d1%>%ggvis(~t,~c)%>%
#   layer_points(fill:=input_select(c("red","blue","green")))%>%
#   layer_points(x=~t, y=~c.hat.dp,fill:=input_select(c("red","blue","green")))%>%
#   layer_points(x=~t, y=~y,fill:=input_select(c("red","blue","green")))%>%
#   add_axis("x",title="Time")%>%
#   add_axis("y",title="Concentrations")
# 
# #julia errors in Y are DP
# thetajuliadp1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmcdp1.csv")
# apply(thetajuliadp1, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
# plot(thetajuliadp1[,3],type="l")
# acf(thetajuliadp1[,2], lag.max=1000)
# xjuliadp1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmcdp1.csv")
# pred.juliadp1<-apply(xjuliadp1, 2, mean)
# quartz()
# plot(exp(logyt),col="green")
# points(c2,col="blue")
# points(pred.juliadp1,col="red")

#DP all
#jags
cat("model{
    for (i in 1:N){ 
    #omega[i+1]~dnorm(0,bomega[z2[i+1]])
    omega[i+1]~dnorm(aomega[z1[i+1]],bomega[z2[i+1]])T(-((1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V),)
    c[i+1]<-((1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V)+(omega[i+1])
    et[i+1]~dnorm(muv2[z4[i+1]],sigmav[z3[i+1]])T(-c[i+1],)
    yt[i+1]~dnorm((c[i+1])+muv2[z4[i+1]],sigmav[z3[i+1]])
    yhat[i+1]~dnorm((c[i+1])+muv2[z4[i+1]],sigmav[z3[i+1]])T(0,)
    #yhat[i+1]<-ifelse(yhatu[i+1]<0,0,yhatu[i+1])
    loglik[i+1]<-log(dnorm(yt[i+1],(c[i+1])+muv2[z4[i+1]],sigmav[z3[i+1]]))
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    z4[i+1]~dcat(p4[])
    }
    loglik[1]<-log(dnorm(yt[1],c[1],sigmav[z3[2]]))
    yhat[1]~dnorm((c[1])+muv2[z4[2]],sigmav[z3[2]])T(0,)
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    p4[1]<-r4[1]
    for(j in 2:(m-1)){p1[j]<-r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1]+0.0001
    p2[j]<-r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1]+0.0001
    p3[j]<-r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1]+0.0001
    p4[j]<-r4[j]*(1-r4[j-1])*p4[j-1]/r4[j-1]+0.0001}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    r4[l]~dbeta(1,alpha+0.0001)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    ps4<-sum(p4[1:(m-1)])
    p4[m]<-1-ps4
    for(l in 1:m){
    sigmav[l]~dgamma(0.01,0.01)
    muv[l]~dnorm(0,10)
    wmuv[l]<-muv[l]*p4[l]
    aomega[l]~dnorm(basemu,basetau)T(-1,)
    #waomega[l]<-aomega[l]*p1[l]
    bomega[l]~dgamma(0.01,0.01)
    #aomega[l]~dunif(0.01,1)
    #bomega[l]~dunif(0.01,1)
    }
    for(l in 1:m){
    muv2[l]<-muv[l]-sum(wmuv)
    #aomega2[l]<-aomega[l]-sum(waomega)
    #muv2[l]<-muv[l]-mean(muv)
    }
    basemu<-0
    basetau<-1
    omega[1]~dnorm(0,bomega[z2[2]])T(-yt[1],)
    c[1]<-yt[1]+omega[1]
    alpha~dgamma(0.01,0.01)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsdp2.jag")

jags.dp2<-jags.model("modeljagsdp2.jag",data=list("yt"=exp(logyt), N=99, 
                                                  "h"=0.01,"V"=V,"m"=5))
update(jags.dp2, 2000)
mcmcjags.dp2<-jags.samples(jags.dp2,
                           c('Qprime','Gprime','kl','sigmav','c','loglik','et','yhat'),
                           1000)
mcmcsamplesjags.dp2<-coda.samples(jags.dp2,
                                  c('Qprime','Gprime','kl','sigmav','c','loglik','et','yhat'),
                                  1000)
m1.mcmc.dp2<-(as.mcmc(mcmcsamplesjags.dp2))
m1.mat.dp2<-as.matrix(mcmcsamplesjags.dp2)
m1.dat.dp2<-as.data.frame(m1.mat.dp2)
qprime.post1zone.dp2<-m1.dat.dp2$Qprime
gprime.post1zone.dp2<-m1.dat.dp2$Gprime
kl.post1zone.dp2<-m1.dat.dp2$kl
plot((apply(mcmcjags.dp2$et, 1, mean)))
loglik.post1zone.dp2<-t(mcmcjags.dp2$loglik[,,1])
yhat.post1zone.dp2<-t(mcmcjags.dp2$yhat[,,1])
quantile(qprime.post1zone.dp2,c(0.025,0.5,0.975))
quantile(gprime.post1zone.dp2,c(0.025,0.5,0.975))
quantile(kl.post1zone.dp2,c(0.025,0.5,0.975))
waic(loglik.post1zone.dp2)

# Compute CRPS of simulated sample
mean(crps_sample(y = exp(logyt), dat = mcmcjags.dp2$yhat[,,1]))
c.hat.dp2 <- apply(mcmcjags.dp2$c, 1, mean)#+apply(mcmcjags.dp2$et, 1, mean)
c.hatlower50<- apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.25))
c.hatupper50<-apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.75))
yhat.dp2<-apply((mcmcjags.dp2$yhat[,,1]), 1, mean)
yhatlower50<-apply((mcmcjags.dp2$yhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50<-apply((mcmcjags.dp2$yhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90<-apply((mcmcjags.dp2$yhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90<-apply((mcmcjags.dp2$yhat[,,1]), 1, function(x) quantile(x,0.975))
c.hatlower90<- apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.05))
c.hatupper90<-apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.95))
mean(yhatupper90-yhatlower90)
mean(yhatupper50-yhatlower50)
cov50dp<-cov90dp<-NULL
for(i in 1:length(logyt)){
  cov50dp[i]<-(ifelse(yhatlower50[i]<=exp(logyt[i])&exp(logyt[i])<=yhatupper50[i],1,0))
  cov90dp[i]<-(ifelse(yhatlower90[i]<=exp(logyt[i])&exp(logyt[i])<=yhatupper90[i],1,0))
}
sum(cov50dp)/100
sum(cov90dp)/100
sum((c-c.hat.dp2)^2)/100
list<-apply(mcmcjags.dp2$yhat[,,1],2,function(x) empirical_cdf(x,exp(logyt)))
cdf<-NULL
for(i in 1:100){
  cdf[i]<-mean(sapply(list, function(X) X$CDF[i]))
}
hist(cdf)
#cdfdp1<-apply(mcmcjags.dp2$c[,,1],1,function(x) sum(ifelse(x<c,1,0))/1000)
#yhatdfdp1<-apply(mcmcjags.dp2$yhat[,,1],1,function(x) sum(ifelse(x<exp(logyt),1,0))/1000)
#ghatdp1<-ecdf(c)
ghatdp1<-empirical_cdf(exp(logyt),exp(logyt))
#gyhatdp1<-ecdf(exp(logyt))
diffy<-(cdf-ghatdp1$CDF)
quartz()
par(mfrow=c(2,1))
hist(cdf,breaks=15)
plot(diffy)
abline(h=0)
quartz()
plot(c.hat.dp2, type="l",ylim=c(0,max(yhatupper90)))
lines(yhatlower50,col="blue")
lines(yhatupper50, col="blue")
lines(yhatlower90,col="red")
lines(yhatupper90, col="red")
lines(exp(logyt), col="green")
lines(c, col="purple", lty=3)
plotCI(c, c.hat.dp2[1:100],ui=c.hatupperdp2[1:100], li=c.hatlowerdp2[1:100])
abline(0,1)
quartz()
plot((c2),col="blue")
points(c.hat.dp2,col="red")
points(exp(logyt), col="green")
legend("bottomright",legend=c("True","Predicted","Observed"),cex=1, pch=15,
       col=c("blue","red","green"))
#simulations
jags.dp2sim<-mcmcjags.dp2sim<-mcmcsamplesjags.dp2sim<-vector("list")
for(i in 1:100){
  jags.dp2sim[[i]]<-jags.model("modeljagsdp2.jag",data=list("yt"=exp(ytsim1[[i]]), N=99, 
                                                          "h"=0.01,"V"=V,"m"=5))
  update(jags.dp2sim[[i]], 2000)
  mcmcjags.dp2sim[[i]]<-jags.samples(jags.dp2sim[[i]],
                             c('Qprime','Gprime','kl','sigmav','c','loglik','et','yhat'),
                             1000)
  mcmcsamplesjags.dp2sim[[i]]<-coda.samples(jags.dp2sim[[i]],
                                    c('Qprime','Gprime','kl','sigmav','c','loglik','et','yhat'),
                                    1000)
}
m1.mcmc.dp2sim<-m1.mat.dp2sim<-m1.dat.dp2sim<-vector("list")
for(i in 1:100){
m1.mcmc.dp2sim[[i]]<-(as.mcmc(mcmcsamplesjags.dp2sim[[i]]))
m1.mat.dp2sim[[i]]<-as.matrix(mcmcsamplesjags.dp2sim[[i]])
m1.dat.dp2sim[[i]]<-as.data.frame(m1.mat.dp2sim[[i]])
}
crpssim1<-yhat90sim1<-yhat50sim1<-NULL
yhat.dpsim1<-yhatlower50sim1<-yhatupper50sim1<-yhatlower90sim1<-yhatupper90sim1<-
  chat.dpsim1<-c.hatlower50sim1<-c.hatupper50sim1<-c.hatlower90sim1<-c.hatupper90sim1<-vector("list")
for(i in 1:100){
  crpssim1[i]<-mean(crps_sample(y = exp(ytsim1[[i]]), dat = (mcmcjags.dp2sim[[i]])$yhat[,,1]))
  yhat.dpsim1[[i]]<-apply(((mcmcjags.dp2sim[[i]])$yhat[,,1]), 1, mean)
  yhatlower50sim1[[i]]<-apply(((mcmcjags.dp2sim[[i]])$yhat[,,1]), 1, function(x) quantile(x,0.25))
  yhatupper50sim1[[i]]<-apply(((mcmcjags.dp2sim[[i]])$yhat[,,1]), 1, function(x) quantile(x,0.75))
  yhatlower90sim1[[i]]<-apply(((mcmcjags.dp2sim[[i]])$yhat[,,1]), 1, function(x) quantile(x,0.05))
  yhatupper90sim1[[i]]<-apply(((mcmcjags.dp2sim[[i]])$yhat[,,1]), 1, function(x) quantile(x,0.975))
  chat.dpsim1[[i]]<-apply(((mcmcjags.dp2sim[[i]])$c), 1, mean)
  c.hatlower50sim1[[i]]<- apply(mcmcjags.dp2sim[[i]]$c, 1, function(x) quantile(x,0.25))
  c.hatupper50sim1[[i]]<-apply(mcmcjags.dp2sim[[i]]$c, 1, function(x) quantile(x,0.75))
  c.hatlower90sim1[[i]]<- apply(mcmcjags.dp2sim[[i]]$c, 1, function(x) quantile(x,0.05))
  c.hatupper90sim1[[i]]<-apply(mcmcjags.dp2sim[[i]]$c, 1, function(x) quantile(x,0.95))
  yhat90sim1[i]<-mean(yhatupper90sim1[[i]]-yhatlower90sim1[[i]])
  yhat50sim1[i]<-mean(yhatupper50sim1[[i]]-yhatlower50sim1[[i]])
}
cov50dpsim<-cov90dpsim<-msedpsim<-vector("list")
for(i in 1:100){
  cov50dpsim[[i]]<-(ifelse(c.hatlower50sim1[[i]]<=c2&c2<=c.hatupper50sim1[[i]],1,0))
  cov90dpsim[[i]]<-(ifelse(c.hatlower90sim1[[i]]<=c2&c2<=c.hatupper90sim1[[i]],1,0))
  msedpsim[[i]]<-sum((chat.dpsim1[[i]]-c2)^2)/100
  }
mean(unlist(lapply(cov50dpsim, function(x) sum(x)/100)))
mean(unlist(lapply(cov90dpsim, function(x) sum(x)/100)))
mean(unlist(msedpsim))
signal2noise<-vector("list")
meansignal2noise<-NULL
for(i in 1:100){
  signal2noise[[i]]<-c2/(sigmasim1[i])
  meansignal2noise[i]<-mean(signal2noise[[i]])
}
plot(meansignal2noise,yhat90sim1)
plot(meansignal2noise,crpssim1)
mean(crpssim1)
mean(yhat90sim1)
listsim<-vector("list")
for (i in 1:100){
  listsim[[i]]<-apply(mcmcjags.dp2sim[[i]]$yhat[,,1],2,function(x) empirical_cdf(x,exp(ytsim1[[i]])))
}
cdfsim<-matrix(0,100,100)
for(i in 1:100){
  for(j in 1:100){
  cdfsim[i,j]<-mean(sapply(listsim[[i]], function(X) X$CDF[j]))
}}
hist(cdfsim[1,])
ghatdp1sim<-diffysim<-vector("list")
for(i in 1:100){
  ghatdp1sim[[i]]<-empirical_cdf(exp(ytsim1[[i]]),exp(ytsim1[[i]]))
  diffysim[[i]]<-(cdfsim[i,]-ghatdp1sim[[i]]$CDF)
}

plot(unlist(diffysim[[1]]),type="l",col="red",lty=1,lwd=2 ,ylab="F.forecast-F.observed",ylim=c(-0.1,0.08) )
for( i in 2:100) { 
  lines(unlist(diffysim[[i]]),col=i,lty=1,lwd=2)
} 
abline(h=0)
#data one compartment
data1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/Modeling_Master_v5.csv")
dim(data1)
data1zone=matrix(0,nrow=204, ncol=(81*6))
for(i in 1:(81*6)){
  data1zone[,i]=data1[,4*i]
}
#true: V=11.9, Q=0.04-0.07, 0.23-0.27 and 0.47-0.77 m^3/min, 
#G=39.6, 79.1 and 118.7 for l, med, h
#note that each was measured at 6 locations
cat("model{
    for (i in 1:N){ 
    #omega[i+1]~dnorm(0,bomega[z2[i+1]])
    omega[i+1]~dnorm(aomega2[z1[i+1]],bomega[z2[i+1]])T(-((1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V),)
    c[i+1]<-((1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V)+(omega[i+1])
    et[i+1]~dnorm(muv2[z4[i+1]],sigmav[z3[i+1]])T(-c[i+1],)
    yt[i+1]~dnorm((c[i+1])+muv2[z4[i+1]],sigmav[z3[i+1]])
    yhat[i+1]~dnorm((c[i+1])+muv2[z4[i+1]],sigmav[z3[i+1]])T(0,)
    #yhat[i+1]<-ifelse(yhatu[i+1]<0,0,yhatu[i+1])
    loglik[i+1]<-log(dnorm(yt[i+1],(c[i+1])+muv2[z4[i+1]],sigmav[z3[i+1]]))
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    z4[i+1]~dcat(p4[])
    
    }
    loglik[1]<-log(dnorm(yt[1],(c[1]),sigmav[z3[2]]))
    yhat[1]<-(c[1])
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    p4[1]<-r4[1]
    for(j in 2:(m-1)){
    p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.0001
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.0001
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.0001
    p4[j]<-(r4[j]*(1-r4[j-1])*p4[j-1]/r4[j-1])}
    for(l in 1:(m-1)){
    r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    r4[l]~dbeta(1,alpha)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    ps4<-sum(p4[1:(m-1)])
    p4[m]<-1-ps4
    for(l in 1:m){
    muv[l]~dnorm(basemu, basetau)
    wmuv[l]<-muv[l]*p4[l]
    sigmav[l]~dgamma(0.01,0.01)
    aomega[l]~dnorm(basemu, basetau)
    womega[l]<-aomega[l]*p1[l]
    bomega[l]~dgamma(0.01,0.01)
    muv2[l]<-muv[l]-sum(wmuv)
    aomega2[l]<-aomega[l]-sum(womega)
    }
    basemu<-0
    basetau<-10
    omega[1]~dnorm(0,bomega[z2[2]])T(-yt[1],)
    c[1]<-(yt[1])+omega[1]
    #yt[1]~dnorm((c[1]),sigmav[z3[2]])
    alpha~dgamma(0.01,0.01)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,.001)
    }",file="modeljagsdatadp.jag")
jagslowdp<-jags.model("modeljagsdatadp.jag",data=list("yt"=(data1zone[1:201,2]), 
                                 N=200,"h"=0.01, "V"=11.9,"m"=10))
update(jagslowdp, 1000)
mcmcjagslowdp<-jags.samples(jagslowdp,
                          c('Qprime','Gprime',"kl",'c','loglik','et','yhat'),
                          2000)
mcmcsamplesjagslowdp<-coda.samples(jagslowdp,
                                 c('Qprime','Gprime','kl','c','loglik','et','yhat'),
                                 2000)
#low ventilation Q=0.04-0.07 G=43.18
m1.mcmclowdp<-(as.mcmc(mcmcsamplesjagslowdp[[1]]))
m1.matlowdp<-as.matrix(mcmcsamplesjagslowdp)
m1.datlowdp<-as.data.frame(m1.matlowdp)
aomega.post1zonelowdp<-m1.datlowdp$aomega
bomega.post1zonelowdp<-m1.datlowdp$bomega
qprime.post1zonelowdp<-m1.datlowdp$Qprime
gprime.post1zonelowdp<-m1.datlowdp$Gprime
kl.post1zonelowdp<-m1.datlowdp$kl
plot(apply(mcmcjagslowdp$et, 1, mean))
loglik.post1zoneldp<-t(mcmcjagslowdp$loglik[,,1])
waic(loglik.post1zoneldp)
mean(crps_sample(y = data1zone[1:201,2], dat = mcmcjagslowdp$yhat[,,1]))
quantile((gprime.post1zonelowdp),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonelowdp),c(0.025,0.25,0.5,0.75,0.975))
quantile((kl.post1zonelowdp),c(0.025,0.25,0.5,0.75,0.975))
c.hatlowdp <- apply(mcmcjagslowdp$c, 1, mean)
plot(c.hatlowdp)
points(data1zone[1:201,2],col="red")
sum((data1zone[1:201,2]-c.hatlowdp)^2)/201
c.hatlower50low<- apply(mcmcjagslowdp$c, 1, function(x) quantile(x,0.25))
c.hatupper50low<-apply(mcmcjagslowdp$c, 1, function(x) quantile(x,0.75))
c.hatlower90low<- apply(mcmcjagslowdp$c, 1, function(x) quantile(x,0.05))
c.hatupper90low<-apply(mcmcjagslowdp$c, 1, function(x) quantile(x,0.95))
yhat.dplow<-apply((mcmcjagslowdp$yhat[,,1]), 1, mean)
yhatlower50low<-apply((mcmcjagslowdp$yhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50low<-apply((mcmcjagslowdp$yhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90low<-apply((mcmcjagslowdp$yhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90low<-apply((mcmcjagslowdp$yhat[,,1]), 1, function(x) quantile(x,0.975))
mean(yhatupper50low-yhatlower50low)
mean(yhatupper90low-yhatlower90low)
cov50dpl<-cov90dpl<-NULL
for(i in 1:length(data1zone[1:201,2])){
  cov50dpl[i]<-(ifelse(yhatlower50low[i]<=data1zone[i,2]&data1zone[i,2]<=yhatupper50low[i],1,0))
  cov90dpl[i]<-(ifelse(yhatlower90low[i]<=data1zone[i,2]&data1zone[i,2]<=yhatupper90low[i],1,0))
}
sum(cov50dpl)/201
sum(cov90dpl)/201
cdfdp1low<-NULL
cdfdp1low<-apply(mcmcjagslowdp$c[,,1],1,function(x) sum(ifelse(x<data1zone[1:201,2],1,0))/(2000))
yhatdfdp1low<-NULL
yhatdfdp1low<-apply(mcmcjagslowdp$yhat[,,1],1,function(x) sum(ifelse(x<data1zone[1:201,2],1,0))/(2000))
quartz()
hist((yhatdfdp1low))
llist<-apply(mcmcjagslowdp$yhat[,,1],2,function(x) empirical_cdf(x,data1zone[1:201,2]))
lcdf<-NULL
for(i in 1:201){
  lcdf[i]<-mean(sapply(llist, function(X) X$CDF[i]))
}
hist(lcdf)
ghatdp1low<-ecdf(data1zone[1:201,2])
ghatdp1low<-empirical_cdf(data1zone[1:201,2], ubounds=data1zone[1:201,2])
difflow<-(lcdf-ghatdp1low$CDF)
quartz()
plot(difflow)
abline(h=0)
quartz()
plot(c.hatlowdp , type="l")
lines(yhatlower50low,col="blue")
lines(yhatupper50low, col="blue")
lines(yhatlower90low,col="red")
lines(yhatupper90low, col="red")
lines(data1zone[1:201,2], col="green")
#medium Q Q=0.23-0.27 G=43.18
jagsmeddp<-jags.model("modeljagsdatadp.jag",data=list("yt"=(data1zone[1:81,30]), N=80, "V"=11.9,"h"=0.01,"m"=10))
update(jagsmeddp, 1000)
mcmcjagsmeddp<-jags.samples(jagsmeddp,
                          c('Qprime','Gprime',"kl",'c','loglik','bomega','yhat'),
                          2000)
mcmcsamplesjagsmeddp<-coda.samples(jagsmeddp,
                                 c('Qprime','Gprime','kl','c','loglik','bomega','yhat'),
                                 2000)
m1.mcmcmeddp<-(as.mcmc(mcmcsamplesjagsmeddp))
m1.matmeddp<-as.matrix(mcmcsamplesjagsmeddp)
m1.datmeddp<-as.data.frame(m1.matmeddp)
gprime.post1zonemeddp<-m1.datmeddp$Gprime
qprime.post1zonemeddp<-m1.datmeddp$Qprime
kl.post1zonemeddp<-m1.datmeddp$kl
loglik.post1zonemdp<-t(mcmcjagsmeddp$loglik[,,1])
waic(loglik.post1zonemdp)
mean(crps_sample(y = data1zone[1:81,30], dat = mcmcjagsmeddp$yhat[,,1]))
quantile((gprime.post1zonemeddp),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonemeddp),c(0.025,0.25,0.5,0.75,0.975))
quantile((kl.post1zonemeddp),c(0.025,0.25,0.5,0.75,0.975))
c.hatmeddp <- apply(mcmcjagsmeddp$c, 1, mean)
plot(c.hatmeddp)
points(data1zone[,30],col="red")
sum((data1zone[1:81,30]-c.hatmeddp)^2)/81
c.hatlower50med<- apply(mcmcjagsmeddp$c, 1, function(x) quantile(x,0.25))
c.hatupper50med<-apply(mcmcjagsmeddp$c, 1, function(x) quantile(x,0.75))
c.hatlower90med<- apply(mcmcjagsmeddp$c, 1, function(x) quantile(x,0.05))
c.hatupper90med<-apply(mcmcjagsmeddp$c, 1, function(x) quantile(x,0.95))
yhat.dpm<-apply((mcmcjagsmeddp$yhat[,,1]), 1, mean)
yhatlower50m<-apply((mcmcjagsmeddp$yhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50m<-apply((mcmcjagsmeddp$yhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90m<-apply((mcmcjagsmeddp$yhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90m<-apply((mcmcjagsmeddp$yhat[,,1]), 1, function(x) quantile(x,0.975))
mean(yhatupper50m-yhatlower50m)
mean(yhatupper90m-yhatlower90m)
cov50dpm<-cov90dpm<-NULL
for(i in 1:length(data1zone[1:81,30])){
  cov50dpm[i]<-(ifelse(yhatlower50m[i]<=data1zone[i,30]&data1zone[i,30]<=yhatupper50m[i],1,0))
  cov90dpm[i]<-(ifelse(yhatlower90m[i]<=data1zone[i,30]&data1zone[i,30]<=yhatupper90m[i],1,0))
}
sum(cov50dpm)/81
sum(cov90dpm)/81
cdfdp1med<-NULL
cdfdp1med<-apply(mcmcjagsmeddp$yhat[,,1],1,function(x) sum(ifelse(x<data1zone[1:81,30],1,0))/(2000))
quartz()
hist((cdfdp1med))
mlist<-apply(mcmcjagsmeddp$yhat[,,1],2,function(x) empirical_cdf(x,data1zone[1:81,30]))
mcdf<-NULL
for(i in 1:81){
  mcdf[i]<-mean(sapply(mlist, function(X) X$CDF[i]))
}
hist(mcdf)
ghatdp1med<-ecdf(data1zone[1:81,30])
ghatdp1med<-empirical_cdf(data1zone[1:81,30], ubounds=data1zone[1:81,30])
plot(data1zone[1:81,30],ghatdp1med$CDF)
diffmed<-(mcdf-ghatdp1med$CDF)
quartz()
plot(diffmed)
abline(h=0)
quartz()
plot(c.hatmeddp , type="l")
lines(yhatlower50m,col="blue")
lines(yhatupper50m, col="blue")
lines(yhatlower90m,col="red")
lines(yhatupper90m, col="red")
lines(data1zone[1:81,30], col="green")
#high ventilation Q=0.47-0.77 low G=39.55
jagshdp<-jags.model("modeljagsdatadp.jag",data=list("yt"=(data1zone[1:41,51*4]), N=40, "V"=11.9,"h"=0.01,"m"=5))
update(jagshdp, 1000)
mcmcjagshdp<-jags.samples(jagshdp,
                        c('Qprime','Gprime',"kl",'c','loglik','bomega','yhat'),
                        2000)
mcmcsamplesjagshdp<-coda.samples(jagshdp,
                               c('Qprime','Gprime','kl','c','loglik','bomega','yhat'),
                               2000)
m1.mcmchdp<-(as.mcmc(mcmcsamplesjagshdp))
m1.mathdp<-as.matrix(mcmcsamplesjagshdp)
m1.dathdp<-as.data.frame(m1.mathdp)
gprime.post1zonehdp<-m1.dathdp$Gprime
qprime.post1zonehdp<-m1.dathdp$Qprime
loglik.post1zonehdp<-t(mcmcjagshdp$loglik[,,1])
waic(loglik.post1zonehdp)
mean(crps_sample(y = data1zone[1:41,51*4], dat = mcmcjagshdp$yhat[,,1]))
quantile((gprime.post1zonehdp),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonehdp),c(0.025,0.25,0.5,0.75,0.975))
c.hathdp <- apply(mcmcjagshdp$c, 1, mean)
sum((data1zone[1:41,51*4]-c.hathdp)^2)/41
c.hatlower50h<- apply(mcmcjagshdp$c, 1, function(x) quantile(x,0.25))
c.hatupper50h<-apply(mcmcjagshdp$c, 1, function(x) quantile(x,0.75))
c.hatlower90h<- apply(mcmcjagshdp$c, 1, function(x) quantile(x,0.05))
c.hatupper90h<-apply(mcmcjagshdp$c, 1, function(x) quantile(x,0.95))
yhat.dph<-apply((mcmcjagshdp$yhat[,,1]), 1, mean)
yhatlower50h<-apply((mcmcjagshdp$yhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50h<-apply((mcmcjagshdp$yhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90h<-apply((mcmcjagshdp$yhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90h<-apply((mcmcjagshdp$yhat[,,1]), 1, function(x) quantile(x,0.975))
mean(yhatupper50h-yhatlower50h)
mean(yhatupper90h-yhatlower90h)
cov50dph<-cov90dph<-NULL
for(i in 1:length(data1zone[1:41,51*4])){
  cov50dph[i]<-(ifelse(yhatlower50h[i]<=data1zone[i,51*4]&data1zone[i,51*4]<=yhatupper50h[i],1,0))
  cov90dph[i]<-(ifelse(yhatlower90h[i]<=data1zone[i,51*4]&data1zone[i,51*4]<=yhatupper90h[i],1,0))
}
sum(cov50dph)/41
sum(cov90dph)/41
cdfdp1h<-NULL
cdfdp1h<-apply(mcmcjagshdp$yhat[,,1],1,function(x) sum(ifelse(x<data1zone[1:41,51*4],1,0))/(2000))
hlist<-apply(mcmcjagshdp$yhat[,,1],2,function(x) empirical_cdf(x, data1zone[1:41,51*4]))
hcdf<-NULL
for(i in 1:41){
hcdf[i]<-mean(sapply(hlist, function(X) X$CDF[i]))
}
quartz()
hist((hcdf))
#ghatdp1h<-ecdf(data1zone[1:41,51*4])
ghatdp1h<-empirical_cdf(data1zone[1:41,51*4], ubounds = data1zone[1:41,51*4])
hist(ghatdp1h$CDF)
diffh<-((hcdf)-ghatdp1h$CDF)
quartz()
plot(diffh)
abline(h=0)
quartz()
par(mfrow=c(1,3))
plot(c.hatlowdp , type="l", main="low")
lines(yhatlower50low,col="blue")
lines(yhatupper50low, col="blue")
lines(yhatlower90low,col="red")
lines(yhatupper90low, col="red")
lines(data1zone[1:201,2], col="green")
plot(c.hatmeddp , type="l", main="med")
lines(yhatlower50m,col="blue")
lines(yhatupper50m, col="blue")
lines(yhatlower90m,col="red")
lines(yhatupper90m, col="red")
lines(data1zone[1:81,30], col="green")
plot(c.hathdp , type="l", main="high")
lines(yhatlower50h,col="blue")
lines(yhatupper50h, col="blue")
lines(yhatlower90h,col="red")
lines(yhatupper90h, col="red")
lines(data1zone[1:41,51*4], col="green")

#all PIT
quartz()
par(mfrow=c(3,1))
hist((lcdf), main="low", breaks=15)
hist((mcdf), main="medium", breaks=15)
hist((hcdf),main="high", breaks=15)
#all diff plots
quartz()
par(mfrow=c(3,1))
plot((difflow), main="low")
abline(h=0)
plot((diffmed), main="medium")
abline(h=0)
plot((diffh),main="high")
abline(h=0)

quartz()
par(mfrow=c(1,3))
plot(c.hatlowdp, main="low ventilation")
points(data1zone[1:201,2],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat low","measured"))
plot(c.hatmeddp, main="medium ventilation")
points(data1zone[1:81,30],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat medium","measured"))
plot(c.hathdp,main="High ventilation")
points(data1zone[,51*4],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat high","measured"))

#two compartment model simulation
n=100
Qprime=13.8
Gprime=351.54
VF=3.8
VN=0.03141593
Sigma<-diag(0.01,2)
kl=0.1
beta=5
V=diag(0.1,2)
r=3
set.seed(2017)
wt<-rlnorm.rplus(n,c(0,0),diag(2))
#wt<-rmvnorm(n,c(0,0),V)
h=0.01
A<-matrix(c(-(beta)/VN, (beta)/VN, beta/VF, -(beta+Qprime)/VF+kl),nrow=2, ncol=2,
          byrow=TRUE)
eigen(A)
g<-matrix(c(Gprime/VN,0),nrow=1, ncol=2)

yt<-NULL
c<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
for(i in 1:(n-1)){
  c[1,]<-c(0,0.5)
  c[i+1,]<-(c[i,])%*%t(h*A+diag(1,2))+(g*h)
}
c
plot(c[,1],ylab="mg/m3", xlab="time")
points(c[,2])

#exact 
c2<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
#imp

for(i in 1:(n-1)){
  c2[1,]<-c(0,0.5)
  c2[i+1,]<- expm(h*i*A)%*%c2[1,]+solve(A)%*%(expm(h*i*A)-diag(2))%*%t(g)
}
plot(c2[,1])
points(c2[,2])
plot(c[,1],c2[,1])
abline(0,1)
plot(c[,2],c2[,2])
abline(0,1)

c<-c2

set.seed(123)
vt<-rmvnorm(n,c(0,0),Sigma)
logyt2<-log((c+0.05))+(vt)
exp(logyt2)

#simulations
var2sim<-vtsim<-logyt2sim<-vector("list")
for(i in 1:100){
  var2sim[[i]]<-Sigma*exp(rnorm(1,0,1))
  vtsim[[i]]<-rmvnorm(n,c(0,0),var2sim[[i]])
  logyt2sim[[i]]<-log((c+0.05))+(vtsim[[i]])  
}
#winbugs
library(R2WinBUGS)
cat("model{
    muc[1] <- 0
    muc[2] <- 0
    for (i in 1:N){
    lim1[1] <- (-((c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h)))
    lim1[2] <- (-(c[i,1]*cc*h+c[i,2]*((d*h)+1)))
    lim2[1] <- 10000
    lim2[2] <- 10000
    omega[i+1,1:2] ~ dmnorm(muc[],vari[1:2,1:2,z3[i+1]])I(lim1[],lim2[])
    c[i+1,1] <- (c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h)+omega[i+1,1]
    c[i+1,2] <- (c[i,1]*cc*h+c[i,2]*((d*h)+1))+omega[i+1,2]
    yt[i+1,1:2] ~ dmnorm((c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
    yhat[i+1,1:2] ~ dmnorm((c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
    z1[i+1] ~ dcat(p1[])
    z2[i+1] ~ dcat(p2[])
    z3[i+1] ~ dcat(p3[])
    z4[i+1] ~ dcat(p4[])
    }
    yhat[1,1:2] ~ dmnorm((c[1,1:2])+c(aomega1w[z1[2]],aomega2w[z1[2]]), prec[1:2,1:2,z4[2]])
    p1[1] <- r1[1]
    p2[1] <- r2[1]
    p3[1] <- r3[1]
    p4[1] <- r4[1]
    for(j in 2:(m-1)){
    p1[j] <- (r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.01
    p2[j] <- (r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.01
    p3[j] <- (r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.01
    p4[j] <- (r4[j]*(1-r4[j-1])*p4[j-1]/r4[j-1])+0.01}
    for(l in 1:(m-1)){
    r1[l] ~ dbeta(1,alpha+0.01)
    r2[l] ~ dbeta(1,alpha+0.01)
    r3[l] ~ dbeta(1,alpha+0.01)
    r4[l] ~ dbeta(1,alpha+0.01)
    }
    ps1 <- sum(p1[1:(m-1)])
    p1[m] <- 1-ps1
    ps2 <- sum(p2[1:(m-1)])
    p2[m] <- 1-ps2
    ps3 <- sum(p3[1:(m-1)])
    p3[m] < -1-ps3
    ps4 <- sum(p4[1:(m-1)])
    p4[m] <- 1-ps4
    for(l in 1:m){
    aomega1[l] ~ dnorm(basemu, basetau)
    aomega1w[l] <- aomega1[l]-sum(aomega1[l]*p1[l])
    aomega2[l] ~ dnorm(basemu, basetau)
    aomega2w[l] <- aomega1[l]-sum(aomega2[l]*p1[l])
    prec[1:2,1:2,l] ~ dwish(I,3)
    vari[1:2,1:2,l] ~ dwish(I,3)
    }
    basemu <- 0
    basetau <- 10000
    alpha ~ dgamma(0.01,0.01)
    a <- (-(beta)/VN)
    b <- (beta)/VN
    cc <- beta/VF
    d <- (-(beta+Qprime)/VF+kl)
    g1 <- Gprime/VN
    c[1,1] <- 0
    c[1,2] < -0.5
    yt[1,1:2] ~ dmnorm((c[1,1:2]),I2)
    Qprime ~ dunif(11,17)
    Gprime ~ dunif(281,482)
    kl ~ dunif(0,0.5)
    beta ~ dunif(0,10)
    }",file="twozone.txt")
data=list("yt"=exp(logyt2), N=99,"I"=diag(1,2), "VN"=VN, "VF"=VF,
          "h"=0.01,"m"=10,"I2"=diag(100,2))
parm<-list('Qprime','Gprime','beta',"kl",'c','yhat')
inits=list(Qprime=13,Gprime=350,beta=5,kl=0.1,c=exp(logyt2),yhat=exp(logyt2))
bugs(data, inits, parm, model.file = "twozone.txt", n.chains =6, n.iter =2000, 
     n.burnin =500,bugs.directory="C:/Users/Nada/Downloads/winbugs14 (1)/WinBUGS14",debug=TRUE)
#jags
cat("model{
    for (i in 1:N){
    # mean[1:2,i]=(dmnorm(lim1[1:2,i],c(0,0),vari[1:2,1:2,z3[i+1]]))/
    # (1-exp(logdensity.mnorm(lim1[1:2,i],c(0,0),vari[1:2,1:2,z4[i+1]])))
    # variance[1:2,1:2,i]=vari[1:2,1:2,z3[i+1]]*(1+(lim1[1:2,i]*dmnorm(lim1[1:2,i],c(0,0),vari[1:2,1:2,z3[i+1]]))/
    # (1-exp(logdensity.mnorm(lim1[1:2,i],c(0,0),vari[1:2,1:2,z4[i+1]])))-((dmnorm(lim1[1:2,i],c(0,0),vari[1:2,1:2,z3[i+1]]))/
    # (1-exp(logdensity.mnorm(lim1[1:2,i],c(0,0),vari[1:2,1:2,z4[i+1]]))))^2)
    lim1[1:2,i] <- c(-((c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h)),-(c[i,1]*cc*h+c[i,2]*((d*h)+1)))
    omega[i+1,1:2]~dmnorm(c(0,0),vari[1:2,1:2,z3[i+1]])
    omegac[i+1,1]<-ifelse(omega[i+1,1]<lim1[1,i],lim1[1,i],omega[i+1,1])
    omegac[i+1,2]<-ifelse(omega[i+1,2]<lim1[2,i],lim1[2,i],omega[i+1,2])
     #omega[i+1,1]~dnorm(0,(1/sigma1[z3[i+1]]))T(lim1[1,i],)
     #omega[i+1,2]~dnorm((sigma12[z3[i+1]])*((1/sigma1[z3[i+1]]))*omega[i+1,1],(1/(sigma2[z3[i+1]]-
     #sigma12[z3[i+1]]^2*(1/sigma1[z3[i+1]]))))T(lim1[2,i],)
    c[i+1,1]<-(c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h)+omegac[i+1,1]
    c[i+1,2]<-(c[i,1]*cc*h+c[i,2]*((d*h)+1))+omegac[i+1,2]
    yt[i+1,1:2]~dmnorm((c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
    yhatc[i+1,1:2]~dmnorm((c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
yhat[i+1,1]<-ifelse(yhatc[i+1,1]<0,yt[i+1,1],yhatc[i+1,1])
yhat[i+1,2]<-ifelse(yhatc[i+1,2]<0,yt[i+1,2],yhatc[i+1,2])
    #yhatu[i+1,1:2]~dmnorm(c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
    #yhat[i+1,1]~dnorm(c[i+1,1]+yhatu[i+1,1],10000000)T(0,)
    #yhat[i+1,2]~dnorm(c[i+1,2]+yhatu[i+1,2],10000000)T(0,)
    #loglik[i+1]<-(logdensity.mnorm(yt[i+1,1:2],(c[i+1,1:2]),prec[1:2,1:2,z4[i+1]]))
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    z4[i+1]~dcat(p4[])
    }
    yhat[1,1:2]~dmnorm((c[1,1:2])+c(aomega1w[z1[2]],aomega2w[z1[2]]), prec[1:2,1:2,z4[2]])
    #loglik[1]<-(logdensity.mnorm(yt[1,1:2],(c[1,1:2]),prec[1:2,1:2,z4[2]]))
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    p4[1]<-r4[1]
    for(j in 2:(m-1)){p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.01
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.01
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.01
    p4[j]<-(r4[j]*(1-r4[j-1])*p4[j-1]/r4[j-1])+0.01}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.01)
    r2[l]~dbeta(1,alpha+0.01)
    r3[l]~dbeta(1,alpha+0.01)
    r4[l]~dbeta(1,alpha+0.01)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    ps4<-sum(p4[1:(m-1)])
    p4[m]<-1-ps4
    for(l in 1:m){
    aomega1[l]~dnorm(basemu, basetau)
    aomega1w[l]<-aomega1[l]-sum(aomega1[l]*p1[l])
    aomega2[l]~dnorm(basemu, basetau)
    aomega2w[l]<-aomega1[l]-sum(aomega2[l]*p1[l])
    prec[1:2,1:2,l]~dwish(I,3)
    vari[1:2,1:2,l]~dwish(I,3)
    #varit[1:2,1:2,l]<-inverse(vari[1:2,1:2,l])
    # sigma1[l]<-varit[1,1,l]
    # sigma2[l]<-varit[2,2,l]
    # sigma12[l]<-varit[1,2,l]
    }
    basemu<-0
    basetau<-10000
    alpha~dgamma(0.01,0.01)
    a<- -(beta)/VN
    b<-(beta)/VN
    cc<-beta/VF
    d<- -(beta+Qprime)/VF+kl
    g1<-Gprime/VN
    c[1,1]<-0
    c[1,2]<-0.5
    yt[1,1:2]~dmnorm((c[1,1:2]),I2)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.5)
    beta~dunif(0,10)
    }",file="modeljags2zonedp.jag")
jagsdp2<-jags.model("modeljags2zonedp.jag",data=list("yt"=exp(logyt2), N=99,"I"=diag(1,2), "VN"=VN, "VF"=VF,
                                                  "h"=0.01,"m"=6,"I2"=diag(100,2)))
update(jagsdp2, 1000)
mcmcjags.2zonedp2<-jags.samples(jagsdp2,
                             c('Qprime','Gprime','beta',"kl",'c','loglik','yhat'),
                             1000)
mcmcsamplesjags.2zonedp2<-coda.samples(jagsdp2,
                                    c('Qprime','Gprime','beta',"kl",'c','loglik','yhat'),
                                    1000)

m1.mcmc.2zonedp<-(as.mcmc(mcmcsamplesjags.2zonedp2))
m1.mat.2zonedp<-as.matrix(mcmcsamplesjags.2zonedp2)
m1.dat.2zonedp<-as.data.frame(m1.mat.2zonedp)
qprime.post.2zonedp<-m1.dat.2zonedp$Qprime
gprime.post.2zonedp<-m1.dat.2zonedp$Gprime
beta.post.2zonedp<-m1.dat.2zonedp$beta
kl.post.2zonedp<-m1.dat.2zonedp$kl
#loglik.post2zonedp<-t(mcmcjags.2zonedp2$loglik[,,1])
#waic(loglik.post2zonedp)
# yhat1<-yhat2<-vector("list")
# for(i in 1:100){
# yhat1[[i]]<-subset(mcmcjags.2zonedp2$yhat[i,1,,1],mcmcjags.2zonedp2$yhat[i,1,,1]>0)
# yhat2[[i]]<-subset(mcmcjags.2zonedp2$yhat[i,2,,1],mcmcjags.2zonedp2$yhat[i,2,,1]>0)
# }
mean(crps_sample(y = exp(logyt2[,1]), dat = mcmcjags.2zonedp2$yhat[,1,,1]))+
  mean(crps_sample(y = exp(logyt2[,2]), dat = mcmcjags.2zonedp2$yhat[,2,,1]))
quantile(beta.post.2zonedp,c(0.025,0.25,0.5,0.75,0.975))
quantile(qprime.post.2zonedp,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonedp,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zone <- apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, mean)
c2.hat.2zone <- apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, mean)
yhat1.2zone <- apply(mcmcjags.2zonedp2$yhat[,1,,1], 1, mean)
yhat2.2zone <- apply(mcmcjags.2zonedp2$yhat[,2,,1], 1, mean)

c1.hatlower.2zone<- apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, function(x) quantile(x,0.025))
c1.hatupper.2zone<-apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, function(x) quantile(x,0.975))
c2.hatlower.2zone<- apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, function(x) quantile(x,0.025))
c2.hatupper.2zone<-apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, function(x) quantile(x,0.975))
plotCI(c[,1], c1.hat.2zone,li=c1.hatlower.2zone, ui=c1.hatupper.2zone)
abline(c(0,1))
plotCI(c[,2], c2.hat.2zone,li=c2.hatlower.2zone, ui=c2.hatupper.2zone)
abline(c(0,1))
c.hatlower50n<- apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, function(x) quantile(x,0.25))
c.hatupper50n<-apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, function(x) quantile(x,0.75))
yhatlower50n<- apply(mcmcjags.2zonedp2$yhat[,1,,1], 1, function(x) quantile(x,0.25))
yhatupper50n<-apply(mcmcjags.2zonedp2$yhat[,1,,1], 1, function(x) quantile(x,0.75))
mean(yhatupper50n-yhatlower50n)
c.hatlower50f<- apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, function(x) quantile(x,0.25))
c.hatupper50f<-apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, function(x) quantile(x,0.75))
yhatlower50f<- apply(mcmcjags.2zonedp2$yhat[,2,,1], 1, function(x) quantile(x,0.25))
yhatupper50f<-apply(mcmcjags.2zonedp2$yhat[,2,,1], 1, function(x) quantile(x,0.75))
mean(yhatupper50f-yhatlower50f)
c.hatlower90n<- apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, function(x) quantile(x,0.05))
c.hatupper90n<-apply(mcmcsamplesjags.2zonedp2[[1]][,4:103], 2, function(x) quantile(x,0.95))
yhatlower90n<- apply(mcmcjags.2zonedp2$yhat[,1,,1], 1, function(x) quantile(x,0.05))
yhatupper90n<-apply(mcmcjags.2zonedp2$yhat[,1,,1], 1, function(x) quantile(x,0.95))
mean(yhatupper90n-yhatlower90n)
c.hatlower90f<- apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, function(x) quantile(x,0.05))
c.hatupper90f<-apply(mcmcsamplesjags.2zonedp2[[1]][,104:203], 2, function(x) quantile(x,0.95))
yhatlower90f<- apply(mcmcjags.2zonedp2$yhat[,2,,1], 1, function(x) quantile(x,0.05))
yhatupper90f<-apply(mcmcjags.2zonedp2$yhat[,2,,1], 1, function(x) quantile(x,0.95))
mean(yhatupper90f-yhatlower90f)
cov502ndp<-cov502fdp<-cov902ndp<-cov902fdp<-NULL
for(i in 1:length(logyt2[,1])){
  cov502ndp[i]<-(ifelse(yhatlower50n[i]<=exp(logyt2[i,1])&exp(logyt2[i,1])<=yhatupper50n[i],1,0))
  cov502fdp[i]<-(ifelse(yhatlower50f[i]<=exp(logyt2[i,2])&exp(logyt2[i,2])<=yhatupper50f[i],1,0))
  cov902ndp[i]<-(ifelse(yhatlower90n[i]<=exp(logyt2[i,1])&exp(logyt2[i,1])<=yhatupper90n[i],1,0))
  cov902fdp[i]<-  (ifelse(yhatlower90f[i]<=exp(logyt2[i,2])&exp(logyt2[i,2])<=yhatupper90f[i],1,0))
}
sum(cov502ndp)/200+sum(cov502fdp)/200
sum(cov902ndp)/200+sum(cov902fdp)/200
cdfdp2<-NULL
#cdfdp2n<-apply(mcmcjags.2zonedp2$c[,1,,1],1,function(x) sum(ifelse(x<=c[,1],1,0))/1000)
#cdfdp2f<-apply(mcmcjags.2zonedp2$c[,2,,1],1,function(x) sum(ifelse(x<=c[,2],1,0))/1000)
cdfdp2n<-apply(mcmcjags.2zonedp2$yhat[,1,,1],1,function(x) empirical_cdf(x,exp(logyt2[,1])))
cdfdp2f<-apply(mcmcjags.2zonedp2$yhat[,2,,1],1,function(x) empirical_cdf(x,exp(logyt2[,2])))
cdfn<-cdff<-NULL
for(i in 1:100){
  cdfn[i]<-mean(sapply(cdfdp2n, function(X) X$CDF[i]))
  cdff[i]<-mean(sapply(cdfdp2f, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,1))
hist((cdfn),breaks = 15, main="Near Field")
hist((cdff),breaks = 15, main="Far Field")
ghatdp2n<-empirical_cdf(exp(logyt2[,1]),exp(logyt2[,1]))
ghatdp2f<-empirical_cdf(exp(logyt2[,2]),exp(logyt2[,2]))
diff2n<-((cdfn)-ghatdp2n$CDF)
diff2f<-((cdff)-ghatdp2f$CDF)
quartz()
par(mfrow=c(2,1))
plot(diff2n,main="Near Field")
abline(h=0)
plot(diff2f, main="Far Field")
abline(h=0)
quartz()
par(mfrow=c(2,1))
plot(c1.hat.2zone, type="l",ylim=c(0, max(yhatupper90n)))
lines(yhatlower50n,col="blue")
lines(yhatupper50n, col="blue")
lines(yhatlower90n,col="red")
lines(yhatupper90n, col="red")
lines(exp(logyt2[,1]), col="green")
lines(c[,1], col="purple", lty=3)
plot(c2.hat.2zone, type="l",ylim=c(0, max(yhatupper90f)))
lines(yhatlower50f,col="blue")
lines(yhatupper50f, col="blue")
lines(yhatlower90f,col="red")
lines(yhatupper90f, col="red")
lines(exp(logyt2[,2]), col="green")
lines(c[,2], col="purple", lty=3)
quartz()
par(mfrow=c(2,1))
plot(c[,1],col="red")
points(c1.hat.2zone, col="blue")
points(exp(logyt2[,1]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtN","chatN","ytN"),bg="white", lty=1,cex=0.8)
plot(c[,2],col="red")
points(c2.hat.2zone, col="blue")
points(exp(logyt2[,2]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtF","chatF","ytF"),bg="white", lty=1,cex=0.8)
#sim
jagsdp2sim<-mcmcjags.2zonedp2sim<-mcmcsamplesjags.2zonedp2sim<-vector("list")
for(i in 1:100){
  jagsdp2sim[[i]]<-jags.model("modeljags2zonedp.jag",data=list("yt"=exp(logyt2sim[[i]]), N=99,"I"=diag(1,2), "VN"=VN, "VF"=VF,
                                                       "h"=0.01,"m"=10,"I2"=diag(100,2)))
  update(jagsdp2sim[[i]], 1000)
  mcmcjags.2zonedp2sim[[i]]<-jags.samples(jagsdp2sim[[i]],
                                  c('Qprime','Gprime','beta',"kl",'c','yhat'),
                                  1000)
  mcmcsamplesjags.2zonedp2sim[[i]]<-coda.samples(jagsdp2sim[[i]],
                                         c('Qprime','Gprime','beta',"kl",'c','yhat'),
                                         1000)
}
m1.mcmc.2zonedpsim<-m1.mat.2zonedpsim<-m1.dat.2zonedpsim<-yhatlower50nsim<-
  yhatupper50nsim<-yhatlower50fsim<-yhatupper50fsim<-yhatlower90nsim<-
  yhatupper90nsim<-yhatlower90fsim<-yhatupper90fsim<-yhatnsim<-yhatfsim<-chatlower50nsim<-
  chatupper50nsim<-chatlower50fsim<-chatupper50fsim<-chatlower90nsim<-
  chatupper90nsim<-chatlower90fsim<-chatupper90fsim<-chatnsim<-chatfsim<-vector("list")
for(i in 1:100){
  m1.mcmc.2zonedpsim[[i]]<-(as.mcmc(mcmcsamplesjags.2zonedp2sim[[i]]))
  m1.mat.2zonedpsim[[i]]<-as.matrix(mcmcsamplesjags.2zonedp2sim[[i]])
  m1.dat.2zonedpsim[[i]]<-as.data.frame(m1.mat.2zonedpsim[[i]])
  yhatnsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$yhat[,1,,1], 1, mean)
  yhatfsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$yhat[,2,,1], 1,mean)
  yhatlower50nsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$yhat[,1,,1], 1, function(x) quantile(x,0.25))
  yhatupper50nsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$yhat[,1,,1], 1, function(x) quantile(x,0.75))
yhatlower50fsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$yhat[,2,,1], 1, function(x) quantile(x,0.25))
yhatupper50fsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$yhat[,2,,1], 1, function(x) quantile(x,0.75))
yhatlower90nsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$yhat[,1,,1], 1, function(x) quantile(x,0.05))
yhatupper90nsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$yhat[,1,,1], 1, function(x) quantile(x,0.95))
yhatlower90fsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$yhat[,2,,1], 1, function(x) quantile(x,0.05))
yhatupper90fsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$yhat[,2,,1], 1, function(x) quantile(x,0.95))
chatnsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$c[,1,,1], 1, mean)
chatfsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$c[,2,,1], 1,mean)
chatlower50nsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.25))
chatupper50nsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.75))
chatlower50fsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.25))
chatupper50fsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.75))
chatlower90nsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.05))
chatupper90nsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.95))
chatlower90fsim[[i]]<- apply((mcmcjags.2zonedp2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.05))
chatupper90fsim[[i]]<-apply((mcmcjags.2zonedp2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.95))
}
cov502ndpsim<-cov502fdpsim<-cov902ndpsim<-cov902fdpsim<-mse2dpnsim<-mse2dpfsim<-crps2sim<-vector("list")
for(i in 1:100){
  cov502ndpsim[[i]]<-(ifelse(chatlower50nsim[[i]]<=(c[,1])&(c[,1])<=chatupper50nsim[[i]],1,0))
  cov502fdpsim[[i]]<-(ifelse(chatlower50fsim[[i]]<=(c[,2])&(c[,2])<=chatupper50fsim[[i]],1,0))
  cov902ndpsim[[i]]<-(ifelse(chatlower90nsim[[i]]<=(c[,1])&(c[,1])<=chatupper90nsim[[i]],1,0))
  cov902fdpsim[[i]]<- (ifelse(chatlower90fsim[[i]]<=(c[,2])&(c[,2])<=chatupper90fsim[[i]],1,0))
  mse2dpnsim[[i]]<-sum((chatnsim[[i]]-(c[,1]))^2)/100
  mse2dpfsim[[i]]<-sum((chatfsim[[i]]-(c[,2]))^2)/100
  crps2sim[[i]]<-mean(crps_sample(y = exp(logyt2sim[[i]])[,1], dat = (mcmcjags.2zonedp2sim[[i]])$yhat[,1,,1]))+
    mean(crps_sample(y =exp(logyt2sim[[i]])[,2], dat = (mcmcjags.2zonedp2sim[[i]])$yhat[,2,,1]))
  }
0.5*(mean(unlist(lapply(cov502ndpsim, function(x) sum(x)/100)))+mean(unlist(lapply(cov502fdpsim, function(x) sum(x)/100))))
0.5*(mean(unlist(lapply(cov902ndpsim, function(x) sum(x)/100)))+mean(unlist(lapply(cov902fdpsim, function(x) sum(x)/100))))
mean(unlist(mse2dpnsim))+mean(unlist(mse2dpfsim))
mean(unlist(crps2sim))
qprime.post.2zonedpsim<-gprime.post.2zonedpsim<-beta.post.2zonedpsim<-kl.post.2zonedpsim<-vector("list")
for(i in 1:100){
  qprime.post.2zonedpsim[[i]]<-m1.dat.2zonedpsim[[i]]$Qprime
  gprime.post.2zonedpsim[[i]]<-m1.dat.2zonedpsim[[i]]$Gprime
  beta.post.2zonedpsim[[i]]<-m1.dat.2zonedpsim[[i]]$beta
  kl.post.2zonedpsim[[i]]<-m1.dat.2zonedpsim[[i]]$kl
}
quantile(unlist(qprime.post.2zonedpsim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(gprime.post.2zonedpsim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(beta.post.2zonedpsim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(kl.post.2zonedpsim),c(0.025,0.25,0.5,0.75,0.975))
cdfdp2nsim<-cdfdp2fsim<-vector("list")
for(i in 1:100){
  cdfdp2nsim[[i]]<-apply(mcmcjags.2zonedp2sim[[i]]$yhat[,1,,1],1,function(x) empirical_cdf(x,exp(logyt2sim[[i]][,1])))
  cdfdp2fsim[[i]]<-apply(mcmcjags.2zonedp2sim[[i]]$yhat[,2,,1],1,function(x) empirical_cdf(x,exp(logyt2sim[[i]][,2])))
}

cdfnsim<-cdffsim<-matrix(0,100,100)
for(i in 1:100){
  for(j in 1:100){
    cdfnsim[i,j]<-mean(sapply(cdfdp2nsim[[i]], function(X) X$CDF[j]))
    cdffsim[i,j]<-mean(sapply(cdfdp2fsim[[i]], function(X) X$CDF[j]))
  }
}
hist((cdfnsim[10,]),breaks = 15, main="Near Field")
hist((cdffsim[20,]),breaks = 15, main="Far Field")
ghatdp2nsim<-ghatdp2fsim<-diff2nsim<-diff2fsim<-vector("list")
for(i in 1:100){
  ghatdp2nsim[[i]]<-empirical_cdf(exp(logyt2sim[[i]][,1]),exp(logyt2sim[[i]][,1]))
  ghatdp2fsim[[i]]<-empirical_cdf(exp(logyt2sim[[i]][,2]),exp(logyt2sim[[i]][,2]))
  diff2nsim[[i]]<-((cdfnsim[i,])-ghatdp2nsim[[i]]$CDF)
  diff2fsim[[i]]<-((cdffsim[i,])-ghatdp2fsim[[i]]$CDF)
}
par(mfrow=c(1,2))
plot(unlist(diff2nsim[[1]]),main="Near Field",type="l", lwd=2, ylim=c(-0.13,0.13),ylab=
       "F.forecast-F.observed")
for(i in 2:100){
  lines(unlist(diff2nsim[[i]]),col=i,lwd=2)
}
abline(h=0)
plot(unlist(diff2fsim[[1]]),main="Far Field",type="l", lwd=2, ylim=c(-0.13,0.13),ylab=
       "F.forecast-F.observed")
for(i in 2:100){
  lines(unlist(diff2fsim[[i]]),col=i,lwd=2)
}
abline(h=0)
#data
data2zone=matrix(0,nrow=204, ncol=(81*4*6))
data2zonen=matrix(0,nrow=204, ncol=81)
data2zonef1=matrix(0,nrow=204, ncol=81)
data2zonef2=matrix(0,nrow=204, ncol=81)
data2zonef3=matrix(0,nrow=204, ncol=81)
dim(data1)
for(i in 1:(81*6)){
  data2zone[,i]=data1[,(81*4*6)+(4*i)]
}
for(i in 1:81){
  data2zonen[,i]=data2zone[,1+(6*(i-1))]
  data2zonef1[,i]=data2zone[,4+(6*(i-1))]
  data2zonef2[,i]=data2zone[,5+(6*(i-1))]
  data2zonef3[,i]=data2zone[,6+(6*(i-1))]
}

#true values Q=0.059, 0.258 and 0.595 Vn=0.105 Vf=11.79 beta=0.24-1.24
cat("model{
    for (i in 1:N){
    omega[i+1,1:2]~dmnorm(c(0,0),vari[1:2,1:2,z3[i+1]])
    c[i+1,1]<-(c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h)+(omega[i+1,1])
    c[i+1,2]<-(c[i,1]*cc*h+c[i,2]*((d*h)+1))+(omega[i+1,2])
    yt[i+1,1:2]~dmnorm((c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
    yhat[i+1,1:2]~dmnorm((c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]), prec[1:2,1:2,z4[i+1]])
    loglik[i+1]<-(logdensity.mnorm(yt[i+1,1:2],(c[i+1,1:2])+c(aomega1w[z1[i+1]],aomega2w[z1[i+1]]),prec[1:2,1:2,z4[i+1]]))
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    z4[i+1]~dcat(p4[])
    }
    yhat[1,1:2]~dmnorm((c[1,1:2])+c(aomega1w[z1[2]],aomega2w[z1[2]]), prec[1:2,1:2,z4[2]])
    loglik[1]<-(logdensity.mnorm(yt[1,1:2],(c[1,1:2]),prec[1:2,1:2,z4[2]]))
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    p4[1]<-r4[1]
    for(j in 2:(m-1)){p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.0001
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.0001
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.0001
    p4[j]<-(r4[j]*(1-r4[j-1])*p4[j-1]/r4[j-1])+0.0001}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    r4[l]~dbeta(1,alpha+0.0001)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    ps4<-sum(p4[1:(m-1)])
    p4[m]<-1-ps4
    for(l in 1:m){
    aomega1[l]~dnorm(basemu, basetau)
    aomega2[l]~dnorm(basemu, basetau)
    aomega1w[l]<-aomega1[l]-sum(aomega1[l]*p1[l])
    aomega2w[l]<-aomega2[l]-sum(aomega2[l]*p1[l])
    vari[1:2,1:2,l]~dwish(I,3)
    prec[1:2,1:2,l]~dwish(I,3)
    }
    basemu<-0
    basetau<-100
    alpha~dgamma(0.01,0.01)
    a<- -(beta)/VN
    b<-(beta)/VN
    cc<-beta/VF
    d<- -(beta+Qprime)/VF+kl
    g1<-Gprime/VN
    omega[1,1:2]~dmnorm(c(0,0),100*I)
    c[1,1:2]<-yt[1,1:2]+omega[1,1:2]
    #yt[1,1:2]~dmnorm((c[1,1:2]), I)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,0.001)
    beta~dunif(0,5)
    }",file="modeljags2zonedpd.jag")
#med ventilation Q=0.23-0.27 med G=86.36
jagsmdp<-jags.model("modeljags2zonedpd.jag",data=list("yt"=cbind((data2zonen[1:75,13]),(data2zonef1[1:75,13])),
                                                 N=74,"I"=diag(1,2), "VN"=0.1, "VF"=11.9,"h"=0.01,"m"=10))
update(jagsmdp, 1000)
mcmcjags.2zonemdp<-jags.samples(jagsmdp,
                              c('Qprime','Gprime','beta',"kl",'c',"prec",'loglik','yhat'),
                              1000)
mcmcsamplesjags.2zonemdp<-coda.samples(jagsmdp,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec",'loglik','yhat'),
                                     1000)

m1.mcmc.2zonemdp<-(as.mcmc(mcmcsamplesjags.2zonemdp))
m1.mat.2zonemdp<-as.matrix(mcmcsamplesjags.2zonemdp)
m1.dat.2zonemdp<-as.data.frame(m1.mat.2zonemdp)
qprime.post.2zonemdp<-m1.dat.2zonemdp$Qprime
gprime.post.2zonemdp<-m1.dat.2zonemdp$Gprime
beta.post.2zonemdp<-m1.dat.2zonemdp$beta
loglik.post2zonemdp<-t(mcmcjags.2zonemdp$loglik[,,1])
waic(loglik.post2zonemdp)
mean(crps_sample(y = data2zonen[1:75,13], dat = mcmcjags.2zonemdp$yhat[,1,,1]))+
  mean(crps_sample(y = data2zonef1[1:75,13], dat =mcmcjags.2zonemdp$yhat[,2,,1]))
kl.post.2zonemdp<-m1.dat.2zonemdp$kl
quantile(qprime.post.2zonemdp,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonemdp,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonemdp,c(0.025,0.25,0.5,0.75,0.975))
quantile(kl.post.2zonemdp,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonemdp <- apply(mcmcsamplesjags.2zonemdp[[1]][,4:(78)], 2, mean)
c2.hat.2zonemdp <- apply(mcmcsamplesjags.2zonemdp[[1]][,(79):(79+74)], 2, mean)
plot(c2.hat.2zonemdp,data2zonef1[1:75,13])
abline(0,1)
(sum((data2zonen[1:75,13]-(c1.hat.2zonemdp))^2)+
    sum((data2zonef1[1:75,13]-(c2.hat.2zonemdp))^2))/150
c.hatlowermed50n<- apply(mcmcsamplesjags.2zonemdp[[1]][,4:(78)], 2, function(x) quantile(x,0.25))
c.hatuppermed50n<-apply(mcmcsamplesjags.2zonemdp[[1]][,4:(78)], 2, function(x) quantile(x,0.75))
yhatlowermed50n<- apply(mcmcsamplesjags.2zonemdp[[1]][,270:344], 2, function(x) quantile(x,0.25))
yhatuppermed50n<-apply(mcmcsamplesjags.2zonemdp[[1]][,270:344], 2, function(x) quantile(x,0.75))
mean(yhatuppermed50n-yhatlowermed50n)
c.hatlowermed50f<- apply(mcmcsamplesjags.2zonemdp[[1]][,(79):(79+74)], 2, function(x) quantile(x,0.25))
c.hatuppermed50f<-apply(mcmcsamplesjags.2zonemdp[[1]][,(79):(79+74)], 2, function(x) quantile(x,0.75))
yhatlowermed50f<- apply(mcmcsamplesjags.2zonemdp[[1]][,(345):(419)], 2, function(x) quantile(x,0.25))
yhatuppermed50f<-apply(mcmcsamplesjags.2zonemdp[[1]][,(345):(419)], 2, function(x) quantile(x,0.75))
mean(yhatuppermed50f-yhatlowermed50f)
c.hatlowermed90n<- apply(mcmcsamplesjags.2zonemdp[[1]][,4:(78)], 2, function(x) quantile(x,0.05))
c.hatuppermed90n<-apply(mcmcsamplesjags.2zonemdp[[1]][,4:(78)], 2, function(x) quantile(x,0.95))
yhatlowermed90n<- apply(mcmcsamplesjags.2zonemdp[[1]][,270:344], 2, function(x) quantile(x,0.05))
yhatuppermed90n<-apply(mcmcsamplesjags.2zonemdp[[1]][,270:344], 2, function(x) quantile(x,0.95))
mean(yhatuppermed90n-yhatlowermed90n)
c.hatlowermed90f<- apply(mcmcsamplesjags.2zonemdp[[1]][,(79):(79+74)], 2, function(x) quantile(x,0.05))
c.hatuppermed90f<-apply(mcmcsamplesjags.2zonemdp[[1]][,(79):(79+74)], 2, function(x) quantile(x,0.95))
yhatlowermed90f<- apply(mcmcsamplesjags.2zonemdp[[1]][,(345):(419)], 2, function(x) quantile(x,0.05))
yhatuppermed90f<-apply(mcmcsamplesjags.2zonemdp[[1]][,(345):(419)], 2, function(x) quantile(x,0.95))
mean(yhatuppermed90f-yhatlowermed90f)
cov502nmdp<-cov502fmdp<-cov902nmdp<-cov902fmdp<-NULL
for(i in 1:length(data2zonen[1:75,13])){
  cov502nmdp[i]<-(ifelse(yhatlowermed50n[i]<=data2zonen[i,13]&data2zonen[i,13]<=yhatuppermed50n[i],1,0))
  cov502fmdp[i]<-(ifelse(yhatlowermed50f[i]<=data2zonef1[i,13]&data2zonef1[i,13]<=yhatuppermed50f[i],1,0))
  cov902nmdp[i]<-(ifelse(yhatlowermed90n[i]<=data2zonen[i,13]&data2zonen[i,13]<=yhatuppermed90n[i],1,0))
  cov902fmdp[i]<-  (ifelse(yhatlowermed90f[i]<=data2zonef1[i,13]&data2zonef1[i,13]<=yhatuppermed90f[i],1,0))
}
sum(cov502nmdp)/150+sum(cov502fmdp)/150
sum(cov902nmdp)/150+sum(cov902fmdp)/150
cdfdp2medn<-apply(mcmcjags.2zonemdp$c[,1,,1],1,function(x) empirical_cdf(x,data2zonen[1:75,13]))
cdfdp2medf<-apply(mcmcjags.2zonemdp$c[,2,,1],1,function(x) empirical_cdf(x,data2zonef1[1:75,13]))
cdfnm<-cdffm<-NULL
for(i in 1:75){
  cdfnm[i]<-mean(sapply(cdfdp2medn, function(X) X$CDF[i]))
  cdffm[i]<-mean(sapply(cdfdp2medf, function(X) X$CDF[i]))
}
hist((cdfnm),breaks = 15, main="Near Field")
hist((cdffm),breaks = 15, main="Far Field")
ghatdp2medn<-empirical_cdf(data2zonen[1:75,13],data2zonen[1:75,13])
ghatdp2medf<-empirical_cdf(data2zonef1[1:75,13],data2zonef1[1:75,13])
diff2medn<-((cdfnm)-ghatdp2medn$CDF)
diff2medf<-((cdffm)-ghatdp2medf$CDF)
plot(diff2medn,main="Near Field",type="l")
abline(h=0)
plot(diff2medf, main="Far Field",type="l")
abline(h=0)
plot(c1.hat.2zonemdp, type="l")
lines(yhatlowermed50n,col="blue")
lines(yhatuppermed50n, col="blue")
lines(yhatlowermed90n,col="red")
lines(yhatuppermed90n, col="red")
lines(data2zonen[1:75,13], col="green")
plot(c2.hat.2zonemdp, type="l")
lines(yhatlowermed50f,col="blue")
lines(yhatuppermed50f, col="blue")
lines(yhatlowermed90f,col="red")
lines(yhatuppermed90f, col="red")
lines(data2zonef1[1:75,13], col="green")
#high
jagshdp<-jags.model("modeljags2zonedpd.jag",data=list("yt"=cbind((data2zonen[1:41,81]),(data2zonef1[1:41,81])),
                                                 N=40,"I"=diag(1,2), "VN"=0.1, "VF"=11.9,"h"=0.01,"m"=5))
update(jagshdp, 1000)
mcmcjags.2zonehdp<-jags.samples(jagshdp,
                              c('Qprime','Gprime','beta',"kl",'c',"prec",'loglik','yhat'),
                              1000)
mcmcsamplesjags.2zonehdp<-coda.samples(jagshdp,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec",'loglik','yhat'),
                                     1000)
m1.mcmc.2zonehdp<-(as.mcmc(mcmcsamplesjags.2zonehdp))
m1.mat.2zonehdp<-as.matrix(mcmcsamplesjags.2zonehdp)
m1.dat.2zonehdp<-as.data.frame(m1.mat.2zonehdp)
qprime.post.2zonehdp<-m1.dat.2zonehdp$Qprime
gprime.post.2zonehdp<-m1.dat.2zonehdp$Gprime
beta.post.2zonehdp<-m1.dat.2zonehdp$beta
kl.post.2zonehdp<-m1.dat.2zonehdp$kl
loglik.post2zonehdp<-t(mcmcjags.2zonehdp$loglik[,,1])
waic(loglik.post2zonehdp)
mean(crps_sample(y = data2zonen[1:41,81], dat = mcmcjags.2zonehdp$yhat[,1,,1]))+
  mean(crps_sample(y = data2zonef1[1:41,81], dat =mcmcjags.2zonehdp$yhat[,2,,1]))
quantile(qprime.post.2zonehdp,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonehdp,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonehdp,c(0.025,0.25,0.5,0.75,0.975))
quantile(alpha.post.2zonehdp,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonehdp <- apply(mcmcsamplesjags.2zonehdp[[1]][,4:(44)], 2, mean)
c2.hat.2zonehdp <- apply(mcmcsamplesjags.2zonehdp[[1]][,(45):85], 2, mean)
plot(c2.hat.2zonehdp,data2zonef1[1:41,81])
abline(0,1)
(sum((data2zonen[1:41,81]-(c1.hat.2zonehdp))^2)+
    sum((data2zonef1[1:41,81]-(c2.hat.2zonehdp))^2))/82
c.hatlowerh50n<- apply(mcmcsamplesjags.2zonehdp[[1]][,4:(44)], 2, function(x) quantile(x,0.25))
c.hatupperh50n<-apply(mcmcsamplesjags.2zonehdp[[1]][,4:(44)], 2, function(x) quantile(x,0.75))
yhatlowerh50n<- apply(mcmcsamplesjags.2zonehdp[[1]][,148:188], 2, function(x) quantile(x,0.25))
yhatupperh50n<-apply(mcmcsamplesjags.2zonehdp[[1]][,148:188], 2, function(x) quantile(x,0.75))
mean(yhatupperh50n-yhatlowerh50n)
c.hatlowerh50f<- apply(mcmcsamplesjags.2zonehdp[[1]][,(45):(85)], 2, function(x) quantile(x,0.25))
c.hatupperh50f<-apply(mcmcsamplesjags.2zonehdp[[1]][,(45):(85)], 2, function(x) quantile(x,0.75))
yhatlowerh50f<- apply(mcmcsamplesjags.2zonehdp[[1]][,189:229], 2, function(x) quantile(x,0.25))
yhatupperh50f<-apply(mcmcsamplesjags.2zonehdp[[1]][,189:229], 2, function(x) quantile(x,0.75))
mean(yhatupperh50f-yhatlowerh50f)
c.hatlowerh90n<- apply(mcmcsamplesjags.2zonehdp[[1]][,4:(44)], 2, function(x) quantile(x,0.05))
c.hatupperh90n<-apply(mcmcsamplesjags.2zonehdp[[1]][,4:(44)], 2, function(x) quantile(x,0.95))
yhatlowerh90n<- apply(mcmcsamplesjags.2zonehdp[[1]][,148:188], 2, function(x) quantile(x,0.05))
yhatupperh90n<-apply(mcmcsamplesjags.2zonehdp[[1]][,148:188], 2, function(x) quantile(x,0.95))
mean(yhatupperh90n-yhatlowerh90n)
c.hatlowerh90f<- apply(mcmcsamplesjags.2zonehdp[[1]][,(45):(85)], 2, function(x) quantile(x,0.05))
c.hatupperh90f<-apply(mcmcsamplesjags.2zonehdp[[1]][,(45):(85)], 2, function(x) quantile(x,0.95))
yhatlowerh90f<- apply(mcmcsamplesjags.2zonehdp[[1]][,189:229], 2, function(x) quantile(x,0.05))
yhatupperh90f<-apply(mcmcsamplesjags.2zonehdp[[1]][,189:229], 2, function(x) quantile(x,0.95))
mean(yhatupperh90f-yhatlowerh90f)
cov502nhdp<-cov502fhdp<-cov902nhdp<-cov902fhdp<-NULL
for(i in 1:length(data2zonen[1:41,81])){
  cov502nhdp[i]<-(ifelse(yhatlowerh50n[i]<=data2zonen[i,81]&data2zonen[i,81]<=yhatupperh50n[i],1,0))
  cov502fhdp[i]<-(ifelse(yhatlowerh50f[i]<=data2zonef1[i,81]&data2zonef1[i,81]<=yhatupperh50f[i],1,0))
  cov902nhdp[i]<-(ifelse(yhatlowerh90n[i]<=data2zonen[i,81]&data2zonen[i,81]<=yhatupperh90n[i],1,0))
  cov902fhdp[i]<-  (ifelse(yhatlowerh90f[i]<=data2zonef1[i,81]&data2zonef1[i,81]<=yhatupperh90f[i],1,0))
}
sum(cov502nhdp)/82+sum(cov502fhdp)/82
sum(cov902nhdp)/82+sum(cov902fhdp)/82
cdfdp2hn<-apply(mcmcjags.2zonehdp$yhat[,1,,1],1,function(x) empirical_cdf(x,data2zonen[1:41,81]))
cdfdp2hf<-apply(mcmcjags.2zonehdp$yhat[,2,,1],1,function(x) empirical_cdf(x,data2zonef1[1:41,81]))
cdfnh<-cdffh<-NULL
for(i in 1:41){
  cdfnh[i]<-mean(sapply(cdfdp2hn, function(X) X$CDF[i]))
  cdffh[i]<-mean(sapply(cdfdp2hf, function(X) X$CDF[i]))
}
hist((cdfnh),breaks = 15, main="Near Field")
hist((cdffh),breaks = 15, main="Far Field")
ghatdp2hn<-empirical_cdf(data2zonen[1:41,81],data2zonen[1:41,81])
ghatdp2hf<-empirical_cdf(data2zonef1[1:41,81],data2zonef1[1:41,81])
diff2hn<-((cdfnh)-ghatdp2hn$CDF)
diff2hf<-((cdffh)-ghatdp2hf$CDF)
plot(diff2hn,main="Near Field")
abline(h=0)
plot(diff2hf, main="Far Field")
abline(h=0)
quartz()
par(mfrow=c(2,1))
plot(c1.hat.2zonehdp, type="l")
lines(yhatlowerh50n,col="blue")
lines(yhatupperh50n, col="blue")
lines(yhatlowerh90n,col="red")
lines(yhatupperh90n, col="red")
lines(data2zonen[1:41,81], col="green")
plot(c2.hat.2zonehdp, type="l")
lines(yhatlowerh50f,col="blue")
lines(yhatupperh50f, col="blue")
lines(yhatlowerh90f,col="red")
lines(yhatupperh90f, col="red")
lines(data2zonef1[1:41,81], col="green")
#low Q=0.056
jagsldp<-jags.model("modeljags2zonedpd.jag",data=list("yt"=cbind((data2zonen[1:176,1]),(data2zonef1[1:176,1])),
                                                 N=175,"I"=diag(1,2), "VN"=0.1, "VF"=11.9,"h"=0.01,"m"=10))
update(jagsldp, 1000)
mcmcjags.2zoneldp<-jags.samples(jagsldp,
                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','yhat'),
                              2000)
mcmcsamplesjags.2zoneldp<-coda.samples(jagsldp,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec",'loglik','yhat'),
                                     2000)
m1.mcmc.2zoneldp<-(as.mcmc(mcmcsamplesjags.2zoneldp))
m1.mat.2zoneldp<-as.matrix(mcmcsamplesjags.2zoneldp)
m1.dat.2zoneldp<-as.data.frame(m1.mat.2zoneldp)
qprime.post.2zoneldp<-m1.dat.2zoneldp$Qprime
gprime.post.2zoneldp<-m1.dat.2zoneldp$Gprime
beta.post.2zoneldp<-m1.dat.2zoneldp$beta
kl.post.2zoneldp<-m1.dat.2zoneldp$kl
loglik.post2zoneldp<-t(mcmcjags.2zoneldp$loglik[,,1])
waic(loglik.post2zoneldp)
mean(crps_sample(y = data2zonen[1:176,1], dat = mcmcjags.2zoneldp$yhat[,1,,1]))+
  mean(crps_sample(y = data2zonef1[1:176,1], dat =mcmcjags.2zoneldp$yhat[,2,,1]))
quantile(qprime.post.2zoneldp,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zoneldp,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zoneldp,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zoneldp <- apply(mcmcsamplesjags.2zoneldp[[1]][,4:(103+76)], 2, mean)
c2.hat.2zoneldp <- apply(mcmcsamplesjags.2zoneldp[[1]][,(102+78):355], 2, mean)
quartz()
plot(c2.hat.2zoneldp,data2zonef1[1:176,1])
abline(0,1)
(sum((data2zonen[1:176,1]-(c1.hat.2zoneldp))^2)+
    sum((data2zonef1[1:176,1]-(c2.hat.2zoneldp))^2))/(176*2)
c.hatlowerl50n<- apply(mcmcsamplesjags.2zoneldp[[1]][,4:(103+76)], 2, function(x) quantile(x,0.25))
c.hatupperl50n<-apply(mcmcsamplesjags.2zoneldp[[1]][,4:(103+76)], 2, function(x) quantile(x,0.75))
yhatlowerl50n<- apply(mcmcjags.2zoneldp$yhat[,1,,1], 1, function(x) quantile(x,0.25))
yhatupperl50n<-apply(mcmcjags.2zoneldp$yhat[,1,,1], 1, function(x) quantile(x,0.75))
mean(yhatupperl50n-yhatlowerl50n)
c.hatlowerl50f<- apply(mcmcsamplesjags.2zoneldp[[1]][,(102+78):355], 2, function(x) quantile(x,0.25))
c.hatupperl50f<-apply(mcmcsamplesjags.2zoneldp[[1]][,(102+78):355], 2, function(x) quantile(x,0.75))
yhatlowerl50f<- apply(mcmcjags.2zoneldp$yhat[,2,,1], 1, function(x) quantile(x,0.25))
yhatupperl50f<-apply(mcmcjags.2zoneldp$yhat[,2,,1], 1, function(x) quantile(x,0.75))
mean(yhatupperl50f-yhatlowerl50f)
c.hatlowerl90n<- apply(mcmcsamplesjags.2zoneldp[[1]][,4:(103+76)], 2, function(x) quantile(x,0.05))
c.hatupperl90n<-apply(mcmcsamplesjags.2zoneldp[[1]][,4:(103+76)], 2, function(x) quantile(x,0.95))
yhatlowerl90n<- apply(mcmcjags.2zoneldp$yhat[,1,,1], 1, function(x) quantile(x,0.05))
yhatupperl90n<-apply(mcmcjags.2zoneldp$yhat[,1,,1], 1, function(x) quantile(x,0.95))
mean(yhatupperl90n-yhatlowerl90n)
c.hatlowerl90f<- apply(mcmcsamplesjags.2zoneldp[[1]][,(102+78):355], 2, function(x) quantile(x,0.05))
c.hatupperl90f<-apply(mcmcsamplesjags.2zoneldp[[1]][,(102+78):355], 2, function(x) quantile(x,0.95))
yhatlowerl90f<- apply(mcmcjags.2zoneldp$yhat[,2,,1], 1, function(x) quantile(x,0.05))
yhatupperl90f<-apply(mcmcjags.2zoneldp$yhat[,2,,1], 1, function(x) quantile(x,0.95))
mean(yhatupperl90f-yhatlowerl90f)
cov502nldp<-cov502fldp<-cov902nldp<-cov902fldp<-NULL
for(i in 1:length(data2zonen[1:176,1])){
  cov502nldp[i]<-(ifelse(yhatlowerl50n[i]<=data2zonen[i,1]&data2zonen[i,1]<=yhatupperl50n[i],1,0))
  cov502fldp[i]<-(ifelse(yhatlowerl50f[i]<=data2zonef1[i,1]&data2zonef1[i,1]<=yhatupperl50f[i],1,0))
  cov902nldp[i]<-(ifelse(yhatlowerl90n[i]<=data2zonen[i,1]&data2zonen[i,1]<=yhatupperl90n[i],1,0))
  cov902fldp[i]<-  (ifelse(yhatlowerl90f[i]<=data2zonef1[i,1]&data2zonef1[i,1]<=yhatupperl90f[i],1,0))
}
sum(cov502nldp)/(176*2)+sum(cov502fldp)/(176*2)
sum(cov902nldp)/(176*2)+sum(cov902fldp)/(176*2)
cdfdp2ln<-apply(mcmcjags.2zoneldp$yhat[,1,,1],1,function(x) empirical_cdf(x,data2zonen[1:176,1]))
cdfdp2lf<-apply(mcmcjags.2zoneldp$yhat[,2,,1],1,function(x) empirical_cdf(x,data2zonef1[1:176,1]))
cdfnl<-cdffl<-NULL
for(i in 1:176){
  cdfnl[i]<-mean(sapply(cdfdp2ln, function(X) X$CDF[i]))
  cdffl[i]<-mean(sapply(cdfdp2lf, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,3))
hist(cdfnl,breaks = 15, main="Near Field")
hist(cdfnm,breaks = 15, main="Near Field")
hist(cdfnh,breaks = 15, main="Near Field")
hist(cdffl,breaks = 15, main="Far Field")
hist(cdffm,breaks = 15, main="Far Field")
hist(cdffh,breaks = 15, main="Far Field")
ghatdp2ln<-empirical_cdf(data2zonen[1:176,1],data2zonen[1:176,1])
ghatdp2lf<-empirical_cdf(data2zonef1[1:176,1],data2zonef1[1:176,1])
diff2ln<-(cdfnl-ghatdp2ln$CDF)
diff2lf<-(cdffl-ghatdp2ln$CDF)
quartz()
par(mfrow=c(2,3))
plot(diff2ln,main="Near Field")
abline(h=0)
plot(diff2hn,main="Near Field")
abline(h=0)
plot(diff2medn,main="Near Field")
abline(h=0)
plot(diff2lf, main="Far Field")
abline(h=0)
plot(diff2hf, main="Far Field")
abline(h=0)
plot(diff2medf, main="Far Field")
abline(h=0)
quartz()
par(mfrow=c(2,3))
plot(c1.hat.2zoneldp, type="l",main="Near Field")
lines(yhatlowerl50n,col="blue")
lines(yhatupperl50n, col="blue")
lines(yhatlowerl90n,col="red")
lines(yhatupperl90n, col="red")
lines(data2zonen[1:176,1], col="green")
plot(c1.hat.2zonemdp, type="l",main="Near Field")
lines(yhatlowermed50n,col="blue")
lines(yhatuppermed50n, col="blue")
lines(yhatlowermed90n,col="red")
lines(yhatuppermed90n, col="red")
lines(data2zonen[1:75,13], col="green")
plot(c1.hat.2zonehdp, type="l",main="Near Field")
lines(yhatlowerh50n,col="blue")
lines(yhatupperh50n, col="blue")
lines(yhatlowerh90n,col="red")
lines(yhatupperh90n, col="red")
lines(data2zonen[1:41,81], col="green")
plot(c2.hat.2zoneldp, type="l", main="Far Field")
lines(yhatlowerl50f,col="blue")
lines(yhatupperl50f, col="blue")
lines(yhatlowerl90f,col="red")
lines(yhatupperl90f, col="red")
lines(data2zonef1[1:176,1], col="green")
plot(c2.hat.2zonemdp, type="l", main="Far Field")
lines(yhatlowermed50f,col="blue")
lines(yhatuppermed50f, col="blue")
lines(yhatlowermed90f,col="red")
lines(yhatuppermed90f, col="red")
lines(data2zonef1[1:75,13], col="green")
plot(c2.hat.2zonehdp, type="l", main="Far Field")
lines(yhatlowerh50f,col="blue")
lines(yhatupperh50f, col="blue")
lines(yhatlowerh90f,col="red")
lines(yhatupperh90f, col="red")
lines(data2zonef1[1:41,81], col="green")

quartz()
par(mfrow=c(3,2))
plot(c1.hat.2zoneldp, col="blue",main="Low")
points(data2zonen[1:176,1],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zoneldp, col="blue")
points(data2zonef1[3:176,1],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

plot(c1.hat.2zonem, col="blue",main="Medium")
points(data2zonen[1:75,13],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonem, col="blue")
points(data2zonef1[1:75,13],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

plot(c1.hat.2zonehdp, col="blue",main="High")
points(data2zonen[1:41,81],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonehdp, col="blue")
points(data2zonef1[1:41,81],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

#Turbulent eddy-diffusion model
library(spacetime)
library(gstat)
library(sp)
Dt=1
nloc<-10
ntime<-100
h=1
Gprime=351.4
x<-seq(1,nloc,length.out=ntime)
y<-seq(1,nloc,length.out=ntime)
matdist=NULL
#matdist <- as.matrix(dist(c(x,y),upper=TRUE, diag=TRUE))
for (i in 1:nloc){
  matdist[i]=sqrt(y[i]^2+x[i]^2)
}
#matdist=c(sqrt(y[1]^2+x[100]^2),matdist[1:10],sqrt(y[100]^2+x[1]^2),sqrt(y[100]^2+x[100]^2))
erf<-function(x){
  (2/sqrt(pi))*exp(-x^2)
}
csp<-matrix(0,nloc,ntime)
for(i in 1:ntime){
  for(j in 1:nloc){
    csp[j,i]<-Gprime/(2*pi*Dt*(matdist[j]))*(1-
                                               (integrate(erf,lower=0,upper=(matdist[j])/sqrt(4*Dt*h*i)))$value)
  } 
}
Gprime/(2*pi*Dt*matdist[3])
#quartz()
plot(csp[1,])
points(csp[2,])
points(csp[3,])
points(csp[4,])
points(csp[5,])
csp2<-matrix(0,nloc,ntime)
csp2[,1]<-csp[,1]
for(i in 2:ntime){
  for(j in 1:nloc){
    csp2[,i]<-csp2[,(i-1)]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
      (exp(-(matdist[j])^2/(4*Dt*i)))
  }
}
csp2
plot(csp[1,], csp2[1,])
abline(c(0,1))

#quartz()
plot(csp2[1,], ylim=c(0,max(csp2[1,])))
points(csp2[2,])
points(csp2[3,])
points(csp2[4,])
points(csp2[5,])


#sim from random error
yt3<-matrix(0,nrow=nloc, ncol=ntime)
for(i in 1:ntime){
  for(j in 1:nloc){
    yt3[j,i]<-exp(rnorm(1,log(csp[j,i]),0.1))
  }
}
plot(csp[1,])
points((yt3[1,]),col="red")

plot(csp[2,], ylim=c(min(csp[2,]),max((yt3[2,]))))
points((yt3[2,]),col="red")

plot(csp[3,],ylim=c(min(csp[3,]),max((yt3[3,]))))
points((yt3[3,]),col="red")

plot(csp[4,])
points((yt3[4,]),col="red")

plot(csp[5,])
points((yt3[5,]),col="red")
#sim from additive time space cove
phi<-1
sigma<-0.1
tausq<-0.25

#library(geoR)
Sigma.exp <- function(x,y,phi,sigma) {
  Sigma <- matrix(rep(0, length(x)*length(y)), nrow=length(x))
  Sigma<- as.matrix(sigma*exp(-(dist(data.frame(cbind(x,y)),upper=TRUE, diag=TRUE,method="euclidean")*phi)))
  return(Sigma)
}
sigma.sp<-Sigma.exp(x,y,phi,sigma)
dim(sigma.sp)
set.seed(1111)
et2<-rmvnorm(ntime,rep(0,nloc),((sigma.sp[1:nloc,1:nloc])), method=( "svd"))
wt2<-matrix(0,ntime,nloc)
wt2[1,1:nloc]<-rep(0,nloc)
for(i in 2:ntime){
  wt2[i,]<-wt2[i-1,]+et2[i,]
}
wt2
yt32<-matrix(0,nrow=nloc, ncol=ntime)
for(i in 1:ntime){
  for(j in 1:nloc){
    yt32[j,i]<-exp(rnorm(1,log(csp[j,i]+et2[i,j]),0.1))
  }
}
plot(csp[1,])
points((yt32[1,]),col="red")
plot(csp[2,])
points((yt32[2,]),col="red")
plot(csp[1,])
points((yt32[3,]),col="red")
plot(csp[3,])
points((yt32[4,]),col="red")
plot(csp[5,])
points((yt32[5,]),col="red")

#sim from ns cov
nloc=5
k<-matrix(0,nloc*ntime,nloc*ntime)
a<-b<-1
for(j in 1:nloc){
  for(l in 1:nloc){
    for(r in 1:(ntime)){
      for(s in 1:(ntime)){
        k[(l*r)+((l-1)*((ntime-1)-(r-1))),(j*s)+((j-1)*((ntime-1)-(s-1)))]<-
          (sigma/(a*((r-s))^2+1))*exp(-b*dist[l,j]/(a*((r-s))^2+1))
 }}}}

et3<-rmvnorm(1,rep(0,nloc*ntime),k)
yt33<-matrix(0,nrow=nloc, ncol=ntime)
for(i in 1:ntime){
  for(j in 1:nloc){
    yt33[j,i]<-exp(rnorm(1,log(csp[j,i]+et3[i*j]),0.1))
  }
}
plot(csp[1,])
points((yt33[1,]),col="red")

plot(csp[2,], ylim=c(min(csp[2,]),max((yt3[2,]))))
points((yt33[2,]),col="red")

plot(csp[3,],ylim=c(min(csp[3,]),max((yt3[3,]))))
points((yt33[3,]),col="red")

plot(csp[4,])
points((yt33[4,]),col="red")

plot(csp[5,])
points((yt33[5,]),col="red")
#100 simulations
#sim from random error
vareddysim<-NULL
sigma=0.1
yteddysim<- array(0, dim = c(50, 5, 100))
for(i in 1:50){
  vareddysim[i]<-sigma*exp(rnorm(1,0,0.1))
  for(j in 1:100){
    for(k in 1:5){
      yteddysim[i,k,j]<-exp(rnorm(1,mean=(log(csp[k,j])),sd=(vareddysim[i]) ))
    }
  }
 
}
plot(c(yteddysim[20,5,]),col="red")
points(csp[5,],col="blue")
#sim from additive time space cove
nloc=5
vareddysim2<-vector("list")
yteddysim2<-et2sim<-wt2sim<- array(0, dim = c(50, 5, 100))
for(i in 1:50){
  vareddysim2[[i]]<-sigma.sp[1:nloc,1:nloc]*exp(rnorm(1,0,1))
  wt2sim[i,1,1:5]<-rep(0,5)
  et2sim[i,,]<-t(rmvnorm(ntime,rep(0,nloc),vareddysim2[[i]], method=( "svd")))
  for(j in 1:99){
    for(k in 1:5){
        wt2sim[i,,(j+1)]<-wt2sim[i,,j]+et2sim[i,,j]
        yteddysim2[i,k,j]<-exp(rnorm(1,log(csp[k,j]+et2sim[i,k,j]),0.1))
      }
  }
}

st<-data.frame(x=c(rep(c(x[1:nloc]),
                       ((100)))),y=c(rep(c(y[1:nloc]),(100))),
               z=c(rep(1,nloc*100)))
stsp<-SpatialPoints(st)
times<-rep(c(1:100),rep(nloc,100))
stt<-as.POSIXct(times,origin="2016-01-01")
allval<-c(c(yt3[1,]),c(yt3[2,]),c(yt3[3,]),c(yt3[4,]),c(yt3[5,]),
          c(yt3[6,]),c(yt3[7,]),c(yt3[8,]),c(yt3[9,]),c(yt3[10,]))
stdf<-STIDF(stsp,stt,data.frame(allval=exp(allval)))
head(stdf@data)
quartz()
stplot(stdf)

varst<-variogramST(allval~1,data=stdf)
head(varst)

quartz()
plot(varst,map=T)
quartz()
plot(varst,map=F)
quartz()
plot(varst,wireframe=T)

fit.StVariogram(varst)

distances<-as.matrix(dist(data.frame(cbind(c(x),c(y))),upper=TRUE, diag=TRUE,method="euclidean"))
cat("model{
    for (i in 1:N){
    for(j in 1:nloc){
 omega[i+1,j]~dnorm(aomega[z1[i+1]],bomega[z2[i+1]])T(-(csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))),)
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1,j])
    logyt[j,i+1]~dnorm((csp[j,i+1])+et2[j,z3[i+1]], tausq)
    #logyt[j,i+1]~dnorm((csp[j,i+1])+w[i+1,j], tausq)
    yhat[j,i+1]~dnorm((csp[j,i+1])+et2[j,z3[i+1]], tausq)T(0,)
    #yhat[j,i+1]<-ifelse(yhatu[j,i+1]<0,0,yhatu[j,i+1])
    }
    loglik[i+1]<-logdensity.mnorm(logyt[1:nloc,i+1],(csp[1:nloc,i+1])+et2[1:nloc,z3[i+1]],tausq*I)
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    }
    loglik[1]<-logdensity.mnorm(logyt[1:nloc,1],(csp[1:nloc,1])+w[1,1:nloc],tausq*I)
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-(phi+0.000001)*dist[m,j])
    }
    }
    kinv<-inverse(k)
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    for(j in 2:(m-1)){p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.0001
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.0001
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.0001}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    for(l in 1:m){
    aomega[l]~dnorm(0,10000)
    bomega[l]~dgamma(0.01,0.01)
    aomega2[l]<-aomega[l]-mean(aomega)
    #vari[1:2,1:2,l]~dwish(I,3)
    et[1:nloc,l]~dmnorm(muet,kinv)
    for(k in 1:nloc){
    etw[k,l]<-et[k,l]*p3[l]
    et2[k,l]<-et[k,l]-mean(et[k,1:m])
    }}
    w[1,1:nloc]<-c(0,0,0,0,0)
    for(t in 1:N){
    w[t+1,1:nloc]<-w[t,1:nloc]+et2[1:nloc,z3[t+1]]
    }
    for(j in 1:nloc){
    omega[1,j]~dgamma(0.01,0.01)
    yhat[j,1]~dnorm((csp[j,1])+et2[1,j], tausq)T(0,)
    #yhat[j,1]<-ifelse(yhatu[j,1]<0,0,yhatu[j,1])
    csp[j,1]<-(logyt[j,1])+omega[1,j]
    }
    Gprime~dunif(281,482)
    Dt~dunif(0,3)
    sigma~dgamma(0.01,0.01)
    phi~dgamma(0.01,0.01)
    tausq~dgamma(0.01,0.01)
    alpha~dgamma(0.01,0.01)
    }",file="modeljagseddydp.jag")
jags.eddydp1<-jags.model("modeljagseddydp.jag",data=list("logyt"=(yt3[1:5,]), N=99,nloc=5,"h"=h,"muet"=rep(0,5),
                                                    "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=10,"I"=diag(5)))
update(jags.eddydp1, 1000)
mcmcjags.eddydp1<-jags.samples(jags.eddydp1,
                            c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                            1000)
mcmcsamplesjags.eddydp1<-coda.samples(jags.eddydp1,
                                   c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                   1000)
# all1<-matrix(0,1000,5)
# for(i in 1:5){
# all1[,i]<-(mcmcsamplesjags.eddydp[[1]][,502+(5*i)])
# }
# plot(apply(all1,1,mean))
m1.mcmc.eddydp1<-(as.mcmc(mcmcsamplesjags.eddydp1))
m1.mat.eddydp1<-as.matrix(mcmcsamplesjags.eddydp1)
m1.dat.eddydp1<-as.data.frame(m1.mat.eddydp1)
Dt.post.eddydp1<-m1.dat.eddydp1$Dt
gprime.post.eddydp1<-m1.dat.eddydp1$Gprime
aomega.post.eddydp1<-m1.dat.eddydp1$aomega
bomega.post.eddydp1<-m1.dat.eddydp1$bomega
sigma.post.eddydp1<-m1.dat.eddydp1$sigma
phi.post.eddydp1<-m1.dat.eddydp1$phi
quantile(Dt.post.eddydp1,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddydp1,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddydp1,c(0.025,0.25,0.5,0.75,0.975))
loglik.posteddydp1<-t(mcmcjags.eddydp1$loglik[,,1])
#dim(loglik.posteddydp)
waic(loglik.posteddydp)
mean(crps_sample(y = yt3[1,], dat = mcmcjags.eddydp1$yhat[1,,,1]))+
  mean(crps_sample(y = yt3[2,], dat =mcmcjags.eddydp1$yhat[2,,,1]))+
  mean(crps_sample(y = yt3[3,], dat =mcmcjags.eddydp1$yhat[3,,,1]))+
  mean(crps_sample(y = yt3[4,], dat =mcmcjags.eddydp1$yhat[4,,,1]))+
  mean(crps_sample(y = yt3[5,], dat =mcmcjags.eddydp1$yhat[5,,,1]))
csp.hatdp1<-yhatdp1<-matrix(0, nrow=5, ncol=100)
c.hatlower50e1<-yhatlower50e1<-matrix(0, nrow=5, ncol=100)
c.hatupper50e1<-yhatupper50e1<-matrix(0, nrow=5, ncol=100)
c.hatlower90e1<-yhatlower90e1<-matrix(0, nrow=5, ncol=100)
c.hatupper90e1<-yhatupper90e1<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  for(i in 1:100){
    csp.hatdp1[j,i] <- mean(mcmcsamplesjags.eddydp1[[1]][,2+(5*i+1-(6-j))])
    c.hatlower50e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,2+(5*i+1-(6-j))],0.25)
    c.hatupper50e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,2+(5*i+1-(6-j))],0.75)
    c.hatlower90e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,2+(5*i+1-(6-j))],0.05)
    c.hatupper90e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,2+(5*i+1-(6-j))],0.95)
    yhatdp1[j,i] <- mean(mcmcsamplesjags.eddydp1[[1]][,653+(5*i+1-(6-j))])
    yhatlower50e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,653+(5*i+1-(6-j))],0.25)
    yhatupper50e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,653+(5*i+1-(6-j))],0.75)
    yhatlower90e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,653+(5*i+1-(6-j))],0.05)
    yhatupper90e1[j,i]<- quantile(mcmcsamplesjags.eddydp1[[1]][,653+(5*i+1-(6-j))],0.95)
  }
}
yhatdp1
mean501<-mean901<-NULL
for(i in 1:5){
  mean501[i]<-mean(yhatupper50e1[i,]-yhatlower50e1[i,])
  mean901[i]<-mean(yhatupper90e1[i,]-yhatlower90e1[i,])
}
mean(mean501)
mean(mean901)
cov501A<-cov901A<-matrix(0,100,5)
for(i in 1:length(yt3[1,])){
  for(j in 1:5){
    cov501A[i,j]<-(ifelse(yhatlower50e1[j,i]<=(yt3[j,i])&(yt3[j,i])<=yhatupper50e1[j,i],1,0))
    cov901A[i,j]<-(ifelse(yhatlower90e1[j,i]<=(yt3[j,i])&(yt3[j,i])<=yhatupper90e1[j,i],1,0))
  }}
sum(cov501A)/500
sum(cov901A)/500
cdfdpe11<-apply(mcmcjags.eddydp1$yhat[1,,,1],1,function(x) empirical_cdf(x,(yt3[1,])))
cdfdpe21<-apply(mcmcjags.eddydp1$yhat[2,,,1],1,function(x) empirical_cdf(x,(yt3[2,])))
cdfdpe31<-apply(mcmcjags.eddydp1$yhat[3,,,1],1,function(x) empirical_cdf(x,(yt3[3,])))
cdfdpe41<-apply(mcmcjags.eddydp1$yhat[4,,,1],1,function(x) empirical_cdf(x,(yt3[4,])))
cdfdpe51<-apply(mcmcjags.eddydp1$yhat[5,,,1],1,function(x) empirical_cdf(x,(yt3[5,])))
cdf11<-cdf21<-cdf31<-cdf41<-cdf51<-NULL
for(i in 1:100){
  cdf11[i]<-mean(sapply(cdfdpe11, function(X) X$CDF[i]))
  cdf21[i]<-mean(sapply(cdfdpe21, function(X) X$CDF[i]))
  cdf31[i]<-mean(sapply(cdfdpe31, function(X) X$CDF[i]))
  cdf41[i]<-mean(sapply(cdfdpe41, function(X) X$CDF[i]))
  cdf51[i]<-mean(sapply(cdfdpe51, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,1))
hist((cdf11),breaks = 15, main="L1")
hist((cdf21),breaks = 15, main="L2")
hist((cdf31),breaks = 15, main="L3")
hist((cdf41),breaks = 15, main="L4")
hist((cdf51),breaks = 15, main="L5")
ghatdpl11<-empirical_cdf(yt3[1,],yt3[1,])
ghatdpl21<-empirical_cdf(yt3[2,],yt3[2,])
ghatdpl31<-empirical_cdf(yt3[3,],yt3[3,])
ghatdpl41<-empirical_cdf(yt3[4,],yt3[4,])
ghatdpl51<-empirical_cdf(yt3[5,],yt3[5,])
diffel11<-((cdf11)-ghatdpl11$CDF)
diffel21<-((cdf21)-ghatdpl21$CDF)
diffel31<-((cdf31)-ghatdpl31$CDF)
diffel41<-((cdf41)-ghatdpl41$CDF)
diffel51<-((cdf51)-ghatdpl51$CDF)

plot(diffel1,main="L1")
abline(h=0)
plot(diffel2, main="L2")
abline(h=0)

plot(  csp.hatdp1[1,], type="l")
lines(yhatlower50e1[1,],col="blue")
lines(yhatupper50e1[1,], col="blue")
lines(yhatlower90e1[1,],col="red")
lines(yhatupper90e1[1,], col="red")
lines((yt3[1,]), col="green")
lines(csp[1,], col="purple", lty=3)
plot(  csp.hatdp1[2,], type="l")
lines(yhatlower50e1[2,],col="blue")
lines(yhatupper50e1[2,], col="blue")
lines(yhatlower90e1[2,],col="red")
lines(yhatupper90e1[2,], col="red")
lines((yt3[2,]), col="green")
lines(csp[2,], col="purple", lty=3)

plot(csp[5,],csp.hatdp[5,])
abline(0,1)

plot(csp[1,],col="red")
points(csp.hatdp1[1,],col="blue")
points((yt3[1,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)
plot(csp[2,],col="red")
points(csp.hatdp1[2,],col="blue")
points((yt3[2,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
plot(csp[3,],col="red")
points(csp.hatdp1[3,],col="blue")
points((yt3[3,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct3","chat3","yt3"),bg="white", lty=1,cex=0.8)
plot(csp[4,],col="red")
points(csp.hatdp1[4,],col="blue")
points((yt3[4,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct4","chat4","yt4"),bg="white", lty=1,cex=0.8)
plot(csp[5,],col="red")
points(csp.hatdp1[5,],col="blue")
points((yt3[5,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct5","chat5","yt5"),bg="white", lty=1,cex=0.8)
#all sim
jagseddysim<-mcmcjags.eddysim<-mcmcsamplesjags.eddysim<-vector("list")
for(i in 1:50){
  jagseddysim[[i]]<-jags.model("modeljagseddydp.jag",data=list("logyt"=(yteddysim[i,,]), N=99,nloc=5,"h"=h,"muet"=rep(0,5),
                                                               "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=10,"I"=diag(5)))
  update(jagseddysim[[i]], 1000)
  mcmcjags.eddysim[[i]]<-jags.samples(jagseddysim[[i]],
                                      c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                      1000)
  mcmcsamplesjags.eddysim[[i]]<-coda.samples(jagseddysim[[i]],
                                             c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                             1000)
}
m1.mcmc.eddysim<-m1.mat.eddysim<-m1.dat.eddysim<-vector("list")
yhateddysim<-yhatlower50eddysim<-yhatupper50eddysim<-yhatlower90eddysim<-yhatupper90eddysim<-
  chateddysim<-chatlower50eddysim<-chatupper50eddysim<-chatlower90eddysim<-chatupper90eddysim<-array(0,dim = c(50, 100, 5))
for(i in 1:50){
  m1.mcmc.eddysim[[i]]<-(as.mcmc(mcmcsamplesjags.eddysim[[i]]))
  m1.mat.eddysim[[i]]<-as.matrix(mcmcsamplesjags.eddysim[[i]])
  m1.dat.eddysim[[i]]<-as.data.frame(m1.mat.eddysim[[i]])
  for(j in 1:100){
      yhateddysim[i,j,]<- apply(((mcmcjags.eddysim[[i]]$yhat[,j,,])), 1, mean)
      yhatlower50eddysim[i,j,]<- apply(((mcmcjags.eddysim[[i]]$yhat[,j,,])), 1, function(x) quantile(x,0.25))
      yhatupper50eddysim[i,j,]<-apply((mcmcjags.eddysim[[i]]$yhat[,j,,]), 1, function(x) quantile(x,0.75))
      yhatlower90eddysim[i,j,]<- apply(mcmcjags.eddysim[[i]]$yhat[,j,,], 1, function(x) quantile(x,0.05))
      yhatupper90eddysim[i,j,]<-apply(mcmcjags.eddysim[[i]]$yhat[,j,,], 1, function(x) quantile(x,0.95))
      chateddysim[i,j,]<- apply(mcmcjags.eddysim[[i]]$csp[,j,,], 1, mean)
      chatlower50eddysim[i,j,]<- apply(mcmcjags.eddysim[[i]]$csp[,j,,], 1, function(x) quantile(x,0.25))
      chatupper50eddysim[i,j,]<-apply(mcmcjags.eddysim[[i]]$csp[,j,,], 1, function(x) quantile(x,0.75))
      chatlower90eddysim[i,j,]<- apply(mcmcjags.eddysim[[i]]$csp[,j,,], 1, function(x) quantile(x,0.05))
      chatupper90eddysim[i,j,]<-apply(mcmcjags.eddysim[[i]]$csp[,j,,], 1, function(x) quantile(x,0.95))
  }
}
cov50eddysim<-cov90eddysim<-mseeddysim<-crpseddysim<-array(0, dim=c(50,100,5))
for(i in 1:50){
  for (k in 1:5){
    cov50eddysim[i,,k]<-(ifelse(chatlower50eddysim[i,,k]<=(csp[k,])&(csp[k,])<=chatupper50eddysim[i,,k],1,0))
    cov90eddysim[i,,k]<-(ifelse(chatlower90eddysim[i,,k]<=(csp[k,])&(csp[k,])<=chatupper90eddysim[i,,k],1,0))
    mseeddysim[i,,k]<-sum((chateddysim[i,,k]-(csp[,k]))^2)/100
    crpseddysim[i,,k]<-mean(crps_sample(y = c(yteddysim[i,k,]), dat =mcmcjags.eddysim[[i]]$yhat[k,,,]))
  }
 
}
mean((apply(cov50eddysim,1 ,function(x) sum(x)/500)))
mean((apply(cov90eddysim,1, function(x) sum(x)/500)))
mean((apply(mseeddysim,1,mean)))
mean((apply(crpseddysim,1,mean)))
gprime.post.eddysim<-Dt.post.eddysim<-vector("list")
for(i in 1:50){
  Dt.post.eddysim[[i]]<-(m1.dat.eddysim[[i]])$Dt
  gprime.post.eddysim[[i]]<-(m1.dat.eddysim[[i]])$Gprime
}
quantile(unlist(Dt.post.eddysim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(gprime.post.eddysim),c(0.025,0.25,0.5,0.75,0.975))
plot(chateddysim[i,,5])
points(c(yteddysim[i,5,]),col="red")
points(csp[5,],col="blue")
cdfdpe11sim<-cdfdpe21sim<-cdfdpe31sim<-cdfdpe41sim<-cdfdpe51sim<-vector("list")
for(i in 1:50){
  cdfdpe11sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[1,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,1,])))
  cdfdpe21sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[2,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,2,])))
  cdfdpe31sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[3,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,3,])))
  cdfdpe41sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[4,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,4,])))
  cdfdpe51sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[5,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,5,])))
}

cdf11sim<-cdf21sim<-cdf31sim<-cdf41sim<-cdf51sim<-matrix(0,50,100)
for(i in 1:50){
  for(j in 1:100){
  cdf11sim[i,j]<-mean(sapply(cdfdpe11sim[[i]], function(X) X$CDF[j]))
  cdf21sim[i,j]<-mean(sapply(cdfdpe21sim[[i]], function(X) X$CDF[j]))
  cdf31sim[i,j]<-mean(sapply(cdfdpe31sim[[i]], function(X) X$CDF[j]))
  cdf41sim[i,j]<-mean(sapply(cdfdpe41sim[[i]], function(X) X$CDF[j]))
  cdf51sim[i,j]<-mean(sapply(cdfdpe51sim[[i]], function(X) X$CDF[j]))
}}
hist((cdf11sim[1,]),breaks = 15, main="L1")
hist((cdf21sim[20,]),breaks = 15, main="L2")
hist((cdf31sim[40,]),breaks = 15, main="L3")
hist((cdf41sim[50,]),breaks = 15, main="L4")
hist((cdf51sim[11,]),breaks = 15, main="L5")
ghatdpl11sim<-ghatdpl21sim<-ghatdpl31sim<-ghatdpl41sim<-ghatdpl51sim<-diffel11sim<-diffel21sim<-diffel31sim<-diffel41sim<-diffel51sim<-vector("list")
for(i in 1:50){
  ghatdpl11sim[[i]]<-empirical_cdf(yteddysim[i,1,],yteddysim[i,1,])
  ghatdpl21sim[[i]]<-empirical_cdf(yteddysim[i,2,],yteddysim[i,2,])
  ghatdpl31sim[[i]]<-empirical_cdf(yteddysim[i,3,],yteddysim[i,3,])
  ghatdpl41sim[[i]]<-empirical_cdf(yteddysim[i,4,],yteddysim[i,4,])
  ghatdpl51sim[[i]]<-empirical_cdf(yteddysim[i,5,],yteddysim[i,5,])
  diffel11sim[[i]]<-((cdf11sim[i,])-ghatdpl11sim[[i]]$CDF)
  diffel21sim[[i]]<-((cdf21sim[i,])-ghatdpl21sim[[i]]$CDF)
  diffel31sim[[i]]<-((cdf31sim[i,])-ghatdpl31sim[[i]]$CDF)
  diffel41sim[[i]]<-((cdf41sim[i,])-ghatdpl41sim[[i]]$CDF)
  diffel51sim[[i]]<-((cdf51sim[i,])-ghatdpl51sim[[i]]$CDF)
  
}
par(mfrow=c(1,5))
plot(unlist(diffel11sim[[1]]),main="L1",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel11sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel21sim[[1]]),main="L2",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel21sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel31sim[[1]]),main="L3",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel31sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel41sim[[1]]),main="L4",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel41sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel51sim[[1]]),main="L5",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel51sim[[i]]),col=i, lwd=2)
}
abline(h=0)
#model2
jags.eddydp2<-jags.model("modeljagseddydp.jag",data=list("logyt"=(yt32[1:5,]), N=99,nloc=5,"h"=h,"muet"=rep(0,5),
                                                         "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=10,"I"=diag(5)))
update(jags.eddydp2, 1000)
mcmcjags.eddydp2<-jags.samples(jags.eddydp2,
                               c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                               1000)
mcmcsamplesjags.eddydp2<-coda.samples(jags.eddydp2,
                                      c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                      1000)
m1.mcmc.eddydp2<-(as.mcmc(mcmcsamplesjags.eddydp2))
m1.mat.eddydp2<-as.matrix(mcmcsamplesjags.eddydp2)
m1.dat.eddydp2<-as.data.frame(m1.mat.eddydp2)
Dt.post.eddydp2<-m1.dat.eddydp2$Dt
gprime.post.eddydp2<-m1.dat.eddydp2$Gprime
aomega.post.eddydp2<-m1.dat.eddydp2$aomega
bomega.post.eddydp2<-m1.dat.eddydp2$bomega
sigma.post.eddydp2<-m1.dat.eddydp2$sigma
phi.post.eddydp2<-m1.dat.eddydp2$phi
quantile(Dt.post.eddydp2,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddydp2,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddydp2,c(0.025,0.25,0.5,0.75,0.975))
loglik.posteddydp2<-t(mcmcjags.eddydp2$loglik[,,1])
waic(loglik.posteddydp2)
mean(crps_sample(y = yt32[1,], dat = mcmcjags.eddydp2$yhat[1,,,1]))+
  mean(crps_sample(y = yt32[2,], dat =mcmcjags.eddydp2$yhat[2,,,1]))+
  mean(crps_sample(y = yt32[3,], dat =mcmcjags.eddydp2$yhat[3,,,1]))+
  mean(crps_sample(y = yt32[4,], dat =mcmcjags.eddydp2$yhat[4,,,1]))+
  mean(crps_sample(y = yt32[5,], dat =mcmcjags.eddydp2$yhat[5,,,1]))
csp.hatdp2<-yhatdp2<-matrix(0, nrow=5, ncol=100)
c.hatlower50e2<-yhatlower50e2<-matrix(0, nrow=5, ncol=100)
c.hatupper50e2<-yhatupper50e2<-matrix(0, nrow=5, ncol=100)
c.hatlower90e2<-yhatlower90e2<-matrix(0, nrow=5, ncol=100)
c.hatupper90e2<-yhatupper90e2<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  for(i in 1:100){
    csp.hatdp2[j,i] <- mean(mcmcsamplesjags.eddydp2[[1]][,2+(5*i+1-(6-j))])
    c.hatlower50e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,2+(5*i+1-(6-j))],0.25)
    c.hatupper50e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,2+(5*i+1-(6-j))],0.75)
    c.hatlower90e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,2+(5*i+1-(6-j))],0.05)
    c.hatupper90e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,2+(5*i+1-(6-j))],0.95)
    yhatdp2[j,i] <- mean(mcmcsamplesjags.eddydp2[[1]][,653+(5*i+1-(6-j))])
    yhatlower50e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,653+(5*i+1-(6-j))],0.25)
    yhatupper50e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,653+(5*i+1-(6-j))],0.75)
    yhatlower90e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,653+(5*i+1-(6-j))],0.05)
    yhatupper90e2[j,i]<- quantile(mcmcsamplesjags.eddydp2[[1]][,653+(5*i+1-(6-j))],0.95)
  }
}
yhatdp2
mean502<-mean902<-NULL
for(i in 1:5){
  mean502[i]<-mean(yhatupper50e2[i,]-yhatlower50e2[i,])
  mean902[i]<-mean(yhatupper90e2[i,]-yhatlower90e2[i,])
}
mean(mean502)
mean(mean902)
cov502A<-cov902A<-matrix(0,100,5)
for(i in 1:length(yt3[1,])){
  for(j in 1:5){
    cov502A[i,j]<-(ifelse(yhatlower50e2[j,i]<=(yt32[j,i])&(yt32[j,i])<=yhatupper50e2[j,i],1,0))
    cov902A[i,j]<-(ifelse(yhatlower90e2[j,i]<=(yt32[j,i])&(yt32[j,i])<=yhatupper90e2[j,i],1,0))
  }}
sum(cov502A)/500
sum(cov902A)/500
cdfdpe12<-apply(mcmcjags.eddydp2$yhat[1,,,1],1,function(x) empirical_cdf(x,(yt32[1,])))
cdfdpe22<-apply(mcmcjags.eddydp2$yhat[2,,,1],1,function(x) empirical_cdf(x,(yt32[2,])))
cdfdpe32<-apply(mcmcjags.eddydp2$yhat[3,,,1],1,function(x) empirical_cdf(x,(yt32[3,])))
cdfdpe42<-apply(mcmcjags.eddydp2$yhat[4,,,1],1,function(x) empirical_cdf(x,(yt32[4,])))
cdfdpe52<-apply(mcmcjags.eddydp2$yhat[5,,,1],1,function(x) empirical_cdf(x,(yt32[5,])))
cdf12<-cdf22<-cdf32<-cdf42<-cdf52<-NULL
for(i in 1:100){
  cdf12[i]<-mean(sapply(cdfdpe12, function(X) X$CDF[i]))
  cdf22[i]<-mean(sapply(cdfdpe22, function(X) X$CDF[i]))
  cdf32[i]<-mean(sapply(cdfdpe32, function(X) X$CDF[i]))
  cdf42[i]<-mean(sapply(cdfdpe42, function(X) X$CDF[i]))
  cdf52[i]<-mean(sapply(cdfdpe52, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,1))
hist((cdf12),breaks = 15, main="L1")
hist((cdf22),breaks = 15, main="L2")
hist((cdf32),breaks = 15, main="L3")
hist((cdf42),breaks = 15, main="L4")
hist((cdf52),breaks = 15, main="L5")
ghatdpl12<-empirical_cdf(yt32[1,],yt32[1,])
ghatdpl22<-empirical_cdf(yt32[2,],yt32[2,])
ghatdpl32<-empirical_cdf(yt32[3,],yt32[3,])
ghatdpl42<-empirical_cdf(yt32[4,],yt32[4,])
ghatdpl52<-empirical_cdf(yt32[5,],yt32[5,])
diffel12<-((cdf12)-ghatdpl12$CDF)
diffel22<-((cdf22)-ghatdpl22$CDF)
diffel32<-((cdf32)-ghatdpl32$CDF)
diffel42<-((cdf42)-ghatdpl42$CDF)
diffel52<-((cdf52)-ghatdpl52$CDF)

plot(diffel12,main="L1")
abline(h=0)
plot(diffel22, main="L2")
abline(h=0)

plot(  csp.hatdp2[1,], type="l")
lines(yhatlower50e2[1,],col="blue")
lines(yhatupper50e2[1,], col="blue")
lines(yhatlower90e2[1,],col="red")
lines(yhatupper90e2[1,], col="red")
lines((yt32[1,]), col="green")
lines(csp[1,], col="purple", lty=3)
plot(  csp.hatdp2[2,], type="l")
lines(yhatlower50e2[2,],col="blue")
lines(yhatupper50e2[2,], col="blue")
lines(yhatlower90e2[2,],col="red")
lines(yhatupper90e2[2,], col="red")
lines((yt32[2,]), col="green")
lines(csp[2,], col="purple", lty=3)

plot(csp[5,],csp.hatdp[5,])
abline(0,1)

plot(csp[1,],col="red")
points(csp.hatdp2[1,],col="blue")
points((yt32[1,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)
plot(csp[2,],col="red")
points(csp.hatdp2[2,],col="blue")
points((yt32[2,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
plot(csp[3,],col="red")
points(csp.hatdp2[3,],col="blue")
points((yt32[3,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct3","chat3","yt3"),bg="white", lty=1,cex=0.8)
plot(csp[4,],col="red")
points(csp.hatdp2[4,],col="blue")
points((yt32[4,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct4","chat4","yt4"),bg="white", lty=1,cex=0.8)
plot(csp[5,],col="red")
points(csp.hatdp2[5,],col="blue")
points((yt32[5,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct5","chat5","yt5"),bg="white", lty=1,cex=0.8)
#all sim

#model3
jags.eddydp3<-jags.model("modeljagseddydp.jag",data=list("logyt"=(yt33[1:5,]), N=99,nloc=5,"h"=h,"muet"=rep(0,5),
                                                         "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=10,"I"=diag(5)))
update(jags.eddydp3, 1000)
mcmcjags.eddydp3<-jags.samples(jags.eddydp3,
                               c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                               1000)
mcmcsamplesjags.eddydp3<-coda.samples(jags.eddydp3,
                                      c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                      1000)
m1.mcmc.eddydp3<-(as.mcmc(mcmcsamplesjags.eddydp3))
m1.mat.eddydp3<-as.matrix(mcmcsamplesjags.eddydp3)
m1.dat.eddydp3<-as.data.frame(m1.mat.eddydp3)
Dt.post.eddydp3<-m1.dat.eddydp3$Dt
gprime.post.eddydp3<-m1.dat.eddydp3$Gprime
aomega.post.eddydp3<-m1.dat.eddydp3$aomega
bomega.post.eddydp3<-m1.dat.eddydp3$bomega
sigma.post.eddydp3<-m1.dat.eddydp3$sigma
phi.post.eddydp3<-m1.dat.eddydp3$phi
quantile(Dt.post.eddydp3,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddydp3,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddydp3,c(0.025,0.25,0.5,0.75,0.975))
loglik.posteddydp3<-t(mcmcjags.eddydp3$loglik[,,1])
waic(loglik.posteddydp3)
mean(crps_sample(y = yt33[1,], dat = mcmcjags.eddydp3$yhat[1,,,1]))+
  mean(crps_sample(y = yt33[2,], dat =mcmcjags.eddydp3$yhat[2,,,1]))+
  mean(crps_sample(y = yt33[3,], dat =mcmcjags.eddydp3$yhat[3,,,1]))+
  mean(crps_sample(y = yt33[4,], dat =mcmcjags.eddydp3$yhat[4,,,1]))+
  mean(crps_sample(y = yt33[5,], dat =mcmcjags.eddydp3$yhat[5,,,1]))
csp.hatdp3<-yhatdp3<-matrix(0, nrow=5, ncol=100)
c.hatlower50e3<-yhatlower50e3<-matrix(0, nrow=5, ncol=100)
c.hatupper50e3<-yhatupper50e3<-matrix(0, nrow=5, ncol=100)
c.hatlower90e3<-yhatlower90e3<-matrix(0, nrow=5, ncol=100)
c.hatupper90e3<-yhatupper90e3<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  for(i in 1:100){
    csp.hatdp3[j,i] <- mean(mcmcsamplesjags.eddydp3[[1]][,2+(5*i+1-(6-j))])
    c.hatlower50e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,2+(5*i+1-(6-j))],0.25)
    c.hatupper50e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,2+(5*i+1-(6-j))],0.75)
    c.hatlower90e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,2+(5*i+1-(6-j))],0.05)
    c.hatupper90e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,2+(5*i+1-(6-j))],0.95)
    yhatdp3[j,i] <- mean(mcmcsamplesjags.eddydp3[[1]][,653+(5*i+1-(6-j))])
    yhatlower50e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,653+(5*i+1-(6-j))],0.25)
    yhatupper50e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,653+(5*i+1-(6-j))],0.75)
    yhatlower90e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,653+(5*i+1-(6-j))],0.05)
    yhatupper90e3[j,i]<- quantile(mcmcsamplesjags.eddydp3[[1]][,653+(5*i+1-(6-j))],0.95)
  }
}
yhatdp3
mean503<-mean903<-NULL
for(i in 1:5){
  mean503[i]<-mean(yhatupper50e3[i,]-yhatlower50e3[i,])
  mean903[i]<-mean(yhatupper90e3[i,]-yhatlower90e3[i,])
}
mean(mean503)
mean(mean903)
cov503A<-cov903A<-matrix(0,100,5)
for(i in 1:length(yt3[1,])){
  for(j in 1:5){
    cov503A[i,j]<-(ifelse(yhatlower50e3[j,i]<=(yt33[j,i])&(yt33[j,i])<=yhatupper50e3[j,i],1,0))
    cov903A[i,j]<-(ifelse(yhatlower90e3[j,i]<=(yt33[j,i])&(yt33[j,i])<=yhatupper90e3[j,i],1,0))
  }}
sum(cov503A)/500
sum(cov903A)/500
cdfdpe13<-apply(mcmcjags.eddydp3$yhat[1,,,1],1,function(x) empirical_cdf(x,(yt33[1,])))
cdfdpe23<-apply(mcmcjags.eddydp3$yhat[2,,,1],1,function(x) empirical_cdf(x,(yt33[2,])))
cdfdpe33<-apply(mcmcjags.eddydp3$yhat[3,,,1],1,function(x) empirical_cdf(x,(yt33[3,])))
cdfdpe43<-apply(mcmcjags.eddydp3$yhat[4,,,1],1,function(x) empirical_cdf(x,(yt33[4,])))
cdfdpe53<-apply(mcmcjags.eddydp3$yhat[5,,,1],1,function(x) empirical_cdf(x,(yt33[5,])))
cdf13<-cdf23<-cdf33<-cdf43<-cdf53<-NULL
for(i in 1:100){
  cdf13[i]<-mean(sapply(cdfdpe13, function(X) X$CDF[i]))
  cdf23[i]<-mean(sapply(cdfdpe23, function(X) X$CDF[i]))
  cdf33[i]<-mean(sapply(cdfdpe33, function(X) X$CDF[i]))
  cdf43[i]<-mean(sapply(cdfdpe43, function(X) X$CDF[i]))
  cdf53[i]<-mean(sapply(cdfdpe53, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,1))
hist((cdf13),breaks = 15, main="L1")
hist((cdf23),breaks = 15, main="L2")
hist((cdf33),breaks = 15, main="L3")
hist((cdf43),breaks = 15, main="L4")
hist((cdf53),breaks = 15, main="L5")
ghatdpl13<-empirical_cdf(yt33[1,],yt33[1,])
ghatdpl23<-empirical_cdf(yt33[2,],yt33[2,])
ghatdpl33<-empirical_cdf(yt33[3,],yt33[3,])
ghatdpl43<-empirical_cdf(yt33[4,],yt33[4,])
ghatdpl53<-empirical_cdf(yt33[5,],yt33[5,])
diffel13<-((cdf13)-ghatdpl13$CDF)
diffel23<-((cdf23)-ghatdpl23$CDF)
diffel33<-((cdf33)-ghatdpl33$CDF)
diffel43<-((cdf43)-ghatdpl43$CDF)
diffel53<-((cdf53)-ghatdpl53$CDF)

plot(diffel13,main="L1")
abline(h=0)
plot(diffel23, main="L2")
abline(h=0)

plot(  csp.hatdp3[1,], type="l")
lines(yhatlower50e3[1,],col="blue")
lines(yhatupper50e3[1,], col="blue")
lines(yhatlower90e3[1,],col="red")
lines(yhatupper90e3[1,], col="red")
lines((yt33[1,]), col="green")
lines(csp[1,], col="purple", lty=3)
plot(  csp.hatdp3[2,], type="l")
lines(yhatlower50e3[2,],col="blue")
lines(yhatupper50e3[2,], col="blue")
lines(yhatlower90e3[2,],col="red")
lines(yhatupper90e3[2,], col="red")
lines((yt33[2,]), col="green")
lines(csp[2,], col="purple", lty=3)

plot(csp[5,],csp.hatdp[5,])
abline(0,1)

plot(csp[1,],col="red")
points(csp.hatdp3[1,],col="blue")
points((yt33[1,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)
plot(csp[2,],col="red")
points(csp.hatdp3[2,],col="blue")
points((yt33[2,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
plot(csp[3,],col="red")
points(csp.hatdp3[3,],col="blue")
points((yt33[3,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct3","chat3","yt3"),bg="white", lty=1,cex=0.8)
plot(csp[4,],col="red")
points(csp.hatdp3[4,],col="blue")
points((yt33[4,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct4","chat4","yt4"),bg="white", lty=1,cex=0.8)
plot(csp[5,],col="red")
points(csp.hatdp3[5,],col="blue")
points((yt33[5,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct5","chat5","yt5"),bg="white", lty=1,cex=0.8)

cat("model{
    for (i in 1:N){
    omega[i+1]~dnorm(aomega2[z1[i+1]],bomega[z2[i+1]])T(-(csp[1,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[1])^2/(4*Dt*i)))),)
    for(j in 1:nloc){
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1])
    logyt[j,i+1]~dnorm((csp[j,i+1])+et2[j*(i+1),z3[i+1]], tausq)
    yhat[j,i+1]~dnorm((csp[j,i+1])+et2[j*(i+1),z3[i+1]], tausq)T(0,)
    #yhat[j,i+1]<-ifelse(yhatu[j,i+1]<0,0,yhatu[j,i+1])
    }
    #loglik[i+1]<-logdensity.mnorm(logyt[1:nloc,i+1],(csp[1:nloc,i+1])+et2[(1:nloc)*(i+1),z3[i+1]],tausq*I)
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    }
    #loglik[1]<-logdensity.mnorm(logyt[1:nloc,1],(csp[1:nloc,1])+et2[(1:nloc)*(1),z3[2]],tausq*I)
    for(j in 1:nloc){
    for(l in 1:nloc){
    for(r in 1:(N+1)){
    for(s in 1:(N+1)){
    k[(l*r)+((l-1)*(N-(r-1))),(j*s)+((j-1)*(N-(s-1)))]<-
    (sigma/(a*((r-s))^2+1))*exp(-b*dist[l,j]/(a*((r-s))^2+1))
    }}}}
    kinv<-inverse(k)
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    for(j in 2:(m-1)){p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.0001
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.0001
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.0001}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    for(l in 1:m){
    aomega[l]~dnorm(0,100)
    bomega[l]~dgamma(0.01,0.01)
    aomega2[l]<-aomega[l]-mean(aomega)
    et[1:(nloc*(N+1)),l]~dmnorm(muet,kinv)
     for(k in 1:(nloc*(N+1))){
     etw[k,l]<-et[k,l]*p3[l]
     et2[k,l]<-et[k,l]-mean(et[k,1:m])
     }
    }
    for(j in 1:nloc){
    yhat[j,1]~dnorm((csp[j,1])+et2[j*(1),z3[2]], tausq)T(0,)
    #yhat[j,1]<-ifelse(yhatu[j,1]<0,0,yhatu[j,1])
    csp[j,1]<-(logyt[j,1])+omega[1]
    }
    omega[1]~dgamma(0.01,0.01)
    Gprime~dunif(281,482)
    Dt~dunif(0,3)
    sigma~dgamma(0.01,0.01)
    b<-1
    a<-1
    tausq~dgamma(0.01,0.01)
    alpha~dgamma(0.01,0.01)
    }",file="modeljagseddydpns.jag")
jags.eddydpns<-jags.model("modeljagseddydpns.jag",data=list("logyt"=(yt3[1:5,]), N=99,nloc=5,"h"=h,"muet"=rep(0,100*5),
                                                        "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=5,"I"=diag(5)))
update(jags.eddydpns, 1000)
mcmcjags.eddydpns<-jags.samples(jags.eddydpns,
                              c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                              1000)
mcmcsamplesjags.eddydpns<-coda.samples(jags.eddydpns,
                                     c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                     1000)
all1ns<-matrix(0,1000,5)
for(i in 1:5){
  all1ns[,i]<-(mcmcsamplesjags.eddydpns[[1]][,502+(5*i)])
}
plot(apply(all1,1,mean))
m1.mcmc.eddydpns<-(as.mcmc(mcmcsamplesjags.eddydpns))
m1.mat.eddydpns<-as.matrix(mcmcsamplesjags.eddydpns)
m1.dat.eddydpns<-as.data.frame(m1.mat.eddydpns)
Dt.post.eddydpns<-m1.dat.eddydpns$Dt
gprime.post.eddydpns<-m1.dat.eddydpns$Gprime
aomega.post.eddydpns<-m1.dat.eddydpns$aomega
bomega.post.eddydpns<-m1.dat.eddydpns$bomega
sigma.post.eddydpns<-m1.dat.eddydpns$sigma
phi.post.eddydpns<-m1.dat.eddydpns$phi
quantile(Dt.post.eddydpns,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddydpns,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddydpns,c(0.025,0.25,0.5,0.75,0.975))
#loglik.posteddydp<-t(mcmcjags.eddydp$loglik[,,1])
#dim(loglik.posteddydp)
#waic(loglik.posteddydp[,2:100])
mean(crps_sample(y = yt3[1,], dat = mcmcjags.eddydpns$yhat[1,,,1]))+
  mean(crps_sample(y = yt3[2,], dat =mcmcjags.eddydpns$yhat[2,,,1]))+
  mean(crps_sample(y = yt3[3,], dat =mcmcjags.eddydpns$yhat[3,,,1]))+
  mean(crps_sample(y = yt3[4,], dat =mcmcjags.eddydpns$yhat[4,,,1]))+
  mean(crps_sample(y = yt3[5,], dat =mcmcjags.eddydpns$yhat[5,,,1]))
csp.hatdpns<-yhatdpns<-matrix(0, nrow=5, ncol=100)
c.hatlower50ens<-yhatlower50ens<-matrix(0, nrow=5, ncol=100)
c.hatupper50ens<-yhatupper50ens<-matrix(0, nrow=5, ncol=100)
c.hatlower90ens<-yhatlower90ens<-matrix(0, nrow=5, ncol=100)
c.hatupper90ens<-yhatupper90ens<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  yhatdpns[j,] <- apply( mcmcjags.eddydpns$yhat[j,,,1],1,mean)
  yhatlower50ens[j,]<- apply( mcmcjags.eddydpns$yhat[j,,,1],1,function(x) quantile(x,0.25))
  yhatupper50ens[j,]<- apply( mcmcjags.eddydpns$yhat[j,,,1],1,function(x) quantile(x,0.75))
  yhatlower90ens[j,]<- apply( mcmcjags.eddydpns$yhat[j,,,1],1,function(x) quantile(x,0.05))
  yhatupper90ens[j,]<- apply( mcmcjags.eddydpns$yhat[j,,,1],1,function(x) quantile(x,0.95))
  for(i in 1:100){
    csp.hatdpns[j,i] <- mean(mcmcsamplesjags.eddydpns[[1]][,2+(5*i+1-(6-j))])
    c.hatlower50ens[j,i]<- quantile(mcmcsamplesjags.eddydpns[[1]][,2+(5*i+1-(6-j))],0.25)
    c.hatupper50ens[j,i]<- quantile(mcmcsamplesjags.eddydpns[[1]][,2+(5*i+1-(6-j))],0.75)
    c.hatlower90ens[j,i]<- quantile(mcmcsamplesjags.eddydpns[[1]][,2+(5*i+1-(6-j))],0.05)
    c.hatupper90ens[j,i]<- quantile(mcmcsamplesjags.eddydpns[[1]][,2+(5*i+1-(6-j))],0.95)
  }
}
yhatdp
csp.hatdpns
plot(csp[5,],csp.hatdpns[5,])
abline(0,1)
mean50ns<-mean90ns<-NULL
for(i in 1:5){
  mean50ns[i]<-mean(yhatupper50ens[i,]-yhatlower50ens[i,])
  mean90ns[i]<-mean(yhatupper90ens[i,]-yhatlower90ens[i,])
}
mean(mean50ns)
mean(mean90ns)
cov501b<-cov901b<-matrix(0,100,5)
for(i in 1:length(yt3[1,])){
  for(j in 1:5){
    cov501b[i,j]<-(ifelse(yhatlower50ens[j,i]<=(yt3[j,i])&(yt3[j,i])<=yhatupper50ens[j,i],1,0))
    cov901b[i,j]<-(ifelse(yhatlower90ens[j,i]<=(yt3[j,i])&(yt3[j,i])<=yhatupper90ens[j,i],1,0))
  }}
sum(cov501b)/500
sum(cov901b)/500
cdfdpe1ns<-apply(mcmcjags.eddydpns$yhat[1,,,1],1,function(x) empirical_cdf(x,yt3[1,]))
cdfdpe2ns<-apply(mcmcjags.eddydpns$yhat[2,,,1],1,function(x) empirical_cdf(x,yt3[2,]))
cdfdpe3ns<-apply(mcmcjags.eddydpns$yhat[3,,,1],1,function(x) empirical_cdf(x,yt3[3,]))
cdfdpe4ns<-apply(mcmcjags.eddydpns$yhat[4,,,1],1,function(x) empirical_cdf(x,yt3[4,]))
cdfdpe5ns<-apply(mcmcjags.eddydpns$yhat[5,,,1],1,function(x) empirical_cdf(x,yt3[5,]))
cdf1ns<-cdf2ns<-cdf3ns<-cdf4ns<-cdf5ns<-NULL
for(i in 1:100){
  cdf1ns[i]<-mean(sapply(cdfdpe1ns, function(X) X$CDF[i]))
  cdf2ns[i]<-mean(sapply(cdfdpe2ns, function(X) X$CDF[i]))
  cdf3ns[i]<-mean(sapply(cdfdpe3ns, function(X) X$CDF[i]))
  cdf4ns[i]<-mean(sapply(cdfdpe4ns, function(X) X$CDF[i]))
  cdf5ns[i]<-mean(sapply(cdfdpe5ns, function(X) X$CDF[i]))
}

quartz()
par(mfrow=c(2,5))
hist((cdf1ns),breaks = 15, main="L1")
hist((cdf2ns),breaks = 15, main="L2")
hist((cdf3ns),breaks = 15, main="L3")
hist((cdf4ns),breaks = 15, main="L4")
hist((cdf5ns),breaks = 15, main="L5")
#additive
hist((cdf1),breaks = 15, main="L1")
hist((cdf2),breaks = 15, main="L2")
hist((cdf3),breaks = 15, main="L3")
hist((cdf4),breaks = 15, main="L4")
hist((cdf5),breaks = 15, main="L5")
diffel1ns<-((cdf1ns)-ghatdpl1$CDF)
diffel2ns<-((cdf2ns)-ghatdpl2$CDF)
diffel3ns<-((cdf3ns)-ghatdpl3$CDF)
diffel4ns<-((cdf4ns)-ghatdpl4$CDF)
diffel5ns<-((cdf5ns)-ghatdpl5$CDF)
quartz()
par(mfrow=c(2,5))
plot(diffel1ns,main="L1")
abline(h=0)
plot(diffel2ns, main="L2")
abline(h=0)
plot(diffel3ns,main="L3")
abline(h=0)
plot(diffel4ns, main="L4")
abline(h=0)
plot(diffel5ns,main="L5")
abline(h=0)


plot(diffel1,main="L1")
abline(h=0)
plot(diffel2, main="L2")
abline(h=0)
plot(diffel3,main="L3")
abline(h=0)
plot(diffel4, main="L4")
abline(h=0)
plot(diffel5,main="L5")
abline(h=0)
#sim1
quartz()
par(mfrow=c(2,5))
plot((yt3[1,]), col="green", type="l",main="Simulation1")
lines(  csp.hatdpns[1,], main="L1")
lines(yhatlower50ens[1,],col="blue")
lines(yhatupper50ens[1,], col="blue")
lines(yhatlower90ens[1,],col="red")
lines(yhatupper90ens[1,], col="red")
lines(csp[1,], col="purple", lty=3)
plot((yt3[2,]), col="green", type="l")
lines(  csp.hatdpns[2,], main="L2")
lines(yhatlower50ens[2,],col="blue")
lines(yhatupper50ens[2,], col="blue")
lines(yhatlower90ens[2,],col="red")
lines(yhatupper90ens[2,], col="red")
lines(csp[2,], col="purple", lty=3)
plot( (yt3[3,]), col="green" , type="l", main="L3")
lines(yhatlower50ens[3,],col="blue")
lines(yhatupper50ens[3,], col="blue")
lines(yhatlower90ens[3,],col="red")
lines(yhatupper90ens[3,], col="red")
lines(csp.hatdpns[3,])
lines(csp[3,], col="purple", lty=3)
plot( (yt3[4,]), col="green", type="l", main="L4")
lines(yhatlower50ens[4,],col="blue")
lines(yhatupper50ens[4,], col="blue")
lines(yhatlower90ens[4,],col="red")
lines(yhatupper90ens[4,], col="red")
lines( csp.hatdpns[4,])
lines(csp[4,], col="purple", lty=3)
plot((yt3[5,]), col="green"  , type="l", main="L5")
lines(yhatlower50ens[5,],col="blue")
lines(yhatupper50ens[5,], col="blue")
lines(yhatlower90ens[5,],col="red")
lines(yhatupper90ens[5,], col="red")
lines(csp.hatdpns[5,])
lines(csp[5,], col="purple", lty=3)


plot( (yt3[1,]), col="green" , type="l", main="L1")
lines(yhatlower50e1[1,],col="blue")
lines(yhatupper50e1[1,], col="blue")
lines(yhatlower90e1[1,],col="red")
lines(yhatupper90e1[1,], col="red")
lines(csp.hatdp1[1,])
lines(csp[1,], col="purple", lty=3)
plot((yt3[2,]), col="green" , type="l", main="L2")
lines(yhatlower50e1[2,],col="blue")
lines(yhatupper50e1[2,], col="blue")
lines(yhatlower90e1[2,],col="red")
lines(yhatupper90e1[2,], col="red")
lines( csp.hatdp1[2,])
lines(csp[2,], col="purple", lty=3)
plot((yt3[3,]), col="green" , type="l", main="L3")
lines(yhatlower50e1[3,],col="blue")
lines(yhatupper50e1[3,], col="blue")
lines(yhatlower90e1[3,],col="red")
lines(yhatupper90e1[3,], col="red")
lines(csp.hatdp1[3,])
lines(csp[3,], col="purple", lty=3)
plot((yt3[4,]), col="green" , type="l", main="L4")
lines(yhatlower50e1[4,],col="blue")
lines(yhatupper50e1[4,], col="blue")
lines(yhatlower90e1[4,],col="red")
lines(yhatupper90e1[4,], col="red")
lines(csp.hatdp1[4,])
lines(csp[4,], col="purple", lty=3)
plot((yt3[5,]), col="green"  , type="l", main="L5")
lines(yhatlower50e1[5,],col="blue")
lines(yhatupper50e1[5,], col="blue")
lines(yhatlower90e1[5,],col="red")
lines(yhatupper90e1[5,], col="red")
lines(csp.hatdp1[5,])
lines(csp[5,], col="purple", lty=3)

quartz()
par(mfrow=c(3,2))
plot(csp[1,],col="red")
points(csp.hatdpns[1,],col="blue")
points((yt3[1,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)
plot(csp[2,],col="red")
points(csp.hatdpns[2,],col="blue")
points((yt3[2,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
plot(csp[3,],col="red")
points(csp.hatdpns[3,],col="blue")
points((yt3[3,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct3","chat3","yt3"),bg="white", lty=1,cex=0.8)
plot(csp[4,],col="red")
points(csp.hatdpns[4,],col="blue")
points((yt3[4,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct4","chat4","yt4"),bg="white", lty=1,cex=0.8)
plot(csp[5,],col="red")
points(csp.hatdpns[5,],col="blue")
points((yt3[5,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct5","chat5","yt5"),bg="white", lty=1,cex=0.8)
#allsim from random error
jagseddysimns<-mcmcjags.eddysimns<-mcmcsamplesjags.eddysimns<-vector("list")
for(i in 41:50){
  jagseddysimns[[i]]<- jags.model("modeljagseddydpns.jag",data=list("logyt"=(yteddysim[i,,]), N=99,nloc=5,"h"=h,"muet"=rep(0,100*5),
                                                                  "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=5,"I"=diag(5)))
  update(jagseddysimns[[i]], 1000)
  mcmcjags.eddysimns[[i]]<-jags.samples(jagseddysimns[[i]],
                                      c('Dt','Gprime','sigma','csp','yhat'),
                                      1000)
  # mcmcsamplesjags.eddysimns[[i]]<-coda.samples(jagseddysimns[[i]],
  #                                            c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
  #                                            1000)
}
#save.image(file="temp1to30.RData")
#rm(list=ls())
#load(file="temp.RData")
m1.mcmc.eddysimns<-m1.mat.eddysimns<-m1.dat.eddysimns<-vector("list")
yhateddysimns<-yhatlower50eddysimns<-yhatupper50eddysimns<-yhatlower90eddysimns<-yhatupper90eddysimns<-
  chateddysimns<-chatlower50eddysimns<-chatupper50eddysimns<-chatlower90eddysimns<-chatupper90eddysimns<-array(0,dim = c(49, 100, 5))
for(i in 1:49){
  #m1.mcmc.eddysimns[[i]]<-(as.mcmc(mcmcsamplesjags.eddysimns[[i]]))
  #m1.mat.eddysimns[[i]]<-as.matrix(mcmcsamplesjags.eddysimns[[i]])
  #m1.dat.eddysimns[[i]]<-as.data.frame(m1.mat.eddysimns[[i]])
  for(j in 1:100){
    yhateddysimns[i,j,]<- apply(((mcmcjags.eddysimns[[i]]$yhat[,j,,])), 1, mean)
    yhatlower50eddysimns[i,j,]<- apply(((mcmcjags.eddysimns[[i]]$yhat[,j,,])), 1, function(x) quantile(x,0.25))
    yhatupper50eddysimns[i,j,]<-apply((mcmcjags.eddysimns[[i]]$yhat[,j,,]), 1, function(x) quantile(x,0.75))
    yhatlower90eddysimns[i,j,]<- apply(mcmcjags.eddysimns[[i]]$yhat[,j,,], 1, function(x) quantile(x,0.05))
    yhatupper90eddysimns[i,j,]<-apply(mcmcjags.eddysimns[[i]]$yhat[,j,,], 1, function(x) quantile(x,0.95))
    chateddysimns[i,j,]<- apply(mcmcjags.eddysimns[[i]]$csp[,j,,], 1, mean)
    chatlower50eddysimns[i,j,]<- apply(mcmcjags.eddysimns[[i]]$csp[,j,,], 1, function(x) quantile(x,0.25))
    chatupper50eddysimns[i,j,]<-apply(mcmcjags.eddysimns[[i]]$csp[,j,,], 1, function(x) quantile(x,0.75))
    chatlower90eddysimns[i,j,]<- apply(mcmcjags.eddysimns[[i]]$csp[,j,,], 1, function(x) quantile(x,0.05))
    chatupper90eddysimns[i,j,]<-apply(mcmcjags.eddysimns[[i]]$csp[,j,,], 1, function(x) quantile(x,0.95))
  }
}
cov50eddysimns<-cov90eddysimns<-mseeddysimns<-crpseddysimns<-array(0, dim=c(49,100,5))
for(i in 1:49){
  for (k in 1:5){
    cov50eddysimns[i,,k]<-(ifelse(chatlower50eddysimns[i,,k]<=(csp[k,])&(csp[k,])<=chatupper50eddysimns[i,,k],1,0))
    cov90eddysimns[i,,k]<-(ifelse(chatlower90eddysimns[i,,k]<=(csp[k,])&(csp[k,])<=chatupper90eddysimns[i,,k],1,0))
    mseeddysimns[i,,k]<-sum((chateddysimns[i,,k]-(csp[,k]))^2)/100
    crpseddysimns[i,,k]<-mean(crps_sample(y = c(yteddysim[i,k,]), dat =mcmcjags.eddysimns[[i]]$yhat[k,,,]))
  }
  
}
mean((apply(cov50eddysimns,1 ,function(x) sum(x)/500)))
mean((apply(cov90eddysimns,1, function(x) sum(x)/500)))
mean((apply(mseeddysimns,1,mean)))
mean((apply(crpseddysimns,1,mean)))
gprime.post.eddysim<-Dt.post.eddysim<-vector("list")
for(i in 1:49){
  Dt.post.eddysim[[i]]<-mcmcjags.eddysimns[[i]]$Dt
  gprime.post.eddysim[[i]]<-mcmcjags.eddysimns[[i]]$Gprime
}
quantile(unlist(Dt.post.eddysim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(gprime.post.eddysim),c(0.025,0.25,0.5,0.75,0.975))
plot(chateddysimns[i,,5])
points(c(yteddysim[i,5,]),col="red")
points(csp[5,],col="blue")
cdfdpe11sim<-cdfdpe21sim<-cdfdpe31sim<-cdfdpe41sim<-cdfdpe51sim<-vector("list")
for(i in 1:49){
  cdfdpe11sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[1,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,1,])))
  cdfdpe21sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[2,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,2,])))
  cdfdpe31sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[3,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,3,])))
  cdfdpe41sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[4,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,4,])))
  cdfdpe51sim[[i]]<-apply(mcmcjags.eddysim[[i]]$yhat[5,,,1],1,function(x) empirical_cdf(x,(yteddysim[i,5,])))
}

cdf11sim<-cdf21sim<-cdf31sim<-cdf41sim<-cdf51sim<-matrix(0,50,100)
for(i in 1:49){
  for(j in 1:100){
    cdf11sim[i,j]<-mean(sapply(cdfdpe11sim[[i]], function(X) X$CDF[j]))
    cdf21sim[i,j]<-mean(sapply(cdfdpe21sim[[i]], function(X) X$CDF[j]))
    cdf31sim[i,j]<-mean(sapply(cdfdpe31sim[[i]], function(X) X$CDF[j]))
    cdf41sim[i,j]<-mean(sapply(cdfdpe41sim[[i]], function(X) X$CDF[j]))
    cdf51sim[i,j]<-mean(sapply(cdfdpe51sim[[i]], function(X) X$CDF[j]))
  }}
hist((cdf11sim[1,]),breaks = 15, main="L1")
hist((cdf21sim[20,]),breaks = 15, main="L2")
hist((cdf31sim[40,]),breaks = 15, main="L3")
hist((cdf41sim[49,]),breaks = 15, main="L4")
hist((cdf51sim[11,]),breaks = 15, main="L5")
ghatdpl11sim<-ghatdpl21sim<-ghatdpl31sim<-ghatdpl41sim<-ghatdpl51sim<-diffel11sim<-diffel21sim<-diffel31sim<-diffel41sim<-diffel51sim<-vector("list")
for(i in 1:49){
  ghatdpl11sim[[i]]<-empirical_cdf(yteddysim[i,1,],yteddysim[i,1,])
  ghatdpl21sim[[i]]<-empirical_cdf(yteddysim[i,2,],yteddysim[i,2,])
  ghatdpl31sim[[i]]<-empirical_cdf(yteddysim[i,3,],yteddysim[i,3,])
  ghatdpl41sim[[i]]<-empirical_cdf(yteddysim[i,4,],yteddysim[i,4,])
  ghatdpl51sim[[i]]<-empirical_cdf(yteddysim[i,5,],yteddysim[i,5,])
  diffel11sim[[i]]<-((cdf11sim[i,])-ghatdpl11sim[[i]]$CDF)
  diffel21sim[[i]]<-((cdf21sim[i,])-ghatdpl21sim[[i]]$CDF)
  diffel31sim[[i]]<-((cdf31sim[i,])-ghatdpl31sim[[i]]$CDF)
  diffel41sim[[i]]<-((cdf41sim[i,])-ghatdpl41sim[[i]]$CDF)
  diffel51sim[[i]]<-((cdf51sim[i,])-ghatdpl51sim[[i]]$CDF)
  
}
par(mfrow=c(1,5))
plot(unlist(diffel11sim[[1]]),main="L1",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel11sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel21sim[[1]]),main="L2",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel21sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel31sim[[1]]),main="L3",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel31sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel41sim[[1]]),main="L4",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel41sim[[i]]),col=i, lwd=2)
}
abline(h=0)
plot(unlist(diffel51sim[[1]]),main="L5",type="l",lwd=2,ylim=c(-0.1,0.1),ylab="F.forecast-F.observed")
for(i in 2:50){
  lines(unlist(diffel51sim[[i]]),col=i, lwd=2)
}
abline(h=0)
#model2
jags.eddydpns2<-jags.model("modeljagseddydpns.jag",data=list("logyt"=(yt32[1:5,]), N=99,nloc=5,"h"=h,"muet"=rep(0,100*5),
                                                            "pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=5,"I"=diag(5)))
update(jags.eddydpns2, 1000)
mcmcjags.eddydpns2<-jags.samples(jags.eddydpns2,
                                c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                1000)
mcmcsamplesjags.eddydpns2<-coda.samples(jags.eddydpns2,
                                       c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                       1000)

m1.mcmc.eddydpns2<-(as.mcmc(mcmcsamplesjags.eddydpns2))
m1.mat.eddydpns2<-as.matrix(mcmcsamplesjags.eddydpns2)
m1.dat.eddydpns2<-as.data.frame(m1.mat.eddydpns2)
Dt.post.eddydpns2<-m1.dat.eddydpns2$Dt
gprime.post.eddydpns2<-m1.dat.eddydpns2$Gprime
aomega.post.eddydpns2<-m1.dat.eddydpns2$aomega
bomega.post.eddydpns2<-m1.dat.eddydpns2$bomega
sigma.post.eddydpns2<-m1.dat.eddydpns2$sigma
phi.post.eddydpns2<-m1.dat.eddydpns2$phi
quantile(Dt.post.eddydpns2,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddydpns2,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddydpns2,c(0.025,0.25,0.5,0.75,0.975))
mean(crps_sample(y = yt32[1,], dat = mcmcjags.eddydpns2$yhat[1,,,1]))+
  mean(crps_sample(y = yt32[2,], dat =mcmcjags.eddydpns2$yhat[2,,,1]))+
  mean(crps_sample(y = yt32[3,], dat =mcmcjags.eddydpns2$yhat[3,,,1]))+
  mean(crps_sample(y = yt32[4,], dat =mcmcjags.eddydpns2$yhat[4,,,1]))+
  mean(crps_sample(y = yt32[5,], dat =mcmcjags.eddydpns2$yhat[5,,,1]))
csp.hatdpns2<-yhatdpns2<-matrix(0, nrow=5, ncol=100)
c.hatlower50ens2<-yhatlower50ens2<-matrix(0, nrow=5, ncol=100)
c.hatupper50ens2<-yhatupper50ens2<-matrix(0, nrow=5, ncol=100)
c.hatlower90ens2<-yhatlower90ens2<-matrix(0, nrow=5, ncol=100)
c.hatupper90ens2<-yhatupper90ens2<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  yhatdpns2[j,] <- apply( mcmcjags.eddydpns2$yhat[j,,,1],1,mean)
  yhatlower50ens2[j,]<- apply( mcmcjags.eddydpns2$yhat[j,,,1],1,function(x) quantile(x,0.25))
  yhatupper50ens2[j,]<- apply( mcmcjags.eddydpns2$yhat[j,,,1],1,function(x) quantile(x,0.75))
  yhatlower90ens2[j,]<- apply( mcmcjags.eddydpns2$yhat[j,,,1],1,function(x) quantile(x,0.05))
  yhatupper90ens2[j,]<- apply( mcmcjags.eddydpns2$yhat[j,,,1],1,function(x) quantile(x,0.95))
  for(i in 1:100){
    csp.hatdpns2[j,i] <- mean(mcmcsamplesjags.eddydpns2[[1]][,2+(5*i+1-(6-j))])
    c.hatlower50ens2[j,i]<- quantile(mcmcsamplesjags.eddydpns2[[1]][,2+(5*i+1-(6-j))],0.25)
    c.hatupper50ens2[j,i]<- quantile(mcmcsamplesjags.eddydpns2[[1]][,2+(5*i+1-(6-j))],0.75)
    c.hatlower90ens2[j,i]<- quantile(mcmcsamplesjags.eddydpns2[[1]][,2+(5*i+1-(6-j))],0.05)
    c.hatupper90ens2[j,i]<- quantile(mcmcsamplesjags.eddydpns2[[1]][,2+(5*i+1-(6-j))],0.95)
  }
}
mean50ns2<-mean90ns2<-NULL
for(i in 1:5){
  mean50ns2[i]<-mean(yhatupper50ens2[i,]-yhatlower50ens2[i,])
  mean90ns2[i]<-mean(yhatupper90ens2[i,]-yhatlower90ens2[i,])
}
mean(mean50ns2)
mean(mean90ns2)
cov502b<-cov902b<-matrix(0,100,5)
for(i in 1:length(yt3[1,])){
  for(j in 1:5){
    cov502b[i,j]<-(ifelse(yhatlower50ens2[j,i]<=(yt32[j,i])&(yt32[j,i])<=yhatupper50ens2[j,i],1,0))
    cov902b[i,j]<-(ifelse(yhatlower90ens2[j,i]<=(yt32[j,i])&(yt32[j,i])<=yhatupper90ens2[j,i],1,0))
  }}
sum(cov502b)/500
sum(cov902b)/500
cdfdpe1ns2<-apply(mcmcjags.eddydpns2$yhat[1,,,1],1,function(x) empirical_cdf(x,yt32[1,]))
cdfdpe2ns2<-apply(mcmcjags.eddydpns2$yhat[2,,,1],1,function(x) empirical_cdf(x,yt32[2,]))
cdfdpe3ns2<-apply(mcmcjags.eddydpns2$yhat[3,,,1],1,function(x) empirical_cdf(x,yt32[3,]))
cdfdpe4ns2<-apply(mcmcjags.eddydpns2$yhat[4,,,1],1,function(x) empirical_cdf(x,yt32[4,]))
cdfdpe5ns2<-apply(mcmcjags.eddydpns2$yhat[5,,,1],1,function(x) empirical_cdf(x,yt32[5,]))
cdf1ns2<-cdf2ns2<-cdf3ns2<-cdf4ns2<-cdf5ns2<-NULL
for(i in 1:100){
  cdf1ns2[i]<-mean(sapply(cdfdpe1ns2, function(X) X$CDF[i]))
  cdf2ns2[i]<-mean(sapply(cdfdpe2ns2, function(X) X$CDF[i]))
  cdf3ns2[i]<-mean(sapply(cdfdpe3ns2, function(X) X$CDF[i]))
  cdf4ns2[i]<-mean(sapply(cdfdpe4ns2, function(X) X$CDF[i]))
  cdf5ns2[i]<-mean(sapply(cdfdpe5ns2, function(X) X$CDF[i]))
}

quartz()
par(mfrow=c(2,5))
hist((cdf1ns2),breaks = 15, main="L1")
hist((cdf2ns2),breaks = 15, main="L2")
hist((cdf3ns2),breaks = 15, main="L3")
hist((cdf4ns2),breaks = 15, main="L4")
hist((cdf5ns2),breaks = 15, main="L5")
#additive
hist((cdf12),breaks = 15, main="L1")
hist((cdf22),breaks = 15, main="L2")
hist((cdf32),breaks = 15, main="L3")
hist((cdf42),breaks = 15, main="L4")
hist((cdf52),breaks = 15, main="L5")
diffel1ns2<-((cdf1ns2)-ghatdpl12$CDF)
diffel2ns2<-((cdf2ns2)-ghatdpl22$CDF)
diffel3ns2<-((cdf3ns2)-ghatdpl32$CDF)
diffel4ns2<-((cdf4ns2)-ghatdpl42$CDF)
diffel5ns2<-((cdf5ns2)-ghatdpl52$CDF)
quartz()
par(mfrow=c(2,5))
plot(diffel1ns2,main="L1")
abline(h=0)
plot(diffel2ns2, main="L2")
abline(h=0)
plot(diffel3ns2,main="L3")
abline(h=0)
plot(diffel4ns2, main="L4")
abline(h=0)
plot(diffel5ns2,main="L5")
abline(h=0)


plot(diffel12,main="L1")
abline(h=0)
plot(diffel22, main="L2")
abline(h=0)
plot(diffel32,main="L3")
abline(h=0)
plot(diffel42, main="L4")
abline(h=0)
plot(diffel52,main="L5")
abline(h=0)

quartz()
par(mfrow=c(2,5))
plot((yt32[1,]), col="green", type="l", main="Simulation2")
lines(  csp.hatdpns2[1,], main="L1")
lines(yhatlower50ens2[1,],col="blue")
lines(yhatupper50ens2[1,], col="blue")
lines(yhatlower90ens2[1,],col="red")
lines(yhatupper90ens2[1,], col="red")
lines(csp[1,], col="purple", lty=3)
plot((yt32[2,]), col="green", type="l")
lines(  csp.hatdpns2[2,], main="L2")
lines(yhatlower50ens2[2,],col="blue")
lines(yhatupper50ens2[2,], col="blue")
lines(yhatlower90ens2[2,],col="red")
lines(yhatupper90ens2[2,], col="red")
lines(csp[2,], col="purple", lty=3)
plot( (yt32[3,]), col="green" , type="l", main="L3")
lines(yhatlower50ens2[3,],col="blue")
lines(yhatupper50ens2[3,], col="blue")
lines(yhatlower90ens2[3,],col="red")
lines(yhatupper90ens2[3,], col="red")
lines(csp.hatdpns2[3,])
lines(csp[3,], col="purple", lty=3)
plot( (yt32[4,]), col="green", type="l", main="L4")
lines(yhatlower50ens2[4,],col="blue")
lines(yhatupper50ens2[4,], col="blue")
lines(yhatlower90ens2[4,],col="red")
lines(yhatupper90ens2[4,], col="red")
lines( csp.hatdpns2[4,])
lines(csp[4,], col="purple", lty=3)
plot((yt32[5,]), col="green"  , type="l", main="L5")
lines(yhatlower50ens2[5,],col="blue")
lines(yhatupper50ens2[5,], col="blue")
lines(yhatlower90ens2[5,],col="red")
lines(yhatupper90ens2[5,], col="red")
lines(csp.hatdpns2[5,])
lines(csp[5,], col="purple", lty=3)

plot((yt32[1,]), col="green", type="l")
lines(  csp.hatdp2[1,], main="L1")
lines(yhatlower50e2[1,],col="blue")
lines(yhatupper50e2[1,], col="blue")
lines(yhatlower90e2[1,],col="red")
lines(yhatupper90e2[1,], col="red")
lines(csp[1,], col="purple", lty=3)
plot((yt32[2,]), col="green", type="l")
lines(  csp.hatdp2[2,], main="L2")
lines(yhatlower50e2[2,],col="blue")
lines(yhatupper50e2[2,], col="blue")
lines(yhatlower90e2[2,],col="red")
lines(yhatupper90e2[2,], col="red")
lines(csp[2,], col="purple", lty=3)
plot( (yt32[3,]), col="green" , type="l", main="L3")
lines(yhatlower50e2[3,],col="blue")
lines(yhatupper50e2[3,], col="blue")
lines(yhatlower90e2[3,],col="red")
lines(yhatupper90e2[3,], col="red")
lines(csp.hatdp2[3,])
lines(csp[3,], col="purple", lty=3)
plot( (yt32[4,]), col="green", type="l", main="L4")
lines(yhatlower50e2[4,],col="blue")
lines(yhatupper50e2[4,], col="blue")
lines(yhatlower90e2[4,],col="red")
lines(yhatupper90e2[4,], col="red")
lines( csp.hatdp2[4,])
lines(csp[4,], col="purple", lty=3)
plot((yt32[5,]), col="green"  , type="l", main="L5")
lines(yhatlower50e2[5,],col="blue")
lines(yhatupper50e2[5,], col="blue")
lines(yhatlower90e2[5,],col="red")
lines(yhatupper90e2[5,], col="red")
lines(csp.hatdp2[5,])
lines(csp[5,], col="purple", lty=3)
#model3
jags.eddydpns3<-jags.model("modeljagseddydpns.jag",data=list("logyt"=(yt33[1:5,]), N=99,nloc=5,"h"=h,"muet"=rep(0,100*5),
"pi"=pi, "dist"=(distances), "matdist"=(matdist),"m"=5,"I"=diag(5)))
update(jags.eddydpns3, 1000)
mcmcjags.eddydpns3<-jags.samples(jags.eddydpns3,
                                c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                1000)
mcmcsamplesjags.eddydpns3<-coda.samples(jags.eddydpns3,
                                       c('Dt','Gprime','sigma','csp','loglik','et2','yhat'),
                                       1000)
m1.mcmc.eddydpns3<-(as.mcmc(mcmcsamplesjags.eddydpns3))
m1.mat.eddydpns3<-as.matrix(mcmcsamplesjags.eddydpns3)
m1.dat.eddydpns3<-as.data.frame(m1.mat.eddydpns3)
Dt.post.eddydpns3<-m1.dat.eddydpns3$Dt
gprime.post.eddydpns3<-m1.dat.eddydpns3$Gprime
aomega.post.eddydpns3<-m1.dat.eddydpns3$aomega
bomega.post.eddydpns3<-m1.dat.eddydpns3$bomega
sigma.post.eddydpns3<-m1.dat.eddydpns3$sigma
phi.post.eddydpns3<-m1.dat.eddydpns3$phi
quantile(Dt.post.eddydpns3,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddydpns3,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddydpns3,c(0.025,0.25,0.5,0.75,0.975))
mean(crps_sample(y = yt33[1,], dat = mcmcjags.eddydpns3$yhat[1,,,1]))+
  mean(crps_sample(y = yt33[2,], dat =mcmcjags.eddydpns3$yhat[2,,,1]))+
  mean(crps_sample(y = yt33[3,], dat =mcmcjags.eddydpns3$yhat[3,,,1]))+
  mean(crps_sample(y = yt33[4,], dat =mcmcjags.eddydpns3$yhat[4,,,1]))+
  mean(crps_sample(y = yt33[5,], dat =mcmcjags.eddydpns3$yhat[5,,,1]))
csp.hatdpns3<-yhatdpns3<-matrix(0, nrow=5, ncol=100)
c.hatlower50ens3<-yhatlower50ens3<-matrix(0, nrow=5, ncol=100)
c.hatupper50ens3<-yhatupper50ens3<-matrix(0, nrow=5, ncol=100)
c.hatlower90ens3<-yhatlower90ens3<-matrix(0, nrow=5, ncol=100)
c.hatupper90ens3<-yhatupper90ens3<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  yhatdpns3[j,] <- apply( mcmcjags.eddydpns3$yhat[j,,,1],1,mean)
  yhatlower50ens3[j,]<- apply( mcmcjags.eddydpns3$yhat[j,,,1],1,function(x) quantile(x,0.25))
  yhatupper50ens3[j,]<- apply( mcmcjags.eddydpns3$yhat[j,,,1],1,function(x) quantile(x,0.75))
  yhatlower90ens3[j,]<- apply( mcmcjags.eddydpns3$yhat[j,,,1],1,function(x) quantile(x,0.05))
  yhatupper90ens3[j,]<- apply( mcmcjags.eddydpns3$yhat[j,,,1],1,function(x) quantile(x,0.95))
  for(i in 1:100){
    csp.hatdpns3[j,i] <- mean(mcmcsamplesjags.eddydpns3[[1]][,2+(5*i+1-(6-j))])
    c.hatlower50ens3[j,i]<- quantile(mcmcsamplesjags.eddydpns3[[1]][,2+(5*i+1-(6-j))],0.25)
    c.hatupper50ens3[j,i]<- quantile(mcmcsamplesjags.eddydpns3[[1]][,2+(5*i+1-(6-j))],0.75)
    c.hatlower90ens3[j,i]<- quantile(mcmcsamplesjags.eddydpns3[[1]][,2+(5*i+1-(6-j))],0.05)
    c.hatupper90ens3[j,i]<- quantile(mcmcsamplesjags.eddydpns3[[1]][,2+(5*i+1-(6-j))],0.95)
  }
}

mean50ns3<-mean90ns3<-NULL
for(i in 1:5){
  mean50ns3[i]<-mean(yhatupper50ens3[i,]-yhatlower50ens3[i,])
  mean90ns3[i]<-mean(yhatupper90ens3[i,]-yhatlower90ens3[i,])
}
mean(mean50ns3)
mean(mean90ns3)
cov503b<-cov903b<-matrix(0,100,5)
for(i in 1:length(yt3[1,])){
  for(j in 1:5){
    cov503b[i,j]<-(ifelse(yhatlower50ens3[j,i]<=(yt33[j,i])&(yt33[j,i])<=yhatupper50ens3[j,i],1,0))
    cov903b[i,j]<-(ifelse(yhatlower90ens3[j,i]<=(yt33[j,i])&(yt33[j,i])<=yhatupper90ens3[j,i],1,0))
  }}
sum(cov503b)/500
sum(cov903b)/500
cdfdpe1ns3<-apply(mcmcjags.eddydpns3$yhat[1,,,1],1,function(x) empirical_cdf(x,yt33[1,]))
cdfdpe2ns3<-apply(mcmcjags.eddydpns3$yhat[2,,,1],1,function(x) empirical_cdf(x,yt33[2,]))
cdfdpe3ns3<-apply(mcmcjags.eddydpns3$yhat[3,,,1],1,function(x) empirical_cdf(x,yt33[3,]))
cdfdpe4ns3<-apply(mcmcjags.eddydpns3$yhat[4,,,1],1,function(x) empirical_cdf(x,yt33[4,]))
cdfdpe5ns3<-apply(mcmcjags.eddydpns3$yhat[5,,,1],1,function(x) empirical_cdf(x,yt33[5,]))
cdf1ns3<-cdf2ns3<-cdf3ns3<-cdf4ns3<-cdf5ns3<-NULL
for(i in 1:100){
  cdf1ns3[i]<-mean(sapply(cdfdpe1ns3, function(X) X$CDF[i]))
  cdf2ns3[i]<-mean(sapply(cdfdpe2ns3, function(X) X$CDF[i]))
  cdf3ns3[i]<-mean(sapply(cdfdpe3ns3, function(X) X$CDF[i]))
  cdf4ns3[i]<-mean(sapply(cdfdpe4ns3, function(X) X$CDF[i]))
  cdf5ns3[i]<-mean(sapply(cdfdpe5ns3, function(X) X$CDF[i]))
}

quartz()
par(mfrow=c(2,5))
hist((cdf1ns3),breaks = 15, main="L1")
hist((cdf2ns3),breaks = 15, main="L2")
hist((cdf3ns3),breaks = 15, main="L3")
hist((cdf4ns3),breaks = 15, main="L4")
hist((cdf5ns3),breaks = 15, main="L5")
#additive
hist((cdf13),breaks = 15, main="L1")
hist((cdf23),breaks = 15, main="L2")
hist((cdf33),breaks = 15, main="L3")
hist((cdf43),breaks = 15, main="L4")
hist((cdf53),breaks = 15, main="L5")
diffel1ns3<-((cdf1ns3)-ghatdpl13$CDF)
diffel2ns3<-((cdf2ns3)-ghatdpl23$CDF)
diffel3ns3<-((cdf3ns3)-ghatdpl33$CDF)
diffel4ns3<-((cdf4ns3)-ghatdpl43$CDF)
diffel5ns3<-((cdf5ns3)-ghatdpl53$CDF)
quartz()
par(mfrow=c(2,5))
plot(diffel1ns3,main="L1")
abline(h=0)
plot(diffel2ns3, main="L2")
abline(h=0)
plot(diffel3ns3,main="L3")
abline(h=0)
plot(diffel4ns3, main="L4")
abline(h=0)
plot(diffel5ns3,main="L5")
abline(h=0)


plot(diffel13,main="L1")
abline(h=0)
plot(diffel23, main="L2")
abline(h=0)
plot(diffel33,main="L3")
abline(h=0)
plot(diffel43, main="L4")
abline(h=0)
plot(diffel53,main="L5")
abline(h=0)

quartz()
par(mfrow=c(2,5))
plot((yt33[1,]), col="green", type="l", main="Simulation3")
lines(  csp.hatdpns3[1,], main="L1")
lines(yhatlower50ens3[1,],col="blue")
lines(yhatupper50ens3[1,], col="blue")
lines(yhatlower90ens3[1,],col="red")
lines(yhatupper90ens3[1,], col="red")
lines(csp[1,], col="purple", lty=3)
plot((yt33[2,]), col="green", type="l")
lines(  csp.hatdpns3[2,], main="L2")
lines(yhatlower50ens3[2,],col="blue")
lines(yhatupper50ens3[2,], col="blue")
lines(yhatlower90ens3[2,],col="red")
lines(yhatupper90ens3[2,], col="red")
lines(csp[2,], col="purple", lty=3)
plot( (yt33[3,]), col="green" , type="l", main="L3")
lines(yhatlower50ens3[3,],col="blue")
lines(yhatupper50ens3[3,], col="blue")
lines(yhatlower90ens3[3,],col="red")
lines(yhatupper90ens3[3,], col="red")
lines(csp.hatdpns3[3,])
lines(csp[3,], col="purple", lty=3)
plot( (yt33[4,]), col="green", type="l", main="L4")
lines(yhatlower50ens3[4,],col="blue")
lines(yhatupper50ens3[4,], col="blue")
lines(yhatlower90ens3[4,],col="red")
lines(yhatupper90ens3[4,], col="red")
lines( csp.hatdpns3[4,])
lines(csp[4,], col="purple", lty=3)
plot((yt33[5,]), col="green"  , type="l", main="L5")
lines(yhatlower50ens3[5,],col="blue")
lines(yhatupper50ens3[5,], col="blue")
lines(yhatlower90ens3[5,],col="red")
lines(yhatupper90ens3[5,], col="red")
lines(csp.hatdpns3[5,])
lines(csp[5,], col="purple", lty=3)

plot((yt33[1,]), col="green", type="l", main="Simulation3")
lines(  csp.hatdp3[1,], main="L1")
lines(yhatlower50e3[1,],col="blue")
lines(yhatupper50e3[1,], col="blue")
lines(yhatlower90e3[1,],col="red")
lines(yhatupper90e3[1,], col="red")
lines(csp[1,], col="purple", lty=3)
plot((yt33[2,]), col="green", type="l")
lines(  csp.hatdp3[2,], main="L2")
lines(yhatlower50e3[2,],col="blue")
lines(yhatupper50e3[2,], col="blue")
lines(yhatlower90e3[2,],col="red")
lines(yhatupper90e3[2,], col="red")
lines(csp[2,], col="purple", lty=3)
plot( (yt33[3,]), col="green" , type="l", main="L3")
lines(yhatlower50e3[3,],col="blue")
lines(yhatupper50e3[3,], col="blue")
lines(yhatlower90e3[3,],col="red")
lines(yhatupper90e3[3,], col="red")
lines(csp.hatdp3[3,])
lines(csp[3,], col="purple", lty=3)
plot( (yt33[4,]), col="green", type="l", main="L4")
lines(yhatlower50e3[4,],col="blue")
lines(yhatupper50e3[4,], col="blue")
lines(yhatlower90e3[4,],col="red")
lines(yhatupper90e3[4,], col="red")
lines( csp.hatdp3[4,])
lines(csp[4,], col="purple", lty=3)
plot((yt33[5,]), col="green"  , type="l", main="L5")
lines(yhatlower50e3[5,],col="blue")
lines(yhatupper50e3[5,], col="blue")
lines(yhatlower90e3[5,],col="red")
lines(yhatupper90e3[5,], col="red")
lines(csp.hatdp3[5,])
lines(csp[5,], col="purple", lty=3)
#data
setwd("/Users/n_a_abdallah/Desktop/spatial/Project2/acetone/")
acetone <- list.files(pattern = ".csv$")
acetonedata<-lapply(acetone, read.csv)
plot(acetonedata[[2]][,2])

# stdata<-data.frame(x=c(0,(1.07-0.41)),y=c((1.07-0.41),0),
#                z=c(rep(1,2*89)))
# stspdata<-SpatialPoints(stdata)
# timesdata<-rep(c(1:89),rep(2,89))
# sttdata<-as.POSIXct(timesdata,origin="2016-01-01")
# allvaldata<-c(exp(c(logacetone00484$loc1.measured,logacetone00484$loc2.measured)))
# stdfdata<-STIDF(stspdata,sttdata,data.frame(allvaldata=(allvaldata)))
# head(stdfdata@data)
# quartz()
# stplot(stdfdata)
# 
# varstdata<-variogramST(allvaldata~1,data=stdfdata)
# head(varstdata)
# 
# quartz()
# plot(varstdata,map=T)
# quartz()
# plot(varstdata,map=F)
# quartz()
# plot(varstdata,wireframe=T)

matdist<-c(0.41,1.07)
distance<-matrix(c(0,(1.07-0.41),(1.07-0.41),0),nrow=2, ncol=2, byrow=TRUE)

#k=0.00484
h=1
acetone18<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/acetone/acetone18.csv")
N<-length(acetone18[,2])-1
logacetone00484<-matrix(0,N+1,2)
logacetone00484<-log(acetone18[,3:4]*(58)/(24.45)+1)
cat("model{
    for (i in 1:N){
    omega[i+1]~dnorm(aomega2[z1[i+1]],bomega[z2[i+1]])T(-(csp[1,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[1])^2/(4*Dt*i)))),)
    for(j in 1:nloc){
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1])
    logacetone00484[i+1,j]~dnorm((csp[j,i+1])+et[j*(i+1),z3[i+1]], tausq)
    yhat[i+1,j]~dnorm((csp[j,i+1])+et[j*(i+1),z3[i+1]], tausq)T(0,)
    }
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    }
    for(j in 1:nloc){
    for(l in 1:nloc){
    for(r in 1:(N+1)){
    for(s in 1:(N+1)){
    k[(l*r)+((l-1)*(N-(r-1))),(j*s)+((j-1)*(N-(s-1)))]<-
    (sigma/(a*((r-s))^2+1))*exp(-b*dist[l,j]/(a*((r-s))^2+1))
    }}}}
    kinv<-inverse(k)
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    for(j in 2:(m-1)){p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.0001
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.0001
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.0001}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    for(l in 1:m){
    aomega[l]~dnorm(0,100)
    bomega[l]~dgamma(0.01,0.01)
    aomega2[l]<-aomega[l]-mean(aomega)
    et[1:(nloc*(N+1)),l]~dmnorm(muet,kinv)
    }
    for(j in 1:nloc){
    yhat[1,j]~dnorm((csp[j,1])+et[j*(1),z3[2]], tausq)T(0,)
    csp[j,1]<-(logacetone00484[1,j])+omega[1]
    }
    omega[1]~dnorm(0,100)
    sigma~dgamma(0.01,0.01)
    b~dunif(0.5,3)
    a~dunif(0.5,3)
    tausq~dgamma(0.01,0.01)
    alpha~dgamma(0.01,0.01)
    Gprime~dunif(1104,1650)
    Dt~dunif(0,1)
    }",file="modeljagseddydatadp.jag")
jags.eddyddp<-jags.model("modeljagseddydatadp.jag",data=list("logacetone00484"=exp(logacetone00484), N=N,nloc=2,"h"=h,"muet"=rep(0,2*(N+1)),
                                                         "pi"=pi, "dist"=distance, "matdist"=matdist,"prec"=diag(10,2),"m"=5))
update(jags.eddyddp, 1000)
mcmcjags.eddyddp<-jags.samples(jags.eddyddp,
                             c('Dt','Gprime','sigma','csp','loglik','yhat'),
                             1000)
mcmcsamplesjags.eddyddp<-coda.samples(jags.eddyddp,
                                    c('Dt','Gprime','sigma','csp','loglik','yhat'),
                                    1000)

m1.mcmc.eddyddp<-(as.mcmc(mcmcsamplesjags.eddyddp))
m1.mat.eddyddp<-as.matrix(mcmcsamplesjags.eddyddp)
m1.dat.eddyddp<-as.data.frame(m1.mat.eddyddp)
Dt.post.eddyddp<-m1.dat.eddyddp$Dt
gprime.post.eddyddp<-m1.dat.eddyddp$Gprime
sigma.post.eddyddp<-m1.dat.eddyddp$sigma
phi.post.eddyddp<-m1.dat.eddyddp$phi
quantile(Dt.post.eddyddp,c(0.025,0.25,0.5,0.75,0.975))/120
quantile(gprime.post.eddyddp,c(0.025,0.25,0.5,0.75,0.975))
yhateddp<-matrix(0, nrow=2, ncol=89)
yhatlower50ddp<-matrix(0, nrow=2, ncol=89)
yhatupper50ddp<-matrix(0, nrow=2, ncol=89)
yhatlower90ddp<-matrix(0, nrow=2, ncol=89)
yhatupper90ddp<-matrix(0, nrow=2, ncol=89)
for(j in 1:2){
  yhateddp[j,] <- apply( (mcmcjags.eddyddp$yhat[,j,,1]),1,mean)
  yhatlower50ddp[j,]<- apply((mcmcjags.eddyddp$yhat[,j,,1]),1,function(x) quantile(x,0.25))
  yhatupper50ddp[j,]<- apply((mcmcjags.eddyddp$yhat[,j,,1]),1,function(x) quantile(x,0.75))
  yhatlower90ddp[j,]<- apply( (mcmcjags.eddyddp$yhat[,j,,1]),1,function(x) quantile(x,0.05))
  yhatupper90ddp[j,]<- apply((mcmcjags.eddyddp$yhat[,j,,1]),1,function(x) quantile(x,0.95))
}
csp.hatddp<-matrix(0, nrow=(N+1), ncol=2)
c.hatlower50ddp<-matrix(0, nrow=(N+1), ncol=2)
c.hatupper50ddp<-matrix(0, nrow=(N+1), ncol=2)
c.hatlower90ddp<-matrix(0,nrow=(N+1), ncol=2)
c.hatupper90ddp<-matrix(0, nrow=(N+1), ncol=2)
for(i in 1:(N+1)){
  csp.hatddp[i,1] <- mean(mcmcsamplesjags.eddyddp[[1]][,2+1+(2*(i-1))])
  csp.hatddp[i,2] <- mean(mcmcsamplesjags.eddyddp[[1]][,2+(2*(i))])
  c.hatlower50ddp[i,1] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+1+(2*(i-1))],0.25)
  c.hatupper50ddp[i,1] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+1+(2*(i-1))],0.75)
  c.hatlower90ddp[i,1] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+1+(2*(i-1))],0.05)
  c.hatupper90ddp[i,1] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+1+(2*(i-1))],0.95)
  c.hatlower50ddp[i,2] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+(2*(i))],0.25)
  c.hatupper50ddp[i,2] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+(2*(i))],0.75)
  c.hatlower90ddp[i,2] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+(2*(i))],0.05)
  c.hatupper90ddp[i,2] <- quantile(mcmcsamplesjags.eddyddp[[1]][,2+(2*(i))],0.95)
  }
cov50ednp<-cov90ednp<-matrix(0,89,2)
for(i in 1:89){
  for(j in 1:2){
    cov50ednp[i,j]<-(ifelse(yhatlower50ddp[j,i]<=(exp(logacetone00484[i,j]))&exp(logacetone00484[i,j])<=yhatupper50ddp[j,i],1,0))
    cov90ednp[i,j]<-(ifelse(yhatlower90ddp[j,i]<=exp(logacetone00484[i,j])&exp(logacetone00484[i,j])<=yhatupper90ddp[j,i],1,0))
  }}
sum(cov50ednp)/(89*2)
sum(cov90ednp)/(89*2)
mean(crps_sample(y = exp(logacetone00484[,1]), dat = (mcmcjags.eddyddp$yhat[,1,,1])))+
  mean(crps_sample(y = exp(logacetone00484[,2]), dat =(mcmcjags.eddyddp$yhat[,2,,1])))
mean50ddp<-mean90ddp<-NULL
cdfdpe1nsd<-cdfdpe2nsd<-vector("list")
for(i in 1:2){
  mean50ddp[i]<-mean(c.hatupper50ddp[i,]-c.hatlower50ddp[i,])
  mean90ddp[i]<-mean(c.hatupper90ddp[i,]-c.hatlower90ddp[i,])
}
cdfdpe1nsd<-apply(mcmcjags.eddyddp$yhat[,1,,1],1,function(x) empirical_cdf(x,exp(logacetone00484[,1])))
cdfdpe2nsd<-apply(mcmcjags.eddyddp$yhat[,2,,1],1,function(x) empirical_cdf(x,exp(logacetone00484[,2])))

cdf1nsd<-cdf2nsd<-NULL
for(i in 1:89){
  cdf1nsd[i]<-mean(sapply(cdfdpe1nsd, function(X) X$CDF[i]))
  cdf2nsd[i]<-mean(sapply(cdfdpe2nsd, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,1))
hist((cdf1nsd),breaks = 15, main="L1")
hist((cdf2nsd),breaks = 15, main="L2")
ghatddpl1<-empirical_cdf(exp(logacetone00484)[,1],exp(logacetone00484)[,1])
ghatddpl2<-empirical_cdf(exp(logacetone00484)[,2],exp(logacetone00484)[,2])
diffel1ddp<-((cdf1nsd)-ghatddpl1$CDF)
diffel2ddp<-(cdf2nsd-ghatddpl2$CDF)
quartz()
par(mfrow=c(2,1))
plot(diffel1ddp,main="L1")
abline(h=0)
plot(diffel2ddp, main="L2")
abline(h=0)
quartz()
par(mfrow=c(2,1))
plot(  csp.hatddp[,1], type="l")
lines(yhatlower50ddp[1,],col="blue")
lines(yhatupper50ddp[1,], col="blue")
lines(yhatlower90ddp[1,],col="red")
lines(yhatupper90ddp[1,], col="red")
lines(exp(logacetone00484)[,1], col="green")
plot(  csp.hatddp[,2], type="l")
lines(yhatlower50ddp[2,],col="blue")
lines(yhatupper50ddp[2,], col="blue")
lines(yhatlower90ddp[2,],col="red")
lines(yhatupper90ddp[2,], col="red")
lines(exp(logacetone00484)[,2], col="green")

quartz()
plot(acetone18[,3]*54/24,csp.hatddp[,1])
abline(0,1)
quartz()
par(mfrow=c(1,2))
plot(acetone18[,3]*54/24,col="green")
points(csp.hatddp[,1],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
plot(acetone18[,4]*54/24,col="green")
points(csp.hatddp[,2],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)

#additive cov
cat("model{
    for (i in 1:N){
    omega[i+1]~dnorm(aomega2[z1[i+1]],bomega[z2[i+1]])T(-(csp[1,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[1])^2/(4*Dt*i)))),)
    for(j in 1:nloc){
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1])
    logacetone00484[i+1,j]~dnorm((csp[j,i+1])+w[i+1,j], tausq)
    yhat[i+1,j]~dnorm((csp[j,i+1])+w[i+1,j], tausq)T(0,)
    }
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    }
    for(j in 1:nloc){
    yhat[1,j]~dnorm((csp[j,1])+w[1,j], tausq)T(0,)
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    }
    }
    kinv<-inverse(k)
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    for(j in 2:(m-1)){p1[j]<-(r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1])+0.0001
    p2[j]<-(r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1])+0.0001
    p3[j]<-(r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1])+0.0001}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha+0.0001)
    r2[l]~dbeta(1,alpha+0.0001)
    r3[l]~dbeta(1,alpha+0.0001)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    for(l in 1:m){
    aomega[l]~dnorm(0,10000)
    bomega[l]~dgamma(0.01,0.01)
    aomega2[l]<-aomega[l]-mean(aomega)
    #vari[1:2,1:2,l]~dwish(I,3)
    et[1:nloc,l]~dmnorm(muet,kinv)
    for(k in 1:nloc){
    etw[k,l]<-et[k,l]*p3[l]
    et2[k,l]<-et[k,l]-mean(et[k,1:m])
    }}
    w[1,1:nloc]<-c(0,0)
    for(t in 1:N){
    w[t+1,1:nloc]<-w[t,1:nloc]+et2[1:nloc,z3[t+1]]
    }
    csp[1:2,1]<-(logacetone00484[1,1:2])+omega[1]
    omega[1]~dnorm(0,100)
    sigma~dgamma(0.01,0.01)
    phi~dgamma(0.01,0.01)
    tausq~dgamma(0.01,0.01)
    alpha~dgamma(0.01,0.01)
    Gprime~dunif(1104,1650)
    Dt~dunif(0,1)
    }",file="modeljagseddydatadp2.jag")
jags.eddyddp2<-jags.model("modeljagseddydatadp2.jag",data=list("logacetone00484"=exp(logacetone00484), N=N,nloc=2,"h"=h,"muet"=rep(0,2),
                                                             "pi"=pi, "dist"=distance, "matdist"=matdist,"prec"=diag(10,2),"m"=5))
update(jags.eddyddp2, 1000)
mcmcjags.eddyddp2<-jags.samples(jags.eddyddp2,
                               c('Dt','Gprime','sigma','csp','loglik','yhat'),
                               1000)
mcmcsamplesjags.eddyddp2<-coda.samples(jags.eddyddp2,
                                      c('Dt','Gprime','sigma','csp','loglik','yhat'),
                                      1000)

m1.mcmc.eddyddp2<-(as.mcmc(mcmcsamplesjags.eddyddp2))
m1.mat.eddyddp2<-as.matrix(mcmcsamplesjags.eddyddp2)
m1.dat.eddyddp2<-as.data.frame(m1.mat.eddyddp2)
Dt.post.eddyddp2<-m1.dat.eddyddp2$Dt
gprime.post.eddyddp2<-m1.dat.eddyddp2$Gprime
sigma.post.eddyddp2<-m1.dat.eddyddp2$sigma
phi.post.eddyddp2<-m1.dat.eddyddp2$phi
quantile(Dt.post.eddyddp2,c(0.025,0.25,0.5,0.75,0.975))/120
quantile(gprime.post.eddyddp2,c(0.025,0.25,0.5,0.75,0.975))

csp.hatddp2<-matrix(0, nrow=(N+1), ncol=2)
c.hatlower50ddp2<-matrix(0, nrow=(N+1), ncol=2)
c.hatupper50ddp2<-matrix(0, nrow=(N+1), ncol=2)
c.hatlower90ddp2<-matrix(0,nrow=(N+1), ncol=2)
c.hatupper90ddp2<-matrix(0, nrow=(N+1), ncol=2)
for(i in 1:(N+1)){
  csp.hatddp2[i,1] <- mean(mcmcsamplesjags.eddyddp2[[1]][,2+1+(2*(i-1))])
  csp.hatddp2[i,2] <- mean(mcmcsamplesjags.eddyddp2[[1]][,2+(2*(i))])
  c.hatlower50ddp2[i,1] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+1+(2*(i-1))],0.25)
  c.hatupper50ddp2[i,1] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+1+(2*(i-1))],0.75)
  c.hatlower90ddp2[i,1] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+1+(2*(i-1))],0.05)
  c.hatupper90ddp2[i,1] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+1+(2*(i-1))],0.95)
  c.hatlower50ddp2[i,2] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+(2*(i))],0.25)
  c.hatupper50ddp2[i,2] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+(2*(i))],0.75)
  c.hatlower90ddp2[i,2] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+(2*(i))],0.05)
  c.hatupper90ddp2[i,2] <- quantile(mcmcsamplesjags.eddyddp2[[1]][,2+(2*(i))],0.95)
}
yhatddp2<-matrix(0, nrow=2, ncol=89)
yhatlower50ddp2<-matrix(0, nrow=2, ncol=89)
yhatupper50ddp2<-matrix(0, nrow=2, ncol=89)
yhatlower90ddp2<-matrix(0, nrow=2, ncol=89)
yhatupper90ddp2<-matrix(0, nrow=2, ncol=89)
for(j in 1:2){
  yhatddp2[j,] <- apply( (mcmcjags.eddyddp2$yhat[,j,,1]),1,mean)
  yhatlower50ddp2[j,]<- apply((mcmcjags.eddyddp2$yhat[,j,,1]),1,function(x) quantile(x,0.25))
  yhatupper50ddp2[j,]<- apply((mcmcjags.eddyddp2$yhat[,j,,1]),1,function(x) quantile(x,0.75))
  yhatlower90ddp2[j,]<- apply( (mcmcjags.eddyddp2$yhat[,j,,1]),1,function(x) quantile(x,0.05))
  yhatupper90ddp2[j,]<- apply((mcmcjags.eddyddp2$yhat[,j,,1]),1,function(x) quantile(x,0.95))
}
cov50ednp2<-cov90ednp2<-matrix(0,89,2)
for(i in 1:89){
  for(j in 1:2){
    cov50ednp2[i,j]<-(ifelse(yhatlower50ddp2[j,i]<=(exp(logacetone00484[i,j]))&exp(logacetone00484[i,j])<=yhatupper50ddp2[j,i],1,0))
    cov90ednp2[i,j]<-(ifelse(yhatlower90ddp2[j,i]<=exp(logacetone00484[i,j])&exp(logacetone00484[i,j])<=yhatupper90ddp2[j,i],1,0))
  }}
sum(cov50ednp2)/(89*2)
sum(cov90ednp2)/(89*2)
mean(crps_sample(y = exp(logacetone00484[,1]), dat = (mcmcjags.eddyddp2$yhat[,1,,1])))+
  mean(crps_sample(y = exp(logacetone00484[,2]), dat =(mcmcjags.eddyddp2$yhat[,2,,1])))
mean50ddp2<-mean90ddp2<-NULL
cdfdpeddp2<-matrix(0,nrow=2, ncol=N+1)
for(i in 1:2){
  mean50ddp2[i]<-mean(c.hatupper50ddp2[i,]-c.hatlower50ddp2[i,])
  mean90ddp2[i]<-mean(c.hatupper90ddp2[i,]-c.hatlower90ddp2[i,])
}
cdfdpe1d<-cdfdpe2d<-vector("list")

cdfdpe1d<-apply(mcmcjags.eddyddp2$yhat[,1,,1],1,function(x) empirical_cdf(x,exp(logacetone00484[,1])))
cdfdpe2d<-apply(mcmcjags.eddyddp2$yhat[,2,,1],1,function(x) empirical_cdf(x,exp(logacetone00484[,2])))

cdf1d<-cdf2d<-NULL
for(i in 1:89){
  cdf1d[i]<-mean(sapply(cdfdpe1d, function(X) X$CDF[i]))
  cdf2d[i]<-mean(sapply(cdfdpe2d, function(X) X$CDF[i]))
}
quartz()
par(mfrow=c(2,2))
hist((cdf1d),breaks = 15, main="L1")
hist((cdf2d),breaks = 15, main="L2")
hist((cdf1nsd),breaks = 15, main="L1 ns")
hist((cdf2nsd),breaks = 15, main="L2 ns")

diffel1ddp2<-(cdf1d-ghatddpl1$CDF)
diffel2ddp2<-(cdf2d-ghatddpl2$CDF)
quartz()
par(mfrow=c(2,2))
plot(diffel1ddp2,main="L1")
abline(h=0)
plot(diffel2ddp2, main="L2")
abline(h=0)
plot(diffel1ddp,main="L1 ns")
abline(h=0)
plot(diffel2ddp, main="L2 ns")
abline(h=0)
quartz()
par(mfrow=c(2,2))
plot(exp(logacetone00484)[,1], col="green"  , type="l",ylab="Acetone concentrations L1", xlab="Time")
lines(yhatlower50ddp2[1,],col="blue")
lines(yhatupper50ddp2[1,], col="blue")
lines(yhatlower90ddp2[1,],col="red")
lines(yhatupper90ddp2[1,], col="red")
lines(csp.hatddp2[,1])
legend("bottomright",col=c("green","blue","red","black"),text.font=4,legend=
         c("measured","50 % PI","90% PI","predicted"),bg="white", lty=1,cex=0.8)
plot(exp(logacetone00484)[,2], col="green"  , type="l",ylab="Acetone concentrations L2", xlab="Time")
lines(yhatlower50ddp2[2,],col="blue")
lines(yhatupper50ddp2[2,], col="blue")
lines(yhatlower90ddp2[2,],col="red")
lines(yhatupper90ddp2[2,], col="red")
lines(csp.hatddp2[,2])
legend("bottomright",col=c("green","blue","red","black"),text.font=4,legend=
         c("measured","50 % PI","90% PI","predicted"),bg="white", lty=1,cex=0.8)
plot( exp(logacetone00484)[,1], col="green", main="NS" , type="l",ylab="Acetone concentrations L1", xlab="Time")
lines(yhatlower50ddp[1,],col="blue")
lines(yhatupper50ddp[1,], col="blue")
lines(yhatlower90ddp[1,],col="red")
lines(yhatupper90ddp[1,], col="red")
lines(csp.hatddp[,1])
legend("bottomright",col=c("green","blue","red","black"),text.font=4,legend=
         c("measured","50 % PI","90% PI","predicted"),bg="white", lty=1,cex=0.8)
plot( exp(logacetone00484)[,2], col="green"  ,main="NS", type="l",ylab="Acetone concentrations L2", xlab="Time")
lines(yhatlower50ddp[2,],col="blue")
lines(yhatupper50ddp[2,], col="blue")
lines(yhatlower90ddp[2,],col="red")
lines(yhatupper90ddp[2,], col="red")
lines(csp.hatddp[,2])
legend("bottomright",col=c("green","blue","red","black"),text.font=4,legend=
         c("measured","50 % PI","90% PI","predicted"),bg="white", lty=1,cex=0.8)
quartz()
plot(acetone18[,3]*54/24,csp.hatddp2[,1])
abline(0,1)
quartz()
par(mfrow=c(1,2))
plot(acetone18[,3]*54/24,col="green")
points(csp.hatddp2[,1],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
plot(acetone18[,4]*54/24,col="green")
points(csp.hatddp2[,2],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
