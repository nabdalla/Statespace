#Bayesian dynamic model
#general model: measurement equation: Y(t)=C(t)+vt vt iid Pv
#transition equation: C(t+1)=(1-Q'+K.lV/V)C(t)+G'/V+wt, wt iid Pw
#Q'=Q+erf QR+Ql, G'=(1-el)G
#Simulation
library(nimble)
library(rjags)
library(plotrix)
library(compositions)
library(mvtnorm)
library(expm)
library(DPpackage)
library(matrixStats)
library(loo)
library(msm)
library(scoringRules)
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
#values of c very huge or very small and hence yt so taking log is problem
for(i in 1:(n-1)){
  c[1]<-1
  c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V
}
plot(c)
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

#jags
cat("model{
    for (i in 1:N){ 
    omega[i+1]~dgamma(aomega,bomega)
    c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+(omega[i+1])
    logyt[i+1]~dnorm(log(c[i+1]), prec)
    logyhat[i+1]~dnorm(log(c[i+1]),prec)
    loglik[i+1]<-log(dnorm(logyt[i+1],log(c[i+1]),prec))
    }
    loglik[1]<-log(dnorm(logyt[1],log(c[1]),prec))
    logyhat[1]~dnorm(log(c[1]),prec)
    omega[1]~dgamma(aomega,bomega)
    c[1]<-1
    logyt[1]~dnorm(log(c[1])+omega[1], prec)
    sigma<-1/prec
    aomega~dunif(0,1)
    bomega~dunif(0,1)
    prec~dgamma(1,0.001)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljags.jag")
jags<-jags.model("modeljags.jag",data=list("logyt"=logyt[1:100], N=99, 
                                           "V"=V,"h"=h))
update(jags, 1000)

mcmcjags<-jags.samples(jags,
                       c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','loglik','logyhat'),
                       2000)
mcmcsamplesjags<-coda.samples(jags,
                              c('Qprime','Gprime','kl','sigma','c','aomega','bomega','loglik','logyhat'),
                              2000)

m1.mcmc<-(as.mcmc(mcmcsamplesjags))
m1.mat<-as.matrix(mcmcsamplesjags)
m1.dat<-as.data.frame(m1.mat)
aomega.post1zone<-m1.dat$aomega
bomega.post1zone<-m1.dat$bomega
qprime.post1zone<-m1.dat$Qprime
gprime.post1zone<-m1.dat$Gprime
kl.post1zone<-m1.dat$kl
sigma.post1zone<-sqrt(m1.dat$sigma)
loglik.post1zone<-t(mcmcjags$loglik[,,1])
dim(loglik.post1zone)
quantile(qprime.post1zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post1zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(kl.post1zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post1zone,c(0.025,0.25,0.5,0.75,0.975))
waic(loglik.post1zone)
mean(crps_sample(y = exp(logyt), dat = exp(mcmcjags$logyhat[,,1])))
c.hat <- apply(mcmcjags$c, 1, mean)
yhatp<-apply(exp(mcmcjags$logyhat[,,1]), 1, mean)
yhatlower50p<-apply(exp(mcmcjags$logyhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50p<-apply(exp(mcmcjags$logyhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90p<-apply(exp(mcmcjags$logyhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90p<-apply(exp(mcmcjags$logyhat[,,1]), 1, function(x) quantile(x,0.975))
#nominal coverage
plot(exp(logyt))
lines(yhatlower50p,col="blue")
lines(yhatupper50p,col="blue")
cov50<-cov90<-NULL
for(i in 1:length(logyt)){
  cov50[i]<-(ifelse(yhatlower50p[i]<=exp(logyt[i])&exp(logyt[i])<=yhatupper50p[i],1,0))
  cov90[i]<-(ifelse(yhatlower90p[i]<=exp(logyt[i])&exp(logyt[i])<=yhatupper90p[i],1,0))
  }
sum(cov50)/100
sum(cov90)/100
#simulations
jags.sim<-mcmcjags.sim<-mcmcsamplesjags.sim<-vector("list")
for(i in 1:100){
  jags.sim[[i]]<-jags.model("modeljags.jag",data=list("logyt"=ytsim1[[i]], N=99, 
                                                      "V"=V,"h"=h))
  update(jags.sim[[i]], 2000)
  mcmcjags.sim[[i]]<-jags.samples(jags.sim[[i]],
                                  c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','loglik','logyhat'),
                                     1000)
  mcmcsamplesjags.sim[[i]]<-coda.samples(jags.sim[[i]],
                                         c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','loglik','logyhat'),
                                            1000)
}
m1.mcmc.sim<-m1.mat.sim<-m1.dat.sim<-vector("list")
for(i in 1:100){
  m1.mcmc.sim[[i]]<-(as.mcmc(mcmcsamplesjags.sim[[i]]))
  m1.mat.sim[[i]]<-as.matrix(mcmcsamplesjags.sim[[i]])
  m1.dat.sim[[i]]<-as.data.frame(m1.mat.sim[[i]])
}
crpssim1par<-yhat90sim1par<-yhat50sim1par<-NULL
yhat.sim1par<-yhatlower50sim1par<-yhatupper50sim1par<-yhatlower90sim1par<-yhatupper90sim1par<-
  chat.sim1par<-chatlower50sim1par<-chatupper50sim1par<-chatlower90sim1par<-chatupper90sim1par<-vector("list")
for(i in 1:100){
  crpssim1par[i]<-mean(crps_sample(y = exp(ytsim1[[i]]), dat = exp((mcmcjags.sim[[i]])$logyhat[,,1])))
  yhat.sim1par[[i]]<-apply((exp((mcmcjags.sim[[i]])$logyhat[,,1])), 1, mean)
  yhatlower50sim1par[[i]]<-apply((exp((mcmcjags.sim[[i]])$logyhat[,,1])), 1, function(x) quantile(x,0.25))
  yhatupper50sim1par[[i]]<-apply((exp((mcmcjags.sim[[i]])$logyhat[,,1])), 1, function(x) quantile(x,0.75))
  yhatlower90sim1par[[i]]<-apply((exp((mcmcjags.sim[[i]])$logyhat[,,1])), 1, function(x) quantile(x,0.05))
  yhatupper90sim1par[[i]]<-apply((exp((mcmcjags.sim[[i]])$logyhat[,,1])), 1, function(x) quantile(x,0.975))
  chat.sim1par[[i]]<-apply((((mcmcjags.sim[[i]])$c)), 1, mean)
  chatlower50sim1par[[i]]<-apply((((mcmcjags.sim[[i]])$c)), 1, function(x) quantile(x,0.25))
  chatupper50sim1par[[i]]<-apply((((mcmcjags.sim[[i]])$c)), 1, function(x) quantile(x,0.75))
  chatlower90sim1par[[i]]<-apply((((mcmcjags.sim[[i]])$c)), 1, function(x) quantile(x,0.05))
  chatupper90sim1par[[i]]<-apply((((mcmcjags.sim[[i]])$c)), 1, function(x) quantile(x,0.975))
  yhat90sim1par[i]<-mean(yhatupper90sim1par[[i]]-yhatlower90sim1par[[i]])
  yhat50sim1par[i]<-mean(yhatupper50sim1par[[i]]-yhatlower50sim1par[[i]])
}
signal2noisepar<-vector("list")
meansignal2noisepar<-NULL
for(i in 1:100){
  signal2noisepar[[i]]<-c2/(sigmasim1[i])
  meansignal2noisepar[i]<-mean(signal2noisepar[[i]])
}
plot(meansignal2noisepar,yhat90sim1par)
plot(meansignal2noisepar,crpssim1par)
mean(crpssim1par)
cov50sim<-cov90sim<-msesim<-vector("list")
for(i in 1:100){
  cov50sim[[i]]<-(ifelse(chatlower50sim1par[[i]]<=c2&c2<=chatupper50sim1par[[i]],1,0))
  cov90sim[[i]]<-(ifelse(chatlower90sim1par[[i]]<=c2&c2<=chatupper90sim1par[[i]],1,0))
  msesim[[i]]<-sum((chat.sim1par[[i]]-c2)^2)/100
  }
mean(unlist(lapply(cov50sim, function(x) sum(x)/100)))
mean(unlist(lapply(cov90sim, function(x) sum(x)/100)))
mean(unlist(msesim))

#D=G+P
zrep<-matrix(0,100,1000)
for(j in 1:1000){
  
    zrep[,j]<-exp(rnorm(100,log(mcmcjags$c[,j,1]),(sigma.post1zone[j])))
}
G1<- (c-apply((zrep),1,mean))^2
P1<-apply((zrep),1,var)
print(D1<-sum(G1)+sum(P1))
cs1<-matrix(0,2000,100)
cs1[,100]<-mcmcsamplesjags[[1]][,(4+100)]
omega.post1zone=rgamma(aomega.post1zone,bomega.post1zone)
N=n-1
for(j in 2:(N)){
  cs1[,(j+N-2*(j-1))]<-(mcmcsamplesjags[[1]][,4+(j+N-2*(j-1)+1)]-h*(gprime.post1zone)/V-(h*(omega.post1zone)))/(1-h*((qprime.post1zone)+
                                                                                    (kl.post1zone)*V)/V)
} 
cs1[,1]<-(mcmcsamplesjags[[1]][,6]-h*(gprime.post1zone)/V+h*(omega.post1zone))/(1-h*((qprime.post1zone)+(kl.post1zone)*V)/V)
cs.hat <- apply(cs1, 2, mean)
cs.hat
plot(cs.hat)
points(c2,col="red")
points(c.hat,col="blue")

quartz()
plot(exp(logyt),col="green")
points(c2,col="blue")
points(c.hat,col="red")
points(cs.hat,col="black")
print(mse1sim<-sum(((c[1:100])-(c.hat[1:100]))^2)/100)
#dic
Lst1<-log(dnorm(logyt,log(c.hat),mean(sigma.post1zone)))
Pst1<-matrix(0,1000,100)
pst1sum<-NULL
for(i in 1:100){
  for(j in 1:1000){
    Pst1[j,i]<-log(dnorm(logyt[i],log(mcmcjags$c[i,j,1]),sigma.post1zone[(j)]))
    pst1sum[j]<-sum(Pst1[j,])
    }
  }

-2*(sum(Lst1[1:100])-(2*(sum(Lst1[1:100])-mean(pst1sum))))
#Pf import from julia output
#thetajulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc22.csv")
thetajuliacorr<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmccorr.csv")
apply(thetajuliacorr, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
plot(thetajuliacorr[,6],type="l")
acf(thetajuliacorr[,2], lag.max=1000)
xjuliacorr<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmccorr.csv")
ypf<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/logytpf.csv")
pred.juliacorr<-apply(xjuliacorr, 2, mean)
# #dic
# Lpf1<-log(dnorm(logyt,log(pred.julia),0.08))
# Ppf1<-matrix(0,1000,100)
# ppf1sum<-NULL
# for(i in 1:100){
#   for(j in 1:1000){
#     Ppf1[j,i]<-log(dnorm(logyt[i],log(xjulia[1000+j,i]),thetajulia$x5[(1000+j)]))
#     ppf1sum[j]<-sum(Ppf1[j,])
#   }
# }
# 
# -2*(sum(Lpf1)-(2*(sum(Lpf1)-mean(ppf1sum))))
zrepph<-matrix(0,100,1000)
for(j in 1:1000){
  for(i in 1:100){
  zrepph[i,j]<-exp(rnorm(1,log(xjuliacorr[j,i]),sqrt(thetajuliacorr$x6[j])))
}}
G1ph<- (exp(ypf[,1])-apply((zrepph),1,mean))^2
P1ph<-apply((zrepph),1,var)
print(D1ph<-sum(G1ph)+sum(P1ph))
print(mse1simpf<-sum(((c[1:100])-(pred.juliacorr[1:100]))^2)/100)
print(mse1simkf<-sum(((c[1:100])-(predkf[1:100]))^2)/100)

loglikpf1<-matrix(0,1000,100)
for(i in 1:100){
  for(j in 1:1000){
  loglikpf1[j,i]<-log(dnorm(logyt[i],log(abs(xjulia[j+2000,i])),thetajulia$x6[j+2000]))
  }}
# for(i in 1:100){
#   lpdhat[i]<-max(loglikpf1[,i])+(log(mean(exp(loglikpf1[,i]-max(loglikpf1[,i])))))
#   pwaic[i]<-(max(loglikpf1[,i])+(log(sum(exp(loglikpf1[,i]-max(loglikpf1[,i]))-
#                                            mean(exp(loglikpf1[,i]-max(loglikpf1[,i])))^2))))/1000
# }
# -2*(sum(lpdhat)-sum(pwaic))
# for(i in 1:100){
#   lpdhat2[i]<-(log(mean(exp(loglik.post1zone[,i]))))
#   pwaic2[i]<-(var(loglik.post1zone[,i]))
#   }
# -2*(sum(lpdhat2)-sum(pwaic2))
# dim(loglikpf1)
waic(loglikpf1[,1:100])
waic(loglik.post1zone)
quartz()
plot(exp(logyt),col="green")
points(c2,col="blue")
points(pred.juliacorr,col="red")
#nonlinear regression
cat("model{
    for (i in 1:N){ 
    c[i+1]<-exp(-i*h*(Qprime+kl*V)/V)*c[1]+1/(-(Qprime+kl*V)/V)*
    (exp(-i*h*(Qprime+kl*V)/V)-1)*Gprime/V
    logyt[i+1]~dnorm(log(c[i+1]), prec)
    loglik[i+1]<-log(dnorm(logyt[i+1],log(c[i+1]),prec))
    }
    loglik[1]<-log(dnorm(logyt[1],log(c[1]),prec))
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), prec)
    sigma<-1/prec
    prec~dgamma(2,0.01)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsnl.jag")
jagsnl<-jags.model("modeljagsnl.jag",data=list("logyt"=logyt[1:100], N=99, "V"=V,"h"=0.01))
update(jagsnl, 1000)
mcmcjagsnl<-jags.samples(jagsnl,
                       c('Qprime','Gprime',"kl",'sigma','c','loglik'),
                       1000)
mcmcsamplesjagsnl<-coda.samples(jagsnl,
                              c('Qprime','Gprime','kl','sigma','c','loglik'),
                              1000)

summary(mcmcsamplesjagsnl)
m1.mcmcnl<-(as.mcmc(mcmcsamplesjagsnl))
m1.matnl<-as.matrix(mcmcsamplesjagsnl)
m1.datnl<-as.data.frame(m1.matnl)
c.hatnl <- apply(mcmcjagsnl$c, 1, mean)
loglik.post1zonel<-t(mcmcjagsnl$loglik[,,1])
dim(loglik.post1zonel)
waic(loglik.post1zonel)
plot(exp(logyt),col="green")
points(c2,col="blue")
points(c.hatnl,col="red")

#all simulation plots
quartz()
#layout(matrix(c(1,1,1,1,1,0,2,2,2,2,2,3,3,3,3,3,0,4,4,4,4,4,0,0,0,5,5,5,5,5,
#                0,0,0), 3, 11, byrow = TRUE))
layout(matrix(c(1,1,1,1,1,0,2,2,2,2,2,0,0,0,3,3,3,3,3,
                     0,0,0), 2, 11, byrow = TRUE))
# plot(exp(logyt),col="green", main="a", ylab="Concentrations", xlab="Time")
# points(c.hat[1:100],col="red")
# lines(c2,col="blue", lwd=3)
plot(exp(logyt),col="green", main="a", ylab="Concentrations", xlab="Time",type="l",lwd=3)
lines(c.hat[1:100],col="red",lwd=3,lty=6)
lines(c2,col="blue", lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6) ,legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

# plot(exp(logyt),col="green", main="GEnKS",ylab="Measurements", xlab="Time")
# points(predenks.julia, col="red")
# points(c2,col="blue")
# legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
#                                                             "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(exp(logyt),col="green", main="b",ylab="Concentrations", xlab="Time")
# points(pred.julia,col="red")
# lines(c2,col="blue",lwd=3)
# plot(exp(logyt),col="green", main="b",ylab="Concentrations", xlab="Time",type="l",lwd=3)
# lines(pred.julia,col="red",lwd=3)
# lines(c2,col="blue",lwd=2)
# legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
#                                                             "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(exp(logyt), col="green",main="c",ylab="Concentrations", xlab="Time")
# points(predkf,col="red")
# lines(c2,col="blue",lwd=3)
plot(exp(logyt), col="green",main="b",ylab="Concentrations", xlab="Time",type="l",lwd=3)
lines(predkf,col="red",lwd=3,lty=6)
lines(c2,col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6) ,legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

# plot(exp(logyt),col="green", main="Smoothing", ylab="Concentrations", xlab="Time")
# points(cs.hat[1:100],col="purple")
# lines(c2,col="blue", lwd=3)
plot(exp(logyt),col="green", main="Smoothing", ylab="Concentrations", xlab="Time",type="l",lwd=3)
lines(cs.hat[1:100],col="purple",lwd=3,lty=6)
lines(c2,col="blue", lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "purple"), lty=c(1,3,6),legend=c("Measurements","True",
                                                            "Smoothed"),pch=15,cex=0.75,box.lwd = 0)
#data one compartment
data1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/Modeling_Master_v5.csv")
dim(data1)
data1zone=matrix(0,nrow=204, ncol=(81*6))
for(i in 1:(81*6)){
  data1zone[,i]=data1[,4*i]
}
data1zone[,204]
write.csv(data1zone,"/Users/n_a_abdallah/Desktop/spatial/Project2/acetone1zone.csv")
#true: V=11.9, Q=0.04-0.07, 0.23-0.27 and 0.47-0.77 m^3/min, 
#G=39.6, 79.1 and 118.7 for l, med, h
#note that each was measured at 6 locations
quartz()
hist(log(data1zone[,2]))
plot(data1zone[,30])
plot(data1zone[,51*4])
points(data1zone[,400])
points(data1zone[,55])
points(data1zone[,6])

cat("model{
    for (i in 1:N){ 
    omega[i+1]~dgamma(aomega,bomega)
    omega2[i+1]~dgamma(aomega2,bomega2)
    c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+omega[i+1]
    logyt[i+1]~dnorm(log((c[i+1])), prec)
    logyhat[i+1]~dnorm(log(c[i+1]),prec)
    loglik[i+1]<-log(dnorm(logyt[i+1],log(c[i+1]),prec))
    }
    logyhat[1]~dnorm(log(c[1]),prec)
    loglik[1]<-log(dnorm(logyt[1],log(c[1]),prec))
    omega[1]~dgamma(aomega,bomega)
    c[1]<-0+omega[1]
    logyt[1]~dnorm(log(c[1]), prec)
    sigma<-1/prec
    aomega~dunif(0.5,3)
    bomega~dunif(0.5,3)
aomega2~dunif(0.5,3)
    bomega2~dunif(0.5,3)
    prec~dgamma(1,0.001)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,0.8)
    }",file="modeljags.jag")
jagslow<-jags.model("modeljags.jag",data=list("logyt"=log(data1zone[1:201,2]), N=200, "V"=11.9,"h"=0.01))
update(jagslow, 1000)
mcmcjagslow<-jags.samples(jagslow,
                       c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','loglik','logyhat'),
                       1000)
mcmcsamplesjagslow<-coda.samples(jagslow,
                              c('Qprime','Gprime','kl','sigma','c','aomega','bomega','loglik','logyhat'),
                              1000)
#low ventilation Q=0.04-0.07 G=43.18
m1.mcmclow<-(as.mcmc(mcmcsamplesjagslow))
m1.matlow<-as.matrix(mcmcsamplesjagslow)
m1.datlow<-as.data.frame(m1.matlow)
aomega.post1zonelow<-m1.datlow$aomega
bomega.post1zonelow<-m1.datlow$bomega
qprime.post1zonelow<-m1.datlow$Qprime
gprime.post1zonelow<-m1.datlow$Gprime
sigma.post1zonelow<-sqrt(m1.datlow$sigma)
kl.post1zonelow<-m1.datlow$kl
loglik.post1zonel<-t(mcmcjagslow$loglik[,,1])
dim(loglik.post1zonel)
waic(loglik.post1zonel[,1:201])
mean(crps_sample(y = data1zone[1:201,2], dat = exp(mcmcjagslow$logyhat[,,1])))
c.hatl <- apply(mcmcjagslow$c, 1, mean)
yhatpl<-apply(exp(mcmcjagslow$logyhat[,,1]), 1, mean)
yhatlower50pl<-apply(exp(mcmcjagslow$logyhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50pl<-apply(exp(mcmcjagslow$logyhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90pl<-apply(exp(mcmcjagslow$logyhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90pl<-apply(exp(mcmcjagslow$logyhat[,,1]), 1, function(x) quantile(x,0.975))
#nominal coverage
plot(data1zone[1:201,2])
lines(yhatlower50pl,col="blue")
lines(yhatupper50pl,col="blue")
cov50l<-cov90l<-NULL
for(i in 1:length(data1zone[1:201,2])){
  cov50l[i]<-(ifelse(yhatlower50pl[i]<=data1zone[i,2]&data1zone[i,2]<=yhatupper50pl[i],1,0))
  cov90l[i]<-(ifelse(yhatlower90pl[i]<=data1zone[i,2]&data1zone[i,2]<=yhatupper90pl[i],1,0))
}
sum(cov50l)/201
sum(cov90l)/201
quantile((gprime.post1zonelow),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonelow),c(0.025,0.25,0.5,0.75,0.975))
quantile((kl.post1zonelow),c(0.025,0.25,0.5,0.75,0.975))
c.hatlow <- apply(mcmcjagslow$c, 1, mean)
yhatlow<-apply(mcmcjagslow$logyt, 1, mean)
zrepl<-matrix(0,201,1000)
for(j in 1:1000){
  set.seed(123)
  zrepl[,j]<-exp(rnorm(201,log(mcmcjagslow$c[,j,1]),
                   (sigma.post1zonelow[j])))
}
G1l<-(data1zone[1:201,2]-apply((zrepl),1,mean))^2
P1l<-apply((zrepl),1,var)
print(D1l<-sum(G1l)+sum(P1l))
cs1l<-matrix(0,1000,201)
cs1l[,201]<-mcmcsamplesjagslow[[1]][,(4+201)]
omega.post1zonelow=rgamma(aomega.post1zonelow,bomega.post1zonelow)
N=201-1
for(j in 2:(N)){
  cs1l[,(j+N-2*(j-1))]<-(mcmcsamplesjagslow[[1]][,4+(j+N-2*(j-1)+1)]-h*(gprime.post1zonelow)/V-(h*(omega.post1zonelow)))/(1-h*((qprime.post1zonelow)+
                                                                                                                       (kl.post1zonelow)*V)/V)
} 
cs1l[,1]<-(mcmcsamplesjagslow[[1]][,6]-h*(gprime.post1zonelow)/V+h*(omega.post1zonelow))/(1-h*((qprime.post1zonelow)+(kl.post1zonelow)*V)/V)
cs.hatl <- apply(cs1l, 2, mean)
cs.hatl
#medium Q Q=0.23-0.27 G=43.18
jagsmed<-jags.model("modeljags.jag",data=list("logyt"=log(data1zone[3:81,30]), N=78, "V"=11.9,"h"=0.01))
update(jagsmed, 1000)
mcmcjagsmed<-jags.samples(jagsmed,
                       c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','loglik','logyhat'),
                       1000)
mcmcsamplesjagsmed<-coda.samples(jagsmed,
                              c('Qprime','Gprime','kl','sigma','c','aomega','bomega','loglik','logyhat'),
                              1000)
m1.mcmcmed<-(as.mcmc(mcmcsamplesjagsmed))
m1.matmed<-as.matrix(mcmcsamplesjagsmed)
m1.datmed<-as.data.frame(m1.matmed)
gprime.post1zonemed<-m1.datmed$Gprime
qprime.post1zonemed<-m1.datmed$Qprime
aomega.post1zonemed<-m1.datmed$aomega
bomega.post1zonemed<-m1.datmed$bomega
sigma.post1zonemed<-sqrt(m1.datmed$sigma)
kl.post1zonemed<-m1.datmed$kl
quantile((gprime.post1zonemed),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonemed),c(0.025,0.25,0.5,0.75,0.975))
quantile((sigma.post1zonemed),c(0.025,0.25,0.5,0.75,0.975))
loglik.post1zonem<-t(mcmcjagsmed$loglik[,,1])
dim(loglik.post1zonem)
waic(loglik.post1zonem[,1:81])
mean(crps_sample(y = data1zone[3:81,30], dat = exp(mcmcjagsmed$logyhat[,,1])))
yhatpm<-apply(exp(mcmcjagsmed$logyhat[,,1]), 1, mean)
yhatlower50pm<-apply(exp(mcmcjagsmed$logyhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50pm<-apply(exp(mcmcjagsmed$logyhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90pm<-apply(exp(mcmcjagsmed$logyhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90pm<-apply(exp(mcmcjagsmed$logyhat[,,1]), 1, function(x) quantile(x,0.975))
#nominal coverage
plot(data1zone[3:81,30],type="l")
lines(yhatlower50pm,col="blue")
lines(yhatupper50pm,col="blue")
cov50m<-cov90m<-NULL
for(i in 1:length(data1zone[3:81,30])){
  cov50m[i]<-(ifelse(yhatlower50pm[i]<=data1zone[i+2,30]&data1zone[i+2,30]<=yhatupper50pm[i],1,0))
  cov90m[i]<-(ifelse(yhatlower90pm[i]<=data1zone[i+2,30]&data1zone[i+2,30]<=yhatupper90pm[i],1,0))
}
sum(cov50m)/79
sum(cov90m)/79
c.hatmed <- apply(mcmcjagsmed$c, 1, mean)
zrepm<-matrix(0,79,1000)
for(j in 1:1000){
  zrepm[,j]<-exp(rnorm(79,log(mcmcjagsmed$c[,j,1]),(sigma.post1zonemed[j])))
}
G1m<-(data1zone[3:81,30]-apply((zrepm),1,mean))^2
P1m<-apply((zrepm),1,var)
print(D1m<-sum(G1m)+sum(P1m))
cs1m<-matrix(0,1000,81)
cs1m[,81]<-mcmcsamplesjagsmed[[1]][,(4+81)]
omega.post1zonem=rgamma(aomega.post1zonemed,bomega.post1zonemed)
N=81-1
for(j in 2:(N)){
  cs1m[,(j+N-2*(j-1))]<-(mcmcsamplesjagsmed[[1]][,4+(j+N-2*(j-1)+1)]-h*(gprime.post1zonemed)/V-(h*(omega.post1zonem)))/(1-h*((qprime.post1zonemed)+
                                                                                                                                 (kl.post1zonemed)*V)/V)
} 
cs1m[,1]<-(mcmcsamplesjagsmed[[1]][,6]-h*(gprime.post1zonemed)/V+h*(omega.post1zonem))/(1-h*((qprime.post1zonemed)+(kl.post1zonemed)*V)/V)
cs.hatm <- apply(cs1m, 2, mean)
cs.hatm
#high ventilation Q=0.47-0.77 low G=39.55
jagsh<-jags.model("modeljags.jag",data=list("logyt"=log(data1zone[2:41,51*4]), N=39, "V"=11.9,"h"=0.01))
update(jagsh, 1000)
mcmcjagsh<-jags.samples(jagsh,
                       c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','loglik','logyhat'),
                       1000)
mcmcsamplesjagsh<-coda.samples(jagsh,
                              c('Qprime','Gprime','kl','sigma','c','aomega','bomega','loglik','logyhat'),
                              1000)
summary(mcmcsamplesjagsh)
m1.mcmch<-(as.mcmc(mcmcsamplesjagsh))
m1.math<-as.matrix(mcmcsamplesjagsh)
m1.dath<-as.data.frame(m1.math)
gprime.post1zoneh<-m1.dath$Gprime
qprime.post1zoneh<-m1.dath$Qprime
aomega.post1zoneh<-m1.dath$aomega
bomega.post1zoneh<-m1.dath$bomega
sigma.post1zoneh<-sqrt(m1.dath$sigma)
kl.post1zoneh<-m1.dath$kl
quantile((gprime.post1zoneh),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zoneh),c(0.025,0.25,0.5,0.75,0.975))
loglik.post1zoneh<-t(mcmcjagsh$loglik[,,1])
dim(loglik.post1zoneh)
waic(loglik.post1zoneh[,1:41])
mean(crps_sample(y =data1zone[2:41,51*4], dat = exp(mcmcjagsh$logyhat[,,1])))
yhatph<-apply(exp(mcmcjagsh$logyhat[,,1]), 1, mean)
yhatlower50ph<-apply(exp(mcmcjagsh$logyhat[,,1]), 1, function(x) quantile(x,0.25))
yhatupper50ph<-apply(exp(mcmcjagsh$logyhat[,,1]), 1, function(x) quantile(x,0.75))
yhatlower90ph<-apply(exp(mcmcjagsh$logyhat[,,1]), 1, function(x) quantile(x,0.05))
yhatupper90ph<-apply(exp(mcmcjagsh$logyhat[,,1]), 1, function(x) quantile(x,0.975))
#nominal coverage
plot(data1zone[2:41,51*4],type="l")
lines(yhatlower50ph,col="blue")
lines(yhatupper50ph,col="blue")
cov50h<-cov90h<-NULL
for(i in 1:length(data1zone[2:41,51*4])){
  cov50h[i]<-(ifelse(yhatlower50ph[i]<=data1zone[i+1,51*4]&data1zone[i+1,51*4]<=yhatupper50ph[i],1,0))
  cov90h[i]<-(ifelse(yhatlower90ph[i]<=data1zone[i+1,51*4]&data1zone[i+1,51*4]<=yhatupper90ph[i],1,0))
}
sum(cov50h)/40
sum(cov90h)/40
c.hath <- apply(mcmcjagsh$c, 1, mean)
zreph<-matrix(0,40,1000)
for(j in 1:1000){
  zreph[,j]<-rnorm(40,log(mcmcjagsh$c[,j,1]),(sigma.post1zoneh[j]))
}
G1h<-(data1zone[2:41,51*4]-apply(exp(zreph),1,mean))^2
P1h<-apply(exp(zreph),1,var)
print(D1h<-sum(G1h)+sum(P1h))
cs1h<-matrix(0,1000,41)
cs1h[,41]<-mcmcsamplesjagsh[[1]][,(4+41)]
omega.post1zoneh=rgamma(aomega.post1zoneh,bomega.post1zoneh)
N=41-1
for(j in 2:(N)){
  cs1h[,(j+N-2*(j-1))]<-(mcmcsamplesjagsh[[1]][,4+(j+N-2*(j-1)+1)]-h*(gprime.post1zoneh)/V-(h*(omega.post1zoneh)))/(1-h*((qprime.post1zoneh)+
                                                                                                                               (kl.post1zoneh)*V)/V)
} 
cs1h[,1]<-(mcmcsamplesjagsh[[1]][,6]-h*(gprime.post1zoneh)/V+h*(omega.post1zoneh))/(1-h*((qprime.post1zoneh)+(kl.post1zoneh)*V)/V)
cs.hath <- apply(cs1h, 2, mean)
cs.hath
quartz()
par(mfrow=c(1,3))
plot(c.hatlow, main="low ventilation")
points(data1zone[1:201,2],col="red")
points(cs.hatl, col="blue")
legend("bottomright", col=c("black","red"),legend=c("chat low","measured"))
plot(c.hatmed, main="medium ventilation")
points(data1zone[1:81,30],col="red")
points(cs.hatm,col="blue")
legend("bottomright", col=c("black","red"),legend=c("chat medium","measured"))
plot(c.hath,main="High ventilation")
points(data1zone[,51*4],col="red")
points(cs.hath,col="blue")
legend("bottomright", col=c("black","red"),legend=c("chat high","measured"))
#Pf data one zone
#Pf import from julia output
#low
#thetajuliadl<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmcl.csv")
thetajuliadcorrl<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmccorr.csv")
apply(thetajuliadcorrl, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
xjuliadcorrl<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmclcorr.csv")
pred.juliadcorrl<-apply(xjuliadcorrl, 2, mean)
zrepphl<-matrix(0,201,1000)
for(j in 1:1000){
  for(i in 1:201){
    zrepphl[i,j]<-exp(rnorm(1,log(xjuliadcorrl[j,i]),(thetajuliadcorrl$x6[j])))
  }}
G1phl<- (data1zone[1:201,2]-apply((zrepphl),1,mean))^2
P1phl<-apply((zrepphl),1,var)
print(D1phl<-sum(G1phl)+sum(P1phl))
# loglikpfl<-matrix(0,1000,201)
# for(j in 1:1000){
# for(i in 1:201){
# loglikpfl[j,i]<-log(dnorm(log(data1zone[i,2]),log(xjuliadl[j,i]),thetajuliadl$x6[j]))
# }}
# waic(loglikpfl[,1:201])
# hist(c(loglikpfl))
#medium
#thetajuliadm<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmcm.csv")
thetajuliadcorrm<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmcmcorr.csv")
apply(thetajuliadcorrm, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
xjuliadcorrm<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmcmcorr.csv")
pred.juliadcorrm<-apply(xjuliadcorrm, 2, mean)
zrepphm<-matrix(0,81,1000)
for(j in 1:1000){
  for(i in 1:81){
    zrepphm[i,j]<-exp(rnorm(1,log(xjuliadcorrm[j,i]),sqrt(thetajuliadcorrm$x6[j])))
  }}
G1phm<- (data1zone[1:81,30]-apply((zrepphm),1,mean))^2
P1phm<-apply((zrepphm),1,var)
print(D1phm<-sum(G1phm)+sum(P1phm))
# loglikpfm<-matrix(0,1000,81)
# for(j in 1:1000){
# for(i in 1:81){
#   loglikpfm[,i]<-(dnorm(log(data1zone[i,30]+0.05),log(xjuliadm[j+1000,i]+0.05),thetajuliadm$x6[j+1000]))
# }}
# waic(loglikpfm[,1:81])
# compare(waic(loglik.post1zonem[,2:81]),waic(loglikpfm[,2:81]))
#high
#thetajuliadh<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmch.csv")
thetajuliadhcorr<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmchcorr.csv")
apply(thetajuliadhcorr, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
xjuliadhcorr<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmchcorr.csv")
pred.juliadhcorr<-apply(xjuliadhcorr, 2, mean)
zrepphh<-matrix(0,41,1000)
for(j in 1:1000){
  for(i in 2:41){
    zrepphh[i,j]<-exp(rnorm(1,log(xjuliadhcorr[j,i]),sqrt(thetajuliadhcorr$x6[j])))
  }}
G1phh<- (data1zone[1:41,51*4]-apply((zrepphh),1,mean))^2
P1phh<-apply((zrepphh),1,var)
print(D1phh<-sum(G1phh)+sum(P1phh))
loglikpfh<-matrix(0,1000,40)
for(i in 1:40){
  for(j in 1:1000){
  loglikpfh[j,i]<-log(dnorm(log(data1zone[i+1,51*4]),log(xjuliadh[j+1000,i+1]),thetajuliadh$x6[j+1000]))
}}
waic(loglikpfh)

#mse
print(mse1pfl<-sum((data1zone[1:201,2]-(pred.juliadcorrl))^2)/201)
print(mse1stl<-sum((data1zone[1:201,2]-(c.hatlow))^2)/201)
print(mse1pfm<-sum((data1zone[1:81,30]-(pred.juliadcorrm))^2)/81)
print(mse1stm<-sum((data1zone[1:81,30]-(c.hatmed))^2)/81)
print(mse1pfh<-sum((data1zone[1:41,51*4]-(pred.juliadhcorr))^2)/41)
print(mse1sth<-sum((data1zone[1:41,51*4]-(c.hath))^2)/41)
#plot all data
quartz()
par(mfrow=c(3,1))
# plot(pred.juliadl,col="red", main="Low", ylab="Toluene Concentrations",
#      xlab="Time")
# points(c.hatlow, col="blue")
# points(data1zone[1:201,2], col="green",cex=0.8)
# points(predkfdl,col="chocolate",cex=0.5)
# points(cs.hatl,col="purple",cex=0.5)
plot(data1zone[1:201,2], col="green",lwd=3, main="Low", ylab="Toluene Concentrations",
     xlab="Time",type="l")
#lines(pred.juliadl,col="red",lwd=3)
lines(c.hatlow, col="blue",lty=6,lwd=3)
lines(predkfdl,col="chocolate",lwd=3,lty=2)
lines(cs.hatl,col="purple",lwd=3,lty=3)
legend("bottomright", col=c("green", "blue", "chocolate","purple"),lty=c(1,6,2,3) ,legend=c("Measurements",
                                                              "Model a", "Model b","Smoothed"),pch=15,cex=1,box.lwd = 0)
# plot(pred.juliadm,col="red", main="Medium", ylab="Toluene Concentrations",
#      xlab="Time")
# points(c.hatmed, col="blue")
# points(data1zone[1:81,30], col="green",cex=0.8)
# points(predkfdm,col="chocolate",cex=0.5)
# points(cs.hatm,col="purple",cex=0.5)
plot(data1zone[1:81,30], col="green", main="Medium", ylab="Toluene Concentrations",
     xlab="Time",type="l",lwd=3)
#lines(pred.juliadm,col="red",lwd=3)
lines(c.hatmed, col="blue",lwd=3,lty=6)
lines(predkfdm,col="chocolate",lwd=3,lty=2)
lines(cs.hatm,col="purple",lwd=3,lty=3)
legend("bottomright", col=c("green", "blue", "chocolate","purple"), lty=c(1,6,2,3),legend=c("Measurements",
                                                              "Model a", "Model b","Smoothed"),pch=15,cex=1,box.lwd = 0)
# plot(pred.juliadh,col="red", main="High", ylab="Acetone Concentrations",
#        xlab="Time")
# points(c.hath, col="blue")
# points(data1zone[1:41,203], col="green",cex=0.8)
# points(predkfdh,col="chocolate",cex=0.5)
# points(cs.hath,col="purple",cex=0.5)
plot(data1zone[1:41,203], col="green", main="High", ylab="Acetone Concentrations",
     xlab="Time",type="l",lwd=3)
#lines(pred.juliadh,col="red",lwd=3)
lines(c.hath , col="blue",lwd=3,lty=6)
lines(predkfdh,col="chocolate",lwd=3,lty=2)
lines(cs.hath,col="purple",lwd=3,lty=3)
legend("bottomright", col=c("green", "blue", "chocolate","purple"),lty=c(1,6,2,3),legend=c("Measurements",
                                                              "Model a", "Model b","Smoothed"),pch=15,cex=1,box.lwd = 0)
#nonlinear one zone data
cat("model{
    for (i in 1:N){ 
    c[i+1]<-exp(-i*(Qprime+kl*V)/V)*c[1]+1/(-(Qprime+kl*V)/V)*
    (exp(-i*(Qprime+kl*V)/V)-1)*Gprime/V
    logyt[i+1]~dnorm(log((c[i+1])), prec)
    }
    c[1]<-2
    logyt[1]~dnorm(log(c[1]), prec)
    sigma<-1/prec
    prec~dgamma(1,0.01)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,0.8)
    }",file="modeljagsnld.jag")
jagslownl<-jags.model("modeljagsnld.jag",data=list("logyt"=log(data1zone[1:201,2]), N=200, "V"=11.9))
update(jagslownl, 1000)
mcmcjagslownl<-jags.samples(jagslownl,
                          c('Qprime','Gprime',"kl",'sigma','c'),
                          1000)
mcmcsamplesjagslownl<-coda.samples(jagslownl,
                                 c('Qprime','Gprime','kl','sigma','c'),
                                 1000)
#low ventilation Q=0.04-0.07 G=43.18
m1.mcmclownl<-(as.mcmc(mcmcsamplesjagslownl))
m1.matlownl<-as.matrix(mcmcsamplesjagslownl)
m1.datlownl<-as.data.frame(m1.matlownl)
qprime.post1zonelownl<-m1.datlownl$Qprime
gprime.post1zonelownl<-m1.datlownl$Gprime
kl.post1zonelownl<-m1.datlownl$kl
quantile((gprime.post1zonelownl),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonelownl),c(0.025,0.25,0.5,0.75,0.975))
quantile((kl.post1zonelownl),c(0.025,0.25,0.5,0.75,0.975))
c.hatlownl <- apply(mcmcjagslownl$c, 1, mean)

#medium Q Q=0.23-0.27 G=43.18
jagsmednl<-jags.model("modeljagsnld.jag",data=list("logyt"=log(data1zone[,30]+0.05), N=80, "V"=11.9))
update(jagsmednl, 1000)
mcmcjagsmednl<-jags.samples(jagsmednl,
                          c('Qprime','Gprime',"kl",'sigma','c'),
                          1000)
mcmcsamplesjagsmednl<-coda.samples(jagsmednl,
                                 c('Qprime','Gprime','kl','sigma','c'),
                                 1000)
m1.mcmcmednl<-(as.mcmc(mcmcsamplesjagsmednl))
m1.matmednl<-as.matrix(mcmcsamplesjagsmednl)
m1.datmednl<-as.data.frame(m1.matmednl)
gprime.post1zonemednl<-m1.datmednl$Gprime
qprime.post1zonemednl<-m1.datmednl$Qprime
quantile((gprime.post1zonemednl),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonemednl),c(0.025,0.25,0.5,0.75,0.975))
c.hatmednl <- apply(mcmcjagsmednl$c, 1, mean)

#high ventilation Q=0.47-0.77 low G=39.55
jagshnl<-jags.model("modeljagsnld.jag",data=list("logyt"=log(data1zone[,51*4]+0.05), N=40, "V"=11.9))
update(jagshnl, 1000)
mcmcjagshnl<-jags.samples(jagshnl,
                        c('Qprime','Gprime',"kl",'sigma','c'),
                        1000)
mcmcsamplesjagshnl<-coda.samples(jagshnl,
                               c('Qprime','Gprime','kl','sigma','c'),
                               1000)
m1.mcmchnl<-(as.mcmc(mcmcsamplesjagshnl))
m1.mathnl<-as.matrix(mcmcsamplesjagshnl)
m1.dathnl<-as.data.frame(m1.mathnl)
gprime.post1zonehnl<-m1.dathnl$Gprime
qprime.post1zonehnl<-m1.dathnl$Qprime
quantile((gprime.post1zonehnl),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonehnl),c(0.025,0.25,0.5,0.75,0.975))
c.hathnl <- apply(mcmcjagshnl$c, 1, mean)
quartz()
par(mfrow=c(1,3))
plot(c.hatlownl, main="low ventilation")
points(data1zone[1:201,2],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat low","measured"))
plot(c.hatmednl, main="medium ventilation")
points(data1zone[1:81,30],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat medium","measured"))
plot(c.hathnl,main="High ventilation")
points(data1zone[,51*4],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat high","measured"))
#two compartment model simulation
n=100
Qprime=13.8
Gprime=351.4
#VF=11.9
#VN=0.1
VF=3.8
VN=pi*10^-2
Sigma<-diag(0.1,2)
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
#values of c very huge or very small and hence yt so taking log is problem
for(i in 1:(n-1)){
  c[1,]<-c(0,0.5)
  c[i+1,]<-(c[i,])%*%t(h*(A)+diag(1,2))+(g*h)
  #c[i+1,1]<-c[i,1]*((A[1,1]*VN*h)+1)+c[i,2]*A[1,2]*VN*h+VN*g[1]*h
  #c[i+1,2]<-c[i,1]*A[2,1]*VF*h+c[i,2]*((A[2,2]*VF*h)+1)
}
c
plot(c[,1],ylab="mg/m3", xlab="time")
points(c[,2])
compute_CNCF <- function(Beta, Q, G, VN, VF, t){
  lambda1 <- 0.5*(-(Beta*VF+(Beta+Q)*VN)/(VN*VF) + sqrt(((Beta*VF+(Beta+Q)*VN)/(VN*VF))^2 - 4*(Beta*Q/(VN*VF))))
  lambda2 <- 0.5*(-(Beta*VF+(Beta+Q)*VN)/(VN*VF) - sqrt(((Beta*VF+(Beta+Q)*VN)/(VN*VF))^2 - 4*(Beta*Q/(VN*VF))))
  elambda1 <- exp(lambda1*t)
  elambda2 <- exp(lambda2*t)
  CN <- G/Q + G/Beta + G*((Beta*Q + lambda2*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda1 - 
    G*((Beta*Q + lambda1*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda2
  
  CF <- G/Q + G*((lambda1*VN+Beta)/Beta)*((Beta*Q + lambda2*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda1 - 
    G*((lambda2*VN+Beta)/Beta)*((Beta*Q + lambda1*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda2 
  
  return(cbind(CN,CF))
}
lambda1 <- 0.5*(-(Beta*VF+(Beta+Q)*VN)/(VN*VF) + sqrt(((Beta*VF+(Beta+Q)*VN)/(VN*VF))^2 - 4*(Beta*Q/(VN*VF))))
lambda2 <- 0.5*(-(Beta*VF+(Beta+Q)*VN)/(VN*VF) - sqrt(((Beta*VF+(Beta+Q)*VN)/(VN*VF))^2 - 4*(Beta*Q/(VN*VF))))

compute_CNCF(5,13.8,352,1.1,240,2)
hist(c22[,2]-c[,2])

#exact 
c22<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
#imp

for(i in 1:(n-1)){
  c22[1,]<-c(0.01,0.5)
  c22[i+1,]<- expm(h*i*A)%*%c22[1,]+solve(A)%*%(expm(h*i*A)-diag(2))%*%t(g)
}
plot(c22[,1])
points(c22[,2])
plot(c[,1],c22[,1])
abline(0,1)
plot(c[,2],c22[,2])
abline(0,1)

c<-c2

set.seed(123)
vt<-rmvnorm(n,c(0,0),Sigma)
logyt2<-log((c22))+(vt)
exp(logyt2)

#jags
I<-diag(1,2)

cat("model{
    for (i in 1:N){
    logomega[i+1,1:2]~dmnorm(c(0,0),vari)
    c[i+1,1]<-c[i,1]*((A[1,1]*h)+1)+c[i,2]*A[1,2]*h+g[1]*h+h*exp(logomega[i+1,1])
    c[i+1,2]<-c[i,1]*A[2,1]*h+c[i,2]*((A[2,2]*h)+1)+h*exp(logomega[i+1,2])
    logyt[i+1,1:2]~dmnorm(log(abs(c[i+1,1:2])), prec)
    logyhat[i+1,1:2]~dmnorm(log(abs(c[i+1,1:2])), prec)
    loglik[i+1]<-(logdensity.mnorm(logyt[i+1,1:2],log(c[i+1,1:2]),prec))
    }
    logyhat[1,1:2]~dmnorm(log(abs(c[1,1:2])), prec)
    A[1,1]<- -(beta)/VN
    A[1,2]<-(beta)/VN
    A[2,1]<- beta/VF
    A[2,2]<- -(beta+Qprime)/VF+kl
    g[1]<-Gprime/VN
    g[2]<-0
    c[1,1]<-0.01
    c[1,2]<-0.5
    logyt[1,1:2]~dmnorm(log(c[1,1:2]), prec)
    prec~dwish(I,3)
    vari~dwish(I,3)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    beta~dunif(0,10)
    }",file="modeljags2zone.jag")
jags<-jags.model("modeljags2zone.jag",data=list("logyt"=logyt2, N=99,"I"=diag(2), "VN"=pi*10^-2, "VF"=3.8,"h"=h))
update(jags, 1000)
mcmcjags.2zone<-jags.samples(jags,
                             c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                             1000)
mcmcsamplesjags.2zone<-coda.samples(jags,
                                    c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                                    1000)
quartz()
plot(mcmcsamplesjags.2zone[[1]][,1])
summary(mcmcsamplesjags.2zone[[1]])
summary(1/mcmcsamplesjags.2zone[[1]][,205:212])
m1.mcmc.2zone<-(as.mcmc(mcmcsamplesjags.2zone))
m1.mat.2zone<-as.matrix(mcmcsamplesjags.2zone)
m1.dat.2zone<-as.data.frame(m1.mat.2zone)
qprime.post.2zone<-m1.dat.2zone$Qprime
gprime.post.2zone<-m1.dat.2zone$Gprime
beta.post.2zone<-m1.dat.2zone$beta
kl.post.2zone<-m1.dat.2zone$kl
vari11.post.2zone<-m1.dat.2zone$`prec[1,1]`
vari12.post.2zone<-m1.dat.2zone$`prec[1,2]`
vari21.post.2zone<-m1.dat.2zone$`prec[2,1]`
vari22.post.2zone<-m1.dat.2zone$`prec[2,2]`
loglik.post2zone<-t(mcmcjags.2zone$loglik[,,1])
dim(loglik.post2zone)
waic(loglik.post2zone[,2:100])
quantile(kl.post.2zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(qprime.post.2zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zone,c(0.025,0.25,0.5,0.75,0.975))
mean(crps_sample(y = exp(logyt2[,1]), dat = exp(mcmcjags.2zone$logyhat[,1,,1])))+
  mean(crps_sample(y = exp(logyt2[,2]), dat =exp(mcmcjags.2zone$logyhat[,2,,1])))
yhatlower50np<- apply( exp(mcmcjags.2zone$logyhat[,1,,1]), 1, function(x) quantile(x,0.25))
yhatupper50np<-apply( exp(mcmcjags.2zone$logyhat[,1,,1]), 1, function(x) quantile(x,0.75))
yhatlower90np<- apply( exp(mcmcjags.2zone$logyhat[,1,,1]), 1, function(x) quantile(x,0.05))
yhatupper90np<-apply( exp(mcmcjags.2zone$logyhat[,1,,1]), 1, function(x) quantile(x,0.95))
yhatlower50fp<- apply( exp(mcmcjags.2zone$logyhat[,2,,1]), 1, function(x) quantile(x,0.25))
yhatupper50fp<-apply( exp(mcmcjags.2zone$logyhat[,2,,1]), 1, function(x) quantile(x,0.75))
yhatlower90fp<- apply( exp(mcmcjags.2zone$logyhat[,2,,1]), 1, function(x) quantile(x,0.05))
yhatupper90fp<-apply( exp(mcmcjags.2zone$logyhat[,2,,1]), 1, function(x) quantile(x,0.95))
#nominal coverage
plot(exp(logyt2[,2]))
lines(yhatlower50fp,col="blue")
lines(yhatupper50fp,col="blue")
cov502n<-cov502f<-cov902n<-cov902f<-NULL
for(i in 1:length(logyt2[,1])){
  cov502n[i]<-(ifelse(yhatlower50np[i]<=exp(logyt2[i,1])&exp(logyt2[i,1])<=yhatupper50np[i],1,0))
  cov502f[i]<-(ifelse(yhatlower50fp[i]<=exp(logyt2[i,2])&exp(logyt2[i,2])<=yhatupper50fp[i],1,0))
  cov902n[i]<-(ifelse(yhatlower90np[i]<=exp(logyt2[i,1])&exp(logyt2[i,1])<=yhatupper90np[i],1,0))
  cov902f[i]<-  (ifelse(yhatlower90fp[i]<=exp(logyt2[i,2])&exp(logyt2[i,2])<=yhatupper90fp[i],1,0))
}
sum(cov502n)/200+sum(cov502f)/200
sum(cov902n)/200+sum(cov902f)/200
c1.hat.2zone <- apply(mcmcsamplesjags.2zone[[1]][,4:103], 2, mean)
c2.hat.2zone <- apply(mcmcsamplesjags.2zone[[1]][,104:203], 2, mean)
zrep21<-matrix(0,100,1000)
zrep22<-matrix(0,100,1000)
for(j in 1:1000){
  for(i in 1:100){
  zrep21[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zone[[1]][j,3+i]),log(mcmcsamplesjags.2zone[[1]][j,103+i])),
                      sigma=solve(matrix(c(vari11.post.2zone[j],vari12.post.2zone[j],vari21.post.2zone[j],
                                           vari22.post.2zone[j]),2,2,byrow=TRUE)),method="svd")[1]
  zrep22[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zone[[1]][j,3+i]),log(mcmcsamplesjags.2zone[[1]][j,103+i])),
                       sigma=solve(matrix(c(vari11.post.2zone[j],vari12.post.2zone[j],vari21.post.2zone[j],
                                            vari22.post.2zone[j]),2,2,byrow=TRUE)),method="svd")[2]
}}
G2<-(c22[,1]-apply(exp(zrep22),1,mean))^2+(c22[,2]-apply(exp(zrep21),1,mean))^2
P2<-apply(exp(zrep22),1,var)+apply(exp(zrep21),1,var)
print(D2<-sum(G2)+sum(P2))
logomega<-matrix(0,1000,2)
for(i in 1000){
logomega[,1:2]<-rlnorm.rplus(1,c(0,0),matrix(c(m1.dat.2zone[i,308],m1.dat.2zone[i,309],
                                         m1.dat.2zone[i,310],m1.dat.2zone[i,311]),2,2,byrow=TRUE))[1,]
}
N=99
cs12<-matrix(0,1000,100)
cs22<-matrix(0,1000,100)
cs12[,N+1]<-mcmcsamplesjags.2zone[[1]][,103]
cs22[,N+1]<-mcmcsamplesjags.2zone[[1]][,203]
As<-array(rep(0,4000),dim=c(2,2,1000))
Bs<-array(rep(0,4000),dim=c(2,2,1000))
gs<-array(rep(0,2000),dim=c(1,2,1000))
for(i in 1:1000){
for(j in 2:N){
  As[1,1,i]<- -(beta.post.2zone[i])/VN
  As[1,2,i]<-(beta.post.2zone[i])/VN
  As[2,1,i]<- beta.post.2zone[i]/VF
  As[2,2,i]<- -(beta.post.2zone[i]+qprime.post.2zone[i])/VF+kl
  gs[1,1,i]<-gprime.post.2zone[i]/VN
  gs[1,2,i]<-0
  Bs[1,1,i]<-(-h*(beta.post.2zone[i]+qprime.post.2zone[i])/VF+1)/(((h*As[1,1,i]+1)*(h*As[2,2,i]+1))-((h*As[2,1,i]*h*As[1,2,i])))
  Bs[1,2,i]<-(-h*beta.post.2zone[i]/VN)/(((h*As[1,1,i]+1)*(h*As[2,2,i]+1))-((h*As[2,1,i]*h*As[1,2,i])))
  Bs[2,1,i]<- -(h*beta.post.2zone[i]/VF)/(((h*As[1,1,i]+1)*(h*As[2,2,i]+1))-((h*As[2,1,i]*h*As[1,2,i])))
  Bs[2,2,i]<- (-h*beta.post.2zone[i]/VN+1)/(((h*As[1,1,i]+1)*(h*As[2,2,i]+1))-((h*As[2,1,i]*h*As[1,2,i])))
  cs12[i,(j+N-2*(j-1))]<- ((Bs[,,i])%*%t(c(mcmcsamplesjags.2zone[[1]][i,(j+N-2*(j-1)+1+3)],mcmcsamplesjags.2zone[[1]][i,(j+N-2*(j-1)+1+103)])-t(h*gs[1,1:2,i])-
                                   t(c((logomega[(j+N-2*(j-1)+1),1]),(logomega[(j+N-2*(j-1)+1),2])))))[1,1]
  cs22[i,(j+N-2*(j-1))]<- ((Bs[,,i])%*%t(c(mcmcsamplesjags.2zone[[1]][i,(j+N-2*(j-1)+1+3)],mcmcsamplesjags.2zone[[1]][i,(j+N-2*(j-1)+1+103)])-t(h*gs[1,1:2,i])-
                                          t(c((logomega[(j+N-2*(j-1)+1),1]),(logomega[(j+N-2*(j-1)+1),2])))))[2,1]
} 
  cs12[i,1]<-((Bs[,,i])%*%t(c(mcmcsamplesjags.2zone[[1]][i,5],mcmcsamplesjags.2zone[[1]][i,105])-t(h*gs[1,1:2,i])-
                              t(c((logomega[2,1]),(logomega[2,2])))))[1,1]
  cs22[i,1]<-((Bs[,,i])%*%t(c(mcmcsamplesjags.2zone[[1]][i,5],mcmcsamplesjags.2zone[[1]][i,105])-t(h*gs[1,1:2,i])-
                              t(c((logomega[2,1]),(logomega[2,2])))))[1,1]
  }

cs12hat<-apply(cs12,2,mean)
cs22hat<-apply(cs22,2,mean)
plot(cs12hat)
points(cs22hat)

c1.hatlower.2zone<- apply(mcmcsamplesjags.2zone[[1]][,4:103], 2, function(x) quantile(x,0.025))
c1.hatupper.2zone<-apply(mcmcsamplesjags.2zone[[1]][,4:103], 2, function(x) quantile(x,0.975))
c2.hatlower.2zone<- apply(mcmcsamplesjags.2zone[[1]][,104:203], 2, function(x) quantile(x,0.025))
c2.hatupper.2zone<-apply(mcmcsamplesjags.2zone[[1]][,104:203], 2, function(x) quantile(x,0.975))
quartz()
plotCI(c2[,1], c1.hat.2zone,li=c1.hatlower.2zone, ui=c1.hatupper.2zone)
abline(c(0,1))
plotCI(c2[,2], c2.hat.2zone,li=c2.hatlower.2zone, ui=c2.hatupper.2zone)
abline(c(0,1))

quartz()
par(mfrow=c(1,2))
plot(c22[,1],col="red",type="l")
points(c1.hat.2zone, col="blue")
points(exp(logyt2[,1]),col="green")
points(c1s.hat.2zone)
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtN","chatN","ytN"),bg="white", lty=1,cex=0.8)
plot(c22[,2],col="red")
points(c2.hat.2zone, col="blue",type="l")
points(exp(logyt2[,2]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtF","chatF","ytF"),bg="white", lty=1,cex=0.8)
#sims
jags2sim<-mcmcjags.2zone2sim<-mcmcsamplesjags.2zone2sim<-vector("list")
for(i in 1:100){
  jags2sim[[i]]<-jags.model("modeljags2zone.jag",data=list("logyt"=(logyt2sim[[i]]), N=99,"I"=diag(2), "VN"=pi*10^-2, "VF"=3.8,"h"=h))
  update(jags2sim[[i]], 1000)
  mcmcjags.2zone2sim[[i]]<-jags.samples(jags2sim[[i]],
                                        c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                                        1000)
  mcmcsamplesjags.2zone2sim[[i]]<-coda.samples(jags2sim[[i]],
                                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                                              1000)
}
m1.mcmc.2zonesim<-m1.mat.2zonesim<-m1.dat.2zonesim<-yhatlower50nsimpar<-
  yhatupper50nsimpar<-yhatlower50fsimpar<-yhatupper50fsimpar<-yhatlower90nsimpar<-
  yhatupper90nsimpar<-yhatlower90fsimpar<-yhatupper90fsimpar<-yhatnsimpar<-yhatfsimpar<-chatlower50nsimpar<-
  chatupper50nsimpar<-chatlower50fsimpar<-chatupper50fsimpar<-chatlower90nsimpar<-
  chatupper90nsimpar<-chatlower90fsimpar<-chatupper90fsimpar<-chatnsimpar<-chatfsimpar<-vector("list")
for(i in 1:100){
  m1.mcmc.2zonesim[[i]]<-(as.mcmc(mcmcsamplesjags.2zone2sim[[i]]))
  m1.mat.2zonesim[[i]]<-as.matrix(mcmcsamplesjags.2zone2sim[[i]])
  m1.dat.2zonesim[[i]]<-as.data.frame(m1.mat.2zonesim[[i]])
  yhatnsimpar[[i]]<- apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,1,,1])), 1, mean)
  yhatfsimpar[[i]]<- apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,2,,1])), 1,mean)
  yhatlower50nsimpar[[i]]<- apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,1,,1])), 1, function(x) quantile(x,0.25))
  yhatupper50nsimpar[[i]]<-apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,1,,1])), 1, function(x) quantile(x,0.75))
  yhatlower50fsimpar[[i]]<- apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,2,,1])), 1, function(x) quantile(x,0.25))
  yhatupper50fsimpar[[i]]<-apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,2,,1])), 1, function(x) quantile(x,0.75))
  yhatlower90nsimpar[[i]]<- apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,1,,1])), 1, function(x) quantile(x,0.05))
  yhatupper90nsimpar[[i]]<-apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,1,,1])), 1, function(x) quantile(x,0.95))
  yhatlower90fsimpar[[i]]<- apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,2,,1])), 1, function(x) quantile(x,0.05))
  yhatupper90fsimpar[[i]]<-apply((exp(mcmcjags.2zone2sim[[i]]$logyhat[,2,,1])), 1, function(x) quantile(x,0.95))
  chatnsimpar[[i]]<- apply((mcmcjags.2zone2sim[[i]])$c[,1,,1], 1, mean)
  chatfsimpar[[i]]<- apply((mcmcjags.2zone2sim[[i]])$c[,2,,1], 1,mean)
  chatlower50nsimpar[[i]]<- apply((mcmcjags.2zone2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.25))
  chatupper50nsimpar[[i]]<-apply((mcmcjags.2zone2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.75))
  chatlower50fsimpar[[i]]<- apply((mcmcjags.2zone2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.25))
  chatupper50fsimpar[[i]]<-apply((mcmcjags.2zone2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.75))
  chatlower90nsimpar[[i]]<- apply((mcmcjags.2zone2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.05))
  chatupper90nsimpar[[i]]<-apply((mcmcjags.2zone2sim[[i]])$c[,1,,1], 1, function(x) quantile(x,0.95))
  chatlower90fsimpar[[i]]<- apply((mcmcjags.2zone2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.05))
  chatupper90fsimpar[[i]]<-apply((mcmcjags.2zone2sim[[i]])$c[,2,,1], 1, function(x) quantile(x,0.95))
}
cov502nsim<-cov502fsim<-cov902nsim<-cov902fsim<-mse2nsim<-mse2fsim<-crps2simpar<-vector("list")
for(i in 1:100){
  cov502nsim[[i]]<-(ifelse(chatlower50nsimpar[[i]]<=(c[,1])&(c[,1])<=chatupper50nsimpar[[i]],1,0))
  cov502fsim[[i]]<-(ifelse(chatlower50fsimpar[[i]]<=(c[,2])&(c[,2])<=chatupper50fsimpar[[i]],1,0))
  cov902nsim[[i]]<-(ifelse(chatlower90nsimpar[[i]]<=(c[,1])&(c[,1])<=chatupper90nsimpar[[i]],1,0))
  cov902fsim[[i]]<- (ifelse(chatlower90fsimpar[[i]]<=(c[,2])&(c[,2])<=chatupper90fsimpar[[i]],1,0))
  mse2nsim[[i]]<-sum((chatnsimpar[[i]]-(c[,1]))^2)/100
  mse2fsim[[i]]<-sum((chatfsimpar[[i]]-(c[,2]))^2)/100
  crps2simpar[[i]]<-mean(crps_sample(y = exp(logyt2sim[[i]])[,1], dat = exp(mcmcjags.2zone2sim[[i]]$logyhat[,1,,1])))+
    mean(crps_sample(y =exp(logyt2sim[[i]])[,2], dat = exp(mcmcjags.2zone2sim[[i]]$logyhat[,2,,1])))
}
0.5*(mean(unlist(lapply(cov502nsim, function(x) sum(x)/100)))+mean(unlist(lapply(cov502fsim, function(x) sum(x)/100))))
0.5*(mean(unlist(lapply(cov902nsim, function(x) sum(x)/100)))+mean(unlist(lapply(cov902fsim, function(x) sum(x)/100))))
mean(unlist(mse2nsim))+mean(unlist(mse2fsim))
mean(unlist(crps2simpar))
qprime.post.2zonesim<-gprime.post.2zonesim<-beta.post.2zonesim<-kl.post.2zonesim<-vector("list")
for(i in 1:100){
  qprime.post.2zonesim[[i]]<-m1.dat.2zonesim[[i]]$Qprime
  gprime.post.2zonesim[[i]]<-m1.dat.2zonesim[[i]]$Gprime
  beta.post.2zonesim[[i]]<-m1.dat.2zonesim[[i]]$beta
  kl.post.2zonesim[[i]]<-m1.dat.2zonesim[[i]]$kl
}
quantile(unlist(qprime.post.2zonesim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(gprime.post.2zonesim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(beta.post.2zonesim),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(kl.post.2zonesim),c(0.025,0.25,0.5,0.75,0.975))
#theta2julia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2v.csv")
logyt2pf<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/logyt2pf.csv")
theta2juliacorr<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2corr.csv")
apply(theta2juliacorr, 2, function(x){quantile(x)})
xjulia21<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21corr.csv")
xjulia22<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22corr.csv")
pred21.julia<-apply(xjulia21, 2, mean)
pred22.julia<-apply(xjulia22, 2, mean)
print(mse2pf<-(sum(((c22[,1])-(pred21.julia))^2)+sum(((c22[,2])-(pred22.julia))^2))/200)
print(mse2st<-(sum(((c22[,1])-(c1s.hat.2zone))^2)+sum(((c22[,2])-(c2s.hat.2zone))^2))/200)
print(mse2kf<-(sum(((c22[,1])-(predkf1))^2)+sum(((c22[,2])-(predkf2))^2))/200)
zrep21pf<-matrix(0,100,1000)
zrep22pf<-matrix(0,100,1000)
for(j in 1:1000){
  for(i in 1:100){
    zrep21pf[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21[j,i]),log(xjulia22[j,i])),
                         sigma=diag(c(theta2juliacorr$x9[j],theta2juliacorr$x12[j])),method="svd")[1]
    zrep22pf[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21[j,i]),log(xjulia22[j,i])),
                         sigma=diag(c(theta2juliacorr$x9[j],theta2juliacorr$x12[j])),method="svd")[2]
  }}
G2pf<-(exp(logyt2pf[,2])-apply(exp(zrep22pf),1,mean))^2+(exp(logyt2pf[,1])-apply(exp(zrep21pf),1,mean))^2
P2pf<-apply(exp(zrep22pf),1,var)
print(D2pf<-sum(G2pf)+sum(P2pf))
loglikpf2<-matrix(0,1000,100)
for(i in 1:100){
  for(j in 1:1000)
  loglikpf2[j,i]<-log(dmvnorm(logyt2[i,],cbind(log(xjulia21[j,i]+0.05),log(xjulia22[j,i])),
                            diag(c(theta2julia$x9[j],theta2julia$x12[j]))))
}
waic(loglikpf2[,2:100])
compare(waic(loglikpf2[,2:100]),waic(loglik.post2zone[,2:100]))
quartz()
par(mfrow=c(1,2))
plot(c22[,1],col="blue", main="State-by-state Gibbs Near", ylab="Concentrations", xlab="Time")
points(pred21.julia, col="red")
points(exp(logyt2[,1]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
plot(c22[,2],col="blue", main="State-by-state Gibbs Far", ylab="Concentrations", xlab="Time")
points(pred22.julia, col="red")
points(exp(logyt2[,2]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)
plot(cs12hat,c1.hat.2zone)
abline(c(0,1))
#all sim plots
quartz()
par(mfrow=c(4,2))
# plot(c1.hat.2zone, main="a Near",col="red",  ylab="Concentrations",
#      xlab="Time")
# lines(c22[,1],col="blue",type="l",lwd=3)
# points(exp(logyt2[,1]),col="green")
plot(exp(logyt2[,1]),col="green" ,main="a Near",  ylab="Concentrations",
     xlab="Time",type="l",lwd=3)
lines(c1.hat.2zone,col="red",lwd=3,lty=6)
lines(c22[,1],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"), lty=c(1,3,6),legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
# plot(c2.hat.2zone, col="red", main="a Far", ylab="Concentrations",
#      xlab="Time")
# lines(c22[,2],col="blue",type="l",lwd=3)
# points(exp(logyt2[,2]),col="green")
plot(exp(logyt2[,2]),col="green",main="a Far", ylab="Concentrations",
     xlab="Time",type="l",lwd=3)
lines(c2.hat.2zone, col="red",lwd=3,lty=6)
lines(c22[,2],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"), lty=c(1,3,6),legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)
# plot(pred21.julia, col="red", main="b Near",  ylab="Concentrations", 
#      xlab="Time")
# lines(c22[,1],col="blue",type="l",lwd=3)
# points(exp(logyt2[,1]),col="green")
plot(exp(logyt2[,1]),col="green", main="b Near",  ylab="Concentrations",
     xlab="Time",type="l",lwd=3)
lines(predkf1, col="red",lwd=3,lty=6)
lines(c22[,1],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
# plot(pred22.julia[2:100], col="red", main="b Far",  ylab="Concentrations", 
#      xlab="Time")
# lines(c22[2:100,2],col="blue",type="l",lwd=3)
# points(exp(logyt2[2:100,2]),col="green")
plot(exp(logyt2[,2]),col="green", main="b Far",  ylab="Concentrations",
     xlab="Time",type="l",lwd=3)
lines(predkf2, col="red",lwd=3,lty=6)
lines(c22[2:100,2],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)

# plot(exp(logyt2[,1]), col="green", main="c Near",  ylab="Concentrations", xlab="Time")
# points(predkf1,col="red")
# lines(c22[,1],col="blue",type="l",lwd=3)
plot(exp(logyt2[,1]), col="green", main="BNLR",  ylab="Concentrations", xlab="Time",type="l",lwd=3)
#lines(predkf1,col="red",lwd=3)
lines(yhat1p,col="red",lwd=3,lty=6)
lines(c22[,1],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
# plot(exp(logyt2[,2]), col="green", main="c Far", ylab="Concentrations", xlab="Time")
# points(predkf2,col="red")
# lines(c22[,2],col="blue",type="l",lwd=3)
plot(exp(logyt2[,2]), col="green", main="BNLR", ylab="Concentrations", xlab="Time",type="l",lwd=3)
#lines(predkf2,col="red",lwd=3)
lines(yhat2p,col="red",lwd=3,lty=6)
lines(c22[,2],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"), lty=c(1,3,6),legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)
# plot(cs12hat, col="purple", main="Smoothing Near",  ylab="Concentrations", 
#      xlab="Time")
# lines(c22[,1],col="blue",type="l",lwd=3)
# points(exp(logyt2[,1]),col="green")
plot(exp(logyt2[,1]),col="green", main="Smoothing Near",  ylab="Concentrations", 
     xlab="Time",type="l",lwd=3)
lines(c1.hat.2zone, col="purple",lwd=3,lty=6)
lines(c22[,1],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "purple"),lty=c(1,3,6), legend=c("Measurements N","True N",
                                                            "Smoothed N"),pch=15,cex=0.75,box.lwd = 0)
# plot(cs22hat, col="purple", main="Smoothing Far",  ylab="Concentrations", 
#      xlab="Time")
# lines(c22[,2],col="blue",type="l",lwd=3)
# points(exp(logyt2[,2]),col="green")
plot(exp(logyt2[,2]),col="green", main="Smoothing Far",  ylab="Concentrations", 
     xlab="Time",type="l",lwd=3)
lines(c2.hat.2zone, col="purple",lwd=3,lty=6)
lines(c22[,2],col="blue",type="l",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "purple"),lty=c(1,3,6), legend=c("Measurements F","True F",
                                                            "Smoothed F"),pch=15,cex=0.75,box.lwd = 0)

#nonlinear
load.module("msm")
cat("model{
    for (i in 1:N){
    c[i+1,1:2]<-mexp(h*i*A)%*%c[1,1:2]+B%*%(mexp(h*i*A)-I)%*%(g)
    logyt[i+1,1:2]~dmnorm(log(abs(c[i+1,1:2])), prec)
    loglik[i+1]<-(logdensity.mnorm(logyt[i+1,1:2],log(abs(c[i+1,1:2])),prec))
    }
    loglik[1]<-(logdensity.mnorm(logyt[1,1:2],log(c[1,1:2]),prec))
    A[1,1]<- -(beta)/VN
    A[1,2]<- beta/VN
    A[2,1]<- beta/VF
    A[2,2]<- -(beta+Qprime)/VF+kl
    B[1,1]<- (-(beta+Qprime)/VF +kl)/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    B[1,2]<- (-beta/VN)/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    B[2,1]<- -beta*VF/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    B[2,2]<- -((beta)/VN)/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    g[1,1]<- Gprime/VN
    g[2,1]<-0
    c[1,1]<-0.01
    c[1,2]<-0.5
    logyt[1,1:2]~dmnorm(log(c[1,1:2]), prec)
    prec~dwish(I,3)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    beta~dunif(0,10)
    }",file="modeljags2zonenl.jag")
#load.module("msm")
jagsnl2<-jags.model("modeljags2zonenl.jag",data=list("logyt"=logyt2[1:100,], N=99,"I"=diag(0.1,2), "VN"=pi*10^-3, "VF"=3.8,"h"=h))
update(jagsnl2, 1000)
mcmcjags.2zonenl<-jags.samples(jagsnl2,
                             c('Qprime','Gprime','beta',"kl",'c','loglik'),
                             1000)
mcmcsamplesjags.2zonenl<-coda.samples(jagsnl2,
                                    c('Qprime','Gprime','beta',"kl",'c','loglik'),
                                    1000)
m1.mcmc.2zonenl<-(as.mcmc(mcmcsamplesjags.2zonenl))
m1.mat.2zonenl<-as.matrix(mcmcsamplesjags.2zonenl)
m1.dat.2zonenl<-as.data.frame(m1.mat.2zonenl)
qprime.post.2zonenl<-m1.dat.2zonenl$Qprime
gprime.post.2zonenl<-m1.dat.2zonenl$Gprime
beta.post.2zonenl<-m1.dat.2zonenl$beta
kl.post.2zonenl<-m1.dat.2zonenl$kl
loglik.post2zonenl<-t(mcmcjags.2zonenl$loglik[,,1])
waic(loglik.post2zonenl)
quantile(qprime.post.2zonenl,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonenl <- apply(mcmcsamplesjags.2zonenl[[1]][,4:102], 2, mean)
c2.hat.2zonenl <- apply(mcmcsamplesjags.2zonenl[[1]][,103:202], 2, mean)

quartz()
par(mfrow=c(1,2))
plot(c2[,1],col="red")
points(c1.hat.2zonenl, col="blue")
points(exp(logyt2[,1]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtN","chatN","ytN"),bg="white", lty=1,cex=0.8)
plot(c2[,2],col="red")
points(c2.hat.2zonenl, col="blue")
points(exp(logyt2[,2]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtF","chatF","ytF"),bg="white", lty=1,cex=0.8)
#B2z
library(B2Z)
#sim
fit.2z <- B2ZM(data = cbind(seq(1,100,1),exp(logyt2[1:100,1]),exp(logyt2[1:100,2])), priorBeta = "unif(0,10)",
               indep.model = FALSE, priorQ = "unif(11,17)", priorG = "unif(281,482)",
               S=diag(2),v=3, VN =pi*10^-2, VF = 3.8, sampler = "MCMC",
               mcmc.control = list(NUpd = 10000, burnin = 1000, lag = 1, m = 5000))
summary(fit.2z)
attributes(fit.2z)
quartz()
plot(fit.2z)

c1p<-c2p<-y1p<-y2p<-matrix(0,100,1000)
Ap<-gp<-vector('list')
for(j in 1:1000){
  Ap[[j]]<-matrix(c(-(fit.2z$Beta[j+1000])/VN, (fit.2z$Beta[j+1000])/VN, fit.2z$Beta[j+1000]/VF, -(fit.2z$Beta[j+1000]+fit.2z$Q[j+1000])/VF),nrow=2, ncol=2,
                  byrow=TRUE)
  gp[[j]]<-matrix(c(fit.2z$G[j+1000]/VN,0),nrow=1, ncol=2)
  for(i in 1:(100-1)){
    c1p[1,]<-0.01
    c2p[1,]<-0.5
    y1p[1,]<-rmvnorm(n=1,mean=c(log(c1p[1,j]),log(c2p[1,j])),
                     sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[1]
    y2p[1,]<-rmvnorm(n=1,mean=c(log(c1p[1,j]),log(c2p[1,j])),
                     sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[2]
    c1p[i+1,j]<- (expm(i*Ap[[j]])%*%c(c1p[1,j],c2p[1,j])+solve(Ap[[j]])%*%(expm(i*Ap[[j]])-diag(2))%*%t(gp[[j]]))[1]
    c2p[i+1,j]<- (expm(i*Ap[[j]])%*%c(c1p[1,j],c2p[1,j])+solve(Ap[[j]])%*%(expm(i*Ap[[j]])-diag(2))%*%t(gp[[j]]))[2]
    y1p[i+1,j]<- rmvnorm(n=1,mean=c(log(c1p[i+1,j]),log(c2p[i+1,j])),
            sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[1]
    y2p[i+1,j]<- rmvnorm(n=1,mean=c(log(c1p[i+1,j]),log(c2p[i+1,j])),
                         sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[2]
  }}
yhat1p<-exp(apply(y1p,1,mean))
yhat2p<-exp(apply(y2p,1,mean))

print(mse2b2z<-(sum(((c22[,1])-(yhat1p))^2)+sum(((c22[,2])-(yhat2p))^2))/200)
G2b2z<-(c22[,1]-apply(exp(y2p),1,mean))^2+(c22[,2]-apply(exp(y1p),1,mean))^2
P2b2z<-apply(exp(y1p),1,var)+apply(exp(y2p),1,var)
print(D2b2z<-sum(G2b2z)+sum(P2b2z))

# c1p<-c2p<-matrix(0,100,1000)
# Ap<-gp<-vector('list')
# for(j in 1:1000){
#   Ap[[j]]<-matrix(c(-(fit.2z$Beta[j])/VN, (fit.2z$Beta[j])/VN, fit.2z$Beta[j]/VF, -(fit.2z$Beta[j]+fit.2z$Q[j])/VF),nrow=2, ncol=2,
#             byrow=TRUE)
#   gp[[j]]<-matrix(c(fit.2z$G[j]/VN,0),nrow=1, ncol=2)
#   for(i in 1:(100-1)){
#     c1p[1,]<-0.5
#     c2p[1,]<-0.5
#     c1p[i+1,j]<- (expm(h*i*Ap[[j]])%*%c(c1p[1,j],c2p[1,j])+solve(Ap[[j]])%*%(expm(h*i*Ap[[j]])-diag(2))%*%t(gp[[j]]))[1]
#     c2p[i+1,j]<- (expm(h*i*Ap[[j]])%*%c(c1p[1,j],c2p[1,j])+solve(Ap[[j]])%*%(expm(h*i*Ap[[j]])-diag(2))%*%t(gp[[j]]))[2]
#   }
#     zrep21b2z[,j]<-rnorm(n=100,log(c1p[,j]),sqrt(fit.2z$tauN[j]))
#     zrep22b2z[,j]<-rnorm(n=100,log(c2p[,j]),sqrt(fit.2z$tauF[j]))
# }
# zrep21b2z<-apply(zrep21b2z,1,function(x) ifelse(is.na(x),0,x))
# zrep22b2z<-apply(zrep22b2z,1,function(x) ifelse(is.na(x),0.5,x))
# G2b2z<-(exp(logyt2[,2])-apply(exp(zrep22b2z),2,mean))^2+(exp(logyt2[,1])-apply(exp(zrep21b2z),2,mean))^2
# P2b2z<-apply(exp(zrep22b2z),2,var)
# print(D2b2z<-sum(G2b2z)+sum(P2b2z))

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
write.csv(data2zonen,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonen.csv")
write.csv(data2zonef1,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonef1.csv")
write.csv(data2zonef2,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonef2.csv")
write.csv(data2zonef3,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonef3.csv")

#true values Q=0.059, 0.258 and 0.595 Vn=0.105 Vf=11.79 beta=0.24-1.24
quartz()
plot(data2zonef1[,10])
points(data2zonef1[,81])
cat("model{
    for (i in 1:N){
    logomega[i+1,1:2]~dmnorm(c(0,0),vari)
    c[i+1,1]<-c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h+exp(logomega[i+1,1])
    c[i+1,2]<-c[i,1]*cc*h+c[i,2]*((d)+1)+exp(logomega[i+1,2])
    logyt[i+1,1:2]~dmnorm(log(c[i+1,1:2]), prec)
    logyhat[i+1,1:2]~dmnorm(log(c[i+1,1:2]), prec)
    loglik[i+1]<-(logdensity.mnorm(logyt[i+1,1:2],log(c[i+1,1:2]),prec))
    }
    logyhat[1,1:2]~dmnorm(log(c[1,1:2]), prec)
    a<- -(beta)/VN
    b<-(beta)/VN
    cc<-beta/VF
    d<- -(beta+Qprime)/VF+kl
    g1<-Gprime/VN
    c[1,1]<-8
    c[1,2]<-5.8
    logyt[1,1:2]~dmnorm(log(c[1,1:2]), prec)
    prec~dwish(I,3)
    vari~dwish(I,3)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,0.8)
    beta~dunif(0,5)
    }",file="modeljags2zone.jag")
#med ventilation Q=0.23-0.27 med G=86.36
jagsm<-jags.model("modeljags2zone.jag",data=list("logyt"=cbind(log(data2zonen[4:75,13]),log(data2zonef1[4:75,13])),
                                                N=71,"I"=diag(2), "VN"=0.1, "VF"=11.9,"h"=0.01),
                 n.adapt =1000)
update(jagsm, 1000)
mcmcjags.2zonem<-jags.samples(jagsm,
                             c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                             1000)
mcmcsamplesjags.2zonem<-coda.samples(jagsm,
                                    c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                                    1000)

m1.mcmc.2zonem<-(as.mcmc(mcmcsamplesjags.2zonem))
m1.mat.2zonem<-as.matrix(mcmcsamplesjags.2zonem)
m1.dat.2zonem<-as.data.frame(m1.mat.2zonem)
qprime.post.2zonem<-m1.dat.2zonem$Qprime
gprime.post.2zonem<-m1.dat.2zonem$Gprime
beta.post.2zonem<-m1.dat.2zonem$beta
kl.post.2zonem<-m1.dat.2zonem$kl
vari11.post.2zonem<-m1.dat.2zonem$`prec[1,1]`
vari12.post.2zonem<-m1.dat.2zonem$`prec[1,2]`
vari21.post.2zonem<-m1.dat.2zonem$`prec[2,1]`
vari22.post.2zonem<-m1.dat.2zonem$`prec[2,2]`
quantile(qprime.post.2zonem,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonem,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonem,c(0.025,0.25,0.5,0.75,0.975))
loglik.post2zonem<-t(mcmcjags.2zonem$loglik[,,1])
dim(loglik.post2zonem)
waic(loglik.post2zonem[,2:75])
mean(crps_sample(y =data2zonen[4:75,13], dat = exp(mcmcjags.2zonem$logyhat[,1,,1])))+
  mean(crps_sample(y =data2zonef1[4:75,13], dat =exp(mcmcjags.2zonem$logyhat[,2,,1])))
yhatlower50npm<- apply( exp(mcmcjags.2zonem$logyhat[,1,,1]), 1, function(x) quantile(x,0.25))
yhatupper50npm<-apply( exp(mcmcjags.2zonem$logyhat[,1,,1]), 1, function(x) quantile(x,0.75))
yhatlower90npm<- apply( exp(mcmcjags.2zonem$logyhat[,1,,1]), 1, function(x) quantile(x,0.05))
yhatupper90npm<-apply( exp(mcmcjags.2zonem$logyhat[,1,,1]), 1, function(x) quantile(x,0.95))
yhatlower50fpm<- apply( exp(mcmcjags.2zonem$logyhat[,2,,1]), 1, function(x) quantile(x,0.25))
yhatupper50fpm<-apply( exp(mcmcjags.2zonem$logyhat[,2,,1]), 1, function(x) quantile(x,0.75))
yhatlower90fpm<- apply( exp(mcmcjags.2zonem$logyhat[,2,,1]), 1, function(x) quantile(x,0.05))
yhatupper90fpm<-apply( exp(mcmcjags.2zonem$logyhat[,2,,1]), 1, function(x) quantile(x,0.95))
#nominal coverage
cov502nm<-cov502fm<-cov902nm<-cov902fm<-NULL
for(i in 1:length(data2zonen[4:75,13])){
  cov502nm[i]<-(ifelse(yhatlower50npm[i]<=data2zonen[i+3,13]&data2zonen[i+3,13]<=yhatupper50npm[i],1,0))
  cov502fm[i]<-(ifelse(yhatlower50fpm[i]<=data2zonef1[i+3,13]&data2zonef1[i+3,13]<=yhatupper50fpm[i],1,0))
  cov902nm[i]<-(ifelse(yhatlower90npm[i]<=data2zonen[i+3,13]&data2zonen[i+3,13]<=yhatupper90npm[i],1,0))
  cov902fm[i]<-  (ifelse(yhatlower90fpm[i]<=data2zonef1[i+3,13]&data2zonef1[i+3,13]<=yhatupper90fpm[i],1,0))
}
sum(cov502nm)/144+sum(cov502fm)/144
sum(cov902nm)/144+sum(cov902fm)/144
c1.hat.2zonem <- apply(mcmcsamplesjags.2zonem[[1]][,4:(75)], 2, mean)
c2.hat.2zonem <- apply(mcmcsamplesjags.2zonem[[1]][,(76):(75+72)], 2, mean)
zrep21m<-matrix(0,72,1000)
zrep22m<-matrix(0,72,1000)
for(j in 1:1000){
  for(i in 1:72){
    zrep21m[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zonem[[1]][j,3+i]),log(mcmcsamplesjags.2zonem[[1]][j,74+i])),
                         sigma=solve(matrix(c(vari11.post.2zonem[j],vari12.post.2zonem[j],vari21.post.2zonem[j],
                                                vari22.post.2zonem[j]),2,2,byrow=TRUE)),method="svd")[1]
    zrep22m[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zonem[[1]][j,3+i]),log(mcmcsamplesjags.2zonem[[1]][j,74+i])),
                         sigma=solve(matrix(c(vari11.post.2zonem[j],vari12.post.2zonem[j],vari21.post.2zonem[j],
                                                vari22.post.2zonem[j]),2,2,byrow=TRUE)),method="svd")[2]
  }}
G2m<-(data2zonen[4:75,13]-apply(exp(zrep21m),1,mean))^2+(data2zonef1[4:75,13]-apply(exp(zrep22m),1,mean))^2
P2m<-apply(exp(zrep22m),1,var)+apply(exp(zrep21m),1,var)
print(D2m<-sum(G2m)+sum(P2m))
quartz()
plot(c2.hat.2zonem,data2zonef1[4:75,13])
abline(0,1)

quartz()
plot(data2zonen[,3])
points(data2zonef1[,3])
logomegam<-matrix(0,1000,2)
for(i in 1000){
  logomegam[,1:2]<-rlnorm.rplus(1,c(0,0),matrix(c(m1.dat.2zonem[i,233],m1.dat.2zonem[i,234],
                                                 m1.dat.2zonem[i,235],m1.dat.2zonem[i,236]),2,2,byrow=TRUE))[1,]
}
N=71
VN=0.105
VF=11.9
cs12m<-matrix(0,1000,N+1)
cs22m<-matrix(0,1000,N+1)
cs12m[,N+1]<-mcmcsamplesjags.2zonem[[1]][,75]
cs22m[,N+1]<-mcmcsamplesjags.2zonem[[1]][,147]
Asm<-array(rep(0,4000),dim=c(2,2,1000))
Bsm<-array(rep(0,4000),dim=c(2,2,1000))
gsm<-array(rep(0,2000),dim=c(1,2,1000))
for(i in 1:1000){
  for(j in 2:N){
    Asm[1,1,i]<- -(beta.post.2zonem[i])/VN
    Asm[1,2,i]<-(beta.post.2zonem[i])/VN
    Asm[2,1,i]<- beta.post.2zonem[i]/VF
    Asm[2,2,i]<- -(beta.post.2zonem[i]+qprime.post.2zonem[i])/VF
    gsm[1,1,i]<-gprime.post.2zonem[i]/VN
    gsm[1,2,i]<-0
    Bsm[1,1,i]<-(-h*(beta.post.2zonem[i]+qprime.post.2zonem[i])/VF+1)/(((h*Asm[1,1,i]+1)*(h*Asm[2,2,i]+1))-((h*Asm[2,1,i]*h*Asm[1,2,i])))
    Bsm[1,2,i]<-(-h*beta.post.2zonem[i]/VN)/(((h*Asm[1,1,i]+1)*(h*Asm[2,2,i]+1))-((h*Asm[2,1,i]*h*Asm[1,2,i])))
    Bsm[2,1,i]<- -(h*beta.post.2zonem[i]/VF)/(((h*Asm[1,1,i]+1)*(h*Asm[2,2,i]+1))-((h*Asm[2,1,i]*h*Asm[1,2,i])))
    Bsm[2,2,i]<- (-h*beta.post.2zonem[i]/VN+1)/(((h*Asm[1,1,i]+1)*(h*Asm[2,2,i]+1))-((h*Asm[2,1,i]*h*Asm[1,2,i])))
    cs12m[i,(j+N-2*(j-1))]<- ((Bsm[,,i])%*%t(c(mcmcsamplesjags.2zonem[[1]][i,(j+N-2*(j-1)+1+3)],mcmcsamplesjags.2zonem[[1]][i,(j+N-2*(j-1)+1+75)])-t(h*gsm[1,1:2,i])-
                                             t(c((h*logomegam[(j+N-2*(j-1)+1),1]),(h*logomegam[(j+N-2*(j-1)+1),2])))))[1,1]
    cs22m[i,(j+N-2*(j-1))]<- ((Bsm[,,i])%*%t(c(mcmcsamplesjags.2zonem[[1]][i,(j+N-2*(j-1)+1+3)],mcmcsamplesjags.2zonem[[1]][i,(j+N-2*(j-1)+1+75)])-t(h*gsm[1,1:2,i])-
                                             t(c((h*logomegam[(j+N-2*(j-1)+1),1]),(h*logomegam[(j+N-2*(j-1)+1),2])))))[2,1]
  } 
  cs12m[i,2]<-((Bsm[,,i])%*%t(c(mcmcsamplesjags.2zonem[[1]][i,5],mcmcsamplesjags.2zonem[[1]][i,77])-t(h*gsm[1,1:2,i])-
                              t(c((logomegam[2,1]),(logomegam[2,2])))))[1,1]
  cs22m[i,2]<-((Bsm[,,i])%*%t(c(mcmcsamplesjags.2zonem[[1]][i,5],mcmcsamplesjags.2zonem[[1]][i,77])-t(h*gsm[1,1:2,i])-
                              t(c((logomegam[2,1]),(logomegam[2,2])))))[1,1]
}

cs12hatm<-apply(cs12m,2,mean)
cs22hatm<-apply(cs22m,2,mean)
plot(cs12hatm)
points(cs22hatm)
#high
jagsh<-jags.model("modeljags2zone.jag",data=list("logyt"=cbind(log(data2zonen[,81]),log(data2zonef1[,81])),
                                                 N=40,"I"=diag(2), "VN"=0.1, "VF"=11.9,"h"=0.01),
                  n.adapt =1000)
update(jagsh, 2000)
mcmcjags.2zoneh<-jags.samples(jagsh,
                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                              1000)
mcmcsamplesjags.2zoneh<-coda.samples(jagsh,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                                     1000)

m1.mcmc.2zoneh<-(as.mcmc(mcmcsamplesjags.2zoneh))
m1.mat.2zoneh<-as.matrix(mcmcsamplesjags.2zoneh)
m1.dat.2zoneh<-as.data.frame(m1.mat.2zoneh)
qprime.post.2zoneh<-m1.dat.2zoneh$Qprime
gprime.post.2zoneh<-m1.dat.2zoneh$Gprime
beta.post.2zoneh<-m1.dat.2zoneh$beta
kl.post.2zoneh<-m1.dat.2zoneh$kl
vari11.post.2zoneh<-m1.dat.2zoneh$`prec[1,1]`
vari12.post.2zoneh<-m1.dat.2zoneh$`prec[1,2]`
vari21.post.2zoneh<-m1.dat.2zoneh$`prec[2,1]`
vari22.post.2zoneh<-m1.dat.2zoneh$`prec[2,2]`
quantile(qprime.post.2zoneh,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zoneh,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zoneh,c(0.025,0.25,0.5,0.75,0.975))
loglik.post2zoneh<-t(mcmcjags.2zoneh$loglik[,,1])
dim(loglik.post2zoneh)
waic(loglik.post2zoneh[,2:41])
mean(crps_sample(y =data2zonen[1:41,81], dat = exp(mcmcjags.2zoneh$logyhat[,1,,1])))+
  mean(crps_sample(y =data2zonef1[1:41,81], dat =exp(mcmcjags.2zoneh$logyhat[,2,,1])))
yhatlower50nph<- apply( exp(mcmcjags.2zoneh$logyhat[,1,,1]), 1, function(x) quantile(x,0.25))
yhatupper50nph<-apply( exp(mcmcjags.2zoneh$logyhat[,1,,1]), 1, function(x) quantile(x,0.75))
yhatlower90nph<- apply( exp(mcmcjags.2zoneh$logyhat[,1,,1]), 1, function(x) quantile(x,0.05))
yhatupper90nph<-apply( exp(mcmcjags.2zoneh$logyhat[,1,,1]), 1, function(x) quantile(x,0.95))
yhatlower50fph<- apply( exp(mcmcjags.2zoneh$logyhat[,2,,1]), 1, function(x) quantile(x,0.25))
yhatupper50fph<-apply( exp(mcmcjags.2zoneh$logyhat[,2,,1]), 1, function(x) quantile(x,0.75))
yhatlower90fph<- apply( exp(mcmcjags.2zoneh$logyhat[,2,,1]), 1, function(x) quantile(x,0.05))
yhatupper90fph<-apply( exp(mcmcjags.2zoneh$logyhat[,2,,1]), 1, function(x) quantile(x,0.95))
#nominal coverage
cov502nh<-cov502fh<-cov902nh<-cov902fh<-NULL
for(i in 1:length(data2zonen[1:41,81])){
  cov502nh[i]<-(ifelse(yhatlower50nph[i]<=data2zonen[i,81]&data2zonen[i,81]<=yhatupper50nph[i],1,0))
  cov502fh[i]<-(ifelse(yhatlower50fph[i]<=data2zonef1[i,81]&data2zonef1[i,81]<=yhatupper50fph[i],1,0))
  cov902nh[i]<-(ifelse(yhatlower90nph[i]<=data2zonen[i,81]&data2zonen[i,81]<=yhatupper90nph[i],1,0))
  cov902fh[i]<-  (ifelse(yhatlower90fph[i]<=data2zonef1[i,81]&data2zonef1[i,81]<=yhatupper90fph[i],1,0))
}
sum(cov502nh)/82+sum(cov502fh)/82
sum(cov902nh)/82+sum(cov902fh)/82
c1.hat.2zoneh <- apply(mcmcsamplesjags.2zoneh[[1]][,4:(44)], 2, mean)
c2.hat.2zoneh <- apply(mcmcsamplesjags.2zoneh[[1]][,(45):85], 2, mean)
zrep21h<-matrix(0,41,1000)
zrep22h<-matrix(0,41,1000)
for(j in 1:1000){
  for(i in 1:41){
    zrep21h[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zoneh[[1]][j,3+i]),log(mcmcsamplesjags.2zoneh[[1]][j,44+i])),
                          sigma=solve(matrix(c(vari11.post.2zoneh[j],vari12.post.2zoneh[j],vari21.post.2zoneh[j],
                                                       vari22.post.2zoneh[j]),2,2,byrow=TRUE)),method="svd")[1]
    zrep22h[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zoneh[[1]][j,3+i]),log(mcmcsamplesjags.2zoneh[[1]][j,44+i])),
                          sigma=solve(matrix(c(vari11.post.2zoneh[j],vari12.post.2zoneh[j],vari21.post.2zoneh[j],
                                                       vari22.post.2zoneh[j]),2,2,byrow=TRUE)),method="svd")[2]
  }}
G2h<-(data2zonen[1:41,81]-apply(exp(zrep21h),1,mean))^2+(data2zonef1[1:41,81]-apply(exp(zrep22h),1,mean))^2
P2h<-apply(exp(zrep22h),1,var)+apply(exp(zrep21h),1,var)
print(D2h<-sum(G2h)+sum(P2h))
logomegah<-matrix(0,1000,2)
for(i in 1000){
  logomegah[,1:2]<-rlnorm.rplus(1,c(0,0),matrix(c(m1.dat.2zoneh[i,124],m1.dat.2zoneh[i,125],
                                                  m1.dat.2zoneh[i,126],m1.dat.2zoneh[i,127]),2,2,byrow=TRUE))[1,]
}
N=40
VN=0.105
VF=11.9
cs12h<-matrix(0,1000,N+1)
cs22h<-matrix(0,1000,N+1)
cs12h[,N+1]<-mcmcsamplesjags.2zoneh[[1]][,44]
cs22h[,N+1]<-mcmcsamplesjags.2zoneh[[1]][,85]
Ash<-array(rep(0,4000),dim=c(2,2,1000))
Bsh<-array(rep(0,4000),dim=c(2,2,1000))
gsh<-array(rep(0,2000),dim=c(1,2,1000))
for(i in 1:1000){
  for(j in 2:N){
    Ash[1,1,i]<- -(beta.post.2zoneh[i])/VN
    Ash[1,2,i]<-(beta.post.2zoneh[i])/VN
    Ash[2,1,i]<- beta.post.2zoneh[i]/VF
    Ash[2,2,i]<- -(beta.post.2zoneh[i]+qprime.post.2zoneh[i])/VF
    gsh[1,1,i]<-gprime.post.2zoneh[i]/VN
    gsh[1,2,i]<-0
    Bsh[1,1,i]<-(-h*(beta.post.2zoneh[i]+qprime.post.2zoneh[i])/VF+1)/(((h*Ash[1,1,i]+1)*(h*Ash[2,2,i]+1))-((h*Ash[2,1,i]*h*Ash[1,2,i])))
    Bsh[1,2,i]<-(-h*beta.post.2zoneh[i]/VN)/(((h*Ash[1,1,i]+1)*(h*Ash[2,2,i]+1))-((h*Ash[2,1,i]*h*Ash[1,2,i])))
    Bsh[2,1,i]<- -(h*beta.post.2zoneh[i]/VF)/(((h*Ash[1,1,i]+1)*(h*Ash[2,2,i]+1))-((h*Ash[2,1,i]*h*Ash[1,2,i])))
    Bsh[2,2,i]<- (-h*beta.post.2zoneh[i]/VN+1)/(((h*Ash[1,1,i]+1)*(h*Ash[2,2,i]+1))-((h*Asm[2,1,i]*h*Ash[1,2,i])))
    cs12h[i,(j+N-2*(j-1))]<- ((Bsh[,,i])%*%t(c(mcmcsamplesjags.2zoneh[[1]][i,(j+N-2*(j-1)+1+3)],mcmcsamplesjags.2zoneh[[1]][i,(j+N-2*(j-1)+1+44)])-t(h*gsh[1,1:2,i])-
                                               t(c((h*logomegam[(j+N-2*(j-1)+1),1]),(h*logomegah[(j+N-2*(j-1)+1),2])))))[1,1]
    cs22h[i,(j+N-2*(j-1))]<- ((Bsh[,,i])%*%t(c(mcmcsamplesjags.2zoneh[[1]][i,(j+N-2*(j-1)+1+3)],mcmcsamplesjags.2zoneh[[1]][i,(j+N-2*(j-1)+1+44)])-t(h*gsh[1,1:2,i])-
                                               t(c((h*logomegah[(j+N-2*(j-1)+1),1]),(h*logomegah[(j+N-2*(j-1)+1),2])))))[2,1]
  } 
  cs12h[i,2]<-((Bsh[,,i])%*%t(c(mcmcsamplesjags.2zoneh[[1]][i,5],mcmcsamplesjags.2zoneh[[1]][i,46])-t(h*gsh[1,1:2,i])-
                                t(c((logomegah[2,1]),(logomegah[2,2])))))[1,1]
  cs22h[i,2]<-((Bsh[,,i])%*%t(c(mcmcsamplesjags.2zoneh[[1]][i,5],mcmcsamplesjags.2zoneh[[1]][i,46])-t(h*gsh[1,1:2,i])-
                                t(c((logomegah[2,1]),(logomegah[2,2])))))[1,1]
}

cs12hath<-apply(cs12h,2,mean)
cs22hath<-apply(cs22h,2,mean)
plot(cs12hath)
points(cs22hath)
#low Q=0.056
jagsl<-jags.model("modeljags2zone.jag",data=list("logyt"=cbind(log(data2zonen[1:176,1]),log(data2zonef1[1:176,1]+0.05)),
                                                 N=175,"I"=diag(2), "VN"=0.1, "VF"=11.9,"h"=0.01))
update(jagsl, 1000)
mcmcjags.2zonel<-jags.samples(jagsl,
                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                              1000)
mcmcsamplesjags.2zonel<-coda.samples(jagsl,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec","vari",'loglik','logyhat'),
                                     1000)

quartz()
plot(mcmcsamplesjags.2zonel)
summary(mcmcsamplesjags.2zonel[[1]])
m1.mcmc.2zonel<-(as.mcmc(mcmcsamplesjags.2zonel))
m1.mat.2zonel<-as.matrix(mcmcsamplesjags.2zonel)
m1.dat.2zonel<-as.data.frame(m1.mat.2zonel)
qprime.post.2zonel<-m1.dat.2zonel$Qprime
gprime.post.2zonel<-m1.dat.2zonel$Gprime
beta.post.2zonel<-m1.dat.2zonel$beta
kl.post.2zonel<-m1.dat.2zonel$kl
vari11.post.2zonel<-m1.dat.2zonel$`prec[1,1]`
vari12.post.2zonel<-m1.dat.2zonel$`prec[1,2]`
vari21.post.2zonel<-m1.dat.2zonel$`prec[2,1]`
vari22.post.2zonel<-m1.dat.2zonel$`prec[2,2]`
quantile(qprime.post.2zonel,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonel,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonel,c(0.025,0.25,0.5,0.75,0.975))
loglik.post2zonel<-t(mcmcjags.2zonel$loglik[,,1])
dim(loglik.post2zonel)
waic(loglik.post2zonel[,2:176])
mean(crps_sample(y =data2zonen[1:176,1], dat = exp(mcmcjags.2zonel$logyhat[,1,,1])))+
  mean(crps_sample(y =data2zonef1[1:176,1], dat =exp(mcmcjags.2zonel$logyhat[,2,,1])))
yhatlower50npl<- apply( exp(mcmcjags.2zonel$logyhat[,1,,1]), 1, function(x) quantile(x,0.25))
yhatupper50npl<-apply( exp(mcmcjags.2zonel$logyhat[,1,,1]), 1, function(x) quantile(x,0.75))
yhatlower90npl<- apply( exp(mcmcjags.2zonel$logyhat[,1,,1]), 1, function(x) quantile(x,0.05))
yhatupper90npl<-apply( exp(mcmcjags.2zonel$logyhat[,1,,1]), 1, function(x) quantile(x,0.95))
yhatlower50fpl<- apply( exp(mcmcjags.2zonel$logyhat[,2,,1]), 1, function(x) quantile(x,0.25))
yhatupper50fpl<-apply( exp(mcmcjags.2zonel$logyhat[,2,,1]), 1, function(x) quantile(x,0.75))
yhatlower90fpl<- apply( exp(mcmcjags.2zonel$logyhat[,2,,1]), 1, function(x) quantile(x,0.05))
yhatupper90fpl<-apply( exp(mcmcjags.2zonel$logyhat[,2,,1]), 1, function(x) quantile(x,0.95))
#nominal coverage
cov502nl<-cov502fl<-cov902nl<-cov902fl<-NULL
for(i in 1:length(data2zonen[1:176,1])){
  cov502nl[i]<-(ifelse(yhatlower50npl[i]<=data2zonen[i,1]&data2zonen[i,1]<=yhatupper50npl[i],1,0))
  cov502fl[i]<-(ifelse(yhatlower50fpl[i]<=data2zonef1[i,1]&data2zonef1[i,1]<=yhatupper50fpl[i],1,0))
  cov902nl[i]<-(ifelse(yhatlower90npl[i]<=data2zonen[i,1]&data2zonen[i,1]<=yhatupper90npl[i],1,0))
  cov902fl[i]<-  (ifelse(yhatlower90fpl[i]<=data2zonef1[i,1]&data2zonef1[i,1]<=yhatupper90fpl[i],1,0))
}
sum(cov502nl)/(176*2)+sum(cov502fl)/(176*2)
sum(cov902nl)/(176*2)+sum(cov902fl)/(176*2)
c1.hat.2zonel <- apply(mcmcsamplesjags.2zonel[[1]][,4:(103+76)], 2, mean)
c2.hat.2zonel <- apply(mcmcsamplesjags.2zonel[[1]][,(104+76):355], 2, mean)
zrep21l<-matrix(0,176,1000)
zrep22l<-matrix(0,176,1000)
for(j in 1:1000){
  for(i in 1:176){
    zrep21l[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zonel[[1]][j,3+i]),log(mcmcsamplesjags.2zonel[[1]][j,179+i])),
                          sigma=solve(matrix(c(vari11.post.2zonel[j],vari12.post.2zonel[j],vari21.post.2zonel[j],
                                               vari22.post.2zonel[j]),2,2,byrow=TRUE)),method="svd")[1]
    zrep22l[i,j]<-rmvnorm(n=1,mean=c(log(mcmcsamplesjags.2zonel[[1]][j,3+i]),log(mcmcsamplesjags.2zonel[[1]][j,179+i])),
                          sigma=solve(matrix(c(vari11.post.2zonel[j],vari12.post.2zonel[j],vari21.post.2zonel[j],
                                               vari22.post.2zonel[j]),2,2,byrow=TRUE)),method="svd")[2]
  }}
G2l<-(data2zonen[1:176,1]-apply(exp(zrep21l),1,mean))^2+(data2zonef1[1:176,1]-apply(exp(zrep22l),1,mean))^2
P2l<-apply(exp(zrep22l),1,var)+apply(exp(zrep21l),1,var)
print(D2l<-sum(G2l)+sum(P2l))
quartz()
plot(c2.hat.2zonel,data2zonef1[1:176,1])
abline(0,1)
logomegal<-matrix(0,1000,2)
for(i in 1000){
  logomegal[,1:2]<-rlnorm.rplus(1,c(0,0),matrix(c(m1.dat.2zonel[i,394],m1.dat.2zonel[i,395],
                                                  m1.dat.2zonel[i,396],m1.dat.2zonel[i,397]),2,2,byrow=TRUE))[1,]
}
N=175
VN=0.105
VF=11.9
cs12l<-matrix(0,1000,N+1)
cs22l<-matrix(0,1000,N+1)
cs12l[,N+1]<-mcmcsamplesjags.2zonel[[1]][1:1000,180]
cs22l[,N+1]<-mcmcsamplesjags.2zonel[[1]][1:1000,356]
Asl<-array(rep(0,4000),dim=c(2,2,1000))
Bsl<-array(rep(0,4000),dim=c(2,2,1000))
gsl<-array(rep(0,2000),dim=c(1,2,1000))
for(i in 1:1000){
  for(j in 2:N){
    Asl[1,1,i]<- -(beta.post.2zonel[i])/VN
    Asl[1,2,i]<-(beta.post.2zonel[i])/VN
    Asl[2,1,i]<- beta.post.2zonel[i]/VF
    Asl[2,2,i]<- -(beta.post.2zonel[i]+qprime.post.2zonel[i])/VF
    gsl[1,1,i]<-gprime.post.2zonel[i]/VN
    gsl[1,2,i]<-0
    Bsl[1,1,i]<-(-h*(beta.post.2zonel[i]+qprime.post.2zonel[i])/VF+1)/(((h*Asl[1,1,i]+1)*(h*Asl[2,2,i]+1))-((h*Asl[2,1,i]*h*Asl[1,2,i])))
    Bsl[1,2,i]<-(-h*beta.post.2zonel[i]/VN)/(((h*Asl[1,1,i]+1)*(h*Asl[2,2,i]+1))-((h*Asl[2,1,i]*h*Asl[1,2,i])))
    Bsl[2,1,i]<- -(h*beta.post.2zonel[i]/VF)/(((h*Asl[1,1,i]+1)*(h*Asl[2,2,i]+1))-((h*Asl[2,1,i]*h*Asl[1,2,i])))
    Bsl[2,2,i]<- (-h*beta.post.2zonel[i]/VN+1)/(((h*Asl[1,1,i]+1)*(h*Asl[2,2,i]+1))-((h*Asl[2,1,i]*h*Asl[1,2,i])))
    cs12l[i,(j+N-2*(j-1))]<- ((Bsl[,,i])%*%t(c(mcmcsamplesjags.2zonel[[1]][i,(j+N-2*(j-1)+1+4)],mcmcsamplesjags.2zonel[[1]][i,(j+N-2*(j-1)+1+180)])-t(h*gsl[1,1:2,i])-
                                               t(c((h*logomegal[(j+N-2*(j-1)+1),1]),(h*logomegal[(j+N-2*(j-1)+1),2])))))[1,1]
    cs22l[i,(j+N-2*(j-1))]<- ((Bsh[,,i])%*%t(c(mcmcsamplesjags.2zonel[[1]][i,(j+N-2*(j-1)+1+4)],mcmcsamplesjags.2zonel[[1]][i,(j+N-2*(j-1)+1+180)])-t(h*gsl[1,1:2,i])-
                                               t(c((h*logomegal[(j+N-2*(j-1)+1),1]),(h*logomegal[(j+N-2*(j-1)+1),2])))))[2,1]
  } 
  cs12l[i,2]<-((Bsl[,,i])%*%t(c(mcmcsamplesjags.2zonel[[1]][i,6],mcmcsamplesjags.2zonel[[1]][i,182])-t(h*gsl[1,1:2,i])-
                                t(c((logomegal[2,1]),(logomegal[2,2])))))[1,1]
  cs22l[i,2]<-((Bsl[,,i])%*%t(c(mcmcsamplesjags.2zonel[[1]][i,6],mcmcsamplesjags.2zonel[[1]][i,182])-t(h*gsl[1,1:2,i])-
                                t(c((logomegal[2,1]),(logomegal[2,2])))))[1,1]
}

cs12hatl<-apply(cs12l,2,mean)
cs22hatl<-apply(cs22l,2,mean)
plot(cs12hatl)
points(cs22hatl)
quartz()
par(mfrow=c(3,2))
plot(c1.hat.2zonel, col="blue",main="Low")
points(data2zonen[1:176,1],col="red")
points(cs12hatl,col="purple",cex=0.5)
legend("bottomright",col=c("blue","red","purple"),legend=c("chat near","measured near","smoothed near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonel, col="blue")
points(data2zonef1[1:176,1],col="red")
points(cs22hatl, col="purple",cex=0.5)
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far","smoothed far"),
       bg="white", lty=1,cex=0.8)

plot(c1.hat.2zonem, col="blue",main="Medium")
points(data2zonen[1:75,13],col="red")
points(cs12hatm, col="purple",cex=0.5)
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near","smoothed near" ),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonem, col="blue")
points(data2zonef1[1:75,13],col="red")
points(cs22hatm, col="purple",cex=0.5)
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far","smoothed far"),
       bg="white", lty=1,cex=0.8)

plot(c1.hat.2zoneh, col="blue",main="High")
points(data2zonen[1:41,81],col="red")
points(cs12hath, col="purple",cex=0.5)
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near","smoothed near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zoneh, col="blue")
points(data2zonef1[1:41,81],col="red")
points(cs22hath, col="purple",cex=0.5)
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far","smoothed far"),
       bg="white", lty=1,cex=0.8)

#PF
theta2julial<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2l.csv")
apply(theta2julial, 2, function(x){quantile(x)})
xjulia21l<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21l.csv")
xjulia22l<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22l.csv")
pred21.julial<-apply(xjulia21l[1000:5000,], 2, mean)
pred22.julial<-apply(xjulia22l[1000:5000,], 2, mean)
zrep21pfl<-matrix(0,176,1000)
zrep22pfl<-matrix(0,176,1000)
for(j in 1:1000){
  for(i in 1:176){
    zrep21pfl[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21l[j+1000,i]),log(xjulia22l[j+1000,i])),
                           sigma=diag(c(theta2julial$x9[j],theta2julial$x12[j])),method="svd")[1]
    zrep22pfl[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21l[j+1000,i]),log(xjulia22l[j+1000,i])),
                           sigma=diag(c(theta2julial$x9[j],theta2julial$x12[j])),method="svd")[2]
  }}
G2pfl<-(data2zonef1[1:176,1]-apply(exp(zrep22pfl),1,mean))^2+(data2zonen[1:176,1]-apply(exp(zrep21pfl),1,mean))^2
P2pfl<-apply(exp(zrep21pfl),1,var)+apply(exp(zrep22pfl),1,var)
print(D2pfl<-sum(G2pfl)+sum(P2pfl))
loglikpf2l<-matrix(0,1000,176)
for(i in 1:176){
  for(j in 1:1000)
    loglikpf2l[j,i]<-log(dmvnorm(cbind(log(data2zonen[i,1]),log(data2zonef1[i,1])),cbind(log(xjulia21l[j+1000,i]+0.05),log(xjulia22l[j+1000,i])),
                                diag(c(theta2julial$x9[j
                                                       +1000],theta2julial$x12[j+1000])))+0.05)
}
waic(loglikpf2l[,2:176])
compare(waic(loglikpf2l[,2:176]),waic(loglik.post2zonel[,2:176]))
#pf better
theta2juliah<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2m.csv")
apply(theta2juliah, 2, function(x){quantile(x)})
xjulia21h<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21m.csv")
xjulia22h<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22m.csv")
pred21.juliah<-apply(xjulia21h, 2, mean)
pred22.juliah<-apply(xjulia22h, 2, mean)
zrep21pfh<-matrix(0,41,1000)
zrep22pfh<-matrix(0,41,1000)
for(j in 1:1000){
  for(i in 1:41){
    zrep21pfh[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21h[j,i]),log(xjulia22h[j,i])),
                            sigma=diag(c(theta2juliah$x9[j],theta2juliah$x12[j])),method="svd")[1]
    zrep22pfh[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21h[j,i]),log(xjulia22h[j,i])),
                            sigma=diag(c(theta2juliah$x9[j],theta2juliah$x12[j])),method="svd")[2]
  }}
G2pfh<-(data2zonef1[1:41,81]-apply(exp(zrep22pfh),1,mean))^2+(data2zonen[1:41,81]-apply(exp(zrep21pfh),1,mean))^2
P2pfh<-apply(exp(zrep21pfh),1,var)+apply(exp(zrep22pfh),1,var)
print(D2pfh<-sum(G2pfh)+sum(P2pfh))
loglikpf2h<-matrix(0,1000,41)
for(i in 1:41){
  for(j in 1:1000)
    loglikpf2h[j,i]<-log(dmvnorm(cbind(log(data2zonen[i,81]),log(data2zonef1[i,81])),cbind(log(xjulia21h[j+1000,i]+0.05),log(xjulia22h[j+1000,i])),
                                 diag(c(theta2juliah$x9[j+1000],theta2juliah$x12[j+1000])))+0.05)
}
waic(loglikpf2h[,1:41])
#stbs better predictive accuracy
theta2juliam<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2h.csv")
apply(theta2juliam, 2, function(x){quantile(x)})
xjulia21m<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21h.csv")
xjulia22m<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22h.csv")
pred21.juliam<-apply(xjulia21m, 2, mean)
pred22.juliam<-apply(xjulia22m, 2, mean)
zrep21pfm<-matrix(0,75,1000)
zrep22pfm<-matrix(0,75,1000)
for(j in 1:1000){
  for(i in 1:75){
    zrep21pfm[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21m[j,i]),log(xjulia22m[j,i])),
                            sigma=diag(c(theta2juliam$x9[j],theta2juliam$x12[j])),method="svd")[1]
    zrep22pfm[i,j]<-rmvnorm(n=1,mean=c(log(xjulia21m[j,i]),log(xjulia22m[j,i])),
                            sigma=diag(c(theta2juliam$x9[j],theta2juliam$x12[j])),method="svd")[2]
  }}
G2pfm<-(data2zonef1[1:75,13]-apply(exp(zrep22pfm),1,mean))^2+(data2zonen[1:75,13]-apply(exp(zrep21pfm),1,mean))^2
P2pfm<-apply(exp(zrep21pfm),1,var)+apply(exp(zrep22pfm),1,var)
print(D2pfm<-sum(G2pfm)+sum(P2pfm))
print(mse2pfl<-(sum(((data2zonen[1:176,1])-(pred21.julial))^2)+sum((data2zonef1[1:176,1]-(pred22.julial))^2))/302)
print(mse2stl<-(sum((data2zonen[1:176,1]-(c1.hat.2zonel))^2)+sum((data2zonef1[1:176,1]-(c2.hat.2zonel))^2))/302)
print(mse2pfh<-(sum(((data2zonen[1:41,81])-(pred21.juliah))^2)+sum((data2zonef1[1:41,81]-(pred22.juliah))^2))/82)
print(mse2sth<-(sum((data2zonen[1:41,81]-(c1.hat.2zoneh))^2)+sum((data2zonef1[1:41,81]-(c2.hat.2zoneh))^2))/82)
print(mse2pfm<-(sum(((data2zonen[1:75,13])-(pred21.juliam))^2)+sum((data2zonef1[1:75,13]-(pred22.juliam))^2))/150)
print(mse2stm<-(sum((data2zonen[4:75,13]-(c1.hat.2zonem))^2)+sum((data2zonef1[4:75,13]-(c2.hat.2zonem))^2))/144)

loglikpf2m<-matrix(0,1000,75)
for(i in 1:75){
  for(j in 1:1000)
    loglikpf2m[j,i]<-log(dmvnorm(cbind(log(data2zonen[i,13]),log(data2zonef1[i,13])),cbind(log(xjulia21m[j+1000,i]),log(xjulia22m[j+1000,i])),
                                 diag(c(theta2juliam$x9[j+1000],theta2juliam$x12[j+1000])))+0.05)
}
waic(loglikpf2m[,3:75])
compare(waic(loglikpf2m[,3:75]),waic(loglik.post2zonem[,3:75]))
#pf has higher predictive accuracy
quartz()
par(mfrow=c(3,2))
plot(pred21.julial, col="blue",main="Low")
points(data2zonen[1:176,1],col="red")
plot(pred22.julial, col="blue")
points(data2zonef1[1:176,1],col="red")
plot(pred21.juliam, col="blue",main="Medium")
points(data2zonen[1:41,81],col="red")
plot(pred22.juliam, col="blue")
points(data2zonef1[1:41,81],col="red")
plot(pred21.juliah, col="blue",main="High")
points(data2zonen[1:75,13],col="red")
plot(pred22.juliah, col="blue")
points(data2zonef1[1:75,13],col="red")
#plot all data
quartz()
par(mfrow=c(3,2))
# plot(c1.hat.2zonel, col="blue",main="Low Near",ylab="Toluene Concentrations N",xlab="Time")
# points(pred21.julial, col="red")
# points(data2zonen[1:176,1],col="green",cex=0.8)
# points(cs12hatl,col="purple",cex=0.5)

plot(data2zonen[1:176,1],col="green",lwd=3,main="Low Near",ylab="Toluene Concentrations N",xlab="Time",type="l")
#lines(pred21.julial, col="red",lwd=3)
lines(yhat1pl, col="red",lwd=3,lty=5)
lines(c1.hat.2zonel[1:175], col="blue",lwd=3,    lty=6)
lines(predkf1l,col="chocolate",lwd=3,lty=2)
lines(cs12hatl[1:173],col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue","chocolate","red","purple"),lty=c(1,6,2,5,3), legend=c("Measurements N","Model a", "Model b","BNLR","Smoothed"),pch=15,cex=0.75,box.lwd = 0)

# plot(c2.hat.2zonel, col="blue",main="Low Far",ylab="Toluene Concentrations F",xlab="Time")
# points(pred22.julial, col="red")
# points(data2zonef1[1:176,1],col="green",cex=0.8)
# points(cs22hatl,col="purple",cex=0.5)
plot(data2zonef1[1:176,1],col="green",lwd=3,main="Low Far",ylab="Toluene Concentrations F",xlab="Time",type="l")
#lines(pred22.julial, col="red",lwd=3)
lines(yhat2pl, col="red",lwd=3,lty=5)
lines(c2.hat.2zonel[1:175],lty=6, col="blue",lwd=3)
lines(predkf2l,col="chocolate",lwd=3,lty=2)
lines(cs22hatl[1:173],col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue", "chocolate","red","purple"), lty=c(1,6,2,5,3),legend=c("Measurements F","Model a", "Model b","BNLR","Smoothed"),pch=15,cex=0.75,box.lwd = 0)


# plot(data2zonen[1:75,13],col="green",main="Medium Near",ylab="Toluene Concentrations N",xlab="Time",cex=0.8)
# points(c1.hat.2zonem, col="blue")
# points(pred21.juliam, col="red")
# points(cs12hatm,col="purple",cex=0.5)
plot(data2zonen[4:75,13],col="green",main="Medium Near",ylab="Toluene Concentrations N",xlab="Time",type="l",lwd=3)
lines(c1.hat.2zonem, col="blue",lwd=3,lty=6)
lines(yhat1pm, col="red",lwd=3,lty=5)
#lines(pred21.juliam, col="red",lwd=3)
lines(predkf1h,col="chocolate",lwd=3,lty=2)
lines(cs12hatm,col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue", "chocolate","red","purple"), lty=c(1,6,2,5,3),legend=c("Measurements N","Model a", "Model b","BNLR","Smoothed"),pch=15,cex=0.75,box.lwd = 0)

# plot(c2.hat.2zonem, col="blue",main="Medium Far",ylab="Toluene Concentrations F",xlab="Time")
# points(pred22.juliam, col="red")
# points(data2zonef1[1:75,13],col="green",cex=0.8)
# points(cs22hatm,col="purple",cex=0.5)
plot(c2.hat.2zonem,lty=6, col="blue",main="Medium Far",ylab="Toluene Concentrations F",xlab="Time",type="l",lwd=3)
#lines(pred22.juliam, col="red",lwd=3)
lines(yhat2pm, col="red",lwd=3,lty=5)
lines(data2zonef1[1:75,13],col="green",lwd=3)
lines(predkf2h,col="chocolate",lwd=3,lty=2)
lines(cs22hatm,col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue","chocolate","red","purple"), lty=c(1,6,2,5,3),legend=c("Measurements F","Model a", "Model b","BNLR","Smoothed"),pch=15,cex=0.75,box.lwd = 0)

# plot(c1.hat.2zoneh, col="blue",main="High Near",ylab="2-butanone Concentrations N",xlab="Time")
# points(pred21.juliah, col="red")
# points(data2zonen[1:41,81],col="green",cex=0.8)
# points(cs12hath,col="purple",cex=0.5)
plot(c1.hat.2zoneh, lty=6,col="blue",main="High Near",ylab="2-butanone Concentrations N",xlab="Time",type="l",lwd=3)
#lines(pred21.juliah, col="red",lwd=3)
lines(yhat1ph, col="red",lwd=3,lty=5)
lines(data2zonen[1:41,81],col="green",lwd=3)
lines(predkf1m,col="chocolate",lwd=3,lty=2)
lines(cs12hath,col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue","chocolate","red","purple"), lty=c(1,6,2,5,3),legend=c("Measurements N","Model a", "Model b","BNLR","Smoothed"),pch=15,cex=0.75,box.lwd = 0)

# plot(c2.hat.2zoneh, col="blue",main="High Far",ylab="2-butanone Concentrations F",xlab="Time")
# points(pred22.juliah, col="red")
# points(data2zonef1[1:41,81],col="green",cex=0.8)
# points(cs22hath,col="purple",cex=0.5)
plot(c2.hat.2zoneh, lty=6,col="blue",main="High Far",ylab="2-butanone Concentrations F",xlab="Time",type="l",lwd=3)
#lines(pred22.juliah, col="red",lwd=3)
lines(yhat2ph, col="red",lwd=3,lty=5)
lines(data2zonef1[1:41,81],col="green",lwd=3)
lines(predkf2m,col="chocolate",lwd=3,lty=2)
lines(cs22hath,col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue","chocolate","red","purple"),lty=c(1,6,2,5,3), legend=c("Measurements F","Model a", "Model b","BNLR","Smoothed"),pch=15,cex=0.75,box.lwd = 0)

#nonlinear regression data
cat("model{
    for (i in 1:N){
    c[i+1,1:2]<-mexp(i*A)%*%c[1,1:2]+B%*%(mexp(i*A)-I)%*%(g)
    logyt[i+1,1:2]~dmnorm(log(abs(c[i+1,1:2])), prec)
    }
    A[1,1]<- -(beta)/VN
    A[1,2]<- beta/VN
    A[2,1]<- beta/VF
    A[2,2]<- -(beta+Qprime)/VF+kl
    B[1,1]<- (-(beta+Qprime)/VF +kl)/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    B[1,2]<- (-beta/VN)/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    B[2,1]<- -beta*VF/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    B[2,2]<- -((beta)/VN)/((beta*Qprime)/(VF*VN)-(beta*kl)/VN)
    g[1,1]<- Gprime/VN
    g[2,1]<-0
    c[1,1]<-1
    c[1,2]<-1
    logyt[1,1:2]~dmnorm(log(abs(c[1,1:2])), prec)
    prec=I
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,1)
    beta~dunif(0,5)
    }",file="modeljags2zonenld.jag")
jagsmnl<-jags.model("modeljags2zonenld.jag",data=list("logyt"=cbind(log(data2zonen[1:75,13]+0.5),log(data2zonef1[1:75,13]+0.5)),
                                                 N=74,"I"=diag(0.01,2), "VN"=0.1, "VF"=11.9),
                  n.adapt =1000)
update(jagsmnl, 2000)
mcmcjags.2zonemnl<-jags.samples(jagsmnl,
                              c('Qprime','Gprime','beta',"kl",'c'),
                              1000)
mcmcsamplesjags.2zonemnl<-coda.samples(jagsmnl,
                                     c('Qprime','Gprime','beta',"kl",'c'),
                                     1000)

m1.mcmc.2zonemnl<-(as.mcmc(mcmcsamplesjags.2zonemnl))
m1.mat.2zonemnl<-as.matrix(mcmcsamplesjags.2zonemnl)
m1.dat.2zonemnl<-as.data.frame(m1.mat.2zonemnl)
qprime.post.2zonemnl<-m1.dat.2zonemnl$Qprime
gprime.post.2zonemnl<-m1.dat.2zonemnl$Gprime
beta.post.2zonemnl<-m1.dat.2zonemnl$beta
kl.post.2zonemnl<-m1.dat.2zonemnl$kl
quantile(qprime.post.2zonemnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonemnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonemnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(kl.post.2zonemnl,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonemnl <- apply(mcmcsamplesjags.2zonemnl[[1]][,4:(78)], 2, mean)
c2.hat.2zonemnl <- apply(mcmcsamplesjags.2zonemnl[[1]][,(79):(78+75)], 2, mean)

#high
jagshnl<-jags.model("modeljags2zonenld.jag",data=list("logyt"=cbind(log(data2zonen[,81]),log(data2zonef1[,81])),
                                                 N=40,"I"=diag(0.1,2), "VN"=0.1, "VF"=11.9),
                  n.adapt =1000)
update(jagshnl, 2000)
mcmcjags.2zonehnl<-jags.samples(jagshnl,
                              c('Qprime','Gprime','beta',"kl",'c'),
                              1000)
mcmcsamplesjags.2zonehnl<-coda.samples(jagshnl,
                                     c('Qprime','Gprime','beta',"kl",'c'),
                                     1000)
m1.mcmc.2zonehnl<-(as.mcmc(mcmcsamplesjags.2zonehnl))
m1.mat.2zonehnl<-as.matrix(mcmcsamplesjags.2zonehnl)
m1.dat.2zonehnl<-as.data.frame(m1.mat.2zonehnl)
qprime.post.2zonehnl<-m1.dat.2zonehnl$Qprime
gprime.post.2zonehnl<-m1.dat.2zonehnl$Gprime
beta.post.2zonehnl<-m1.dat.2zonehnl$beta
kl.post.2zonehnl<-m1.dat.2zonehnl$kl
quantile(qprime.post.2zonehnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonehnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonehnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(kl.post.2zonehnl,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonehnl <- apply(mcmcsamplesjags.2zonehnl[[1]][,4:(44)], 2, mean)
c2.hat.2zonehnl <- apply(mcmcsamplesjags.2zonehnl[[1]][,(45):85], 2, mean)
#low Q=0.056
jagslnl<-jags.model("modeljags2zonenld.jag",data=list("logyt"=cbind(log(data2zonen[1:176,1]),log(data2zonef1[1:176,1]+0.05)),
                                                 N=175,"I"=diag(2), "VN"=0.1, "VF"=11.9),
                  n.adapt =1000)
update(jagslnl, 2000)
mcmcjags.2zonelnl<-jags.samples(jagslnl,
                              c('Qprime','Gprime','beta',"kl",'c'),
                              1000)
mcmcsamplesjags.2zonelnl<-coda.samples(jagslnl,
                                     c('Qprime','Gprime','beta',"kl",'c'),
                                     1000)

m1.mcmc.2zonelnl<-(as.mcmc(mcmcsamplesjags.2zonelnl))
m1.mat.2zonelnl<-as.matrix(mcmcsamplesjags.2zonelnl)
m1.dat.2zonelnl<-as.data.frame(m1.mat.2zonelnl)
qprime.post.2zonelnl<-m1.dat.2zonelnl$Qprime
gprime.post.2zonelnl<-m1.dat.2zonelnl$Gprime
beta.post.2zonelnl<-m1.dat.2zonelnl$beta
kl.post.2zonelnl<-m1.dat.2zonelnl$kl
quantile(qprime.post.2zonelnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonelnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonelnl,c(0.025,0.25,0.5,0.75,0.975))
quantile(kl.post.2zonelnl,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonelnl <- apply(mcmcsamplesjags.2zonelnl[[1]][,4:(103+76)], 2, mean)
c2.hat.2zonelnl <- apply(mcmcsamplesjags.2zonelnl[[1]][,(104+76):355], 2, mean)

quartz()
par(mfrow=c(3,2))
plot(c1.hat.2zonelnl, col="blue",main="Low")
points(data2zonen[1:176,1],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonelnl, col="blue")
points(data2zonef1[1:176,1],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)
plot(c1.hat.2zonemnl, col="blue",main="Medium")
points(data2zonen[1:75,13],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonemnl, col="blue")
points(data2zonef1[1:75,13],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)
plot(c1.hat.2zonehnl, col="blue",main="High")
points(data2zonen[1:41,81],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonehnl, col="blue")
points(data2zonef1[1:41,81],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

#try B2Z
library(B2Z)
#sim
fit.2z <- B2ZM(data = cbind(seq(h,1,h),exp(logyt2[,1]),exp(logyt2[,2])), priorBeta = "unif(0,10)",
                    indep.model = FALSE, priorQ = "unif(11,17)", priorG = "unif(281,482)",
                    S =diag(2), v = 3, VN =pi*10^-2, VF = 3.8, sampler = "MCMC",
                    mcmc.control = list(NUpd = 10000, burnin = 1000, lag = 1, m = 5000))
summary(fit.2z)
quartz()
plot(fit.2z)
c1p<-c2p<-y1p<-y2p<-matrix(0,100,1000)
Ap<-gp<-vector('list')
for(j in 1:1000){
  Ap[[j]]<-matrix(c(-(fit.2z$Beta[j+1000])/VN, (fit.2z$Beta[j+1000])/VN, fit.2z$Beta[j+1000]/VF, -(fit.2z$Beta[j+1000]+fit.2z$Q[j+1000])/VF),nrow=2, ncol=2,
                  byrow=TRUE)
  gp[[j]]<-matrix(c(fit.2z$G[j+1000]/VN,0),nrow=1, ncol=2)
  for(i in 1:(100-1)){
    c1p[1,]<-0.01
    c2p[1,]<-0.5
    y1p[1,]<-rmvnorm(n=1,mean=c(log(c1p[1,j]),log(c2p[1,j])),
                     sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[1]
    y2p[1,]<-rmvnorm(n=1,mean=c(log(c1p[1,j]),log(c2p[1,j])),
                     sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[2]
    c1p[i+1,j]<- (expm(i*h*Ap[[j]])%*%c(c1p[1,j],c2p[1,j])+solve(Ap[[j]])%*%(expm(i*h*Ap[[j]])-diag(2))%*%t(gp[[j]]))[1]
    c2p[i+1,j]<- (expm(i*h*Ap[[j]])%*%c(c1p[1,j],c2p[1,j])+solve(Ap[[j]])%*%(expm(i*h*Ap[[j]])-diag(2))%*%t(gp[[j]]))[2]
    y1p[i+1,j]<- rmvnorm(n=1,mean=c(log(c1p[i+1,j]),log(c2p[i+1,j])),
                         sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[1]
    y2p[i+1,j]<- rmvnorm(n=1,mean=c(log(c1p[i+1,j]),log(c2p[i+1,j])),
                         sigma=matrix(c(fit.2z$tauN[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauNF[j+1000],fit.2z$tauF[j+1000]),2,2,byrow=TRUE),method="svd")[2]
  }}
yhat1p<-exp(apply(y1p,1,mean))
yhat2p<-exp(apply(y2p,1,mean))
plot(apply(exp(y1p),1,mean))
points(c22[,1])
print(mse2b2z<-(sum(((c22[,1])-(yhat1p))^2)+sum(((c22[,2])-(yhat2p))^2))/200)
G2b2z<-(exp(logyt2[,2])-apply(exp(y2p),1,mean))^2+(exp(logyt2[,1])-apply(exp(y1p),1,mean))^2
P2b2z<-apply(exp(y1p),1,var)+apply(exp(y2p),1,var)
print(D2b2z<-sum(G2b2z)+sum(P2b2z))
#high
fit.2zh <- B2ZM(data = cbind(seq(1,41,1),data2zonen[1:41,81],data2zonef1[1:41,81]), priorBeta = "unif(0,5)",
                    indep.model = FALSE, priorQ = "unif(0,1)", priorG = "unif(28,150)",
                    S = diag(2),v= 3, VN =0.1, VF = 11.9, sampler = "MCMC",
                    mcmc.control = list(NUpd = 10000, burnin = 1000, lag = 1, m = 5000))
summary(fit.dependh)
quartz()
plot(fit.dependh)
c1ph<-c2ph<-y1ph<-y2ph<-matrix(0,41,1000)
Ap<-gp<-vector('list')
for(j in 1:1000){
  Aph[[j]]<-matrix(c(-(fit.2zh$Beta[j+1000])/VN, (fit.2zh$Beta[j+1000])/VN, fit.2zh$Beta[j+1000]/VF, -(fit.2zh$Beta[j+1000]+fit.2zh$Q[j+1000])/VF),nrow=2, ncol=2,
                  byrow=TRUE)
  gph[[j]]<-matrix(c(fit.2zh$G[j+1000]/VN,0),nrow=1, ncol=2)
  for(i in 1:(41-1)){
    c1ph[1,]<-8
    c2ph[1,]<-5
    y1ph[1,]<-rmvnorm(n=1,mean=c(log(c1ph[1,j]),log(c2ph[1,j])),
                     sigma=matrix(c(fit.2zh$tauN[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauF[j+1000]),2,2,byrow=T),method="svd")[1]
    y2ph[1,]<-rmvnorm(n=1,mean=c(log(c1ph[1,j]),log(c2ph[1,j])),
                     sigma=matrix(c(fit.2zh$tauN[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauF[j+1000]),2,2,byrow=T),method="svd")[2]
    c1ph[i+1,j]<- (expm(i*Aph[[j]])%*%c(c1ph[1,j],c2ph[1,j])+solve(Aph[[j]])%*%(expm(i*Aph[[j]])-diag(2))%*%t(gph[[j]]))[1]
    c2ph[i+1,j]<- (expm(i*Aph[[j]])%*%c(c1ph[1,j],c2ph[1,j])+solve(Aph[[j]])%*%(expm(i*Aph[[j]])-diag(2))%*%t(gph[[j]]))[2]
    y1ph[i+1,j]<- rmvnorm(n=1,mean=c(log(c1ph[i+1,j]),log(c2ph[i+1,j])),
                         sigma=matrix(c(fit.2zh$tauN[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauF[j+1000]),2,2,byrow=T),method="svd")[1]
    y2ph[i+1,j]<- rmvnorm(n=1,mean=c(log(c1ph[i+1,j]),log(c2ph[i+1,j])),
                         sigma=matrix(c(fit.2zh$tauN[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauNF[j+1000],fit.2zh$tauF[j+1000]),2,2,byrow=T),method="svd")[2]
  }}
yhat1ph<-exp(apply(y1ph,1,mean))
yhat2ph<-exp(apply(y2ph,1,mean))
print(mse2b2z<-(sum(((data2zonen[1:41,81])-(yhat1ph))^2)+sum(((data2zonef1[1:41,81])-(yhat2ph))^2))/82)
# zrep21b2zh<-matrix(0,41,1000)
# zrep22b2zh<-matrix(0,41,1000)
# c1ph<-c2ph<-matrix(0,41,1000)
# Aph<-gph<-vector('list')
# for(j in 1:1000){
#   Aph[[j]]<-matrix(c(-(fit.2zh$Beta[j])/VN, (fit.2zh$Beta[j])/VN, fit.2zh$Beta[j]/VF, -(fit.2zh$Beta[j]+fit.2zh$Q[j])/VF),nrow=2, ncol=2,
#                   byrow=TRUE)
#   gph[[j]]<-matrix(c(fit.2zh$G[j]/VN,0),nrow=1, ncol=2)
#   for(i in 1:(41-1)){
#     c1ph[1,]<-8
#     c2ph[1,]<-5
#     c1ph[i+1,j]<- (expm(h*i*Aph[[j]])%*%c(c1ph[1,j],c2ph[1,j])+solve(Aph[[j]])%*%(expm(i*Aph[[j]])-diag(2))%*%t(gph[[j]]))[1]
#     c2ph[i+1,j]<- (expm(h*i*Aph[[j]])%*%c(c1ph[1,j],c2ph[1,j])+solve(Aph[[j]])%*%(expm(i*Aph[[j]])-diag(2))%*%t(gph[[j]]))[2]
#   }
#   zrep21b2zh[,j]<-rnorm(n=41,log(c1ph[,j]),sqrt(fit.2zh$tauN[j]))
#   zrep22b2zh[,j]<-rnorm(n=41,log(c2ph[,j]),sqrt(fit.2zh$tauF[j]))
# }
# zrep21b2zh<-apply(zrep21b2zh,1,function(x) ifelse(is.na(x),log(8),x))
# zrep22b2zh<-apply(zrep22b2zh,1,function(x) ifelse(is.na(x),log(5),x))
G2b2zh<-(data2zonen[1:41,81]-apply(exp(y1ph),1,mean))^2+(data2zonef1[1:41,81]-apply(exp(y2ph),1,mean))^2
P2b2zh<-apply(exp(y1ph),1,var)+apply(exp(y2ph),1,var)
print(D2b2zh<-sum(G2b2zh)+sum(P2b2zh))
#med
fit.2zm <- B2ZM(data = cbind(seq(1,71,1),data2zonen[4:75,13],data2zonef1[4:75,13]), priorBeta = "unif(0,5)",
                   indep.model = TRUE, priorQ = "unif(0,1)", priorG = "unif(28,150)",
                   tauN.sh = 2,tauN.sc = 1, tauF.sh = 2,tauF.sc = 1, VN =0.1, VF = 11.9, sampler = "MCMC",
                   mcmc.control = list(NUpd = 10000, burnin = 1000, lag = 1, m = 5000))
summary(fit.2zm)
quartz()
plot(fit.2zm)
c1pm<-c2pm<-y1pm<-y2pm<-matrix(0,71,1000)
Apm<-gpm<-vector('list')
for(j in 1:1000){
  Apm[[j]]<-matrix(c(-(fit.2zm$Beta[j+1000])/VN, (fit.2zm$Beta[j+1000])/VN, fit.2zm$Beta[j+1000]/VF, -(fit.2zm$Beta[j+1000]+fit.2zm$Q[j+1000])/VF),nrow=2, ncol=2,
                   byrow=TRUE)
  gpm[[j]]<-matrix(c(fit.2zm$G[j+1000]/VN,0),nrow=1, ncol=2)
  for(i in 1:(71-1)){
    c1pm[1,]<-7
    c2pm[1,]<-2
    y1pm[1,]<-rmvnorm(n=1,mean=c(log(c1pm[1,j]),log(c2pm[1,j])),
                      sigma=diag(c(fit.2zm$tauN[j+1000],fit.2zm$tauF[j+1000])),method="svd")[1]
    y2pm[1,]<-rmvnorm(n=1,mean=c(log(c1pm[1,j]),log(c2pm[1,j])),
                      sigma=diag(c(fit.2zm$tauN[j+1000],fit.2zm$tauF[j+1000])),method="svd")[2]
    c1pm[i+1,j]<- (expm(i*Apm[[j]])%*%c(c1pm[1,j],c2pm[1,j])+solve(Apm[[j]])%*%(expm(i*Apm[[j]])-diag(2))%*%t(gpm[[j]]))[1]
    c2pm[i+1,j]<- (expm(i*Apm[[j]])%*%c(c1pm[1,j],c2pm[1,j])+solve(Apm[[j]])%*%(expm(i*Apm[[j]])-diag(2))%*%t(gpm[[j]]))[2]
    y1pm[i+1,j]<- rmvnorm(n=1,mean=c(log(c1pm[i+1,j]),log(c2pm[i+1,j])),
                          sigma=diag(c(fit.2zm$tauN[j+1000],fit.2zm$tauF[j+1000])),method="svd")[1]
    y2pm[i+1,j]<- rmvnorm(n=1,mean=c(log(c1pm[i+1,j]),log(c2pm[i+1,j])),
                          sigma=diag(c(fit.2zm$tauN[j+1000],fit.2zm$tauF[j+1000])),method="svd")[2]
  }}
yhat1pm<-exp(apply((y1pm),1,mean))
yhat2pm<-exp(apply((y2pm),1,mean))
plot(yhat1pm)
points(data2zonen[4:75,13],col="red")
print(mse2b2zm<-(sum(((data2zonen[4:75,13])-(yhat1pm))^2)+sum(((data2zonef1[4:75,13])-(yhat2pm))^2))/144)
G2b2zm<-(data2zonen[4:75,13]-apply(exp(y1pm),1,mean))^2+(data2zonef1[4:75,13]-apply(exp(y2pm),1,mean))^2
P2b2zm<-apply(exp(y1pm),1,var)+apply(exp(y2pm),1,var)
print(D2b2zm<-sum(G2b2zm)+sum(P2b2zm))
# zrep21b2zm<-matrix(0,75,1000)
# zrep22b2zm<-matrix(0,75,1000)
# c1pm<-c2pm<-matrix(0,75,1000)
# Apm<-gpm<-vector('list')
# for(j in 1:1000){
#   Apm[[j]]<-matrix(c(-(fit.2zm$Beta[j+1000])/VN, (fit.2zm$Beta[j+1000])/VN, fit.2zm$Beta[j+1000]/VF, -(fit.2zm$Beta[j+1000]+fit.2zm$Q[j+1000])/VF),nrow=2, ncol=2,
#                    byrow=TRUE)
#   gpm[[j]]<-matrix(c(fit.2zm$G[j+1000]/VN,0),nrow=1, ncol=2)
#   for(i in 1:(75-1)){
#     c1pm[1,]<-8
#     c2pm[1,]<-5
#     c1pm[i+1,j]<- (expm(i*Apm[[j]])%*%c(c1pm[1,j],c2pm[1,j])+solve(Apm[[j]])%*%(expm(i*Apm[[j]])-diag(2))%*%t(gpm[[j]]))[1]
#     c2pm[i+1,j]<- (expm(i*Apm[[j]])%*%c(c1pm[1,j],c2pm[1,j])+solve(Apm[[j]])%*%(expm(i*Apm[[j]])-diag(2))%*%t(gpm[[j]]))[2]
#   }
#   zrep21b2zm[,j]<-rnorm(n=75,log(c1pm[,j]),sqrt(fit.2zm$tauN[j+1000]))
#   zrep22b2zm[,j]<-rnorm(n=75,log(c2pm[,j]),sqrt(fit.2zm$tauF[j+1000]))
# }
# zrep21b2zm<-apply(zrep21b2zm,1,function(x) ifelse(is.na(x),log(0.01),x))
# zrep22b2zm<-apply(zrep22b2zm,1,function(x) ifelse(is.na(x),log(2),x))


#low
fit.2zl <- B2ZM(data = cbind(seq(1,176,1),data2zonen[1:176,1]+0.05,data2zonef1[1:176,1]+0.05), priorBeta = "unif(0,5)",
                   indep.model = FALSE, priorQ = "unif(0,1)", priorG = "unif(28,150)",
                   S = diag(2),v = 3, VN =0.1, VF = 11.9, sampler = "MCMC",
                   mcmc.control = list(NUpd = 10000, burnin = 1000, lag = 1, m = 5000))
summary(fit.2zl)
quartz()
plot(fit.2zl)
c1pl<-c2pl<-y1pl<-y2pl<-matrix(0,176,1000)
Apl<-gpl<-vector('list')
for(j in 1:1000){
  Apl[[j]]<-matrix(c(-(fit.2zl$Beta[j+1000])/VN, (fit.2zl$Beta[j+1000])/VN, fit.2zl$Beta[j+1000]/VF, -(fit.2zl$Beta[j+1000]+fit.2zl$Q[j+1000])/VF),nrow=2, ncol=2,
                   byrow=TRUE)
  gpl[[j]]<-matrix(c(fit.2zl$G[j+1000]/VN,0),nrow=1, ncol=2)
  for(i in 1:(176-1)){
    c1pl[1,]<-5
    c2pl[1,]<-0.05
    y1pl[1,]<-rmvnorm(n=1,mean=c(log(c1pl[1,j]),log(c2pl[1,j])),
                      sigma=matrix(c(fit.2zl$tauN[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauF[j+1000]),2,2,byrow=T),method="svd")[1]
    y2pl[1,]<-rmvnorm(n=1,mean=c(log(c1pl[1,j]),log(c2pl[1,j])),
                      sigma=matrix(c(fit.2zl$tauN[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauF[j+1000]),2,2,byrow=T),method="svd")[2]
    c1pl[i+1,j]<- (expm(i*Apl[[j]])%*%c(c1pl[1,j],c2pl[1,j])+solve(Apl[[j]])%*%(expm(i*Apl[[j]])-diag(2))%*%t(gpl[[j]]))[1]
    c2pl[i+1,j]<- (expm(i*Apl[[j]])%*%c(c1pl[1,j],c2pl[1,j])+solve(Apl[[j]])%*%(expm(i*Apl[[j]])-diag(2))%*%t(gpl[[j]]))[2]
    y1pl[i+1,j]<- rmvnorm(n=1,mean=c(log(c1pl[i+1,j]),log(c2pl[i+1,j])),
                          sigma=matrix(c(fit.2zl$tauN[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauF[j+1000]),2,2,byrow=T),method="svd")[1]
    y2pl[i+1,j]<- rmvnorm(n=1,mean=c(log(c1pl[i+1,j]),log(c2pl[i+1,j])),
                          sigma=matrix(c(fit.2zl$tauN[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauNF[j+1000],fit.2zl$tauF[j+1000]),2,2,byrow=T),method="svd")[2]
  }}
yhat1pl<-exp(apply(y1pl,1,mean))
yhat2pl<-exp(apply(y2pl,1,mean))
G2b2zl<-(data2zonen[1:176,1]-apply(exp(y1pl),1,mean))^2+(data2zonef1[1:176,1]-apply(exp(y2pl),1,mean))^2
P2b2zl<-apply(exp(y1pl),1,var)+apply(exp(y2pl),1,var)
print(D2b2zl<-sum(G2b2zl)+sum(P2b2zl))
print(mse2b2zl<-(sum(((data2zonen[1:176,1])-(yhat1pl))^2)+sum(((data2zonef1[1:176,1])-(yhat2pl))^2))/352)
# zrep21b2zl<-matrix(0,176,1000)
# zrep22b2zl<-matrix(0,176,1000)
# c1pl<-c2pl<-matrix(0,176,1000)
# Apl<-gpl<-vector('list')
# for(j in 1:1000){
#   Apl[[j]]<-matrix(c(-(fit.2zl$Beta[j])/VN, (fit.2zl$Beta[j])/VN, fit.2zl$Beta[j]/VF, -(fit.2zl$Beta[j]+fit.2zl$Q[j])/VF),nrow=2, ncol=2,
#                    byrow=TRUE)
#   gpl[[j]]<-matrix(c(fit.2zl$G[j]/VN,0),nrow=1, ncol=2)
#   for(i in 1:(176-1)){
#     c1pl[1,]<-5
#     c2pl[1,]<-0
#     c1pl[i+1,j]<- (expm(h*i*Apl[[j]])%*%c(c1pl[1,j],c2pl[1,j])+solve(Apl[[j]])%*%(expm(i*Apl[[j]])-diag(2))%*%t(gpl[[j]]))[1]
#     c2pl[i+1,j]<- (expm(h*i*Apl[[j]])%*%c(c1pl[1,j],c2pl[1,j])+solve(Apl[[j]])%*%(expm(i*Apl[[j]])-diag(2))%*%t(gpl[[j]]))[2]
#   }
#   zrep21b2zl[,j]<-rnorm(n=176,log(c1pl[,j]),sqrt(fit.2zl$tauN[j]))
#   zrep22b2zl[,j]<-rnorm(n=176,log(c2pl[,j]),sqrt(fit.2zl$tauF[j]))
# }
# zrep21b2zl<-apply(zrep21b2zl,1,function(x) ifelse(is.na(x),log(5),x))
# zrep22b2zl<-apply(zrep22b2zl,1,function(x) ifelse(is.na(x),log(0.01),x))

apply((c1pl),1,mean)+apply((c2pl),1,mean)
plot(data2zonef1[1:176,1])
points(apply(exp(zrep22b2zl),2,mean),col="red")
#Turbulent eddy-diffusion model
Dt=1
nloc<-5
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

phi<-1
sigma<-1
tausq<-0.25
set.seed(123)
eta<-rnorm(ntime,0,0.1)
vt<-NULL
vt[1]<-0
for(i in 2:ntime){
  vt[i]<-vt[i-1]+eta[i]
}
vt
#library(geoR)
Sigma.exp <- function(x,y,phi,sigma) {
  Sigma <- matrix(rep(0, length(x)*length(y)), nrow=length(x))
  Sigma<- as.matrix(sigma*exp(-(dist(data.frame(cbind(x,y)),upper=TRUE, diag=TRUE,method="euclidean")*phi)))
  return(Sigma)
}
sigma.sp<-Sigma.exp(x,y,1,1)
dim(sigma.sp)
set.seed(111)
et<-rmvnorm(1,rep(0,nloc),((diag(nloc)+sigma.sp[1:nloc,1:nloc])))
logyt3<-matrix(0,nrow=nloc, ncol=ntime)
for(i in 1:ntime){
  for(j in 1:nloc){
    logyt3[j,i]<-log(csp[j,i]+et[j])+eta[i]
  }
}
plot(csp[1,])
points(exp(logyt3[1,]),col="red")

plot(csp[2,], ylim=c(min(csp[2,]),max(exp(logyt3[2,]))))
points(exp(logyt3[2,]),col="red")

plot(csp[3,],ylim=c(min(csp[3,]),max(exp(logyt3[3,]))))
points(exp(logyt3[3,]),col="red")

plot(csp[4,])
points(exp(logyt3[4,]),col="red")

plot(csp[5,])
points(exp(logyt3[5,]),col="red")

distances<-as.matrix(dist(data.frame(cbind(c(x),c(y))),upper=TRUE, diag=TRUE,method="euclidean"))

cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(0.001,0.001)
    loglik[i+1]<-logdensity.mnorm(logyt[,i+1],log(csp[,i+1]+abs(et)),prec)
    for(j in 1:nloc){
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1])
    logyt[j,i+1]~dnorm(log(csp[j,i+1]+et[j]), 10)
    logyhat[j,i+1]~dnorm(log(csp[j,i+1]+et[j]), 10)
    }
    }
#    for(i in 1:(N+1)){
 #   eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
  #  }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-I[m,j]+sigma*exp(-dist[m,j])
    }
    }
    kinv<-inverse(k)
    et~dmnorm(muet,kinv)
    #logyt[j,1]~dnorm(log(csp[j,1]+et[j])+eta[1], tausq)
for(j in 1:nloc){
    logyhat[j,1]~dnorm(log(csp[j,1]+et[j]), 10)
    csp[j,1]<-exp(logyt[j,1])+omega[1]
    }
    omega[1]~dgamma(0.001,0.001)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,3)
    sigma~dgamma(0.01,0.01)
    }",file="modeljagseddy.jag")
jags.eddy<-jags.model("modeljagseddy.jag",data=list("logyt"=logyt3, N=99,nloc=5,"h"=h,"muet"=rep(0,nloc),
                                                    "pi"=pi, "dist"=distances, "matdist"=(matdist),"prec"=diag(10,nloc),"I"=diag(nloc)))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','csp','loglik','et','logyhat'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','csp','loglik','et','logyhat'),
                                   1000)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
mean(crps_sample(y = (yt33[1,]), dat = exp(mcmcjags.eddy$logyhat[1,,,1])))+
  mean(crps_sample(y = (yt33[2,]), dat =exp(mcmcjags.eddy$logyhat[2,,,1])))+
  mean(crps_sample(y = (yt33[3,]), dat =exp(mcmcjags.eddy$logyhat[3,,,1])))+
  mean(crps_sample(y = (yt33[4,]), dat =exp(mcmcjags.eddy$logyhat[4,,,1])))+
  mean(crps_sample(y = (yt33[5,]), dat =exp(mcmcjags.eddy$logyhat[5,,,1])))
yhat3<-matrix(0, nrow=5, ncol=100)
yhatlower503<-matrix(0, nrow=5, ncol=100)
yhatupper503<-matrix(0, nrow=5, ncol=100)
yhatlower903<-matrix(0, nrow=5, ncol=100)
yhatupper903<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  yhat3[j,] <- apply( exp(mcmcjags.eddy$logyhat[j,,,1]),1,mean)
  yhatlower503[j,]<- apply(exp(mcmcjags.eddy$logyhat[j,,,1]),1,function(x) quantile(x,0.25))
  yhatupper503[j,]<- apply( exp(mcmcjags.eddy$logyhat[j,,,1]),1,function(x) quantile(x,0.75))
  yhatlower903[j,]<- apply( exp(mcmcjags.eddy$logyhat[j,,,1]),1,function(x) quantile(x,0.05))
  yhatupper903[j,]<- apply(exp(mcmcjags.eddy$logyhat[j,,,1]),1,function(x) quantile(x,0.95))
}
plot((yt3[1,]))
lines(yhat3[2,])
lines( yhatlower503[1,])
lines(yhatupper503[1,])
cov503<-cov903<-matrix(0,100,5)
for(i in 1:length(logyt2[,1])){
  for(j in 1:5){
  cov503[i,j]<-(ifelse(yhatlower503[j,i]<=(yt33[j,i])&(yt33[j,i])<=yhatupper503[j,i],1,0))
  cov903[i,j]<-(ifelse(yhatlower903[j,i]<=(yt33[j,i])&(yt33[j,i])<=yhatupper903[j,i],1,0))
}}
sum(cov503)/500
sum(cov903)/500
zreped1<-matrix(0,100,1000)
zreped2<-matrix(0,100,1000)
zreped3<-matrix(0,100,1000)
zreped4<-matrix(0,100,1000)
zreped5<-matrix(0,100,1000)
for(j in 1:1000){
  for(i in 1:100){
    zreped1[i,j]<-rnorm(1,log(mcmcsamplesjags.eddy[[1]][,2+(5*i+1-(6-1))][j]+et.posteddy[j]),
                          0.1)
    zreped2[i,j]<-rnorm(1,log(mcmcsamplesjags.eddy[[1]][,2+(5*i+1-(6-2))][j]+et.posteddy[j]),
                       0.1)
    zreped3[i,j]<-rnorm(1,log(mcmcsamplesjags.eddy[[1]][,2+(5*i+1-(6-3))][j]+et.posteddy[j]),
                       0.1)
    zreped4[i,j]<-rnorm(1,log(mcmcsamplesjags.eddy[[1]][,2+(5*i+1-(6-4))][j]+et.posteddy[j]),
                       0.1)
    zreped5[i,j]<-rnorm(1,log(mcmcsamplesjags.eddy[[1]][,2+(5*i+1-(6-5))][j]+et.posteddy[j]),
                       0.1)
  }}
Ged<-((csp[1,])-apply(exp(zreped1),1,mean))^2+((csp[2,])-apply(exp(zreped2),1,mean))^2+
  (csp[3,])-apply(exp(zreped3),1,mean)+(csp[4,])-apply((csp),1,mean)+
  (csp[5,])-apply(exp(zreped5),1,mean)
Ped<-apply(exp(zreped1),1,var)+apply(exp(zreped2),1,var)+apply(exp(zreped3),1,var)+apply(exp(zreped4),1,var)+apply(exp(zreped5),1,var)
print(Ded<-sum(Ged)+sum(Ped))
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.eddy,c(0.025,0.25,0.5,0.75,0.975))
quantile(sigma.post.eddy,c(0.025,0.25,0.5,0.75,0.975))
loglik.posteddy<-t(mcmcjags.eddy$loglik[,,1])
dim(loglik.posteddy)
waic(loglik.posteddy[,2:100])
et.posteddy<-matrix(0,nloc,1000)
for(i in 1:nloc){
et.posteddy[i,]<-t(mcmcjags.eddy$et[i,,1])
}

col.br=colorRampPalette(c("blue", "cyan", "yellow", "red"))
col.pal <- col.br(5)
coords <- cbind(c(x[10],x[1:10],x[1]),c(y[1],y[1:10],y[10]))
x.res <- 100
y.res <- 100
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
surf <- mba.surf(cbind(coords,c(1.65, apply(exp(et.posteddy),1,mean),1.65)), no.X = x.res, no.Y = y.res, h = 5, m = 2, extend = FALSE)$xyz.est
quartz()
jpeg("/Users/n_a_abdallah/Desktop/spatial/Project2/spatialplot.jpg",width = 8, height = 4, units = 'in',res = 250)
image.plot(surf, xaxs = "r", yaxs = "r", col = col.br(25))
#contour(surf,add=TRUE) 
dev.off()

csp.hat<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  for(i in 1:100){
    csp.hat[j,i] <- mean(mcmcsamplesjags.eddy[[1]][,2+(5*i+1-(6-j))])
  }
}
csp.hat
plot(csp[5,],csp.hat[5,])
abline(0,1)
cs1ed<-matrix(0,1000,100)
cs2ed<-matrix(0,1000,100)
cs3ed<-matrix(0,1000,100)
cs4ed<-matrix(0,1000,100)
cs5ed<-matrix(0,1000,100)

cs1ed[,100]<-mcmcsamplesjags.eddy[[1]][,(2+(5*100+1-(6-1)))]
cs2ed[,100]<-mcmcsamplesjags.eddy[[1]][,(2+(5*100+1-(6-2)))]
cs3ed[,100]<-mcmcsamplesjags.eddy[[1]][,(2+(5*100+1-(6-3)))]
cs4ed[,100]<-mcmcsamplesjags.eddy[[1]][,(2+(5*100+1-(6-4)))]
cs5ed[,100]<-mcmcsamplesjags.eddy[[1]][,(2+(5*100+1-(6-5)))]

omega.post1zone=rgamma(aomega.post1zone,bomega.post1zone)
N=n-1
for(j in 2:(N)){
  cs1ed[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddy[[1]][,2+(5*(j+N-2*(j-1)+1)+1-(6-1))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(j+N-2*(j-1)+1))^1.5))*
                          (exp(-(matdist[1])^2/(4*Dt.post.eddy*(j+N-2*(j-1)+1)))))
  cs2ed[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddy[[1]][,2+(5*(j+N-2*(j-1)+1)+1-(6-2))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(j+N-2*(j-1)+1))^1.5))*
                            (exp(-(matdist[2])^2/(4*Dt.post.eddy*(j+N-2*(j-1)+1)))))
  cs3ed[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddy[[1]][,2+(5*(j+N-2*(j-1)+1)+1-(6-3))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(j+N-2*(j-1)+1))^1.5))*
                            (exp(-(matdist[3])^2/(4*Dt.post.eddy*(j+N-2*(j-1)+1)))))
  cs4ed[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddy[[1]][,2+(5*(j+N-2*(j-1)+1)+1-(6-4))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(j+N-2*(j-1)+1))^1.5))*
                            (exp(-(matdist[4])^2/(4*Dt.post.eddy*(j+N-2*(j-1)+1)))))
  cs5ed[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddy[[1]][,2+(5*(j+N-2*(j-1)+1)+1-(6-5))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(j+N-2*(j-1)+1))^1.5))*
                            (exp(-(matdist[5])^2/(4*Dt.post.eddy*(j+N-2*(j-1)+1)))))
} 
cs1ed[,1]<-(mcmcsamplesjags.eddy[[1]][,(2+(5*2+1-(6-1)))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(2))^1.5))*
              (exp(-(matdist[1])^2/(4*Dt.post.eddy*(2)))))
cs2ed[,1]<-(mcmcsamplesjags.eddy[[1]][,(2+(5*2+1-(6-2)))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(2))^1.5))*
              (exp(-(matdist[2])^2/(4*Dt.post.eddy*(2)))))
cs3ed[,1]<-(mcmcsamplesjags.eddy[[1]][,(2+(5*2+1-(6-3)))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(2))^1.5))*
              (exp(-(matdist[3])^2/(4*Dt.post.eddy*(2)))))
cs4ed[,1]<-(mcmcsamplesjags.eddy[[1]][,(2+(5*2+1-(6-4)))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(2))^1.5))*
              (exp(-(matdist[4])^2/(4*Dt.post.eddy*(2)))))
cs5ed[,1]<-(mcmcsamplesjags.eddy[[1]][,(2+(5*2+1-(6-5)))]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(2))^1.5))*
              (exp(-(matdist[5])^2/(4*Dt.post.eddy*(2)))))
cs1ed.hat <- apply(cs1ed, 2, mean)
cs2ed.hat <- apply(cs2ed, 2, mean)
cs3ed.hat <- apply(cs3ed, 2, mean)
cs4ed.hat <- apply(cs4ed, 2, mean)
cs5ed.hat <- apply(cs5ed, 2, mean)

quartz()
par(mfrow=c(3,2))
plot(csp[1,],col="red")
points(csp.hat[1,],col="blue")
points(exp(logyt3[1,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)
plot(csp[2,],col="red")
points(csp.hat[2,],col="blue")
points(exp(logyt3[2,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
plot(csp[3,],col="red")
points(csp.hat[3,],col="blue")
points(exp(logyt3[3,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct3","chat3","yt3"),bg="white", lty=1,cex=0.8)
plot(csp[4,],col="red")
points(csp.hat[4,],col="blue")
points(exp(logyt3[4,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct4","chat4","yt4"),bg="white", lty=1,cex=0.8)
plot(csp[5,],col="red")
points(csp.hat[5,],col="blue")
points(exp(logyt3[5,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct5","chat5","yt5"),bg="white", lty=1,cex=0.8)
#sims
#all sim
jagseddysimpar<-mcmcjags.eddysimpar<-mcmcsamplesjags.eddysimpar<-vector("list")
for(i in 1:50){
  jagseddysimpar[[i]]<-jags.model("modeljagseddy.jag",data=list("logyt"=log(yteddysim[i,,]), N=99,nloc=5,"h"=h,"muet"=rep(0,nloc),
                                                                "pi"=pi, "dist"=distances, "matdist"=(matdist),"prec"=diag(10,nloc),"I"=diag(nloc)))
  update(jagseddysimpar[[i]], 1000)
  mcmcjags.eddysimpar[[i]]<-  jags.samples(jagseddysimpar[[i]],
                                           c('Dt','Gprime','sigma','csp','loglik','et','logyhat'),
                                           1000)
  mcmcsamplesjags.eddysimpar[[i]]<-coda.samples(jagseddysimpar[[i]],
                                                c('Dt','Gprime','sigma','csp','loglik','et','logyhat'),
                                                1000)
}
m1.mcmc.eddysimpar<-m1.mat.eddysimpar<-m1.dat.eddysimpar<-vector("list")
yhateddysimpar<-yhatlower50eddysimpar<-yhatupper50eddysimpar<-yhatlower90eddysimpar<-yhatupper90eddysimpar<-
  chateddysimpar<-chatlower50eddysimpar<-chatupper50eddysimpar<-chatlower90eddysimpar<-chatupper90eddysimpar<-array(0,dim = c(50, 100, 5))
for(i in 1:50){
  m1.mcmc.eddysimpar[[i]]<-(as.mcmc(mcmcsamplesjags.eddysimpar[[i]]))
  m1.mat.eddysimpar[[i]]<-as.matrix(mcmcsamplesjags.eddysimpar[[i]])
  m1.dat.eddysimpar[[i]]<-as.data.frame(m1.mat.eddysimpar[[i]])
  for(j in 1:100){
    yhateddysimpar[i,j,]<- apply((exp(mcmcjags.eddysimpar[[i]]$logyhat[,j,,])), 1, mean)
    yhatlower50eddysimpar[i,j,]<- apply((exp(mcmcjags.eddysimpar[[i]]$logyhat[,j,,])), 1, function(x) quantile(x,0.25))
    yhatupper50eddysimpar[i,j,]<-apply(exp(mcmcjags.eddysimpar[[i]]$logyhat[,j,,]), 1, function(x) quantile(x,0.75))
    yhatlower90eddysimpar[i,j,]<- apply(exp(mcmcjags.eddysimpar[[i]]$logyhat[,j,,]), 1, function(x) quantile(x,0.05))
    yhatupper90eddysimpar[i,j,]<-apply(exp(mcmcjags.eddysimpar[[i]]$logyhat[,j,,]), 1, function(x) quantile(x,0.95))
    chateddysimpar[i,j,]<- apply(mcmcjags.eddysimpar[[i]]$csp[,j,,], 1, mean)
    chatlower50eddysimpar[i,j,]<- apply(mcmcjags.eddysimpar[[i]]$csp[,j,,], 1, function(x) quantile(x,0.25))
    chatupper50eddysimpar[i,j,]<-apply(mcmcjags.eddysimpar[[i]]$csp[,j,,], 1, function(x) quantile(x,0.75))
    chatlower90eddysimpar[i,j,]<- apply(mcmcjags.eddysimpar[[i]]$csp[,j,,], 1, function(x) quantile(x,0.05))
    chatupper90eddysimpar[i,j,]<-apply(mcmcjags.eddysimpar[[i]]$csp[,j,,], 1, function(x) quantile(x,0.95))
  }
}
cov50eddysimpar<-cov90eddysimpar<-mseeddysimpar<-crpseddysimpar<-array(0, dim=c(50,100,5))
for(i in 1:50){
  for (k in 1:5){
    cov50eddysimpar[i,,k]<-(ifelse(chatlower50eddysimpar[i,,k]<=(csp[k,])&(csp[k,])<=chatupper50eddysimpar[i,,k],1,0))
    cov90eddysimpar[i,,k]<-(ifelse(chatlower90eddysimpar[i,,k]<=(csp[k,])&(csp[k,])<=chatupper90eddysimpar[i,,k],1,0))
    mseeddysimpar[i,,k]<-sum((chateddysimpar[i,,k]-(csp[,k]))^2)/100
    crpseddysimpar[i,,k]<-mean(crps_sample(y = c(yteddysim[i,k,]), dat =exp(mcmcjags.eddysimpar[[i]]$logyhat[k,,,])))
  }
  
}
mean((apply(cov50eddysimpar,1 ,function(x) sum(x)/500)))
mean((apply(cov90eddysimpar,1, function(x) sum(x)/500)))
mean((apply(mseeddysimpar,1,mean)))
mean((apply(crpseddysimpar,1,mean)))
gprime.post.eddysimpar<-Dt.post.eddysimpar<-vector("list")
for(i in 1:50){
  Dt.post.eddysimpar[[i]]<-(m1.dat.eddysimpar[[i]])$Dt
  gprime.post.eddysimpar[[i]]<-(m1.dat.eddysimpar[[i]])$Gprime
}
quantile(unlist(Dt.post.eddysimpar),c(0.025,0.25,0.5,0.75,0.975))
quantile(unlist(gprime.post.eddysimpar),c(0.025,0.25,0.5,0.75,0.975))
plot(chateddysimpar[i,,4])
points(c(yteddysim[i,4,]),col="red")
points(csp[4,],col="blue")
points(chatlower90eddysimpar[i,,4])
points(chatupper90eddysimpar[i,,4])
#pf results
logyteddy<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/logyteddyb2.csv")
cspeddy<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/cspeddyb2.csv")
thetaeddyjulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmceddyb2.csv")
apply(thetaeddyjulia, 2, function(x){quantile(x)})
xjuliaeddy1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy1b2.csv")
xjuliaeddy2<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy2b2.csv")
xjuliaeddy3<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy3b2.csv")
xjuliaeddy4<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy4b2.csv")
xjuliaeddy5<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy5b2.csv")

prededdy1.julia<-apply(xjuliaeddy1, 2, mean)
prededdy2.julia<-apply(xjuliaeddy2, 2, mean)
prededdy3.julia<-apply(xjuliaeddy3, 2, mean)
prededdy4.julia<-apply(xjuliaeddy4, 2, mean)
prededdy5.julia<-apply(xjuliaeddy5, 2, mean)
print(mseedpf<-(sum(((csp[1,])-(prededdy1.julia))^2)+
                  sum((csp[2,]-(prededdy2.julia))^2)+sum((csp[3,]-(prededdy3.julia))^2)+
                  sum((csp[4,]-(prededdy5.julia))^2)+sum((csp[5,]-(prededdy5.julia))^2))/500)
print(mseedst<-(sum((csp[1,]-(csp.hat[1,]))^2)+sum((csp[2,]-(csp.hat[2,]))^2)+sum((csp[3,]-(csp.hat[3,]))^2)+
                            sum((csp[4,]-(csp.hat[4,]))^2)+sum((csp[5,]-(csp.hat[5,]))^2))/500)
loglikpfeddy<-matrix(0,1000,100)
for(i in 1:100){
  for(j in 1:1000)
    loglikpfeddy[j,i]<-log(dmvnorm(logyt3[,i],cbind(log(xjuliaeddy1[j,i]),log(xjuliaeddy2[j,i]),
                                                    log(xjuliaeddy3[j,i]),log(xjuliaeddy4[j,i]),log(xjuliaeddy5[j,i])),
                                 diag(0.1,5)))
}
waic(loglikpfeddy[,1:100])
#all sim
quartz()
par(mfrow=c(3,3))
# plot(exp(logyt3[1,]),col="green",ylab="Concentrations",xlab="Time", main="a L1")
# points( csp.hat[1,],col="red")
# lines(csp[1,],col="blue",lwd=3)
plot(exp(logyt3[1,]),col="green",ylab="Concentrations",xlab="Time", main="a L1",type="l",lwd=3)
lines( csp.hat[1,],col="red",lwd=3,lty=3)
lines(csp[1,],col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"), lty=c(1,3,6),legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

# plot(unlist(exp(logyteddy[1,])),col="green", ylab="Concentrations",xlab="Time", main="b L1")
# points(prededdy1.julia, col="red")
# lines(unlist(cspeddy[1,]),col="blue",lwd=3)
# plot(unlist(exp(logyt3[1,])),col="green", ylab="Concentrations",xlab="Time", main="b L1",type="l",lwd=3)
# lines(prededdy1.julia, col="red",lwd=3)
# lines(unlist(cspeddy[1,]),col="blue",lwd=2)
# legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
#                                                             "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(unlist(exp(logyt3[1,])),col="green", ylab="Concentrations",xlab="Time", main="b L1",type="l",lwd=3)
lines(predkfeddy1, col="red",lwd=3,lty=6)
lines(unlist(cspeddy[1,]),col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(unlist(exp(logyteddy[1,])),col="green", ylab="Concentrations",xlab="Time", main="Smoothing L1")
# points(cs1ed.hat,col="purple")
# lines(unlist(cspeddy[1,]),col="blue",lwd=3)
plot(unlist(exp(logyt3[1,])),col="green", ylab="Concentrations",xlab="Time", main="Smoothing L1",type="l",lwd=3)
lines(cs1ed.hat,col="purple",lwd=3,lty=6)
lines(unlist(cspeddy[1,]),col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "purple"),lty=c(1,3,6), legend=c("Measurements","True",
                                                            "Smoothed"),pch=15,cex=0.75,box.lwd = 0)

# plot(exp(logyt3[2,]),col="green", ylab="Concentrations",xlab="Time", main="a L2")
# points(csp.hat[2,],col="red")
# lines(csp[2,],col="blue",lwd=3)
plot(exp(logyt3[2,]),col="green", ylab="Concentrations",xlab="Time", main="a L2",type="l",lwd=3)
lines(csp.hat[2,],col="red",lwd=3,lty=6)
lines(csp[2,],col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(unlist(exp(logyteddy[2,])),col="green", ylab="Concentrations",xlab="Time", main="b L2")
# points(prededdy2.julia, col="red")
# lines(unlist(cspeddy[2,]),col="blue",lwd=3)
# plot(unlist(exp(logyt3[2,])),col="green", ylab="Concentrations",xlab="Time", main="b L2",type="l",lwd=3)
# lines(prededdy2.julia, col="red",lwd=3)
# lines(unlist(cspeddy[2,]),col="blue",lwd=2)
# legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
#                                                             "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(unlist(exp(logyt3[2,])),col="green", ylab="Concentrations",xlab="Time", main="b L2",type="l",lwd=3)
lines(predkfeddy2, col="red",lwd=3,lty=6)
lines(unlist(cspeddy[2,]),col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(unlist(exp(logyteddy[2,])),col="green", ylab="Concentrations",xlab="Time", main="Smoothing L2")
# points(cs2ed.hat,col="purple")
# lines(unlist(cspeddy[2,]),col="blue",lwd=3)
plot(unlist(exp(logyt3[2,])),col="green", ylab="Concentrations",xlab="Time", main="Smoothing L2",type="l",lwd=3)
lines(cs2ed.hat,col="purple",lwd=3,lty=6)
lines(unlist(cspeddy[2,]),col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "purple"),lty=c(1,3,6), legend=c("Measurements","True",
                                                               "Smoothed"),pch=15,cex=0.75,box.lwd = 0)

# plot(exp(logyt3[3,]),col="green", ylab="Concentrations",xlab="Time", main="a L3")
# points(csp.hat[3,],col="red")
# lines(csp[3,],col="blue",lwd=3)
plot(exp(logyt3[3,]),col="green", ylab="Concentrations",xlab="Time", main="a L3",type="l",lwd=3)
lines(csp.hat[3,],col="red",lwd=3,lty=6)
lines(csp[3,],col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"),lty=c(1,3,6), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(unlist(exp(logyteddy[3,])),col="green", ylab="Concentrations",xlab="Time", main="b L3")
# points(prededdy3.julia, col="red")
# lines(unlist(cspeddy[3,]),col="blue",lwd=3)
# plot(unlist(exp(logyt3[3,])),col="green", ylab="Concentrations",xlab="Time", main="b L3",type="l",lwd=3)
# lines(prededdy3.julia, col="red",lwd=3)
# lines(unlist(cspeddy[3,]),col="blue",lwd=2)
# legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
#                                                             "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(unlist(exp(logyt3[3,])),col="green", ylab="Concentrations",xlab="Time", main="b L3",type="l",lwd=3)
lines(predkfeddy3, col="red",lwd=3,lty=6)
lines(unlist(cspeddy[3,]),col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "red"), lty=c(1,3,6),legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
# plot(unlist(exp(logyteddy[3,])),col="green", ylab="Concentrations",xlab="Time", main="Smoothing L3")
# points(cs3ed.hat,col="purple")
# lines(unlist(cspeddy[3,]),col="blue",lwd=3)
plot(unlist(exp(logyt3[3,])),col="green", ylab="Concentrations",xlab="Time", main="Smoothing L3",type="l",lwd=3)
lines(cs3ed.hat,col="purple",lwd=3,lty=6)
lines(unlist(cspeddy[3,]),col="blue",lwd=2,lty=3)
legend("bottomright",col=c("green","blue", "purple"),lty=c(1,3,6), legend=c("Measurements","True",
                                                               "Smoothed"),pch=15,cex=0.75,box.lwd = 0)

plot(csp.hat[4,],col="red", ylab="Concentrations",xlab="Time", main="State-by-state update L4")
lines(csp[4,],col="blue",lwd=2)
points(exp(logyt3[4,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(prededdy4.julia, col="red", ylab="Concentrations",xlab="Time", main="PMMH L4")
lines(unlist(cspeddy[4,]),col="blue",lwd=2)
points(unlist(exp(logyteddy[4,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(csp.hat[5,],col="red",ylab="Concentrations",xlab="Time", main="State-by-state update L5")
lines(csp[5,],col="blue", lwd=2)
points(exp(logyt3[5,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(prededdy5.julia, col="red",ylab="Concentrations",xlab="Time", main="PMMH L5")
lines(unlist(cspeddy[5,]),col="blue", lwd=3)
points(unlist(exp(logyteddy[5,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

#data
setwd("/Users/n_a_abdallah/Desktop/spatial/Project2/acetone/")
acetone <- list.files(pattern = ".csv$")
acetonedata<-lapply(acetone, read.csv)
plot(acetonedata[[2]][,2])

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
    omega[i+1]~dgamma(0.001,0.001)
    loglik[i+1]<-logdensity.mnorm(logacetone00484[i+1,],log(csp[i+1,]+et),prec)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1])
    logacetone00484[i+1,j]~dnorm(log(csp[i+1,j]+et[j]), 10)
logyhat[i+1,j]~dnorm(log(csp[i+1,j]+et[j]),10)
    }
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
logyhat[1,j]~dnorm(log(csp[1,j]+et[j]),10)
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00484[1,j])+omega[1]
    }
    omega[1]~dgamma(0.001,0.001)
    #vt[1]<-0
    Gprime~dunif(1104,1650)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
#unstructured cov
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(0.001,0.001)
    loglik[i+1]<-logdensity.mnorm(logacetone00484[i+1,],log(csp[i+1,]+et),prec)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1])
    logacetone00484[i+1,j]~dnorm(log(csp[i+1,j])+et[j], 10)
    logyhat[i+1,j]~dnorm(log(csp[i+1,j])+et[j],10)
    }
    }
    et~dmnorm(c(0,0),kinv)
    for(j in 1:2){
    logyhat[1,j]~dnorm(log(csp[1,j])+et[j],10)
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00484[1,j])+omega[1]
    }
    omega[1]~dgamma(0.001,0.001)
    Gprime~dunif(1104,1650)
    Dt~dunif(0,1)
    kinv~dwish(I,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddyd<-jags.model("modeljagseddydata.jag",data=list("logacetone00484"=logacetone00484, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "pi"=pi, "dist"=distance, "matdist"=matdist,"prec"=diag(10,2),"I"=diag(2)))
update(jags.eddyd, 1000)
mcmcjags.eddyd<-jags.samples(jags.eddyd,
                            c('Dt','Gprime','sigma','csp','loglik','et','logyhat'),
                            1000)
mcmcsamplesjags.eddyd<-coda.samples(jags.eddyd,
                                   c('Dt','Gprime','sigma','csp','loglik','et','logyhat'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddyd<-(as.mcmc(mcmcsamplesjags.eddyd))
m1.mat.eddyd<-as.matrix(mcmcsamplesjags.eddyd)
m1.dat.eddyd<-as.data.frame(m1.mat.eddyd)
aomega.post.eddyd<-m1.dat.eddyd$aomega
bomega.post.eddyd<-m1.dat.eddyd$bomega
Dt.post.eddyd<-m1.dat.eddyd$Dt
gprime.post.eddyd<-m1.dat.eddyd$Gprime
sigma.post.eddyd<-m1.dat.eddyd$sigma
phi.post.eddyd<-m1.dat.eddyd$phi
et.posteddyd<-matrix(0,2,1000)
for(i in 1:2){
  et.posteddyd[i,]<-t(mcmcjags.eddyd$et[i,,1])
}
zreped1d<-matrix(0,88,1000)
zreped2d<-matrix(0,88,1000)
for(j in 1:1000){
  for(i in 1:88){
    zreped1d[i,j]<-rnorm(1,log(mcmcsamplesjags.eddyd[[1]][,2+i][j]+et.posteddyd[j]+abs(min(et.posteddyd))),
                       0.1)
    zreped2d[i,j]<-rnorm(1,log(mcmcsamplesjags.eddyd[[1]][,2+N+i][j]+et.posteddyd[j]+abs(min(et.posteddyd))),
                       0.1)
  }}
Gedd<-(exp(logacetone00484[3:88,1])-apply(exp(zreped1d[3:88,]),1,mean))^2+(exp(logacetone00484[3:88,2])-apply(exp(zreped2d[3:88,]),1,mean))^2
Pedd<-apply(exp(zreped1d[3:88,]),1,var)+apply(exp(zreped2d[3:88,]),1,var)
print(Dedd<-sum(Gedd)+sum(Pedd))
quantile(Dt.post.eddyd,c(0.025,0.25,0.5,0.75,0.975))/120
quantile(gprime.post.eddyd,c(0.025,0.25,0.5,0.75,0.975))
loglik.posteddyd<-t(mcmcjags.eddyd$loglik[,,1])
dim(loglik.posteddyd)
waic(loglik.posteddyd[,2:N])
mean(crps_sample(y = exp(logacetone00484[,1]), dat = exp(mcmcjags.eddyd$logyhat[,1,,1])))+
  mean(crps_sample(y = exp(logacetone00484[,2]), dat =exp(mcmcjags.eddyd$logyhat[,2,,1])))
yhated<-matrix(0, nrow=2, ncol=89)
yhatlower50ed<-matrix(0, nrow=2, ncol=89)
yhatupper50ed<-matrix(0, nrow=2, ncol=89)
yhatlower90ed<-matrix(0, nrow=2, ncol=89)
yhatupper90ed<-matrix(0, nrow=2, ncol=89)
for(j in 1:2){
  yhated[j,] <- apply( exp(mcmcjags.eddyd$logyhat[,j,,1]),1,mean)
  yhatlower50ed[j,]<- apply(exp(mcmcjags.eddyd$logyhat[,j,,1]),1,function(x) quantile(x,0.25))
  yhatupper50ed[j,]<- apply( exp(mcmcjags.eddyd$logyhat[,j,,1]),1,function(x) quantile(x,0.75))
  yhatlower90ed[j,]<- apply( exp(mcmcjags.eddyd$logyhat[,j,,1]),1,function(x) quantile(x,0.05))
  yhatupper90ed[j,]<- apply(exp(mcmcjags.eddyd$logyhat[,j,,1]),1,function(x) quantile(x,0.95))
}
cov50ed<-cov90ed<-matrix(0,89,2)
for(i in 1:89){
  for(j in 1:2){
    cov50ed[i,j]<-(ifelse(yhatlower50ed[j,i]<=(exp(logacetone00484[i,j]))&exp(logacetone00484[i,j])<=yhatupper50ed[j,i],1,0))
    cov90ed[i,j]<-(ifelse(yhatlower90ed[j,i]<=exp(logacetone00484[i,j])&exp(logacetone00484[i,j])<=yhatupper90ed[j,i],1,0))
  }}
sum(cov50ed)/(89*2)
sum(cov90ed)/(89*2)
csp.hatd<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hatd[i,1] <- mean(mcmcsamplesjags.eddyd[[1]][,2+i])
    csp.hatd[i,2] <- mean(mcmcsamplesjags.eddyd[[1]][,2+N+i])
    
  }
}
csp.hatd
quartz()
plot(acetone18[,3]*54/24,csp.hatd[,1])
abline(0,1)

cs1edd<-matrix(0,1000,N+1)
cs2edd<-matrix(0,1000,N+1)

cs1edd[,N+1]<-mcmcsamplesjags.eddyd[[1]][,(2+(N+1))]
cs2edd[,N+1]<-mcmcsamplesjags.eddyd[[1]][,(2+(2*N+1))]

for(j in 2:(N)){
  cs1edd[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddyd[[1]][,2+(j+N-2*(j-1)+1)]-h*(gprime.post.eddyd/(4*(pi*Dt.post.eddyd*(j+N-2*(j-1)+1))^1.5))*
                            (exp(-(matdist[1])^2/(4*Dt.post.eddyd*(j+N-2*(j-1)+1))))-rgamma(1000,0.01,0.01))
  cs2edd[,(j+N-2*(j-1))]<-(mcmcsamplesjags.eddyd[[1]][,2+N+(j+N-2*(j-1)+1)]-h*(gprime.post.eddyd/(4*(pi*Dt.post.eddyd*(j+N-2*(j-1)+1))^1.5))*
                            (exp(-(matdist[2])^2/(4*Dt.post.eddyd*(j+N-2*(j-1)+1))))-rgamma(1000,0.01,0.01))
} 
cs1edd[,1]<-(mcmcsamplesjags.eddyd[[1]][,(2+2)]-h*(gprime.post.eddyd/(4*(pi*Dt.post.eddyd*(2))^1.5))*
              (exp(-(matdist[1])^2/(4*Dt.post.eddyd*(2))))-rgamma(1000,0.01,0.01))
cs2edd[,1]<-(mcmcsamplesjags.eddyd[[1]][,(2+N+2)]-h*(gprime.post.eddy/(4*(pi*Dt.post.eddy*(2))^1.5))*
              (exp(-(matdist[2])^2/(4*Dt.post.eddy*(2))))-rgamma(1000,0.01,0.01))

cs1ed.hatd <- apply(cs1edd, 2, mean)
cs2ed.hatd <- apply(cs2edd, 2, mean)

quartz()
par(mfrow=c(1,2))
plot(acetone18[,3]*54/24,col="green")
points(csp.hatd[,1],col="blue")
points(cs1ed.hatd,col="purple")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
plot(acetone18[,4]*54/24,col="green")
points(csp.hatd[,2],col="blue")
points(cs2ed.hatd,col="purple")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
#PF
logyteddy10048<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy1data22.csv")
logyteddy20048<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy2data22.csv")
thetaeddyjulia0048<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmceddydata22.csv")
apply(thetaeddyjulia0048[1000:2000,], 2, function(x){quantile(x,c(0.025,0.5,0.975))})[,1]/120
apply(thetaeddyjulia0048[1000:2000,], 2, function(x){quantile(x,c(0.025,0.5,0.975))})[,2]
prededdy1.julia0048<-apply(logyteddy10048[1000:2000,], 2, mean)
prededdy2.julia0048<-apply(logyteddy20048[1000:2000,], 2, mean)
print(mseedpf48<-(sum(((acetone18[,3]*54/24)-(prededdy1.julia0048))^2)+
                  sum((acetone18[,4]*54/24-(prededdy2.julia0048))^2))/178)
print(mseedst48<-(sum((acetone18[,3]*54/24-(csp.hatd[,1]))^2)+
                  sum((acetone18[,4]*54/24-(csp.hatd[,2]))^2))/178)
loglikpfeddyd<-matrix(0,1000,N)
for(i in 1:N){
  for(j in 1:1000)
    loglikpfeddyd[j,i]<-log(dmvnorm(logacetone00484[i,],cbind(log((logyteddy10048[j+1000,i])),log((logyteddy20048[j+1000,i]))),
                                   (diag(0.1,2))))
}

waic(loglikpfeddyd[,1:N])
quartz()
par(mfrow=c(1,2))
plot(prededdy1.julia0048)
points(exp(logacetone00484[,1]),col="red")
plot(prededdy2.julia0048)
points(exp(logacetone00484[,2]),col="red")

#both plots
quartz()
par(mfrow=c(1,2))
# plot(acetone18[,3]*54/24,col="green",ylab="Acetone Concentrations L1",xlab="Time", main="Location 1",cex=0.8)
# points(prededdy1.julia0048, col="red")
# points(csp.hatd[,1],col="blue")
# points(cs1ed.hatd,col="purple")
plot(acetone18[,3]*54/24,col="green",ylab="Acetone Concentrations L1",xlab="Time", main="Location 1",cex=0.8,type="l",lwd=3)
#lines(prededdy1.julia0048, col="red",lwd=3)
lines(csp.hatd[,1],col="blue",lwd=3,lty=6)
lines(predkfeddy1d,col="chocolate",lwd=3,lty=2)
lines(cs1ed.hatd,col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue", "chocolate","purple"), lty=c(1,6,2,3), legend=c("Measurements","Model a",
                                                           "Model b","Smoothed"),pch=15,cex=0.75,box.lwd = 0)
# plot(acetone18[,4]*54/24,col="green",ylab="Acetone Concentrations L2",xlab="Time", main="Location 2",cex=0.8)
# points(prededdy2.julia0048,col="red")
# points(csp.hatd[,2],col="blue")
# points(cs2ed.hatd,col="purple")
plot(acetone18[,4]*54/24,col="green",ylab="Acetone Concentrations L2",xlab="Time", main="Location 2",cex=0.8,type="l",lwd=3,
     ylim=c(0,163))
#lines(prededdy2.julia0048,col="red",lwd=3)
lines(csp.hatd[2:89,2],col="blue",lwd=3,lty=6)
lines(predkfeddy2d,col="chocolate",lwd=3,lty=2)
lines(cs2ed.hatd,col="purple",lwd=3,lty=3)
legend("bottomright",col=c("green","blue", "chocolate","purple"), lty=c(1,6,2,3), legend=c("Measurements","Model a",
                                                            "Model b","Smoothed"),pch=15,cex=0.75,box.lwd = 0)

N<-length(acetonedata[[1]][,1])-1
logacetone001<-matrix(0,N+1,2)
logacetone001<-log(acetonedata[[1]][,2:3]+0.5)

cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone001[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone001[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(30,1000)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone001"=logacetone001, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))/120
quantile(gprime.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+91+i])
    
  }
}
csp.hat
quartz()
plot(acetonedata[[1]][,2],csp.hat[,1])
abline(0,1)
quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[1]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[1]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)

#k=0.00098
h=1
N<-length(acetonedata[[2]][,1])-1
logacetone0009<-matrix(0,N+1,2)
logacetone0009<-log(acetonedata[[2]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone0009[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone0009[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone0009"=logacetone0009, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
quartz()
plot(acetonedata[[2]][,2],csp.hat[,1])
abline(0,1)
#quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[2]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[1]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
#k=0.00104
h=1
N<-length(acetonedata[[4]][,2])-1
logacetone00104<-matrix(0,N+1,2)
logacetone00104<-log(acetonedata[[4]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone00104[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00104[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone00104"=logacetone00104, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+85+i])
    
  }
}
csp.hat
plot(acetonedata[[4]][,2],csp.hat[,1])
abline(0,1)
quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[4]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[4]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)

#k=0.0043
h=1
N<-length(acetonedata[[17]][,2])-1
logacetone0043<-matrix(0,N+1,2)
logacetone0043<-log(acetonedata[[17]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone0043[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone0043[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone0043"=logacetone0043, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
plot(acetonedata[[17]][,2],csp.hat[,1])
abline(0,1)
quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[17]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[17]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)

#k=0.00345
h=1
N<-length(acetonedata[[16]][,2])-1
logacetone00345<-matrix(0,N+1,2)
logacetone00345<-log(acetonedata[[16]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone00345[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00345[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone00345"=logacetone00345, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
plot(acetonedata[[16]][,2],csp.hat[,1])
abline(0,1)
quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[16]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[16]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)

#k=0.0028
h=1
N<-length(acetonedata[[15]][,2])-1
logacetone0028<-matrix(0,N+1,2)
logacetone0028<-log(acetonedata[[15]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone0028[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone0028[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone0028"=logacetone0028, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
plot(acetonedata[[15]][,2],csp.hat[,1])
abline(0,1)
#quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[15]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[15]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
#k=0.00272
N<-length(acetonedata[[14]][,2])-1
logacetone00272<-matrix(0,N+1,2)
logacetone00272<-log(acetonedata[[14]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone00272[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00272[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone00272"=logacetone00272, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
plot(acetonedata[[14]][,2],csp.hat[,1])
abline(0,1)
#quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[14]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[14]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
#k=0.00251
N<-length(acetonedata[[13]][,2])-1
logacetone00251<-matrix(0,N+1,2)
logacetone00251<-log(acetonedata[[13]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone00251[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00251[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone00251"=logacetone00251, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
plot(acetonedata[[13]][,2],csp.hat[,1])
abline(0,1)
#quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[13]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[15]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
#k=0.00222
N<-length(acetonedata[[12]][,2])-1
logacetone00222<-matrix(0,N+1,2)
logacetone00222<-log(acetonedata[[12]][,2:3]+0.5)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone00222[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00222[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone00222"=logacetone00222, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}

csp.hat
plot(acetonedata[[12]][,2],csp.hat[,1])
abline(0,1)
#quartz()
#par(mfrow=c(1,2))
plot(acetonedata[[12]][,2],col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)

plot(acetonedata[[15]][,3],col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)



##pf plot
bimodalDistFunc <- function (n,cpct, mu1, mu2, sig1, sig2) {
  y0 <- rnorm(n,mean=mu1, sd = sig1)
  y1 <- rnorm(n,mean=mu2, sd = sig2)
  
  flag <- rbinom(n,size=1,prob=0.5)
  y <- y0*(1 - flag) + y1*flag 
}

bimodalData <- bimodalDistFunc(n=100,0.5,20,-3, 1,1)
quartz()
par(mfrow=c(2,1))
plot(density(bimodalData),xaxt="n",pch=0,yaxt="n",main="",ylab="",bty="n",xlab="")
plot(density(bimodalDistFunc(n=100,0.5,23,-3, 1,1)),xaxt="n",yaxt="n",main="",xlab="",ylab="",bty="n")
