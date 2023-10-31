#####
### product of normals
#####

vg=1
vh=10
lpdf<-function(x) dnorm(x,0,sd = sqrt(vg),log = TRUE)+dnorm(x,0,sd = sqrt(vh),log = TRUE)
xseq=seq(-20,20,by=0.001)


sum(dnorm(xseq,0,sd = sqrt(vg),log = FALSE)*dnorm(xseq,0,sd = sqrt(vh),log = FALSE))*(xseq[2]-xseq[1])
1/(sqrt(2*pi*(vg+vh)))


curve(dnorm(x,sd=sqrt(vh)),xlim=c(-10,10),col="red",ylim=c(0,0.45),ylab="Density")
curve(dnorm(x,sd=sqrt(vg)),xlim=c(-10,10),col="green",add=TRUE)
curve((sqrt(2*pi*(vg+vh)))*dnorm(x,0,sd = sqrt(vg),log = FALSE)*dnorm(x,0,sd = sqrt(vh),log = FALSE),add=TRUE,col="blue")

xp1<-rnorm(1e4,sd=sqrt(vg))
w1<-dnorm(xp1,sd=sqrt(vh))
xp2<-rnorm(1e4,sd=sqrt(vh))
w2<-dnorm(xp2,sd=sqrt(vg))

hist(w1,freq=FALSE,xlim=c(0,0.4),col="green",xlab="w*",main="")
hist(w2,freq=FALSE,add=TRUE,col="red")
plot(density(sample(xp1,1e4,replace = TRUE,prob = w1)))
lines(density(sample(xp2,1e4,replace = TRUE,prob = w2)))

sum(w1)^2/sum(w1^2)
sum(w2)^2/sum(w2^2)

mean(w1)
mean(w2)


sig2g<-1
sig2h<-10

(2*sig2h*sqrt(sig2g)+sig2g*sqrt(sig2g))<
  ((2*sig2g*sqrt(sig2h)+sig2h*sqrt(sig2h)))


x1_mean<-x2_mean<-
  x1_mcse<-x2_mcse<-
  x1_pval<-x2_pval<-
  ess1<-ess2<-0

for(i in 1:1e4){
  
  xx1<-rnorm(1e4,sd=sqrt(sig2g))
  w1<-dnorm(xx1,sd=sqrt(sig2h))
  ess1[i]<-sum(w1)^2/sum(w1^2)
  x1<-sample(xx1,1e4,prob=w1,replace = TRUE)
  # hist(x1,prob=TRUE)
  # lines(density(y),col=2,lwd=2)
  
  
  xx2<-rnorm(1e4,sd=sqrt(sig2h))
  w2<-dnorm(xx2,sd=sqrt(sig2g))
  ess2[i]<-sum(w2)^2/sum(w2^2)
  x2<-sample(xx2,1e4,prob=w2,replace = TRUE)
  # hist(x2,prob=TRUE)
  # lines(density(y),col=2,lwd=2)
  
  a<-0
  x1_mean[i]<-sum(w1*(xx1<a)/sum(w1))
  x2_mean[i]<-sum(w2*(xx2<a)/sum(w2))
  x1_mcse[i]<-sum(((xx1<a)-w1*(xx1<a)/sum(w1))^2*w1^2)/sum(w1)^2
  x2_mcse[i]<-sum(((xx2<a)-w2*(xx2<a)/sum(w2))^2*w2^2)/sum(w2)^2
  
  x1_pval[i]<-ks.test(y,x1)$p.value
  x2_pval[i]<-ks.test(y,x2)$p.value
  
  if(i%%1e3==0) print(i/1e4)
}

boxplot(x1_mean,x2_mean)
boxplot(x1_mcse,x2_mcse)
boxplot(x1_pval,x2_pval)
boxplot(ess1,ess2)


mean(x1)
mean(x2)
mean(y)

var(x1)
var(x2)
var(y)

mean(w1)
mean(w2)

var(w1)
var(w2)

a<-0
sum(w1*(xx1<a)/sum(w1))
sum(w2*(xx2<a)/sum(w2))
sum(((xx1<a)-w1*(xx1<a)/sum(w1))^2*w1^2)/sum(w1)^2
sum(((xx2<a)-w2*(xx2<a)/sum(w2))^2*w2^2)/sum(w2)^2

xseq<-seq(-5,5,by=0.00001)
sum(dnorm(xseq,sd=sqrt(sig2g))*dnorm(xseq,sd=sqrt(sig2h)))*0.00001


ks.test(y,x1)
ks.test(y,x2)


############################
######### HMM
############################

##1:10,1
##2:1,10

n<-1000
nn<-1e2
np<-1e1
sse1<-sse2<-matrix(0,nrow=100,n)
ess1_mean<-ess2_mean<-numeric(100)

for(ii in 1:100){

set.seed(ii+10000)

sig2w<-10
sig2v<-1

beta=0.5
x<-rnorm(1,sd=sqrt(sig2w/(1-beta^2)))
y<-rnorm(1,x,sd=sqrt(sig2v))
for(i in 2:n){
  x[i]<-rnorm(1,beta*x[i-1],sd=sqrt(sig2w))
  y[i]<-rnorm(1,x[i],sd=sqrt(sig2v))
}
plot.ts(y,ylab="y",ylim=c(-15,15))
plot.ts(x,ylab="x",ylim=c(-15,15))


#### single ####
# xp1_save<-matrix(0,nn,n)
# ess1<-NULL
# Vinv=(1-beta^2)/sig2w+1/sig2v
# for(j in 1:nn){
# xp<-rnorm(1,y[1]/sig2v/Vinv,sd = sqrt(1/Vinv))
# for(i in 2:n){
# xx<-rnorm(np,y[i],sqrt(sig2v))
# w<-dnorm(xx,beta*xp[i-1],sd=sqrt(sig2w))
# ess1<-c(ess1,sum(w)^2/sum(w^2))
# xp[i]<-sample(xx,1,prob = w)
# }
# xp1_save[j,]<-xp
# # if(j%%100==0) print(round(j/nn,100))
# }
# 
# 
# xp2_save<-matrix(0,nn,n)
# ess2<-NULL
# Vinv=(1-beta^2)/sig2w+1/sig2v
# for(j in 1:nn){
#   xp<-rnorm(1,y[1]/sig2v/Vinv,sd = sqrt(1/Vinv))
#   for(i in 2:n){
#     xx<-rnorm(np,beta*xp[i-1],sqrt(sig2w))
#     w<-dnorm(y[i],xx,sd=sqrt(sig2v))
#     ess2<-c(ess2,sum(w)^2/sum(w^2))
#     xp[i]<-sample(xx,1,prob = w)
#   }
#   xp2_save[j,]<-xp
#   if(j%%100==0) print(round(j/nn,100))
# }


ess11<-NULL
xp11_save<-matrix(0,nn,n)
Vinv=(1-beta^2)/sig2w+1/sig2v
xp11_save[,1]<-rnorm(nn,y[1]/sig2v/Vinv,sd = sqrt(1/Vinv))
for(i in 2:n){
  xx<-matrix(rnorm(nn,y[i],sqrt(sig2v)))
  w<-dnorm(xx,beta*xp11_save[,i-1],sd=sqrt(sig2w))
  ess11<-c(ess11,apply(w,2,function(x) sum(x)^2/sum(x^2)))
  ind_vals<-apply(w,2, function(x) sample(length(xx),nn,prob=x,replace = TRUE))
  xp11_save[,i]<-xx[ind_vals]
  }

ess22<-NULL
xp22_save<-matrix(0,nn,n)
Vinv=(1-beta^2)/sig2w+1/sig2v
xp22_save[,1]<-rnorm(nn,y[1]/sig2v/Vinv,sd = sqrt(1/Vinv))
for(i in 2:n){
    xx<-rnorm(nn,beta*xp22_save[,i-1],sqrt(sig2w))
    w<-dnorm(y[i],xx,sd=sqrt(sig2v))
    ess22<-c(ess22,sum(w)^2/sum(w^2))
    xp22_save[,i]<-sample(xx,nn,prob = w,replace=TRUE)
  }


# sqrt(mean((apply(xp11_save,2,mean)-x)^2))
# sqrt(mean((apply(xp22_save,2,mean)-x)^2))
# mean(ess11)
# mean(ess22)
# 
# plot.ts(apply(xp1_save,2,mean))
# lines(apply(xp2_save,2,mean),col=2)
# 
# ind<-25
# hist(xp1_save[,ind],prob=TRUE)
# hist(xp2_save[,ind],prob=TRUE,add=TRUE,col=2)
# 
# boxplot(ess22,ess11,names = c("original","reverse"))
# hist(ess2,prob=TRUE,col=3,xlim=c(0,np))
# hist(ess1,prob=TRUE,add=TRUE,col=2)

sse1[ii,]<-(apply(xp11_save,2,mean)-x)^2
sse2[ii,]<-(apply(xp22_save,2,mean)-x)^2

ess1_mean[ii]<-mean(ess11)
ess2_mean[ii]<-mean(ess22)

print(ii/100)
}

mean(sqrt(apply(sse1,2,mean)))
mean(sqrt(apply(sse2,2,mean)))
apply(cbind(ess1_mean,ess2_mean),2,mean)

xp1_smooth<-
xp2_smooth<-matrix(0,nn,n)

for(j in 1:nn){
  
  xp1<-sample(xp1_save[,n],1)
  xp2<-sample(xp2_save[,n],1)
  for(i in 1:(n-1)){
    
    w<-dnorm(xp1[i],beta*xp1_save[,n-i],sqrt(sig2w))
    xp1[i+1]<-sample(xp1_save[,n-i],1,prob=w)
    
    w<-dnorm(xp2[i],beta*xp2_save[,n-i],sqrt(sig2w))
    xp2[i+1]<-sample(xp2_save[,n-i],1,prob=w)
  }
  xp1_smooth[j,]<-rev(xp1)
  xp2_smooth[j,]<-rev(xp2)
  
  if(j%%100==0) print(round(j/nn,100))
}


ffbs = function(y,m0,c0,phi,V,W){
  n = length(y) 
  mf = rep(0,n)
  Cf = rep(0,n)
  m = m0
  C = C0
  for (t in 1:n){
    a = phi*m
    R = phi^2*C+W
    C = 1/(1/R+1/V)
    m = C*(a/R+y[t]/V)
    Cf[t] = C
    mf[t] = m
  }
  x  = rep(0,n)
  x[n] = rnorm(1,mf[n],sqrt(Cf[n]))
  for (t in (n-1):1){
    var = 1/(1/Cf[t]+phi^2/W)
    mean = var*(mf[t]/Cf[t]+phi*x[t+1]/W)
    x[t] = rnorm(1,mean,sqrt(var))
  }
  return(x)
}


m0   =   0.00
C0   = 10.00

ffbs_smooth<-t(replicate(1e5,ffbs(y,m0,C0,beta,V = sig2v,W=sig2w)))

ind<-25
hist(ffbs_smooth[,ind],prob=1)
hist(xp1_smooth[,ind],add=TRUE,col=1,prob=1)
hist(xp2_smooth[,ind],add=TRUE,col=2,prob=1)

ks_xp1<-ks_xp2<-0
for(i in 1:n){
  ind<-i
ks_xp1[i]<-ks.test(ffbs_smooth[,ind],xp1_smooth[,ind])$p.value
ks_xp2[i]<-ks.test(ffbs_smooth[,ind],xp2_smooth[,ind])$p.value
}
boxplot(ks_xp1,ks_xp2)
mean(ks_xp1)
mean(ks_xp2)

sqrt(mean((apply(xp1_smooth,2,mean)-apply(ffbs_smooth,2,mean))^2))
sqrt(mean((apply(xp2_smooth,2,mean)-apply(ffbs_smooth,2,mean))^2))

sqrt(mean((apply(xp1_smooth,2,mean)-x)^2))
sqrt(mean((apply(xp2_smooth,2,mean)-x)^2))

##################################
##################################
### poisson
##################################

n<-1000
nn<-1e2
np<-1e2
sse1<-sse2<-sse3<-matrix(0,nrow=100,n)

ii=2
set.seed(ii)
for(ii in 1:100){
  
  phi<-0.5
  v=1 #sd
  theta=10
  x0<-(rnorm(1,sd=v/sqrt(1-phi^2)))
  x<-rnorm(1,x0,sd=v)
  y<-rpois(1,theta*exp(x))
  for(i in 2:n){
    x[i]<-rnorm(1,phi*x[i-1],sd=v)
    y[i]<-rpois(1,theta*exp(x[i]))
  }
  plot.ts(y)
  lines(x,col=2)
  
  plot.ts(y,ylab="y")
  plot.ts((x),ylab="x")
  
  var(x)
  var(y)
  
  fsamp<-function(x) rlnorm(np,phi*x,sd=v)
  fdens<-function(x) dlnorm(x,phi*x1,sd=v) #log scale
  gdens<-function(x) dpois(x,theta*(x1))
  gsamp<-function(n) rpois(n,theta*(x1))
  
  xf_ess<-NULL
  xf<-matrix(0,nn,n)
  for(j in 1:nn){
    xp<-x[1]
    for(i in 2:n){
      x1<-fsamp(xp[i-1])

      w<-gdens(y[i])
      xf_ess<-c(xf_ess,sum(w)^2/sum(w^2))
      xp[i]<-log(sample(x1,1,replace=TRUE,prob = w))
    }
    xf[j,]<-xp
    if(j%%100==0) print(round(j/nn,100))
  }

  xr_ess<-NULL
  xr<-matrix(0,nn,n)
  for(j in 1:nn){
    xp<-x[1]
    for(i in 2:n){
      xx<-rgamma(np,y[i]+1,theta)
      x1<-xp[i-1]
      w<-fdens(xx)
      xr_ess<-c(xr_ess,sum(w)^2/sum(w^2))
      xp[i]<-log(sample(xx,1,prob=w))
      
    }
    xr[j,]<-xp
    if(j%%100==0) print(round(j/nn,100))
  }
  

  sse1[ii,]<-(apply(xr,2,mean)-x)^2
  sse2[ii,]<-(apply(xf,2,mean)-x)^2
  print(ii/100)
}


##########################################
#### stoch vol  ### 
### code based on:
### https://hedibert.org/r-code-to-our-journal-of-forecasting-review-paper/
##########################################

n<-1000
nn<-1e2
np<-1e3
sse1<-sse2<-matrix(0,nrow=1000,n)
xseq=seq(-50,50,0.1)

ii=104
# while(ii<=100){
for(ii in 100:1000){
set.seed(ii)

    # For t=1,...,n
  #
  #    y[t]|x[t]   ~ N(0,exp(x[t]/2))
  #    x[t]|x[t-1] ~ N(alpha+beta*x[t-1],tau2)
  #
  # and
  #
  #    x[0]       ~ N(m0,C0)
  #    alpha|tau2 ~ N(b0[1],tau2*B0[1])
  #    beta|tau2  ~ N(b0[2],tau2*B0[2])
  #    tau2       ~ IG(nu0/2,nu0*tau20/2)
  
  

  alpha=0.600239
  phi<-0.8722187
  tau2<-0.09059663
  tau<-sqrt(tau2)
  theta=1
  x0<-rnorm(1,alpha/(1-phi),sd=tau/sqrt(1-phi^2))
x<-rnorm(1,alpha+phi*x0,sd=tau)
y<-rnorm(1,0,sd=sqrt(theta)*exp(x/2))
for(i in 2:n){
x[i]<-rnorm(1,alpha+phi*x[i-1],sd=tau)
y[i]<-rnorm(1,0,sd=sqrt(theta)*exp(x[i]/2))
}
plot.ts(y)
plot.ts(x,col=1)

var(x)
var(y)

gdens<-function(y) dnorm(y,0,sd=sqrt(theta)*exp(x1/2))
gsamp<-function(n) rnorm(n,0,sd=sqrt(theta)*exp(x1/2))
fsamp<-function(x) rnorm(np,alpha+phi*x,sd=tau)
fdensx<-function(x) dnorm(x,alpha+phi*x1,sd=tau)
fsampx<-function(n) rnorm(n,alpha+phi*x1,sd=tau)


xseq=seq(-50,50,by=0.0001)
x1=x[i-1]
gf<-function(x) (dnorm(x,alpha+phi*x1,sd=tau)*dnorm(y[i],0,sd=sqrt(theta)*exp(x/2)))
curve(gf(x),xlim=c(-10,10))
sum(gf(xseq)*0.0001)

xf_ess<-numeric(n)
xf<-matrix(0,nn,n)
for(j in 1:nn){
  xp<-x[1]
  for(i in 2:(n)){
    x1<-fsamp(x[i-1])
    w<-gdens(y[i])
    xi=x[i-1]
    xf_ess[i]<-xf_ess[i]+sum(w)^2/sum(w^2)/nn
    xp[i]<-sample(x1,1,replace=TRUE,prob = w)
  xf[j,]<-xp
  if(j%%10==0) print(round(j/nn,100))
}

xr_ess<-numeric(n)
xr<-matrix(0,nn,n)
for(j in 1:nn){
  xp<-x[1]
  for(i in 2:(n)){
    x1<-x[i-1]
    y1<-y[i]
    c=y1^2/2
    xx=-log(rgamma(np,1/2,c))
    w<-fdensx(xx)
    xr_ess[i]<-xr_ess[i]+sum(w)^2/sum(w^2)/nn
    xp[i]<-sample(xx,1,prob=w)
    [j,]<-xp
  if(j%%10==0) print(round(j/nn,100))
}

sse1[ii,]<-(apply(xr,2,mean)-x)^2
sse2[ii,]<-(apply(xf,2,mean)-x)^2
print(ii/100)
ii=ii+1

}

sqrt(mean((apply(xf,2,mean)-x)^2))
sqrt(mean((apply(xr,2,mean)-x)^2))


mean(sqrt(apply(sse1[,,drop=FALSE],2,mean)))
mean(sqrt(apply(sse2[,,drop=FALSE],2,mean)))

plot.ts(xf[1,])
ss<-sample(nn,1e2)
for(i in 1:nn){
lines(xf[ss[i],])
}
lines(x,col=2,lwd=2)

plot.ts(xr[1,])
ss<-sample(nn,1e2)
for(i in 1:nn){
  lines(xr[ss[i],])
}
lines(x,col=2,lwd=2)

plot.ts(apply(xf,2,mean))
lines(apply(xr,2,mean),col=2)
lines(x,col=3,lwd=2)

hist(xf[,2],prob=TRUE);abline(v=x[2],lwd=2)
hist(xr[,2],prob=TRUE);abline(v=x[2],lwd=2)

hist(xf[,50],prob=TRUE);abline(v=x[50],lwd=2)
hist(xr[,50],prob=TRUE);abline(v=x[50],lwd=2)

hist(xf[,n],prob=TRUE);abline(v=x[n],lwd=2)
hist(xr[,n],prob=TRUE);abline(v=x[n],lwd=2)

boxplot(xf_ess,xr_ess,names = c("original","reverse"))
hist(xf_ess,prob=TRUE,col=3,xlim=c(0,np))
hist(xr_ess,prob=TRUE,add=TRUE,col=2)

hist(exp(rnorm(1e6,1,1)),breaks=100,prob=TRUE,xlim=c(0,100))
hist(rlnorm(1e6,1,1),breaks=100,prob=TRUE,xlim=c(0,100))

# Prior hyperparameters
# ---------------------
m0    = 0
C0    = tau2/(1-phi^2)
nu0   = 3
tau20 = 0.01
b0    = c(0,1)
B0    = c(10,10)/tau20
sC0   = sqrt(C0)
sB0   = sqrt(B0)

set.seed(246521)
N        = 1e4
# delta    = 0.975
delta    = 10
xs       = rnorm(N,m0,sC0)
# tau2s    = 1/rgamma(N,nu0/2,nu0*tau20/2)
tau2s=rep(tau2,N)
taus     = sqrt(tau2s)
# alphas   = rnorm(N,b0[1],taus*sB0[1])
alphas   = rep(0.0,N)
# betas    = rnorm(N,b0[2],taus*sB0[2])
betas    = rep(phi,N)
nus      = 2:20
nnu      = length(nus)
pf       = array(0,c(1+nnu,n,4,3))
# like     = matrix(0,n,1+nnu)

LW = function(y,alphas,betas,tau2s,xs,delta){
  n  = length(y)
  N  = length(xs)
  # quants = array(0,c(n,4,3))
  quants = matrix(0,N,n)
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  # like = rep(0,n)
  for (t in 1:n){
    # like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    if (delta<1){
      ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    }else{
      ms1 = ms
    }
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]  	 
    # Storing quantiles
    # quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    # quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    # quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    # quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
    quants[,t]<-xs
  }
  # return(list(like=like,quants=quants))
  return(quants)
}


run      = LW(y,alphas,betas,tau2s,xs,delta)



hist(xf[,2]);abline(v=x[1],lwd=2)
hist(xr[,2]);abline(v=x[1],lwd=2)
hist(run[,2]);abline(v=x[1],lwd=2)

hist(xf[,50]);abline(v=x[1],lwd=2)
hist(xr[,50]);abline(v=x[1],lwd=2)
hist(run[,50]);abline(v=x[50],lwd=2)

hist(xf[,n]);abline(v=x[n],lwd=2)
hist(xr[,n]);abline(v=x[n],lwd=2)
hist(run[,n]);abline(v=x[n],lwd=2)


rmse_xf<-rmse_xr<-rmse_lw<-0
for(i in 1:n){
  rmse_xf[i]<-sqrt(mean((xf[,i]-x[i])^2))
  rmse_xr[i]<-sqrt(mean((xr[,i]-x[i])^2))
  rmse_lw[i]<-sqrt(mean((run[,i]-x[i])^2))
  
  }
mean(rmse_xf)
mean(rmse_xr)
mean(rmse_lw)
mean(rmse_xf<rmse_lw)

plot.ts(apply(xr,2,mean))
lines(apply(run,2,mean),col=2)
lines(x,col=3)
