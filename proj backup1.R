library(MCMCpack)

data1=c(2,1,4,1,6,1,9,1,9,1,9,1,13,1,14,1,18,1,23,1,31,1,32,1,33,1,34,1,43,1,10,0,14,0,14,0,16,0,17,0,18,0,18,0,19,0,20,0,20,0,21,0,21,0,23,0,24,0,29,0,29,0,30,0,30,0,31,0,31,0,31,0,33,0,35,0,37,0,40,0,41,0,42,0,42,0,44,0,46,0,48,0,49,0,51,0,53,0,54,0,54,0,55,0,56,0)
yh=matrix(data = data1,53,2,byrow = T) # data of hormone treatment group
data2=c(1,1,4,1,6,1,7,1,13,1,24,1,25,1,35,1,35,1,39,1,1,0,1,0,3,0,4,0,5,0,8,0,10,0,11,0,13,0,14,0,14,0,15,0,17,0,19,0,20,0,22,0,24,0,24,0,24,0,25,0,26,0,26,0,26,0,28,0,29,0,29,0,32,0,35,0,38,0,39,0,40,0,41,0,44,0,45,0,47,0,47,0,47,0,50,0,50,0,51,0)
yc=matrix(data = data2, ncol=2,byrow = T) # data of control group

del.h=sum(yh[,2]) # sum over all delta-h
del.c=sum(yc[,2]) # sum over all delta-c
x.h= sum(yh[,1])  # sum over all xi-h
x.c= sum(yc[,1])  # sum over all xi-c

log.post <- function(d,c,b,a,v,u) {
  (26+a)*u-(27+a)*log(1+exp(u))+(16+b)*v-(17+b)*log(1+exp(v))-(1233+c)*(exp(u)/(exp(u)+1))-
    (1526+d)*(exp(u)/(exp(u)+1))*(exp(v)/(exp(v)+1))
}

niter = 1*10^6                     # number of iterations
theta = rep(0,niter)               # define (hyper)parameters
u     = rep(0,niter)
v     = rep(0,niter)
tau   = rep(0,niter)
a     = 3
b     = 1 
c     = 60
d     = 120  
j     = 0                         # counter for extra runs to make theta in (0,1)
k     = 0
l     = 0                

# Gibbs Sampler Implementation
set.seed(1)                                  # set random seed

theta[1] = runif(1, min = 0, max = 1)        # Initialization of (hyper)parameters
u[1]     = log(theta[1]/(1-theta[1]))
tau[1]   = runif(1, min = 0, max = 1)
v[1]     = log(tau[1]/(1-tau[1]))


for (i in 2:niter) {                        # for-loop for Gibbs sampler
  
  # update tau in (0,1) and theta in (0,1) together by M-H algo
  temp1=1
  temp2=1
  while (temp1==1 ) {
    v[i] <- v[i-1]+rnorm(1,0,0.05)
    
    R1 <- exp(log.post(d,c,b,a,v[i],u[i-1])
              -log.post(d,c,b,a,v[i-1],u[i-1]))
    if(R1<1){
      if(rbinom(1,1,R1)==0)   {v[i] = v[i-1];  k = k+1}
    }
    temp1=exp(v[i])/(1+exp(v[i]))         # transform v back to tau
  }
  while (temp2==1) {
    u[i] <- u[i-1]+rnorm(1,0,0.05)
    R2 <- exp(log.post(d,c,b,a, v[i],u[i])
              -log.post(d,c,b,a, v[i],u[i-1]))
    if(R2<1){
      if(rbinom(1,1,R2)==0)   { u[i]=u[i-1]; l = l+1}
    }
    temp2=exp(u[i])/(1+exp(u[i]))         # transform u back to theta
    j = j+1
  }
  tau[i]=temp1
  theta[i]=temp2
}
#_______________________________________________________________________
#Plots of the output
burn.in=1:(niter/10)

# sample path plot
par(mfrow=c(2,4))                     # sample path plot by parts
plot(theta[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for theta')
abline(v=100, lty = 3)

plot(u[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for u')
abline(v=100, lty = 3)

plot(tau[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for tau')
abline(v=100, lty = 3)

plot(v[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for v')
abline(v=100, lty = 3)

plot(a[seq(1,niter,by=1000)],type="l",ylab = '',  xlab ='-th 1000 iterations',col = rgb(0,0,1,0.8), main = 'sample path for a')
abline(v=100, lty = 3)

plot(b[seq(1,niter,by=1000)],type="l",ylab = '',  xlab ='-th 1000 iterations',col = rgb(0,1,0,0.8), main = 'sample path for b')
abline(v=100, lty = 3)

plot(c[seq(1,niter,by=1000)],type="l",ylab = '',  xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'sample path for c')
abline(v=100, lty = 3)

plot(d[seq(1,niter,by=1000)],type="l",ylab = '',  xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'sample path for d')
abline(v=100, lty = 3)

# cusum plot for mean of all parameters
theta.temp <- theta[-burn.in]
theta.ntemp <- length(theta.temp)
theta.mtemp <- mean(theta.temp)
theta.temp.1 <- cumsum(theta.temp-theta.mtemp)

u.temp <- u[-burn.in]
u.ntemp <- length(u.temp)
u.mtemp <- mean(u.temp)
u.temp.1 <- cumsum(u.temp-u.mtemp)

tau.temp <- tau[-burn.in]
tau.ntemp <- length(tau.temp)
tau.mtemp <- mean(tau.temp)
tau.temp.1 <- cumsum(tau.temp-tau.mtemp)

v.temp <- v[-burn.in]
v.ntemp <- length(v.temp)
v.mtemp <- mean(v.temp)
v.temp.1 <- cumsum(v.temp-v.mtemp)

a.temp <- a[-burn.in]
a.ntemp <- length(a.temp)
a.mtemp <- mean(a.temp)
a.temp.1 <- cumsum(a.temp-a.mtemp)

b.temp <- b[-burn.in]
b.ntemp <- length(b.temp)
b.mtemp <- mean(b.temp)
b.temp.1 <- cumsum(b.temp-b.mtemp)

c.temp <- c[-burn.in]
c.ntemp <- length(c.temp)
c.mtemp <- mean(c.temp)
c.temp.1 <- cumsum(c.temp-c.mtemp)

d.temp <- d[-burn.in]
d.ntemp <- length(d.temp)
d.mtemp <- mean(d.temp)
d.temp.1 <- cumsum(d.temp-d.mtemp)

par(mfrow=c(2,4))
plot(theta.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for theta')
abline(h=0, lty=3)
plot(u.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for u')
abline(h=0, lty=3)
plot(tau.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(0,0,1,0.8), main = 'cusum plot for tau')
abline(h=0, lty=3)
plot(v.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for v')
abline(h=0, lty=3)
plot(a.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations',col = rgb(0,1,0,0.8), main = 'cusum plot for a')
abline(h=0, lty=3)
plot(b.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'cusum plot for b')
abline(h=0, lty=3)
plot(c.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'cusum plot for c')
abline(h=0, lty=3)
plot(d.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'cusum plot for d')
abline(h=0, lty=3)

# autocorrelation plot
par(mfrow=c(2,4))
acf(theta.temp[seq(1,niter-100000,by=1000)], lag.max=40, type="correlation", xlab="N", main = 'autocorrelation of theta');
acf(u.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of u');
acf(tau.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of tau');
acf(v.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of v');
acf(a.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,0,1,0.8),main = 'autocorrelation of a');
acf(b.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,1,0,0.8),main = 'autocorrelation of b');
acf(c.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,0,1,0.8),main = 'autocorrelation of c');
acf(d.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,0,1,0.8),main = 'autocorrelation of d');

# Histogram
par(mfrow=c(2,4))
hist(theta.temp,nclass = 40,probability = T)
hist(u.temp,nclass = 40, probability = T)
hist(tau.temp,nclass = 40, probability = T)
hist(v.temp,nclass = 40, probability = T)
hist(a.temp,nclass = 40, probability = T)
hist(b.temp,nclass = 40, probability = T)
hist(c.temp,nclass = 40, probability = T)
hist(d.temp,nclass = 40, probability = T)

