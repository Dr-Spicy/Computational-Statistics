niter = 10^6                       # number of iterations
theta = rep(0,niter)               # define (hyper)parameters
tau   = rep(0,niter)
a     = rep(0,niter)
b     = rep(0,niter) 
c     = rep(0,niter)
d     = rep(0,niter)  

# Gibbs Sampler Implementation
set.seed(1)                                  # set random seed

theta[1] = runif(1, min = 0, max = 1)        # Initialization of (hyper)parameters
tau[1]   = runif(1, min = 0, max = 1)
a[1]     = 3
b[1]     = 1
c[1]     = 60
d[1]     = 120

for (i in 2:niter) {                        # for-loop for Gibbs sampler
  # update d
  d[i]<-rexp(1, rate = tau[i-1]*theta[i-1])
  # update c
  c[i]<-rexp(1, rate = theta[i-1])
  # update b
  b[i]<-rexp(1, rate = -log(tau[i-1]))
  # update a
  a[i]<-rexp(1, rate = -log(theta[i-1]))
  # update tau in (0,1)
  temp <- rgamma(1,shape = 16+b[i], rate = theta[i-1]*(1526+d[i]))
  tau[i] <- temp/(temp+1)
  # update theta in (0,1)
  temp <- rgamma(1,shape = 26+a[i], rate = tau[i]*(1526+d[i])+1233+c[i])
  theta[i] <- temp/(temp+1)
}

#_______________________________________________________________________
#Plots of the output
burn.in=1:(niter/10)

# sample path plot
par(mfrow=c(2,3))                     # sample path plot by parts
plot(theta[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for theta')
abline(v=100, lty = 3)

#plot(u[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for u')
#abline(v=100, lty = 3)

plot(tau[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for tau')
abline(v=100, lty = 3)

#plot(v[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for v')
#abline(v=100, lty = 3)

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

#u.temp <- u[-burn.in]
#u.ntemp <- length(u.temp)
#u.mtemp <- mean(u.temp)
#u.temp.1 <- cumsum(u.temp-u.mtemp)

tau.temp <- tau[-burn.in]
tau.ntemp <- length(tau.temp)
tau.mtemp <- mean(tau.temp)
tau.temp.1 <- cumsum(tau.temp-tau.mtemp)

#v.temp <- v[-burn.in]
#v.ntemp <- length(v.temp)
#v.mtemp <- mean(v.temp)
#v.temp.1 <- cumsum(v.temp-v.mtemp)

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

par(mfrow=c(2,3)) 
plot(theta.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for theta')
abline(h=0, lty=3)
#plot(u.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for u')
#abline(h=0, lty=3)
plot(tau.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(0,0,1,0.8), main = 'cusum plot for tau')
abline(h=0, lty=3)
#plot(v.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for v')
#abline(h=0, lty=3)
plot(a.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations',col = rgb(0,1,0,0.8), main = 'cusum plot for a')
abline(h=0, lty=3)
plot(b.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'cusum plot for b')
abline(h=0, lty=3)
plot(c.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'cusum plot for c')
abline(h=0, lty=3)
plot(d.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(1,0,0,0.8), main = 'cusum plot for d')
abline(h=0, lty=3)

# autocorrelation plot
par(mfrow=c(2,3)) 
acf(theta.temp[seq(1,niter-100000,by=1000)], lag.max=40, type="correlation", xlab="N", main = 'autocorrelation of theta');
#acf(u.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of u');
acf(tau.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of tau');
#acf(v.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of v');
acf(a.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,0,1,0.8),main = 'autocorrelation of a');
acf(b.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,1,0,0.8),main = 'autocorrelation of b');
acf(c.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,0,1,0.8),main = 'autocorrelation of c');
acf(d.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N", col = rgb(0,0,1,0.8),main = 'autocorrelation of d');

# Histogram
par(mfrow=c(2,3)) 
hist(theta.temp,nclass = 40,probability = T)
#hist(u.temp,nclass = 40, probability = T)
hist(tau.temp,nclass = 40, probability = T)
#hist(v.temp,nclass = 40, probability = T)
hist(a.temp,nclass = 40, probability = T)
hist(b.temp,nclass = 40, probability = T)
hist(c.temp,nclass = 40, probability = T)
hist(d.temp,nclass = 40, probability = T)

1/theta.mtemp
1/(theta.mtemp*tau.mtemp)
