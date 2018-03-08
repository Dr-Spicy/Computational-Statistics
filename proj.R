library(MCMCpack)

data1=c(2,1,4,1,6,1,9,1,9,1,9,1,13,1,14,1,18,1,23,1,31,1,32,1,33,1,34,1,43,1,10,0,14,0,14,0,16,0,17,0,18,0,18,0,19,0,20,0,20,0,21,0,21,0,23,0,24,0,29,0,29,0,30,0,30,0,31,0,31,0,31,0,33,0,35,0,37,0,40,0,41,0,42,0,42,0,44,0,46,0,48,0,49,0,51,0,53,0,54,0,54,0,55,0,56,0)
yh=matrix(data = data1,53,2,byrow = T) # data of hormone treatment group
data2=c(1,1,4,1,6,1,7,1,13,1,24,1,25,1,35,1,35,1,39,1,1,0,1,0,3,0,4,0,5,0,8,0,10,0,11,0,13,0,14,0,14,0,15,0,17,0,19,0,20,0,22,0,24,0,24,0,24,0,25,0,26,0,26,0,26,0,28,0,29,0,29,0,32,0,35,0,38,0,39,0,40,0,41,0,44,0,45,0,47,0,47,0,47,0,50,0,50,0,51,0)
yc=matrix(data = data2, ncol=2,byrow = T) # data of control group

del.h=sum(yh[,2]) # sum over all delta-h
del.c=sum(yc[,2]) # sum over all delta-c
x.h= sum(yh[,1])  # sum over all xi-h
x.c= sum(yc[,1])  # sum over all xi-c

mu.h.1 = sum(yh[,1] * yh[,2])/sum(yh[,2])          # mean of h group with label 1
mu.h.0 = sum(yh[,1] * (yh[,2]-1)/sum((yh[,2]-1)))  # mean of h group with label 1
mu.c.1 = sum(yc[,1] * yc[,2])/sum(yc[,2])          # mean of c group with label 1
mu.c.0 = sum(yc[,1] * (yc[,2]-1)/sum((yc[,2]-1)))  # mean of c group with label 1


log.post <- function(d,c,b,a,tau,u) {
   (26+a)*u-(27+a)*log(1+exp(u))+(15+b)*log(tau)-(1233+c)*(exp(u)/(exp(u)+1))-
    (1526+d)*(exp(u)/(exp(u)+1))*tau
}

niter = 1*10^6                     # number of iterations
theta = rep(0,niter)               # define (hyper)parameters
u     = rep(0,niter)
tau   = rep(0,niter)
a     = 3
b     = 1 
c     = 60
d     = 120  
j     = 0                         # counter for extra runs to make theta in (0,1) but no (0,1]
l     = 0                

# Gibbs Sampler Implementation
set.seed(1)                                  # set random seed

theta[1] = runif(1, min = 0, max = 1)        # Initialization of (hyper)parameters
u[1]     = log(theta[1]/(1-theta[1]))
tau[1]   = 1



for (i in 2:niter) {                        # for-loop for Gibbs sampler
  # update tau by rgamma 
  tau[i] <- rgamma(1,shape = 16+b, rate = (exp(u[i])/(exp(u[i])+1))*(1526+d))
  # update theta in (0,1)  by M-H algo
  temp2=1
  while (temp2==1) {
    u[i] <- u[i-1]+rnorm(1,0,0.5)
    R2 <- exp(log.post(d,c,b,a,tau[i],u[i])
              -log.post(d,c,b,a,tau[i],u[i-1]))
    if(R2<1){
      if(rbinom(1,1,R2)==0)   { u[i]=u[i-1]; l = l+1}
    }
    temp2=exp(u[i])/(1+exp(u[i]))         # transform u back to theta
    j = j+1
  }
  
  theta[i]=temp2
}

#_______________________________________________________________________
#Plots of the output
burn.in=1:(niter/10)

# sample path plot
par(mfrow=c(1,3))                     # sample path plot by parts
plot(theta[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for theta')
abline(v=100, lty = 3)

plot(u[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for u')
abline(v=100, lty = 3)

plot(tau[seq(1,niter,by=1000)],type="l", ylab = '', xlab ='-th 1000 iterations', main = 'sample path for tau')
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



par(mfrow=c(1,3))
plot(theta.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for theta')
abline(h=0, lty=3)
plot(u.temp.1[seq(1,niter,by=1000)],  ylab ='',type="l",xlab ='-th 1000 iterations', main = 'cusum plot for u')
abline(h=0, lty=3)
plot(tau.temp.1[seq(1,niter,by=1000)], ylab ='', type="l",xlab ='-th 1000 iterations',col = rgb(0,0,1,0.8), main = 'cusum plot for tau')
abline(h=0, lty=3)


# autocorrelation plot
par(mfrow=c(1,3))
acf(theta.temp[seq(1,niter-100000,by=1000)], lag.max=40, type="correlation", xlab="N", main = 'autocorrelation of theta');
acf(u.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of u');
acf(tau.temp[seq(1,niter-100000,by=1000)], lag.max=20, type="correlation", xlab="N",col = rgb(1,0,0,0.8), main = 'autocorrelation of tau');

# Histogram
par(mfrow=c(1,3))
hist(theta.temp,nclass = 40,probability = T)
hist(u.temp,nclass = 40, probability = T)
hist(tau.temp,nclass = 40, probability = T)

