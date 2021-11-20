#import data

load(file=url("https://mfasiolo.github.io/TOI/frogs.RData"))
head(frogs)

#Logistic function

P_l <- function(s,theta){
  e <- theta[1]
  beta <- theta[2]
  phi <- theta[3]
  p <- (exp(e*(phi - s)))/(1 + exp(beta*e*(phi - s)))
  return(p)
}

#Ricker Model

P_r <- function(s,theta){
  a <- theta[1]
  b <- theta[2]
  alpha <- theta[3]
  p <- b*(((s/a)*exp(1-s/a))^alpha)
  return(p)
}

#Frogs ordered by size

frogs <- as.data.frame(frogs)
frogs <- frogs[order(frogs["size"]),]
x <- frogs$size
y <- frogs$killed

#compute negative log-likelihood using the binomial model given

P_l_log <- function(s,y,theta){
  p <- P_l(s,theta)
  loglik <- sum(dbinom(y,10,p,log=TRUE))
  -loglik
}

P_r_log <- function(s,y,theta){
  p <- P_r(s,theta)
  loglik <- sum(dbinom(y,10,p,log=TRUE))
  -loglik
}

#minimise negative log-likelihood

beta0 <- 5
e0 <- -0.25
phi0 <- 13
theta01 <- c(e0, beta0, phi0)
qfit1 <- optim(theta01, P_l_log, method="BFGS", hessian=TRUE, s=x, y=y)
qfit1
print(P_l(qfit1$par[2],qfit1$par))

a0 <- 8
b0 <- 0.25
alpha0 <- 3
theta02 <- c(a0, b0, alpha0)
qfit2 <- optim(theta02, P_r_log, method="BFGS", hessian=TRUE, s=x, y=y)
qfit2

#plot values against data

plot(x,y,xlab="Size",ylab="Number Killed",main="Comparison of models")
lines(x,10*P_l(x, qfit1$par))
lines(x,10*P_r(x,qfit2$par), col="red")
legend("topright",95,legend=c("Logistic function","Ricker Model"),
       col=c("black","red"),lty=1:1,cex=0.8)

#Akaike information criterion

AIC1 <- 2*P_l_log(x,y,qfit1$par)+length(qfit1$par)
AIC2 <- 2*P_r_log(x,y,qfit2$par)+length(qfit2$par)
Comp <- c(AIC1,AIC2)
Comp

#Confidence intervals for model parameters

theta1.sd <- diag(solve(qfit1$hessian))^0.5
ci1 <- c(qfit1$par - 1.96*theta1.sd,qfit1$par + 1.96*theta1.sd)
ci1
conf1 <- c(ci1[1],ci1[4],ci1[2],ci1[5],ci1[3],ci1[6])
conf1
theta2.sd <- diag(solve(qfit2$hessian))^0.5
ci2 <- c(qfit2$par - 1.96*theta2.sd,qfit2$par + 1.96*theta2.sd)
ci2
conf2 <- c(ci2[1],ci2[4],ci2[2],ci2[5],ci2[3],ci2[6])
conf2

#generate sample for a given theta

Samp <- function(theta){
  vals=c()
  for (i in 1:28){
    val <- P_l(i,theta)
    numb <- rbinom(1,10,val)
     vals[i]=numb
  }   
 return(vals)
}

#plot simulation where beta=1 with data 

thetatest <- c(qfit1$par[1],1,qfit1$par[3])
plot(x,Samp(thetatest),type="b",col="red",xlab="Size",ylab="Number Killed",
     main="Data vs Beta=1")
lines(x,y,type="b",pch=22,)
legend("right",95,legend=c("Data","Beta=1"),
       col=c("black","red"),lty=1:1,cex=0.8)

#inverse function

gfunc <- function(theta){
  y <- max(P_l(x,qfit1$par))
  e <- theta[1]
  beta <- theta[2]
  phi <- theta[3]
  s <- phi - log(y) + log(exp(e)-y*exp(e*beta))
}

#s^* estimate

smax <- gfunc(qfit1$par)
smax

#calculate y value for g(theta)

yconst <- max(P_l(x,qfit1$par))
ehat <- qfit1$par[1]
betahat <- qfit1$par[2]
phihat <- qfit1$par[3]

#gradient of g

gradient <- (exp(ehat-yconst*betahat*exp(ehat*betahat)))/(exp(ehat)-yconst*exp(ehat*betahat))
gradbeta <- (yconst*ehat*exp(ehat*betahat))/(yconst*exp(ehat*betahat)-exp(ehat))
gradphi <- 1
gradgfunc <- c(gradient, gradbeta, gradphi)
vargfunc <- t(gradgfunc)%*%qfit1$hessian%*%gradgfunc
stderrorgfunc <- sqrt(vargfunc)

#confidence interval for s^*

upbound <- smax + 1.96*stderrorgfunc
lowbound <- smax - 1.96*stderrorgfunc
print(upbound)
print(lowbound)