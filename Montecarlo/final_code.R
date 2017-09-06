#Partie Importance Sampling
library(MASS)
library(mnormt)

n = 506
p = 13
sigma = 1
lambda = 0.5
kappa = 0.5
Nsim = 100000

data(housing,package = "MASS")
colnames(Boston) = c("Y", "X1", "X2", "X3", "X4", "X5", "x6", "X7", "X8", "X9", "X10", "X11", "X12", "X13")
Boston = as.data.frame(Boston)
for(j in 2:p)#rendre les données centrées réduites
{
  esperance = mean(Boston[, j])
  ecart_type= sqrt(var(Boston[, j]))
  for(i in 1:n)
  {
    Boston[i, j] = (Boston[i, j] - esperance)/ecart_type 
  }
}
esperance=mean(Boston[, 1])
for (i in 1:n)
{Boston[i,1]=Boston[i,1]-esperance}


#Calcul des matrices
XXt = 0
for(i in 1:n)
{
  XXt = XXt + t(as.matrix(Boston[i, 2:14])) %*% as.matrix(Boston[i, 2:14])#somme des matrice des XtX sur toute les observations
}
XXt_inv = solve(XXt)#inverse de la matrice

XY = 0
for(i in 1:n)
{
  XY = XY + t(as.matrix(Boston[i, 2:14])) %*% as.matrix(Boston[i, 1])#somme des vecteurs XY sur toute les observations
}

m= XXt_inv %*% XY #vecteur moyenne de la loi pour quand kappa=0
m = as.vector(m)

Gamma= sigma^2 *XXt_inv
Gamma= as.matrix(Gamma)#matrice de variance covariance  pour quand kappa=0
for(i in 2:p)
{
  for(j in 1:(i-1))
  {
    Gamma[j, i] = Gamma[i, j]
  }
}

g <- function(moy=m,variance=Gamma)#loi a posteriori pour k=0
{
  return(mvrnorm(1, moy, variance))#multivariate distribution, on simule des Beta selon la loi a posteriori de beta si kappa=0
}
# maintenant on simule nos beta suivant g
beta_g= matrix(NA, ncol = Nsim, nrow = p)
for(i in 1:Nsim)
{
  beta_g[,i]= g()
}


w = NULL # on calcul les ponderations
for(i in 1:Nsim)
{
  somme_bi= 0
  for(j in 1:p)
  {
    somme_bi = somme_bi + abs(beta_g[j, i])^kappa
  }
  w[i] = exp(-lambda *somme_bi )
}

est_beta = 0
for(i in 1:Nsim)
{
  est_beta = est_beta + w[i]*beta_g[, i]
}

est_beta=est_beta/sum(w)
est_beta

ess=0#on calcule l'ess
for(i in 1:Nsim)
{
  ess = ess + w[i]^2
}
ess=sum(w)^2/ess
ess

#Partie metropolis


library(tseries)
data(Boston,package = "MASS")
#########Calculate loglikelihood for the conditional law
LogLikelihood<-function(X,Y,beta,lambda,sigma,K){
  log_lik<-0
  for(i in 1:dim(X)[1])
    log_lik=log_lik-((Y[i]-t(X[i,])%*%beta)/sigma)^2*0.5
  for(i in 1:dim(X)[2])
    log_lik=log_lik-lambda*abs(beta[i])^K
  return(log_lik)
}
##############Transform data in order to have zero mean and unit variance
n = 506
p = 13
for(j in 1:p)
{
  esperance = mean(Boston[, j])
  ecart_type= sqrt(var(Boston[, j]))
  for(i in 1:n)
  {
    Boston[i, j] = (Boston[i, j] - esperance)/ecart_type 
  }
}
esperance=mean(Boston[, p+1])
for (i in 1:n)
{Boston[i,p+1]=Boston[i,p+1]-esperance}
#########Get data
X<-Boston[,1:dim(Boston)[2]-1]
Y<-Boston[,dim(Boston)[2]]
X<-data.matrix(X)
#########

lambda<-0.5
K<-0.5
sigma<-1
p<-13
width<-0.05
beta0<-runif(p,-1,1)
beta0=matrix(beta0,nrow=1,ncol=p)
Nsim<-10000
scale_g<-0.0005
#########Jumping width to avoid high autocorrelation
jumping_width<-50
#########Location of first sample
T_0<-100
#########Algorithm with uniform distribution
Metropolis_RW<-function(Nsim,X,Y,K,lambda,sigma,beta0,p,width){
  rho_data<-c()
  beta<-matrix(beta0,1,p)
  for (i in 2:Nsim){
    epsilon<-runif(p,-width,width)
    rho<-exp(LogLikelihood(X,Y,beta[i-1,]+epsilon,lambda,sigma,K)-LogLikelihood(X,Y,beta[i-1,],lambda,sigma,K))
    rho_data<-c(rho_data,min(rho,1))
    print(c(min(rho,1))) #####monitor the ratio of acceptation
    beta<-rbind(beta,beta[i-1,]+epsilon*(runif(1)<rho))
  }
  return(list(beta,rho_data)) ####return ratios of acceptation and simulations
}
###########Calibration of parametre width
###########start time
timestart<-Sys.time()
trajectory<-Metropolis_RW(Nsim,X,Y,K,lambda,sigma,beta0,p,width)
##########Caculate the time
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)
##########Get rate of acceptation
rho<-as.matrix(as.matrix(as.data.frame(trajectory[2])))
mean(rho)
##########Get beta
tra_u<-data.frame(trajectory[1])
################Determine the convergence of simulations
################trace plot
plot(tra_u$X1,type='p',pch='.',col='red',ylab='Value',xlab='Simulation',main='Trace plot of simulations')
################Plot cumulated average
cumulated_aver<-function(data){
  means<-c()
  accu_sum<-0
  for (i in 1:length(data)){
    accu_sum=accu_sum+data[i]
    means=c(means,accu_sum/i)
  }
  plot(means,type='p',pch='.',col='red',ylab='Mean',xlab='Simulations',main='Accumulated mean of simulations')
}
cumulated_aver(tra_u$X1)
###############Determine T_0 by this plot
###############Autocorrelation plot
acf(tra_u$X1,lag=100,col='red',main='Autocorrelation plot')
###########Determine jumping width by this plot
###########Get chosen simulations
new_tra_u=tra_u[seq(T_0,dim(tra_u)[1],by=50),]
###########Check acf
acf(new_tra_u$X1,col='red',main='Autocorrelation plot of new trace')
est_metro_u<-apply(new_tra_u,2,mean)
error_u<-0
for (i in 1:p)
{error_u=error_u+(est_metro_u[i]-real_beta[i])^2}
#############################################
##################Use of gaussian distribution
library(mvtnorm)
Metropolis_RW_Normal<-function(Nsim,X,Y,K,lambda,sigma,beta0,p,width){
  rho_data<-c()
  beta<-matrix(beta0,1,p)
  for (i in 2:Nsim){
    epsilon<-c(rmvnorm(1,mean=rep(0,p),sigma=diag(p)*width))
    rho<-exp(LogLikelihood(X,Y,beta[i-1,]+epsilon,lambda,sigma,K)-LogLikelihood(X,Y,beta[i-1,],lambda,sigma,K))
    rho_data<-c(rho_data,min(rho,1))
    print(min(rho,1))
    beta<-rbind(beta,beta[i-1,]+epsilon*(runif(1)<rho))
  }
  return(list(beta,rho_data))
}
timestart<-Sys.time()
trajectory<-Metropolis_RW_Normal(Nsim,X,Y,K,lambda,sigma,beta0,p,scale_g)
###############Caculate the time
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

rho<-as.matrix(as.matrix(as.data.frame(trajectory[2])))
mean(rho)
tra_g<-data.frame(trajectory[1])
acf(tra_g$X1,lag=100,col='red',main='Autocorrelation plot')
##################
plot(tra_g$X1,type='p',pch='.',col='red',ylab='Value',xlab='Simulation',main='Trace plot of simulations')
##################
cumulated_aver(tra_g$X1)
##############
new_tra_g=tra_g[seq(T_0,dim(tra_g)[1],by=40),]
###########
acf(new_tra_g$X1,col='red',main='Autocorrelation plot of new trace')
###########Estimation
est_metro_g<-apply(new_tra_g,2,mean)
error_g<-0
for (i in 1:p)
{error_g=error_g+(est_metro_g[i]-real_beta[i])^2}
#############Comparason with importance sampling
#############Algorithm of importance sampling

library(MASS)
library(mnormt)

n = 506
p = 13
sigma = 1
lambda = 0.5
kappa = 0.5
Nsim = 10000

data(housing,package = "MASS")
colnames(Boston) = c("X1", "X2", "X3", "X4", "X5", "x6", "X7", "X8", "X9", "X10", "X11", "X12", "X13","Y")
Boston = as.data.frame(Boston)



XXt = 0
for(i in 1:n)
{
  XXt = XXt + t(as.matrix(Boston[i, 1:13])) %*% as.matrix(Boston[i, 1:13])#somme des matrice des XtX sur toute les observations
}
XXt_inv = solve(XXt)#inverse de la matrice

XY = 0
for(i in 1:n)
{
  XY = XY + t(as.matrix(Boston[i, 1:13])) %*% as.matrix(Boston[i, p+1])#somme des vecteurs XY sur toute les observations
}

m= XXt_inv %*% XY #vecteur moyenne de la loi pour quand kappa=0
m = as.vector(m)

Gamma= sigma^2 *XXt_inv
Gamma= as.matrix(Gamma)#matrice de variance covariance  pour quand kappa=0
for(i in 2:p)
{
  for(j in 1:(i-1))
  {
    Gamma[j, i] = Gamma[i, j]
  }
}

m

# nous codons la loi g selon la loi normal(m,Gamma)
g <- function(moy=m,variance=Gamma)
{
  return(mvrnorm(1, moy, variance))#multivariate distribution, on simule des Beta selon la loi a posteriori de beta si kappa=0
}

# maintenant on simule nos beta suivant g
#######start time
timestart<-Sys.time()
beta_g= matrix(NA, ncol = Nsim, nrow = p)
for(i in 1:Nsim)
{
  beta_g[,i]= g()
}

w = NULL # 
for(i in 1:Nsim)
{
  somme_bi= 0
  for(j in 1:p)
  {
    somme_bi = somme_bi + abs(beta_g[j, i])^kappa
  }
  w[i] = exp(-lambda *somme_bi )
}

est_beta = 0
for(i in 1:Nsim)
{
  est_beta = est_beta + w[i]*beta_g[, i]
}

est_beta=est_beta/sum(w)
est_beta
#########end time
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

#Save the "real beta"
#real_beta<-est_beta

ess=0
for(i in 1:Nsim)
{
  ess = ess + w[i]^2
}
ess=sum(w)^2/ess
ess

#######Calcule l'erreur quadratique
error_is<-0
for (i in 1:p)
{error_is=error_is+(est_beta[i]-real_beta[i])^2}
#Partie Gibbs

