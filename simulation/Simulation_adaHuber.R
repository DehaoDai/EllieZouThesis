#-----------Load Packages ------------#
library(adaHuber)
library(glmnet)
#Function adahubert: Generate data and Do model estimation; 
#------------------
#Input: r: ratio, a number between 0-1
#p: dimension
#K number of factors
#theta: true parameter
#s: sparsity level
#----------------
AHubert = function(r,p,K,theta,s){
  gamma0 = 0.5*rep(1,K) #true gamma
  beta0 = c(rep(0.5,s),rep(0,p-s)) #true beta
  B0 = matrix(runif(p*K,-1,1),nrow=p) #factor loading 
  ##-------------- Generate Data---------------------------
  n = ceiling(log(p)/(r/(K+s))^(1+1/theta))#number of observations
  F = matrix(rnorm(n*K),nrow=n) #generate factors
  U = matrix(rnorm(n*p),nrow=n) #generate idiosyncratic component
  X = F%*%t(B0)+U #generate design matrix
  phi_1 = 0.1  #generate noise distribution e
  innovation = rt(n, df  = (2+theta))
  e = arima.sim(model = list(ar = phi_1), n = n, innov = innovation)
  e = as.numeric(e)
  Y = F%*%gamma0+U%*%beta0+e #generate response variable
  
  ##--------------- Estimate Factor Model --------------------------
  Sigma = tcrossprod(X)/n #covariance matrix
  Eig = eigen(Sigma,only.values=TRUE) 
  Eigval = Eig$values #eigenvales
  K.est = which.min(diff(log(Eigval[1:20]))) #estimate number of factors
  Svd = svd(X,nu=K.est,nv=0) 
  Eigvec = Svd$u #eigenvectors
  F.hat = sqrt(n)*Eigvec[,1:K.est] #estimate factors
  mol = lm(X~F.hat-1) #fit model
  U.hat = resid(mol) #estimate idiosyncratic component
  Y_tilde = Y - 1/T*F.hat%*%t(F.hat)%*%Y
  fit1 = glmnet(U.hat, Y_tilde, intercept=FALSE, lambda=cv.glmnet(U.hat, Y_tilde, intercept=FALSE)
                $lambda.1se)
  lambda_beta1 = fit1$lambda
  beta_hat1 = as.vector(fit1$beta)
  ##------------Estimate FARM Model------------------------------
  X.hat = cbind(F.hat, U.hat)
  Huber.est = adaHuber.cv.lasso(X.hat,Y) #Estimate FARM mode using Huber regression 
  Huber.beta = Huber.est$coef 
  Huber.beta = Huber.beta[-(1:(K.est+1))] #estimated beta
  z = sum(abs(Huber.beta-beta0)) #record the difference
  z2 = sum(abs(beta_hat1 - beta0))
  return(c(z, z2)) #return the difference
  return(z)
}

#----------------Main setting---------------------------------
p = 1000; #dimension
K = 2; #number of factors
s = 3; #sparsity level
theta = 1 #parameter
N = 7 #
ratio = seq(0.4,0.7,length=N)
re = 100
L = matrix(0,N,re)#record the results 500 times
L2 = matrix(0, N, re)
#
for(i in 1:N){
  print(i)
  r = ratio[i] #for every r in ratio, use AHubert 500 times to compute the parameter differences.
  #system.time(AHubert(r,p,K,theta,s))
  tmp = replicate(re,AHubert(r,p,K,theta,s))
  #tmp = as.numeric(tmp)
  L[i,] = tmp[1,] #record
  L2[i,]= tmp[2,]
  }
write.csv(L,"p1000_K2_s3_t5_1.csv") #record 500 replication results.
write.csv(L2,"p1000_K2_s3_t10_2.csv") #record 500 replication results.

ratio <- seq(0.4, 0.7, by = 0.05) 


meanline = rowMeans(L)
sd_1 = apply(L[,1:re],1, sd)
lowerbound = meanline - sd_1
upperbound = meanline + sd_1

meanline2 = rowMeans(L2)
sd_2 = apply(L2[,1:re],1, sd)
lowerbound2 = meanline2 - sd_2
upperbound2 = meanline2 + sd_2

## draw 
data <- data.frame(x = rep(ratio, 2),
                   y = c(meanline, meanline2),
                   group = rep(c("Huber","FARM"), each = length(ratio)),
                   lower_interval = c(lowerbound, lowerbound2),
                   upper_interval = c(upperbound, upperbound2))

color_palette <- c("red", "blue")
ggplot(data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = lower_interval, ymax = upper_interval, fill = group), alpha = 0.3) +
  geom_line(data = subset(data, group == "Huber"), aes(color = group)) +
  geom_line(data = subset(data, group == "FARM"), aes(color = group)) +
  scale_color_manual(values = color_palette, guide = guide_legend(override.aes = list(fill = NA))) +
  ylab(TeX("$|\\widehat{\\beta}_\\lambda - \\beta^*|_1$"))+
  xlab(TeX("$S\\sqrt{(\\log p)^{1+2\\nu}/n}, \\nu = 1$"))+
  scale_x_continuous(breaks = seq(0.4, 0.7, by = 0.05))+
  scale_y_continuous(breaks = seq(0, 2, by = 0.2))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = c(0.1,0.7))+
  theme(legend.title = element_blank())
















