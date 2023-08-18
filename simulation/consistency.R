##-------------------------Loading Packages----------------------------------##
library(glmnet)
library(vars)
library(stats)
library(MASS)
library(ggplot2)
library(latex2exp)

##-------------------------Basic Setting to generate data--------------------##
start_time <- proc.time()
p = 100   # dimension
K = 2 # true number of factors
gamma1 = 0.5*rep(1,K) # generate true gamma
S = 3 # sparsity of beta
beta = c(rep(0.5,S),rep(0,p-S)) # generate true beta
B <- matrix(runif(p*K,-1, 1),nrow=p) # generate factor loading
ratio <- seq(0.3, 0.6, by = 0.05) # initialize the x-axis of figure 1.
re = 500 # 500 replications
mat <- matrix(nrow = length(ratio), ncol = 500) # record results 
mat2 <- matrix(nrow= length(ratio), ncol = 500) #record results
d = 2

for (i in 1:length(ratio)){
  dis<-c()
  dis2<-c()
  T = ceiling(S^2*log(p)^d/ratio[i]^2)
  for (j in 1:re){ #500 replications
    #factor
    F<-matrix(rnorm(T*K),nrow=T) 
    
    #idiosyncratic component
    U<-matrix(rnorm(T*p),nrow=T) 
    
    # design matrix
    X = F%*%t(B)+U
    
    phi_1 = 0.1 # 0.3
    innovation = rnorm(T, mean = 0, sd = 0.5)
    e = arima.sim(model = list(ar = phi_1), n = T, innov = innovation)
    e = as.numeric(e)
    
    Y= F%*%gamma1 + U%*%beta + e #Data is generated following Gaussian noise 
    ##---Use eigenvector of covariance matrix to estimate latent factors------##
    cov<-X%*%t(X)/T
    eigvec<-eigen(cov)$vectors #eigenvectors
    eigval<-eigen(cov)$values  #eigenvalues
    K_est<-which.min(diff(log(eigval[1:10]))) #estimate factor numbers using eigenvalue ratios
    hatf<-sqrt(T)*eigvec[,1:K_est] #estimated factors
    hatB<-t(1/T*t(hatf)%*%X) #Estimated Factor Loading
    hatU<-X-hatf%*%t(hatB)   #Estimated Idiosyncratic component
    Y_tilde<-Y-1/T*hatf%*%t(hatf)%*%Y
    fit1 = glmnet(hatU, Y_tilde, intercept=FALSE,
                  lambda=cv.glmnet(hatU,Y_tilde,intercept=FALSE)$lambda.1se)
    fit2= glmnet(X, Y,intercept=FALSE,
                 lambda=cv.glmnet(hatU,Y_tilde,intercept=FALSE)$lambda.1se)
    lambda_beta1<-fit1$lambda    #fitted beta
    beta_hat1<-as.vector(fit1$beta)
    
    lambda_beta2<-fit2$lambda   #fitted beta
    beta_hat2<-as.vector(fit2$beta)
    
    dis<-c(dis,sum(abs(beta_hat1-beta))) #record the difference between estimated FITS beta with true one
    dis2<-c(dis2,sum(abs(beta_hat2-beta))) ##record the difference between estimated LASSO beta with true one
  }
  mat[i,]<-dis    #Repeat 500 times and record as a matrix.
  mat2[i,]<-dis2
  write.csv(mat,"mat.csv")
  write.csv(mat2,"mat2.csv")
}

end_time <- proc.time()
elapsed_time <- end_time - start_time

# Print the elapsed time
print(elapsed_time)

##------------------------------Figure------------------------------##

## draw 
meanline = rowMeans(mat)
sd_1 = apply(mat[,1:re],1, sd)
lowerbound = meanline - sd_1
upperbound = meanline + sd_1

meanline2 = rowMeans(mat2)
sd_2 = apply(mat2[,1:re],1, sd)
lowerbound2 = meanline2 - sd_2
upperbound2 = meanline2 + sd_2

## draw 
data <- data.frame(x = rep(ratio, 2),
                   y = c(meanline, meanline2),
                   group = rep(c("FARM_TS","LASSO"), each = length(ratio)),
                   lower_interval = c(lowerbound, lowerbound2),
                   upper_interval = c(upperbound, upperbound2))

color_palette <- c("red", "blue")
ggplot(data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = lower_interval, ymax = upper_interval, fill = group), alpha = 0.3) +
  geom_line(data = subset(data, group == "FARM_TS"), aes(color = group)) +
  geom_line(data = subset(data, group == "LASSO"), aes(color = group)) +
  scale_color_manual(values = color_palette, guide = guide_legend(override.aes = list(fill = NA))) +
  ylab(TeX("$|\\widehat{\\beta}_\\lambda - \\beta^*|_1$"))+
  xlab(TeX("$S\\sqrt{(\\log p)^{1+2\\nu}/n}, \\nu = 1/2$"))+
  scale_x_continuous(breaks = seq(0.3, 0.6, by = 0.05))+
  scale_y_continuous(breaks = seq(0, 1.6, by = 0.2))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = c(0.8,0.7))+
  theme(legend.title = element_blank())
