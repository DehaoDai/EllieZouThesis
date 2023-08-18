##-------------------------Read Data and Data processing----------------------##
X = read.csv('realdata_processed.csv')
date = X[,1]
X = X[,-1]
X = data.frame(X)
nnn = colnames(X)
##----------------------------------------------------------------------------##

##----------------------Loading required packages-----------------------------##
library(glmnet)
set.seed(914812)
T = 120  ## moving window approach, 30 years as window size                         
M = nrow(X)-T  ## predict sample size 
N = 226 ## number of columns
##----------------------------------------------------------------------------##

## The following values are initialized for storing prediction difference, R^2, 
## predicted values, obtained from different methods like FARMTS, LASSO, RIDGE,
## Sample Mean, PCR, respectively. 
##----------------------Initialize M-by-N Matrix------------------------------##
Pred.FARMTS = matrix(0,M,N)                     ## FARMTS for estimation and prediction
Pred.LASSO = matrix(0,M,N)                      ## Lasso for prediction
Pred.MEAN = matrix(0,M,N)                       ## MEAN for prediction
Pred.PCR = matrix(0,M,N)                        ## PCR for prediction
Pred.RIDGE=matrix(0,M,N)                        ## Ridge Regression for prediction
##----------------------------------------------------------------------------##

## Out-of-sample R^2 
##----------------------------------------------------------------------------##
R2.FARMTS = numeric(N)                          ## FARMTS
R2.LASSO = numeric(N)                           ## LASSO
R2.PCR = numeric(N)                             ## PCR
R2.RIDGE = numeric(N)                           ## Ridge
##----------------------------------------------------------------------------##

## Prediction values
##----------------------------------------------------------------------------##
Prvalue_FARMTS=matrix(0,M,N)                    ## FARMTS
Prvalue_LASSO=matrix(0,M,N)                     ## LASSO
Prvalue_PCR=matrix(0,M,N)                       ## PCR
Prvalue_MEAN=matrix(0,M,N)                      ## Sample MEAN
Prvalue_RIDGE=matrix(0,M,N)                     ## Ridge Regression
##----------------------------------------------------------------------------##

## Predictions are conducted Next
#j = which(nnn %in% 'HOUSTW') #HOUSTW
#j = which(nnn %in% 'IPNCONGD')
for(j in 1:length(nnn)){
D = X[,-j]
Y = X[,j]
for(i in 1:M){
  idx = i:(i+T-1) ## Moving window prediction, every window has length T
  x = D[idx,] ## Training Data Covariate
  y = Y[idx]  ## Training Data Response
  x.pred = D[i+T,] ## Prediction Covariate
  y.pred = Y[i+T] ## Prediction Response
  
  ##----------------------------Do Data Normalization----------------------------##
  x.mu = colMeans(x)
  x.sd = as.numeric(apply(x,2,sd))
  y.mu = mean(y)
  Prvalue_MEAN[i,j]=y.mu 
  y.sd = sd(y)
  x = t((t(x)-x.mu)/x.sd)
  y = (y-y.mu)/y.sd
  x.pred = (x.pred-x.mu)/x.sd
  x.pred = unlist(x.pred)
  y.pred = (y.pred-y.mu)/y.sd
  ##-----------------------------------------------------------------------------##
  
  ## Sample Mean
  Pred.MEAN[i,j] = (y.pred*y.sd)^2
  
  ## Lasso
  cv.fit = cv.glmnet(x,y,intercept=FALSE)
  lambda.fit = cv.fit$lambda.min
  fit_model = glmnet(x,y,intercept=FALSE,lambda=lambda.fit)  ## Lasso Estimation
  beta.hat = as.vector(fit_model$beta)
  Pred.LASSO[i,j] = (y.pred-sum(x.pred*beta.hat))^2*y.sd^2
  Prvalue_LASSO[i,j]=sum(x.pred*beta.hat)*y.sd+y.mu    ## Lasso Prediction
  
  ## Ridge
  cv.fit = cv.glmnet(x,y,intercept=FALSE,alpha=0)	
  lambda.fit = cv.fit$lambda.1se
  fit_model = glmnet(x,y,intercept=FALSE,lambda=lambda.fit,alpha=0)  ## Ridge estimation
  beta.hat = as.vector(fit_model$beta)
  Pred.RIDGE[i,j] = (y.pred-sum(x.pred*beta.hat))^2*y.sd^2
  Prvalue_RIDGE[i,j]=sum(x.pred*beta.hat)*y.sd+y.mu ## Ridge Prediction
  
  ##Factor Estimation
  Sigma.x = tcrossprod(x)/T      #covariance matrix
  eigenx = eigen(Sigma.x)        
  eigvec = eigenx$vectors        #Eigen vector
  eigvalue = eigenx$values       #Eigen values
  K.hat = max(2,which.min(diff(log(eigvalue[1:10]))))      ## Factor estimation
  F.hat = eigvec[,1:K.hat]*sqrt(T)        #Estimated Factor
  B.hat = t(t(F.hat)%*%x)/T          #Estimated Factor Loading
  U.hat = x-F.hat%*%t(B.hat)        #Estimated idiosyncratic component
  
  # FARMTS
  lmY.F = lm(y~F.hat-1)           
  gamma.hat = coef(lmY.F)    
  Y.tilde = resid(lmY.F)                             
  cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	 
  lambda.fit.U = cv.fit.U$lambda.min #tuning parameter selection
  fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U) ## Estimation
  beta.hat.U = as.vector(fit.U$beta) #Obtain the fitted beta
  
  lmx.B = lm(x.pred~B.hat-1)
  f.pred = coef(lmx.B)
  u.pred = resid(lmx.B)
  
  Pred.FARMTS[i,j] = (y.pred-sum(f.pred*gamma.hat)-sum(u.pred*beta.hat.U))^2*y.sd^2 #FARM prediction Difference with the true value
  Prvalue_FARMTS[i,j]=(sum(f.pred*gamma.hat)+sum(u.pred*beta.hat.U))*y.sd+y.mu  #FARM predicted value
  Pred.PCR[i,j] = (y.pred-sum(f.pred*gamma.hat))^2*y.sd^2     #PCR prediction Difference with the true value
  Prvalue_PCR[i,j]=sum(f.pred*gamma.hat)*y.sd+y.mu           #PCR predicted value
  
}

R2.FARMTS[j] = 1-sum(Pred.FARMTS[,j])/sum(Pred.MEAN[,j])	
print(R2.FARMTS[j])                    #Output R^2 value for FARM
R2.PCR[j] = 1-sum(Pred.PCR[,j])/sum(Pred.MEAN[,j])
print(R2.PCR[j])                   #Output R^2 value for PCR
R2.LASSO[j] = 1-sum(Pred.LASSO[,j])/sum(Pred.MEAN[,j])
print(R2.LASSO[j])                 #Output R^2 value for LASSO
R2.RIDGE[j] = 1-sum(Pred.RIDGE[,j])/sum(Pred.MEAN[,j])
print(R2.RIDGE[j])                 #Output R^2 value for RIDGE
}
date <- as.Date(date, "%m/%d/%Y")
plot(date[(T+1):nrow(X)],Y[(T+1):nrow(X)], type = "l",xlab = "",ylab = "", col = 1)
a = 3
lines(date[(T+1):nrow(X)],Prvalue_FARMTS[,j], lty = 2, col = "red")
lines(date[(T+1):nrow(X)],Prvalue_LASSO[,j], lty = a, col = "green")
lines(date[(T+1):nrow(X)],Prvalue_RIDGE[,j], lty = a, col = "purple")
lines(date[(T+1):nrow(X)],Prvalue_MEAN[,j],lty = a,col = "blue")
lines(date[(T+1):nrow(X)],Prvalue_PCR[,j], lty = a, col = "orange")

