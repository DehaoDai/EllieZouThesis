##---------------------------Read Data and Data processing--------------------##
X = read.csv('realdata_processed.csv')
X = X[,-1]
X = data.frame(X)
nnn = colnames(X)
##----------------------------------------------------------------------------##

##---------------------------Loading required packages------------------------##
library(glmnet)
library(elasticnet)
library(randomForest)
library(adaHuber)
set.seed(219812)
#T = 90  ## moving window approach, 90 quarters as window size  
#T = 120  ## moving window approach, 30 years as window size   
T = 150 ##moving window approach, 150 quarters as window size   
M = nrow(X)-T  ## predict sample size 
N = length(nnn) ## number of columns
##----------------------------------------------------------------------------##

## The following values are initialized for storing prediction difference, R^2, 
## predicted values, obtained from different methods like Huber FARM, Lasso, 
## Ridge, PCR, Sample Mean, respectively. 
##----------------------Initialize M-by-N Matrix------------------------------##
Pred.HUBER = matrix(0,M,N)                      ## Huber FARM for estimation and prediction
Pred.LASSO = matrix(0,M,N)                      ## Lasso for prediction
Pred.RIDGE = matrix(0,M,N)                      ## Ridge for prediction
Pred.PCR = matrix(0,M,N)                        ## PCR for prediction
Pred.MEAN = matrix(0,M,N)                       ## MEAN for prediction
##----------------------------------------------------------------------------##

## Store out-of-sample R^2 for several methods mentioned above
##----------------------------------------------------------------------------##
R2.HUBER = numeric(N)                           ## Huber FARM
R2.LASSO = numeric(N)                           ## LASSO
R2.RIDGE = numeric(N)                           ## RIDGE
R2.PCR = numeric(N)                             ## PCR
R2.MEAN = numeric(N)                            ## MEAN
##----------------------------------------------------------------------------##

## Store prediction values for several methods mentioned above
##----------------------------------------------------------------------------##
Prvalue_HUBER=matrix(0,M,N)                     ## Huber FARM
Prvalue_LASSO=matrix(0,M,N)                     ## LASSO
Prvalue_RIDGE=matrix(0,M,N)                     ## RIDGE
Prvalue_PCR=matrix(0,M,N)                       ## PCR
Prvalue_MEAN=matrix(0,M,N)                      ## Sample MEAN
##----------------------------------------------------------------------------##

#colcol =  as.numeric(which( nnn %in% c('PCDGx', "IPDMAT")))
#colcol = as.numeric(which( nnn %in%  "IPDMAT"))
colcol = as.numeric(which( nnn %in%  "PCDGx"))
D = X
for(j in colcol){
  X = D[,-j]
  Y = D[, j]
  for(i in 1:M){
    idx = i:(i+T-1) ## Moving window prediction, every window has length T
    x = X[idx,] ## Training Data Covariate
    y = Y[idx]  ## Training Data Response
    x.pred = X[i+T,] ## Prediction Covariate
    y.pred = Y[i+T] ## Prediction Response
    
    ##----------------------------Do Data Normalization-----------------------##
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
    ##------------------------------------------------------------------------##
    
    ## Sample Mean
    Pred.MEAN[i,j] = (y.pred*y.sd)^2
    
    ## Lasso
    cv.fit = cv.glmnet(x,y,intercept=FALSE)
    lambda.fit = cv.fit$lambda.min
    fit_model = glmnet(x,y,intercept=FALSE,lambda=lambda.fit)  ## Lasso Estimation
    beta.hat = as.vector(fit_model$beta)
    Pred.LASSO[i,j] = (y.pred-sum(x.pred*beta.hat))^2*y.sd^2
    Prvalue_LASSO[i,j]=sum(x.pred*beta.hat)*y.sd+y.mu    ## Lasso Prediction

    ## Ridge estimation
    cv.fit.x = cv.glmnet(x,y,intercept=FALSE,alpha=0)	
    lambda.fit.x = cv.fit.x$lambda.1se
    fit.x = glmnet(x,y,intercept=FALSE,lambda=lambda.fit.x,alpha=0)  ## Ridge estimation
    beta.hat.xr = as.vector(fit.x$beta)
    Pred.RIDGE[i,j] = (y.pred-sum(x.pred*beta.hat.xr))^2*y.sd^2
    Prvalue_RIDGE[i,j]=sum(x.pred*beta.hat.xr)*y.sd+y.mu #Ridge Prediction
    
    ##Factor Estimation
    Sigma.x = tcrossprod(x)/T      #covariance matrix
    eigenx = eigen(Sigma.x)        
    eigvec = eigenx$vectors        #Eigen vector
    eigvalue = eigenx$values       #Eigen values
    K.hat = max(2,which.min(diff(log(eigvalue[1:10]))))      ## Factor estimation
    F.hat = eigvec[,1:K.hat]*sqrt(T)        #Estimated Factor
    B.hat = t(t(F.hat)%*%x)/T          #Estimated Factor Loading
    U.hat = x-F.hat%*%t(B.hat)        #Estimated idiosyncratic component

    ## FARM
    lmY.F = lm(y~F.hat-1)           
    gamma.hat = coef(lmY.F)    #Estimate \gamma vector in the FARM model and PCR
    Y.tilde = resid(lmY.F)                             
    cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	 
    lambda.fit.U = cv.fit.U$lambda.min #tuning parameter selection
    fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U) ##FARM Estimation
    beta.hat.U = as.vector(fit.U$beta) #Obtain the fitted beta for FARM
    lmx.B = lm(x.pred~B.hat-1)
    f.pred = coef(lmx.B)
    u.pred = resid(lmx.B)
    
    ## HuberFARM
    X.hat = cbind(F.hat,U.hat)
    Huber.est = adaHuber.cv.lasso(X.hat,y) #Estimate FARM mode using Huber regression
    Huber.beta = Huber.est$coef
    Huber.const = Huber.beta[1]
    Huber.gamma = Huber.beta[2:(K.hat+1)]
    Huber.beta = Huber.beta[-(1:(K.hat+1))] #estimated beta

    #Huber 
    Pred.HUBER[i,j] = (y.pred-Huber.const-sum(f.pred*Huber.gamma)-sum(u.pred*Huber.beta))^2*y.sd^2 #prediction Difference with the true value
    Prvalue_HUBER[i,j]=(Huber.const + sum(f.pred*Huber.gamma)+sum(u.pred*Huber.beta))*y.sd+y.mu  # predicted value
    
    #PCR
    Pred.PCR[i,j] = (y.pred-sum(f.pred*gamma.hat))^2*y.sd^2     #PCR prediction Difference with the true value
    Prvalue_PCR[i,j]=sum(f.pred*gamma.hat)*y.sd+y.mu           #PCR predicted value
    
  }
  
  R2.HUBER[j] = 1-sum(Pred.HUBER[,j])/sum(Pred.MEAN[,j])	
  print(R2.HUBER[j])                    #Output R^2 value for Huber FARM
  R2.LASSO[j] = 1-sum(Pred.LASSO[,j])/sum(Pred.MEAN[,j])
  print(R2.LASSO[j])                    #Output R^2 value for LASSO
  R2.RIDGE[j] = 1-sum(Pred.RIDGE[,j])/sum(Pred.MEAN[,j])
  print(R2.RIDGE[j])                    #Output R^2 value for RIDGE
  R2.PCR[j] = 1-sum(Pred.PCR[,j])/sum(Pred.MEAN[,j])
  print(R2.PCR[j])                    #Output R^2 value for PCR
  
  }


##----------------------------------------------------------------------------##
##Draw picture of prediction under optimal window size
##----------------------------------------------------------------------------##
c = colcol[1]
Y = D[,c]
a = 3
plot(Y[(T+1):nrow(X)], type = "l",xlab = "",ylab = "", col = 1)
lines(Prvalue_HUBER[,c], lty = 2, col = "red")
lines(Prvalue_LASSO[,c], lty = a, col = "green")
lines(Prvalue_RIDGE[,c], lty = a, col = "purple")
lines(Prvalue_MEAN[,c],lty = a,col = "blue")
lines(Prvalue_PCR[,c], lty = a, col = "orange")




