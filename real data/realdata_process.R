##-------------------Read and Separate Data-----------------------------------##
Z0 = read.csv('2022-12.csv')
nnn = colnames(Z0)
nnn = nnn[-1]
tcode = as.numeric(Z0[1,-1]) # transform code: length 227
date = Z0[-1,1]
Z0 = Z0[-1,]
X0 = Z0[,-1]

# initialize processed matrix
D = matrix(0, 222, 226)
##----------------------------------------------------------------------------##

##-------------------Do Data Transformation-----------------------------------##
for(j in 1:226){
  tmp.code = tcode[j]
  if(tmp.code==1){
    D[,j] = X0[-c(1,2),j]
  }
  
  if(tmp.code==2){
    tmp.v = X0[-1,j]-head(X0[,j],-1)
    D[,j] = tmp.v[-1]
  }
  
  if(tmp.code==3){
    tmp.v = X0[-1,j]-head(X0[,j],-1)
    tmp.v2 = tmp.v[-1]-head(tmp.v,-1)
    D[,j] = tmp.v2
  }
  
  if(tmp.code==4){
    D[,j] = log(X0[-c(1,2),j])
  }
  
  if(tmp.code==5){
    tmp.v = log(X0[-1,j])-log(head(X0[,j],-1))
    D[,j] = tmp.v[-1]
  }
  
  if(tmp.code==6){
    tmp.v = log(X0[-1,j])-log(head(X0[,j],-1))
    tmp.v2 = tmp.v[-1]-head(tmp.v,-1)
    D[,j] = tmp.v2
  }
  
  if(tmp.code==7){
    tmp.v = X0[-1,j]/head(X0[,j],-1)
    D[,j] = tmp.v[-1]-head(tmp.v,-1)
  }
}

# remove the output data, and we will obtain realdata_covariate.csv 
# and realdata_response.csv
date = date[-c(1,2)]
#y = matrix(0, 222,1)
#y = D[, 82] # HOUSTMW # y = D[, 83] # HOUSTNE
#nnn_y = nnn[82]
#D = D[, -82]
#nnn = nnn[-82] # name of every line
##----------------------------------------------------------------------------##

##-------------------Output and Save Data-------------------------------------##
df1 = data.frame(cbind(date, D))
colnames(df1) = c("date", nnn)
write.csv(df1, "realdata_processed.csv", row.names = FALSE)
##----------------------------------------------------------------------------##
