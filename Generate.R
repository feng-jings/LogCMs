library(MCMCpack) # for dirichlet
library(MASS)     # for multivariate normal


# Hyperparameters for simulation data generation ----------
G<-3    # #of individual groups, same as number of canonicals
N<-30   # #of individuals
K<-3    # #of canonical models
P<-4    # length of beta
v<-20   # set dominance: 8 is the least integer to make sure the right clustering of Qall
Ttrain <- 30 #number of training data
Ttest <- 10 #number of testing data

# Membership matrix C: K*N (each column is for one individual) ------
xdoms <- sample(K,N,replace=TRUE)
C <- matrix(0, K, N)
for (i in 1:N) {
  temp_dom <- rep(1,K)
  temp_dom[xdoms[i]] <- v
  C[,i] <- rdirichlet(1,temp_dom) }

# K Canonicals matrix Q: (P+1)*K (each column is one canonical) -----
Q<-matrix(0, nrow=(P+1), ncol=K)
Q[,1]<-c( 0,-1, 2, 2, 1)  #manually
Q[,2]<-c(-5,-1,-1,-1, 0)  #manually
Q[,3]<-c( 5,-2, 0,-4,-2)  #manually

# Real Beta Matrix: (P+1)*K for all individuals ---------------------
Coef_real <- Q%*%C
colnames(Coef_real) <- 1:N

# Similarity matrices L,D,W: N*N ------------------------------------
W <- matrix(0,N,N); rownames(W)<-1:N; colnames(W)<-1:N  #W:N*N
for ( i in 1:N ) {
  for ( j in i:N ) {
    W[i,j] <- exp( -sum( (Coef_real[,i]-Coef_real[,j])^2 ) )
    W[j,i] <- W[i,j]  } }
D<-matrix(0,N,N)
for (i in 1:N) {D[i,i]<-sum(W[i,])}                      #D:N*N
L=D-W                                                    #L:N*N

# Different X distributions for different groups (Maybe not necassary) -----
meanX<-matrix(c(-1, 3,-3,-1,
                -1, 2,-2,-2,
                 1, 0, 1,-3), nrow=3, byrow=TRUE)
varX<-matrix( c(1,1,1,1,1,1,1,1,2,2,2,2), nrow=3, byrow=TRUE)
verror<-0.8
#Data array: 30+10 observations for each individual
DataPool <- rep(0, P+3)
for (i in 1:N) {
  dom <- xdoms[i]
  X <- mvrnorm((Ttrain+Ttest),meanX[dom,],diag(varX[dom,]))
  z <- cbind(1,X) %*% Coef_real[,i]
  Yreal<-ifelse(z>=0,1,0)
  error<-rnorm(1,0,verror)
  Y<-ifelse(z+error>=0,1,0)
  DataPool <- rbind(DataPool, cbind(X,Y,i,Yreal)) }
DataPool <- as.data.frame(DataPool[2:nrow(DataPool),])
colnames(DataPool) <- c('X1','X2','X3','X4','Y','individual','Yreal')
DataPool$Train <- rep(c(rep(1,Ttrain),rep(0,Ttest)),N)

# sampling from the data pool
dsample<-round(runif(N,0.7*Ttrain,Ttrain))     #dense random sampling
ssample<-round(runif(N,0.2*Ttrain,0.4*Ttrain)) #sparse random sampling

index<-0
for (n in 1:N) {
  indexn<-sort((Ttrain+Ttest)*(n-1)+sample(1:Ttrain,dsample[n]))
  indext<-c(((Ttrain+Ttest)*n-Ttest+1):((Ttrain+Ttest)*n))
  index<-c(index,indexn,indext) }
index<-index[-1]

Dense<-DataPool[index,]


write.csv(Dense,file="data_sample.csv",row.names=FALSE)
write.csv(W,file="W.csv",row.names=FALSE)

