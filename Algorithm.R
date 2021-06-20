library(lme4) # for mixed-effect model
library(CVXR) # for convex optimization
library(pROC) # for AUC values

# Basics --------------------------------------------------
posing <- function(x) {return(sum(x[x>0]))}
neging <- function(x) {return(sum(x[x<0]))}
get_Ni <- function(Data) {
  ni <- rep(0, N)
  for (i in 1:N) {
    ni[i] <- nrow(subset(Data, Data$individual==i)) }
  return(ni) }

get_LDW <- function(Coef) {
  N <- ncol(Coef)
  W <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in i:N) {
      W[i,j] <- exp(-sum( (Coef[,i]-Coef[,j])^2 ))
      W[j,i] <- W[i,j]  }  }
  D <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    D[i,i] <- sum(W[i,])  }
  L <- D-W
  Result <- list(L, D, W); names(Result) <- c("L","D","W")
  return(Result)  }
  

get_LD  <- function(Sim) {
  N <- nrow(Sim)
  Sim <- as.matrix(Sim)
  D <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    D[i,i] <- sum(Sim[i,])  }
  L <- D-Sim
  Result <- list(L, D, Sim); names(Result) <- c("L","D","W")
  return(Result) }


# Main Functions for iterative algorithm ------------------
update_Q <- function(Data, C, weight=TRUE) {
  weightedX<-matrix( 0, nrow=1, ncol=((P+1)*K+2) )
  ni <- get_Ni(Data)
  for (i in 1:N) {
    iData <- subset( Data, Data$individual==i )
    if (nrow(iData)==0) { next }
    iX <- cbind( 1, as.matrix( iData[,1:P] ) )
    ic <- C[,i]
    Xtil <- kronecker( t(ic), iX )
    wXnew <- cbind(Xtil, iData$Y, 1/ni[i])
    weightedX<-rbind(weightedX,wXnew) }
  weightedX<-data.frame(weightedX[-1,])
  colnames(weightedX)[((P+1)*K+1)]<-'Y'
  colnames(weightedX)[((P+1)*K+2)]<-'weight'
  
  if (weight==TRUE) { optweight <- weightedX$weight  
    } else { optweight <- rep( 1, nrow(weightedX) ) }
  q <- Variable( rows=(P+1)*K, cols=2 )
  optX <- as.matrix( weightedX[,1:((P+1)*K)] )
  opty <- weightedX$Y
  objective <- Minimize(sum( (log_sum_exp(optX%*%q) - opty*(optX%*%q[,2]))*optweight ) )
  cons <- q[,1]==0
  problem <- Problem(objective,constraints=list(cons))
  result <- solve(problem) #Solver:"ECOS"(default),"SCS"
  qall<-result$getValue(q)[,2] 
  Qnew <- matrix(qall, nrow=(P+1), ncol=K, byrow=FALSE)
  return(Qnew)  }

update_C <- function(Data, Q, C, lambda=0) {
  Cnew <- matrix( 0, nrow=K, ncol=N )
  ni <- get_Ni(Data)
  for (i in 1:N) {
    iData <- subset(Data, Data$individual==i)
    ic <- C[,i]
    if (nrow(iData)==0) {Cnew[,i]<-ic; next}
    iX <- as.matrix(iData[, 1:(P+1)])
    iy <- iData$Y
    for (k in 1:K) {
      xQc <- iX %*% Q %*% ic
      Qxk <- (iX %*% Q)[,k]
      exQc <- exp(xQc)
      exQcf<- exQc/(1+exQc)
      nume <- -neging(exQcf*Qxk) + posing(iy*Qxk) - neging(iy*xQc) + posing(exQcf*xQc)
      deno <-  posing(exQcf*Qxk) - neging(iy*Qxk) + posing(iy*xQc) - neging(exQcf*xQc)
      if (lambda != 0) {
        nume <- nume/(ni[i]) + 2*lambda* ((W%*%t(C))[i,k] + (D%*%t(C))[i,] %*% ic)
        deno <- deno/(ni[i]) + 2*lambda* ((D%*%t(C))[i,k] + (W%*%t(C))[i,] %*% ic) }
      Cnew[k,i]<-C[k,i]*(nume/deno)
      if (is.na(Cnew[k,i])) {Cnew[k,i]<-0.9}
      Cnew[k,i] <- min(Cnew[k,i],0.99)
      Cnew[k,i] <- max(Cnew[k,i],0.01) }
    Cnew[,i] <- Cnew[,i]/(sum(Cnew[,i])) }
  return(Cnew) }


# Evaluation Functions ------------------------------------
get_predict <- function(Data, Beta) {
  Yline <- vector(length=0)
  for (i in 1:N) {
    iData <- subset( Data, Data$individual==i )
    iX <- cbind( 1, as.matrix(iData[,1:P]) )
    yline <- iX %*% Beta[,i]
    Yline <- c(Yline, yline) }
  Yprob <- 1/(1+exp(-Yline))
  Ypred <- ifelse(Yline>=0, 1, 0)
  result <- data.frame(Yreal=Data$Y, Yline=Yline, Yprob=Yprob, Ypred=Ypred, individual=Data$individual)
  return(result) }

calc_objloss <- function(preds, C, lambda=0) {
  loss <- rep(0,N)
  for (i in 1:N) {
    ipreds <- subset( preds, preds$individual==i )
    loss[i] <- mean(log(1+exp(ipreds$Yline))-ipreds$Yreal*ipreds$Yline) }
  objloss <- sum(loss) 
  if (lambda!=0) { objloss <- objloss+lambda*sum(diag(C%*%L%*%t(C))) }
  return( objloss ) }

calc_loss <- function(preds) {
  return(sum(log(1+exp(preds$Yline))-preds$Yreal*preds$Yline)) }

calc_auc <- function(preds) {
  return(roc(preds$Yreal, preds$Yprob)$auc)
}

# Initialization ------------------------------------------
get_MEM <- function(Data, P) {
  mem <- "Y ~"
  for (i in 1:P) { mem <- paste(mem, paste0("X",i), "+") }
  for (i in 1:P) { mem <- paste(mem, paste0("(X",i,"|individual)"), "+") }
  mem <- paste(mem, "(1|individual)")
  mem_model <- glmer(eval(mem), data=Data, family=binomial(link="logit"))
  return( t(as.matrix(coef(mem_model)$individual)))  }


Initiating <- function( Coef, K ) {
  clusters <- kmeans(t(Coef), K)
  Qini <- t(clusters$centers)
  
  Cini <- matrix(0, nrow=K, ncol=N)
  for (i in 1:N) {
    chat <- Variable(K)
    m <- Coef[,i]
    objective<-Minimize(sum((m-Qini%*%chat)^2))
    cons1<- sum(chat)==1
    cons2<- chat>=0.01
    problem<-Problem(objective,constraints=list(cons1,cons2))
    result<-solve(problem)
    Cini[,i]<-result$getValue(chat)  }
  
  Result <- list(Qini, Cini); names(Result) <- c("Qini", "Cini")
  return(Result)  }


# Iteration -----------------------------------------------
Iterating <- function(Data, Qini, Cini, lambda=0,
                      weight=TRUE, IterMax=500, thres.mono=10, thres.conv=0.5) {
  Q <- Qini; C <- Cini
  ni <- get_Ni(Data)
  iter.loss <- calc_objloss( get_predict(Data, Q%*%C), C, lambda=lambda )
  
  tstart <- Sys.time()
  for (iter in 1:IterMax) {
    Cnew <- update_C(Data, Q, C, lambda=lambda)
    iter.loss.new <- calc_objloss( get_predict(Data, Q%*%Cnew), Cnew, lambda=lambda )
    #if ( iter.loss.new > iter.loss[length(iter.loss)]+thres.mono ) { break }
    iter.loss <- c( iter.loss, iter.loss.new )
    C <- Cnew
    
    Qnew <- update_Q(Data, C, weight=weight)
    iter.loss.new <- calc_objloss( get_predict(Data, Qnew%*%C), C, lambda=lambda )
    #if ( iter.loss.new > iter.loss[length(iter.loss)]+thres.mono ) { break }
    iter.loss <- c( iter.loss, iter.loss.new )
    Q <- Qnew
    
    if ( abs(diff(tail(iter.loss, 2))) < thres.conv ) { break }  }
  tend <- Sys.time()
  timediff <- difftime(tend, tstart, units=c("secs"))
  
  Result <- list(Q, C, Q%*%C, timediff, iter.loss)
  names(Result) <- c("Q", "C", "CoefMatrix", "TimeLapse", "IterLoss")
  return(Result) }
