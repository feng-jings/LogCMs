library(CVXR)

# Define auxiliary matrices ---------------------
get_CVXtext <- function(K) {
  CVXtext <- ""
  for ( i in 1:(K*(K-1)/2) ) {
    CVXtext <- paste0(CVXtext, " cvxr_norm(A", i, " %*% q[,2]) +") }
  CVXtext <- substr(CVXtext, 1, nchar(CVXtext)-1)
  return(CVXtext) }
  
get_LOSStext <- function(K) {
  LOSStext <- ""
  for ( i in 1:(K*(K-1)/2) ) {
    LOSStext <- paste0(LOSStext, " sqrt(sum( (A", i, " %*% q)^2 )) +") }
  LOSStext <- substr(LOSStext, 1, nchar(LOSStext)-1)
  return(LOSStext) }


# Main Functions modified for iterative algorithm ---------
update_Qaux <- function(Data, C, mu, weight=TRUE) {
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
  objective <- Minimize(sum( (log_sum_exp(optX%*%q) - opty*(optX%*%q[,2]))*optweight ) +
                          mu*( eval(parse(text=CVXtext)) ))
  cons <- q[,1]==0
  problem <- Problem(objective,constraints=list(cons))
  result <- solve(problem) #Solver:"ECOS"(default),"SCS"
  qall<-result$getValue(q)[,2] 
  Qnew <- matrix(qall, nrow=(P+1), ncol=K, byrow=FALSE)
  return(Qnew)  }


# Evaluation Functions modified ---------------------------
calc_objlossaux <- function(preds, Q, C, mu, lambda=0) {
  loss <- rep(0,N)
  for (i in 1:N) {
    ipreds <- subset( preds, preds$individual==i )
    loss[i] <- mean(log(1+exp(ipreds$Yline))-ipreds$Yreal*preds$Yline) }
  objloss <- sum(loss) 
  if (lambda!=0) { objloss <- objloss+lambda*sum(diag(C%*%L%*%t(C))) }
  q <- c(Q)
  objlossaux <- objloss + mu* (eval(parse(text=LOSStext)))
  return( objlossaux ) }


# Iterating modified --------------------------------------
Iteratingaux <- function(Data, Qini, Cini, mu, lambda=0,
                      weight=TRUE, IterMax=500, thres.mono=10, thres.conv=0.5) {
  Q <- Qini; C <- Cini
  ni <- get_Ni(Data)
  iter.loss <- calc_objlossaux( get_predict(Data, Q%*%C), Q, C, mu, lambda=lambda )
  
  tstart <- Sys.time()
  for (iter in 1:IterMax) {
    Cnew <- update_C(Data, Q, C, lambda=lambda)
    iter.loss.new <- calc_objlossaux( get_predict(Data, Q%*%Cnew), Q, Cnew, mu, lambda=lambda )
    if ( iter.loss.new > iter.loss[length(iter.loss)]+thres.mono ) { break }
    iter.loss <- c( iter.loss, iter.loss.new )
    C <- Cnew
    
    Qnew <- update_Q(Data, C, weight=weight)
    iter.loss.new <- calc_objlossaux( get_predict(Data, Qnew%*%C), Qnew, C, mu, lambda=lambda )
    if ( iter.loss.new > iter.loss[length(iter.loss)]+thres.mono ) { break }
    iter.loss <- c( iter.loss, iter.loss.new )
    Q <- Qnew
    
    if ( abs(diff(tail(iter.loss, 2))) < thres.conv ) { break }  }
  tend <- Sys.time()
  timediff <- difftime(tend, tstart, units=c("secs"))
  
  C <- C[!duplicated(t(round(Q,3))),]; for (n in 1:N) {C[,n]<-C[,n]/sum(C[,n])}
  Q <- Q[,!duplicated(t(round(Q,3)))]
  uniqK <- max(1, ncol(Q))
  
  Result <- list(Q, C, Q%*%C, timediff, iter.loss, uniqK)
  names(Result) <- c("Q", "C", "CoefMatrix", "TimeLapse", "IterLoss", "UniqueK")
  return(Result) }
