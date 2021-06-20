# Define auxiliary matrices ---------------------
pair <- 0
for (i in 1:(K-1)) {
  for (j in (i+1):K) {
    pair <- pair+1
    A <- matrix(0, nrow=(P+1), ncol=(P+1)*K )
    A[, ((i-1)*(P+1)+1):(i*(P+1))] <- diag(P+1)
    A[, ((j-1)*(P+1)+1):(j*(P+1))] <- -diag(P+1)
    assign( paste0("A", pair) ,A ) } }
if (K==2) {A1 <- cbind(diag(P+1), -diag(P+1))}   


CVXtext <- ""
for ( i in 1:(K*(K-1)/2) ) {
  CVXtext <- paste0(CVXtext, " cvxr_norm(A", i, " %*% q[,2]) +") }
CVXtext <- substr(CVXtext, 1, nchar(CVXtext)-1)


LOSStext <- ""
for ( i in 1:(K*(K-1)/2) ) {
  LOSStext <- paste0(LOSStext, " sqrt(sum( (A", i, " %*% q)^2 )) +") }
LOSStext <- substr(LOSStext, 1, nchar(LOSStext)-1)