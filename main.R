# A sample dataset to try the LogCM models ----------------
Data <- read.csv( "data_sample.csv" )
head(Data)
Train <- subset(Data, Data$Train==1)
Test  <- subset(Data, Data$Train!=0)

# Or try to generate a simulation study -------------------
# REQUIRE:
source("Generate.R")

# Now let's get started -----------------------------------
source("Algorithm.R")
# initialization with MEM
# REQUIRE: numbers of individuals, variables, canonicals
N <- 30; P <- 4; K <- 3
MEMcoef <- get_MEM(Train,P)
Initials <- Initiating(MEMcoef, K)

# the LogCM algorithm:
# REQUIRED: N, P, K, Qini, Cini
LogCM <- Iterating(Train, Initials$Qini, Initials$Cini)

# the LogSCM algorithm:
# REQUIRE: N, P, K, Qini, Cini
# REQUIRE: similarity matrix W and L,D
Sims <- get_LDW(MEMcoef) 
# or can use prior knowledge as well
# W <- read.csv("W.csv"); Sims <- get_LD(W)
L <- Sims$L; D <- Sims$D; W <- Sims$W
LogSCM <- Iterating(Train, Initials$Qini, Initials$Cini, lambda=1)


# LogPCM algorithms ---------------------------------------
source("Algorithm_PCM.R")
# REQUIRE: N, P, K, Qini, Cini, LDW
# REQUIRE: a large K0 to start
# single step
K0 <- 5; K <- K0
Initials_aux <- Initiating(get_MEM(Train,P), K)
source("auxs.R")
LogPCM <- Iteratingaux(Train, Initials_aux$Qini, Initials_aux$Cini, mu=2)

# increasing mu's to shrink K -----------------------------
# an adjusting increasing approach: based on whether K shrinked
Qlist_shrink <- vector("list")
mu <- 2; MuList <- c()
K0 <- 6; K <- K0; KList <- c() 
Initials_aux <- Initiating(get_MEM(Train,P), K0)
Q <- Initials_aux$Qini; C <- Initials_aux$Cini
for (iter in 1:500) {
  print(paste0("Now running: mu=", round(mu,3), "; at K=", K))
  MuList <- c(MuList, mu)
  source("auxs.R")
  LogPCM <- Iteratingaux(Train, Q, C, mu, thres.conv = 1)
  newK <- LogPCM$UniqueK; Q <- LogPCM$Q; C <- LogPCM$C
  if (iter==1) { Qlist_shrink[[1]] <- Q }
  if ( newK == K) { mu <- mu*1.2 
  } else {
      Qlist_shrink[[iter]] <- Q
      mu <- mu/(K*(K-1)) * (newK*(newK-1)) }
  K <- newK; KList  <- c(KList, K)
  if (K==1) {break} }

