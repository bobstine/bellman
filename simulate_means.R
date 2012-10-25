###  Code to simulate mean stochastic process

setwd("/Users/bob/C/bellman/sim_details/")

angle <- 210
  n   <- 200
alpha <- 0.05
omega <- 0.50
filename <- paste("bellman.a",angle,".n",n,".o",round(100*omega),".al",round(100*alpha),sep="")

####################################################################
#
#  read simulation details from C++ files
#
####################################################################

WealthArray <- readLines(paste(filename,"wealth",sep="."), n=4)
wealth.desc <- WealthArray[1]
wealth.array<- as.numeric(unlist(strsplit(WealthArray[2]," ")))
wealth.prob <- as.numeric(unlist(strsplit(WealthArray[4]," ")))
wealth.indx <- 1 + as.numeric(unlist(strsplit(WealthArray[3]," ")))

# --- read arrays with optimal process means
mean   <- read.table (paste(filename,"mean",sep="."))
prob   <- read.table (paste(filename,"prob",sep="."))
indx   <- 1 + read.table (paste(filename,"indx",sep="."))
dim(mean)


####################################################################
#
#  run simulation
#
####################################################################

# --- note that R is 1-based indexing
nRounds <- nrow(mean)
iZero <- ncol(mean)-nRounds + 1


# --- simulation fills these vectors
meanProcess <- rep(0, nRounds)
indxProcess <- rep(0, nRounds)


# --- initial conditions
             k   <- iZero 
  meanProcess[1] <-   mean[1,k]   # randomly chosen mean
  indxProcess[1] <-        k      # position
      p.reject   <-   prob[1,k]

# --- need to random per round if reject
unif <- runif(2*nRounds)

# --- run simulation
k <- iZero 
for(round in 2:nRounds) {
	cat("@ k=", k, " mean=", meanProcess[round-1]," p.reject=",p.reject,"\n");
	if (unif[2*round-1] < p.reject) { #reject
		k <- indx[round-1,k] + (wealth.prob[k] < unif[2*round]) }
	else k <- k + 1; # did not reject
	indxProcess[round] <- k
	meanProcess[round] <- mean[round,k]
	p.reject           <- prob[round,k]
	}    
	
	
plot(meanProcess, main=paste("Angle",angle,"  Omega",omega, "  Oracle alpha",alpha,sep=" "),
	xlab="Test Round", ylab="Mean(i)", col=indxProcess)

plot(indxProcess, main=paste("Angle",angle,"  Omega",omega, "  Oracle alpha",alpha,sep=" "),
	xlab="Test Round", ylab="Wealth Index(i)", col=indxProcess)

plot(wealth.array)

#------------------------------------------------------------------
#
#  read risks computed in C++
#
#------------------------------------------------------------------


setwd("/Users/bob/C/bellman/")

n<-200; mu<-3.2; scale<-2;
risk <- read.table(paste("calc_risk_",n,"_",mu,"_",scale,sep="")); dim(risk)
plot(-risk[,1],-risk[,2], xlab="Oracle Risk", ylab="Universal Risk", 
	main="Bellman Risk Calculation",
	sub=paste("n=",n," mu=",mu,", varying p_zero", sep=""))
