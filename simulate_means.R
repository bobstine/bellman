###  Code to simulation mean stochastic process

setwd("/Users/bob/C/bellman/")

filename <- "sim_risk_200_0.5_3.2_2"


# --- read wealth array data
Wealth <- readLines(paste(filename,"wealth",sep="."), n=4)
wealth.desc <- Wealth[1]
wealth.array<- as.numeric(unlist(strsplit(Wealth[2]," ")))
wealth.prob <- as.numeric(unlist(strsplit(Wealth[4]," ")))
wealth.indx <- 1 + as.numeric(unlist(strsplit(Wealth[3]," ")))

# --- read arrays with optimal process means
mean <- read.table (paste(filename,"mean",sep="."))
prob <- read.table (paste(filename,"prob",sep="."))
indx <- 1 + read.table (paste(filename,"indx",sep="."))
dim(mean)


# --- note that R is 1-based indexing
nRounds <- nrow(mean)
iZero <- ncol(mean)-nRounds + 1



# --- simulation fills this vector
meanProcess <- rep(0, nRounds)


# --- initial conditions
             k <- iZero 
meanProcess[1] <- mean[1,k]
      p.reject <- prob[1,k]

# --- need to random per round if reject
unif <- runif(2*nRounds)

# --- run simulation
k <- iZero 
for(round in 2:nRounds) {
	cat("@ k=", k, " mean=", meanProcess[round-1]," p.reject=",p.reject,"\n");
	if (unif[2*round-1] < p.reject) { #reject
		k <- indx[round-1,k] + (wealth.prob[k] < unif[2*round]) }
	else k <- k + 1; # did not reject
	meanProcess[round] <- mean[round,k]
	p.reject <- prob[round,k]
	}    
	
	
plot(meanProcess)


#------------------------------------------------------------------
#
#  read risks computed in C++
#
#------------------------------------------------------------------

calc_risk_200_3.2_2

setwd("/Users/bob/C/bellman/")

n<-200; mu<-3.2; scale<-2;
risk <- read.table(paste("calc_risk_",n,"_",mu,"_",scale,sep="")); dim(risk)
plot(-risk[,1],-risk[,2], xlab="Oracle Risk", ylab="Universal Risk", 
	main="Bellman Risk Calculation",
	sub=paste("n=",n," mu=",mu,", varying p_zero", sep=""))
