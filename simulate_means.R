###  Code to simulate mean stochastic process

setwd("/Users/bob/C/bellman/")

filename <- "bellman.a210.n200"
iZero <- 6

mean <- read.table (paste(filename,"mean",sep=".")); dim(mean)
prob <- read.table (paste(filename,"prob",sep="."))
indx <- 1 + read.table (paste(filename,"indx",sep="."))
nRounds <- nrow(mean)


meanProcess <- rep(0, nRounds)

             k <- iZero 
meanProcess[1] <- mean[1,k]
             p <- prob[1,k]

unif <- runif(nRounds)
for(round in 2:(nRounds-1)) {
	cat("k=", k, " mean = ", meanProcess[round-1]," p=",p,"\n");
	if (unif[round] < p) { #reject
		k <- indx[round-1,k]
		meanProcess[round] <- mean[round,k]
		                 p <- prob[round,k]  }
	else {
		k <- k+1;
		meanProcess[round] <- mean[round,k]
		                 p <- prob[round,k]  }
	}    
		             
row<-1; nCols<-7

prob[row,1:nCols]
mean[row,1:nCols]
indx[row,1:nCols]

