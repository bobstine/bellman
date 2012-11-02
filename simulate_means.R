###  Code to simulate mean stochastic process

set.par <- function (mfrow=c(1,1)) {
	par(mfrow=mfrow, mgp=c(2,1,0),
		fg="white", bg = "slategray4",
		mar=c(3,3,2,1), # bottom left top right
		col="white",col.main="white",col.sub="white",col.lab="white",col.axis="white")
	}
		
reset <- function() { 
	par(mfrow=c(1,1),mar=c(5,4,4,2)+.1,mgp=c(3,1,0), 
	fg="black",bg="white",col="black",	
	col.main="black", col.sub="black", col.lab="black", col.axis="black"); 
	palette("default");  
	}


####################################################################
#
#  Bid as function of wealth
#
####################################################################

sum.of.log.recip <- 3.387735532

bid.int <- function(k,scale) { k<-k+1; scale/(k*(log(k+1))^2) }

bid <- function(w,scale) {
	wealth <- scale * sum.of.log.recip;
	k <- 0; delta <- 0;
	while(wealth > w) { delta<-bid.int(k,scale); wealth <- wealth - delta; k <- k+1; }
	wt <- (w-wealth)/delta;
	wt*delta + (1-wt)*bid.int(k,scale)
	}

# check: graph of bid(wealth) with integer version imposed
n<-10;
x <- 0:n
scale <- 2
wealth <- scale*sum.of.log.recip - c(0,cumsum(sapply(0:(n-1), function(x) bid.int(x,scale))))
bids   <- sapply(x,function(x) bid.int(x,scale))
#    reference points at integers
plot(wealth,bids)
#    add interpolator
w <- seq(min(wealth),max(wealth),length.out=60)
lines(w, sapply(w, function(x) bid(x,scale)))
#    add into plot from below of the tables used in C++
lines(wealth.array, wealth.bids, col="gray")


####################################################################
#
#  Simulation process from C++ files
#
####################################################################

setwd("/Users/bob/C/bellman/sim_details/")

angle <- 334
  n   <- 500
alpha <- 0.05
omega <- 0.50
scale <- 2
filename <- paste("dual_bellman.a",angle,".n",n,".s",round(10*scale),
						".o",round(100*omega),".al",round(100*alpha),sep=""); filename

# --- read simulation details from C++ files

WealthLines <- readLines(paste(filename,"wealth",sep="."), n=7)
wealth.desc <- WealthLines[1]
wealth.array<- as.numeric(unlist(strsplit(WealthLines[2]," ")))
wealth.bids <- as.numeric(unlist(strsplit(WealthLines[3]," ")))
wealth.rIndx<- 1+as.numeric(unlist(strsplit(WealthLines[4]," ")))
wealth.rWts <- as.numeric(unlist(strsplit(WealthLines[5]," ")))
wealth.bIndx<- 1+as.numeric(unlist(strsplit(WealthLines[6]," ")))
wealth.bWts <- as.numeric(unlist(strsplit(WealthLines[7]," ")))
points(wealth.array, wealth.bids)

# --- read arrays with optimal process means
mean   <- read.table (paste(filename,"mean",sep="."))
prob   <- read.table (paste(filename,"prob",sep="."))
dim(mean)


####################################################################
#
#  run simulation
#
####################################################################

# --- note that R is 1-based indexing; need iZero from calculation as omega column
nRounds <- nrow(mean)
iZero <- 80
wealth.array[iZero] # should be omega (or closest that is less)

doit <- function(seed=sample.int(100000,1)) {
	set.seed(seed)
	meanProcess <-indxProcess <- reject.prob <- rep(0, nRounds); 
	riskProcess <-oracleRisk <- rep(0, nRounds); 
	# --- need to random per round if reject
	pval <- runif(nRounds); randomize <- runif(nRounds)
	# --- run simulation
	k <- iZero
	meanProcess[1] <- mean[1,k]   # mean at state 0
	indxProcess[1] <-        k      # wealth state (index)
	reject.prob[1] <- p.reject <- prob[1,k]
	for(round in 2:nRounds) {
		cat("round ", round, "@ k=", k, " mean=", meanProcess[round-1]," p.reject=",p.reject,"\n");
		# oracle risk
		if (pval[round] < 0.05) { oracleRisk[round] <- 1; }
		else                    { oracleRisk[round] <- meanProcess[round-1]^2; }
		# bidder risk
		w <- randomize[round];
		if (pval[round] < p.reject) { k <- wealth.rIndx[k] + (wealth.rWts[k] < w); riskProcess[round] <- 1; }
		else {                        k <- wealth.bIndx[k] + (wealth.bWts[k] < w); riskProcess[round] <- meanProcess[round-1]^2; }
		indxProcess[round] <- k
		meanProcess[round] <- mean[round,k]
		p.reject           <- reject.prob[round] <- prob[round,k]
	}   
	list(bidder.risk=riskProcess, bidder.prob = reject.prob, bidder.state=indxProcess,  
			seed=seed, pval=pval, mean=meanProcess, oracle.risk=oracleRisk,
			bidder.rejects = (pval<reject.prob), oracle.rejects = (pval<0.05))
}
	
simres <- doit()

set.par(mfrow=c(2,2))
	plot(simres$mean, 
		main=paste("Bidder: angle=",angle," al=",alpha," omega=",omega, " scl=",scale,sep=""),
		sub=paste("Dual universal bidder vs Bayesian oracle, seed=",simres$seed,sep=" "),
		xlab="Round", ylab="Mean(i)", pch=".", cex=3, col=c("black","red")[1+simres$bidder.rejects])
	text(nRounds-50,1,paste("risk=",round(sum(simres$bidder.risk),digits=2)))
	lines(max(1,max(simres$mean))*cumsum(simres$bidder.risk)/(max(1,sum(simres$bidder.risk))), col="gray")

	plot(simres$mean, 
		main=paste("Oracle"), col=c("black","red")[1+simres$oracle.rejects],
		xlab="Round", ylab="Mean(i)", pch=".", cex=3,col=c("black","red")[1+simres$oracle.rejects]))
	text(nRounds-50,1,paste("risk=",round(sum(simres$oracle.risk),digits=2)))
	lines(max(1,max(simres$mean))*cumsum(simres$oracle.risk)/(max(1,sum(simres$oracle.risk))), col="gray")

	plot(wealth.array[simres$bidder.state], xlab="Round", ylab="Wealth [i]", pch=20, col=simres$bidder.state)
	
	plot(cumsum(simres$oracle.risk), cumsum(simres$bidder.risk), 
		pch=16,col=c("gray","pink","red")[1+simres$bidder.rejects+simres$oracle.rejects], 
		xlab="Cumulative Oracle Risk", ylab="Bidder Risk")
reset()

# 54057 … one
# 22077 … up for good
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
