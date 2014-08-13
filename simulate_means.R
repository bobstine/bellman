###  Code to simulate mean stochastic process

set.par <- function (mfrow=c(1,1)) {
	fg.col <- "black";
	par(mfrow=mfrow, mgp=c(2,1,0),
		# fg="white", bg = "slategray4",
		fg=fg.col, bg = "white",
		mar=c(3.1,3.1,1,1), # bottom left top right, add space when show several default:   c(5, 4, 4, 2) + 0.1.
		col=fg.col,col.main=fg.col,col.sub=fg.col,col.lab=fg.col,col.axis=fg.col)
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

# directly from paper
alpha.univ <- function(w) { w - log(2)/log(1 + 2^(1/w)) }

# check first bids (check those in paper)
w0 <- 1
w1 <- w0 - alpha.univ(w0); w1
w2 <- w1 - alpha.univ(w1); w2
w3 <- w2 - alpha.univ(w2); w3


########################################################################################
##
##		This code works for the *matrix* utility results; see below for vector
##
########################################################################################

setwd("/Users/bob/C/bellman/sim_details/")

nRounds <- 200

angle   <- 296.565  #165
x.prob  <-   0
x.omega <-   0.25
y.prob  <-   0.001
y.omega <-   0.25

(prefix <- paste("n_",nRounds,"_angle_",angle,"_oracle_",x.prob,"_",x.omega,"_bidder_",y.prob,"_",y.omega,sep=""))

# --- parse wealth details

RowWealth <- readLines(paste(prefix,"row_wealth",sep="."))
(rowWealth.desc  <- RowWealth[1])
rowWealth.array <-   as.numeric(unlist(strsplit(substring(RowWealth[2],first=11)," ")))
rowWealth.bids  <-   as.numeric(unlist(strsplit(substring(RowWealth[3],first= 9)," ")))
rowWealth.rIndx <- 1+as.numeric(unlist(strsplit(substring(RowWealth[4],first=13)," ")))
rowWealth.rWts  <-   as.numeric(unlist(strsplit(substring(RowWealth[5],first=13)," ")))
rowWealth.bIndx <- 1+as.numeric(unlist(strsplit(substring(RowWealth[6],first=13)," ")))
rowWealth.bWts <-    as.numeric(unlist(strsplit(substring(RowWealth[7],first=13)," ")))
plot(rowWealth.array, rowWealth.bids)

ColWealth <- readLines(paste(prefix,"col_wealth",sep="."))
(colWealth.desc  <- ColWealth[1])
colWealth.array <-   as.numeric(unlist(strsplit(substring(ColWealth[2],first=11)," ")))
colWealth.bids  <-   as.numeric(unlist(strsplit(substring(ColWealth[3],first= 9)," ")))
colWealth.rIndx <- 1+as.numeric(unlist(strsplit(substring(ColWealth[4],first=13)," ")))
colWealth.rWts  <-   as.numeric(unlist(strsplit(substring(ColWealth[5],first=13)," ")))
colWealth.bIndx <- 1+as.numeric(unlist(strsplit(substring(ColWealth[6],first=13)," ")))
colWealth.bWts <-    as.numeric(unlist(strsplit(substring(ColWealth[7],first=13)," ")))
plot(colWealth.array, colWealth.bids)


# --- check initial wealth positions (add 1 for 0-based index to 1-based)

iZeroRow <- 60 + 1        # value from first line of RowWealth description
rowWealth.array[iZeroRow] # should be omega (or closest that is less)

iZeroCol <- 60 + 1
colWealth.array[iZeroCol]


# --- simulate

doit <- function(seed=sample.int(100000,1)) {
	set.seed(seed)
	meanProcess 	<- rep(0, nRounds); 
	positions <- risks <- rejectProb		<- matrix(NA, nrow=nRounds, ncol=2)		# row and column
	# --- run simulation
	kR <- iZeroRow; kC <- iZeroCol
	for(round in 1:nRounds) {
		positions[round,] <- c(kR,kC)		# current wealth states
		# --
		mu      <- read.table(paste(prefix,"_",round-1,".mean", sep=""))
		meanProcess[round] <- mu[kR,kC]		# mean at current state
		# --
		rowProb <- read.table(paste(prefix,"_",round-1,".rowProb", sep=""))
		rejectProb[round,1] <- rowProb[kR,kC]
		colProb <- read.table(paste(prefix,"_",round-1,".colProb", sep=""))
		rejectProb[round,2] <- colProb[kR,kC]
		# --
		cat("round ", round, "@ k=", kR,",", kC, 
			" mean=", meanProcess[round]," p(reject)=",rejectProb[round,],"\n");
		# -- random chance
		p <- runif(1)
		w <- runif(1) # linear weights for wealth
		if (p < rejectProb[round,1]) { 
			kR <-rowWealth.rIndx[kR]+(rowWealth.rWts[kR]<w); 
			risks[round,1]<-1; } 
		else {                    
			kR <-rowWealth.bIndx[kR]+(rowWealth.bWts[kR]<w); 
			risks[round,1]<-meanProcess[round]^2; }
		if (p < rejectProb[round,2]) { 
			kC <-colWealth.rIndx[kC]+(colWealth.rWts[kC]<w); 
			risks[round,2]<-1; } 
		else {                    
			kC <-colWealth.bIndx[kC]+(colWealth.bWts[kC]<w); 
			risks[round,2]<-meanProcess[round]^2; 
		}
	}   
	list(risks=risks, positions=positions, means=meanProcess, rejectProbs=rejectProb)
}

simres <- doit(); colSums(simres$risks); 

	plot(simres$means, type="l")

	plot(simres$risks[,1], col="magenta", ylim=range(simres$risks), type="l"); points(simres$risks[,2])
	
	plot(cumsum(simres$risks[,1]),type="l", col="magenta"); lines(cumsum(simres$risks[,2]))
	
	plot(rowWealth.array[simres$positions[,1]], type="l", ylim=c(0,1.5)); lines(simres$means, col="red")
	
	plot (rowWealth.array[simres$positions[,1]], type="l", col="magenta", ylim=c(0,2), ylab="Wealths")
	lines(colWealth.array[simres$positions[,2]], type="l", col="black")
	
	plot (rowWealth.bids[simres$positions[,1]], type="l", col="magenta", log="y")
	lines(colWealth.bids[simres$positions[,2]], type="l", col="black")
	
reset()



########################################################################################
##
##		This code works for the *vector* utility results
##
########################################################################################

####################################################################
#
#  Simulation process from C++ files
#
####################################################################

setwd("/Users/bob/C/bellman/sim_details/")

angle <-  165
  n   <- 1000
alpha <-   1
omega <-   0.5
w0    <-   0.5
filename <- paste("dual_bellman.a",angle,".n",n,".w",round(100*w0),
						".o",round(100*omega),".al",round(100*alpha),sep=""); filename

# --- read simulation details from C++ files

WealthLines <- readLines(paste(filename,"wealth",sep="."), n=7)
wealth.desc <- WealthLines[1]
wealth.array<- as.numeric(unlist(strsplit(substring(WealthLines[2],first=11)," ")))
wealth.bids <- as.numeric(unlist(strsplit(substring(WealthLines[3],first= 9)," ")))
wealth.rIndx<- 1+as.numeric(unlist(strsplit(substring(WealthLines[4],first=13)," ")))
wealth.rWts <-   as.numeric(unlist(strsplit(substring(WealthLines[5],first=13)," ")))
wealth.bIndx<- 1+as.numeric(unlist(strsplit(substring(WealthLines[6],first=13)," ")))
wealth.bWts <-   as.numeric(unlist(strsplit(substring(WealthLines[7],first=13)," ")))
plot(wealth.array, wealth.bids)

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
iZero <- 36
wealth.array[iZero] # should be omega (or closest that is less)

doit <- function(seed=sample.int(100000,1)) {
	set.seed(seed)
	meanProcess <-indxProcess <- reject.prob <- rep(0, nRounds); 
	riskProcess <-oracleRisk <- rep(0, nRounds);
	# --- run simulation
	k <- iZero
	for(round in 1:nRounds) {
		meanProcess[round] <- mean[round,k]   # mean at initial state
		indxProcess[round] <-        k        # wealth state (index)
		p.reject <- reject.prob[round] <- prob[round,k]
		pvals[round] <- 2*(1-pnorm(abs(rnorm(1,mean=meanProcess[round]))))
		cat("round ", round, "@ k=", k, " mean=", meanProcess[round]," p.reject=",p.reject,"\n");
		# oracle risk
		oracleRisk[round] <- min(1,meanProcess[round]^2)
		# bidder risk
		w <- runif(1);
		if (runif(1) < p.reject){ k<-wealth.rIndx[k]+(wealth.rWts[k]<w); riskProcess[round]<-1; }
		else {                    k<-wealth.bIndx[k]+(wealth.bWts[k]<w); riskProcess[round]<-meanProcess[round]^2;}
	}   
	list(bidder.risk=riskProcess, bidder.prob = reject.prob, bidder.state=indxProcess,  
			seed=seed, pvals=pvals, mean=meanProcess, oracle.risk=oracleRisk,
			bidder.rejects = (pvals<reject.prob), oracle.rejects = (pvals<0.05))
}

simres <- doit()

# set.par(mfrow=c(1,2))
	set.par(mfrow=c(1,1))
	plot(simres$mean, 
		# main=paste("Bidder: angle=",angle," al=",alpha," omega=",omega, sep=""),
		# sub=paste("Universal bidder vs Bayesian oracle, seed=",simres$seed,sep=" "),
		xlab="Test", ylab="Mean", pch=".", cex=3) # col=c("gray","red")[1+simres$bidder.rejects])
	text(nRounds-100,1.7,paste("risk=",round(sum(simres$bidder.risk),digits=0)))
	lines(max(1,max(simres$mean))*cumsum(simres$bidder.risk)/(max(1,sum(simres$bidder.risk))), col="gray")

	plot(simres$mean, 
		main=paste("Oracle"), col=c("gray","red")[1+simres$oracle.rejects],
		xlab="Round", ylab="Mean(i)", pch=".", cex=3)
	text(nRounds-150,1,paste("risk=",round(sum(simres$oracle.risk),digits=2)))
	lines(max(1,max(simres$mean))*cumsum(simres$oracle.risk)/(max(1,sum(simres$oracle.risk))), col="gray")

	plot(wealth.array[simres$bidder.state], xlab="Round", ylab="Wealth [i]", pch=20, col=simres$bidder.state)
	
	plot(cumsum(simres$oracle.risk),cumsum(simres$bidder.risk), 
		pch=16,col=c("gray","pink","red")[1+simres$bidder.rejects+simres$oracle.rejects], 
		xlab="Cumulative Oracle Risk", ylab="Bidder Risk")
	abline(a=0,b=atan(2 * pi * (360-angle)/360))
reset()

r <- cbind(simres$mean, simres$pval, simres$bidder.prob, simres$oracle.rejects, simres$bidder.rejects)
colnames(r)<-c("mean","pval","reject.prob", "oracle","bidder")
plot(r[,1],r[,3], xlab="Mean", ylab="reject prob")

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
