# Distributions that may be useful for PRA
#
# Author: Paul Thomas Britton
#

# create a vector of N trianglular samples
# a = min
# b = mode
# c = max
triadist <- function(N,a,b,c) {
	if (!((a<=b)&&(b<=c))) {
		stop("triadist: Invalid Triangle Distribution Parameters")
	}
	U <- runif(N,0,1)
	x <- (b-a)/(c-a)
	Tri <- U <= x
	Tri[Tri==TRUE] <- a+sqrt(U[U<=x]*(c-a)*(b-a))
	Tri[Tri==FALSE] <- c-sqrt((1-U[U>x])*(c-a)*(c-b))
	return(Tri)
}

# create a vector of N lognormal samples
# mean = mean of lognormal
# median = 50th of lognormal
# EF = 95th/50th (of lognormal)
# mu = mean of underlying normal
# sigma = standard dev of underlying normal
#
# user must provide either: mean & EF, median & EF, or mu & sigma
logndist <- function(N,EF=exp(qnorm(0.95)),sigma=log(EF)/qnorm(0.95),
		mean=1,median=exp(log(mean)-(sigma^2)/2),
		mu=log(median)) {
	return(rlnorm(N,mu,sigma))
}

# k dimensional dirichlet probability distribution
# returns an array of k vectors of size N
# each of the k vectors corresponds to the margins of
# the dirichlet distribution
#
# N = number of samples
# ... = alpha vector (as an atomic vector or an argument list)

dirichlet <- function(N,...) {
	args <- list(...)
	len <- length(args)
	if (len == 1) {
		AlphaVector <- as.list(args[[1]])
		k <- length(AlphaVector)
	} else {
		AlphaVector <- args
		k <- len
	}
	if (k < 2) stop(paste("k=",k,"; dirichlet() requires k > 1"))
	Y <- list()
	V <- rep(0,N)
	n <- rep("",k)
	for (i in 1:k) {
		Y[[i]] <- rgamma(N,shape=AlphaVector[[i]],rate=1)
		V <- V + Y[[i]]
		n[i] <- paste("b",i,sep="")
	}
	dividebyV <- function(X) X/V
	diri <- lapply(Y,dividebyV)
	names(diri) <- n
	return(diri)
}
