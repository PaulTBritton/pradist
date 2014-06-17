\name{logndist}
\alias{logndist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create a vector of N lognormal samples
}
\description{
	Create a vector of N lognormal samples from a lognormal distribution with either (mean,EF), (median,EF) or (mu,sigma) as parameters
}
\usage{
logndist(N, EF = exp(qnorm(0.95)), sigma = log(EF)/qnorm(0.95), mean = 1, median = exp(log(mean) - (sigma^2)/2), mu = log(median))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
	Number of samples
}
  \item{EF}{
	EF = error factor = 95th percentile / median
}
  \item{sigma}{
	standard deviation of the underlying normal distribution
}
  \item{mean}{
	mean of the lognormal distribution
}
  \item{median}{
	50th percentile of the lognormal distribution
}
  \item{mu}{
	mean of the underlying normal distribution
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Vector of N numbes sampled from the specified lognormal distribution
}
\references{
\url{http://en.wikipedia.org/wiki/Lognormal_distribution}
}
\author{
Paul Thomas Britton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dirichlet}}, \code{\link{triadist}}
}
\examples{
L <- logndist(5,median=.01,EF=2)
print(L)

## The function is currently defined as
function (N, EF = exp(qnorm(0.95)), sigma = log(EF)/qnorm(0.95), 
    mean = 1, median = exp(log(mean) - (sigma^2)/2), mu = log(median)) 
{
    return(rlnorm(N, mu, sigma))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pradist }
\keyword{ distributions }% __ONLY ONE__ keyword per line
