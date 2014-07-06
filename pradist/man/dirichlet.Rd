\name{dirichlet}
\alias{dirichlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create N samples from a dirichlet distribution
}
\description{
	Returns a list of k vectors, each of size N, where each of the k vectors corresonds to the margins of the dirichlet distribution
}
\usage{
dirichlet(N, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
	Number of samples
}
  \item{\dots}{
	The alpha parameters of the dirichlet distribution
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	List of vectors. The vectors contain N samples each from the margins of the specified dirichlet distribution.
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{b1 }{Marginal distribution in dimension 1}
  \item{\dots }{\dots}
  \item{bk }{Marginal distribution in dimension k}
%% ...
}
\references{
\url{http://en.wikipedia.org/wiki/Dirichlet_distribution}
}
\author{
Paul Thomas Britton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{logndist}}, \code{\link{triadist}}
}
\examples{
diri <- dirichlet(N=20,300,20,1)
print(diri$b1)
## or
E <- dirichlet(N=10,c(400,20))
print(E$b2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pradist }
\keyword{ distributions }% __ONLY ONE__ keyword per line
