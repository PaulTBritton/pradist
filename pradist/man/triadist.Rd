\name{triadist}
\alias{triadist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create a vector of N triangular samples
}
\description{
	Create a vector of N triangular samples from a triangular distribution with parameters a, b and c, where N > 0, a <= b and b <= c.
}
\usage{
triadist(N, a, b, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
	Number of samples
}
  \item{a}{
	a = minimum
}
  \item{b}{
	b = mode
}
  \item{c}{
	c = maximum
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector of N numbers sampled from a triangular distribution with minimum = a, mode = b and maximum = c.
}
\references{
\url{http://en.wikipedia.org/wiki/Triangle_distribution}
}
\author{
Paul Thomas Britton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dirichlet}}, \code{\link{logndist}}
}
\examples{
T <- triadist(10,1,2,3)
print(T)

## The function is currently defined as
function (N, a, b, c) 
{
    if (!((a <= b) && (b <= c))) {
        stop("triadist: Invalid Triangle Distribution Parameters")
    }
    U <- runif(N, 0, 1)
    x <- (b - a)/(c - a)
    Tri <- U <= x
    Tri[Tri == TRUE] <- a + sqrt(U[U <= x] * (c - a) * (b - a))
    Tri[Tri == FALSE] <- c - sqrt((1 - U[U > x]) * (c - a) * 
        (c - b))
    return(Tri)
  }
}
\keyword{ pradist }
\keyword{ distributions }
