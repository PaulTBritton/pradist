\name{pradist-package}
\alias{pradist-package}
\alias{pradist}
\docType{package}
\title{
Addition distributions for PRA
}
\description{
Generate samples from triangle distributions, dirichlet distributions and a reparameterized lognormal distribution.
}
\details{
\tabular{ll}{
Package: \tab pradist\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-17\cr
License: \tab GPL-2\cr
}
}
\author{
Paul Thomas Britton

Maintainer: Paul Thomas Britton <paul.t.britton@nasa.gov>
}
\references{
\url{http://en.wikipedia.org/wiki/Triangle_distribution}
\url{http://en.wikipedia.org/wiki/Lognormal_distribution}
\url{http://en.wikipedia.org/wiki/Dirichlet_distribution}
}
\keyword{ package }
\keyword{ distributions }
\examples{
T <- triadist(10,1,2,3)
L <- logndist(5,median=.01,EF=2)
D <- dirichlet(N=20,300,20,1)
print(T)
print(L)
print(D$b1)
}
