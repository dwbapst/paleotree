\name{probAnc}
\alias{probAnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
probAnc(p, q, R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, q, R) 
{
    Pd <- function(q, T) {
        exp(-q * (T - 1)) - exp(-q * T)
    }
    PN <- function(p, T, N) {
        (exp(-p * T) * (p * T)^N)/factorial(N)
    }
    Qm <- function(p, q, M) {
        x <- (4 * p * q)/((p + q)^2)
        ((p + q)/(2 * p)) * (factorial(2 * M)/((2^(2 * M)) * 
            factorial(M)^2)) * ((x^M)/((2 * M) - 1))
    }
    Pinf <- function(p, q, Pp) {
        res <- numeric()
        for (M in 1:80) {
            res[M] <- Qm(p, q, M) * (1 - (1 - Pp)^M)
        }
        sum(res)
    }
    Pp <- Pp4R(R, q)
    res <- numeric()
    for (t in 1:2000) {
        Nres <- numeric()
        for (N in 1:100) {
            Nres[N] <- PN(p, t, N) * (1 - (Pp)^N)
        }
        res[t] <- (((1 - (1 - R)^t) * Pd(q, t))/Pp) * sum(Nres)
    }
    res <- sum(res)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
