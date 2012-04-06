\name{pqsRate2sProb}
\alias{pqsRate2sProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pqsRate2sProb(r, p, q, dt = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{dt}{
%%     ~~Describe \code{dt} here~~
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
function (r, p, q, dt = 1) 
{
    PDbt <- function(r, dt) {
        1 - exp(-r * dt)
    }
    PDbL <- function(q, r, dt) {
        (((r + (q * exp(-(q + r) * dt)))/(q + r)) - exp(-q * 
            dt))/(1 - exp(-q * dt))
    }
    PDFt <- function(p, r, dt) {
        (((r + (p * exp(-(p + r) * dt)))/(p + r)) - exp(-p * 
            dt))/(1 - exp(-p * dt))
    }
    PDFL <- function(p, q, r, dt) {
        if (p == q) {
            NbNFL <- 1/(exp(-q * dt) + (p * dt) - 1)
            term1 <- (r * dt)/(p + r)
            term2 <- (1 - exp(-p * dt))/p
            term3 <- (p * (1 - exp(-(p + r) * dt)))/((p + r)^2)
            terms <- term1 - term2 + term3
            res <- (NbNFL) * p * terms
        }
        else {
            NbNFL <- 1/(((q * exp((p - q) * dt)) + ((p - q) * 
                exp(-q * dt)) - p)/(p - q))
            term1 <- (p * r * (exp((p - q) * dt) - 1))/((q + 
                r) * (p - q))
            term2 <- (p * q * exp(-(q + r) * dt) * (exp((p + 
                r) * dt) - 1))/((p + r) * (q + r))
            term3 <- exp(-q * dt) * (exp(p * dt) - 1)
            terms <- term1 + term2 - term3
            res <- (NbNFL) * terms
        }
        res
    }
    Pbt <- exp(-q * dt)
    PbL <- (1 - exp(-q * dt))
    PFt <- exp((p - q) * dt) * (1 - exp(-p * dt))
    if (p == q) {
        PFL <- exp(-q * dt) + (p * dt) - 1
    }
    else {
        PFL <- ((q * exp((p - q) * dt)) + ((p - q) * exp(-q * 
            dt)) - p)/(p - q)
    }
    res <- sum(PDbt(r, dt) * Pbt, PDbL(q, r, dt) * PbL, PDFt(p, 
        r, dt) * PFt, PDFL(p, q, r, dt) * PFL)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
