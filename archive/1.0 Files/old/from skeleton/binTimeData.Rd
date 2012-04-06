\name{binTimeData}
\alias{binTimeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
binTimeData(timeData, bin.length = 1, start = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeData}{
%%     ~~Describe \code{timeData} here~~
}
  \item{bin.length}{
%%     ~~Describe \code{bin.length} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
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
function (timeData, bin.length = 1, start = NA) 
{
    if (is.na(start)) {
        start <- max(timeData)
    }
    else {
        if (start < max(timeData)) {
            warning("Start<max(timeData)??")
        }
    }
    end <- start - ceiling((start - min(timeData))/bin.length) * 
        bin.length
    bins <- seq(start, end, by = -bin.length)
    bins <- ifelse(bins < 0, 0, bins)
    fads <- sapply(timeData[, 1], function(x) which(bins < x)[1] - 
        1)
    lads <- sapply(timeData[, 2], function(x) which(bins < x)[1] - 
        1)
    res <- list(int.times = cbind(int.start = bins[1:(length(bins) - 
        1)], int.end = bins[2:length(bins)]), taxon.times = cbind(first.int = fads, 
        last.int = lads))
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
