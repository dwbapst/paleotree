\name{taxicDivCont}
\alias{taxicDivCont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
taxicDivCont(timeData, tblen = 0.1, int.times = NULL, plot = T, output = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeData}{
%%     ~~Describe \code{timeData} here~~
}
  \item{tblen}{
%%     ~~Describe \code{tblen} here~~
}
  \item{int.times}{
%%     ~~Describe \code{int.times} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
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
function (timeData, tblen = 0.1, int.times = NULL, plot = T, 
    output = F) 
{
    FAD <- as.numeric(timeData[, 1])
    LAD <- as.numeric(timeData[, 2])
    if (is.null(int.times)) {
        midtimes <- seq(max(FAD) + 2 * tblen, min(LAD) - 2 * 
            tblen, by = -tblen)
        int.start <- midtimes + (tblen/2)
        int.end <- midtimes - (tblen/2)
    }
    else {
        int.start <- int.times[, 1]
        int.end <- int.times[, 2]
        midtimes <- (int.start + int.end)/2
    }
    div <- sapply(1:length(midtimes), function(x) sum(FAD >= 
        int.end[x]) - sum(LAD > int.start[x]))
    if (plot) {
        plot(midtimes, div, type = "l", xlim = c(max(midtimes), 
            min(midtimes)), xlab = "Time (MY)", ylab = "Taxic Richness")
    }
    if (output) {
        res <- cbind(int.start, int.end, int.div = div)
        return(res)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
