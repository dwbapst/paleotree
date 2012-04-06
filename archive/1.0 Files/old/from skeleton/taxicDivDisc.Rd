\name{taxicDivDisc}
\alias{taxicDivDisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
taxicDivDisc(timeList, plot = T, int.times = NULL, output = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeList}{
%%     ~~Describe \code{timeList} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{int.times}{
%%     ~~Describe \code{int.times} here~~
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
function (timeList, plot = T, int.times = NULL, output = F) 
{
    intMat <- timeList[[1]]
    timeData <- timeList[[2]]
    Fint <- as.numeric(timeData[, 1])
    Lint <- as.numeric(timeData[, 2])
    FAD <- intMat[Fint, 1]
    LAD <- intMat[Lint, 1]
    if (is.null(int.times)) {
        avg_dur <- abs(mean(apply(timeList[[1]], 1, diff)))
        intMat <- rbind(intMat, c(max(intMat) + avg_dur, max(intMat)), 
            c(min(intMat), min(intMat) - avg_dur))
        intMat <- intMat[order(-intMat[, 1]), ]
        int.start <- intMat[, 1]
        int.end <- intMat[, 2]
        midtimes <- apply(intMat, 1, mean)
    }
    else {
        int.start <- int.times[, 1]
        int.end <- int.times[, 2]
        midtimes <- (int.start + int.end)/2
    }
    div <- sapply(1:length(midtimes), function(x) sum(FAD > int.end[x]) - 
        sum(LAD > int.start[x]))
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
