\name{phyloDiv}
\alias{phyloDiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phyloDiv(ttree, tblen = 0.1, int.times = NULL, output = F, plot = T, plotLogRich = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttree}{
%%     ~~Describe \code{ttree} here~~
}
  \item{tblen}{
%%     ~~Describe \code{tblen} here~~
}
  \item{int.times}{
%%     ~~Describe \code{int.times} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{plotLogRich}{
%%     ~~Describe \code{plotLogRich} here~~
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
function (ttree, tblen = 0.1, int.times = NULL, output = F, plot = T, 
    plotLogRich = F) 
{
    require(ape)
    if (is.null(ttree$root.time)) {
        ntime <- dist.nodes(ttree)[, Ntip(ttree) + 1]
        ntime <- max(ntime) - ntime
    }
    else {
        ntime <- dist.nodes(ttree)[, Ntip(ttree) + 1]
        ntime <- ttree$root.time - ntime
    }
    if (is.null(int.times)) {
        midtimes <- seq(max(ntime) + 3 * tblen, min(ntime) - 
            2 * tblen, by = -tblen)
        int.start <- midtimes + (tblen/2)
        int.end <- midtimes - (tblen/2)
    }
    else {
        int.start <- int.times[, 1]
        int.end <- int.times[, 2]
        midtimes <- (int.start + int.end)/2
    }
    LAD <- ntime[1:Ntip(ttree)]
    FAD <- ntime[(Ntip(ttree) + 1):length(ntime)]
    div <- sapply(1:length(midtimes), function(x) 1 + sum(FAD >= 
        int.end[x]) - sum(LAD > int.start[x]))
    if (plot) {
        layout(matrix(1:2, 2, 1))
        par(mar = c(1, 4, 1, 1))
        plot(ladderize(ttree), show.tip.label = F)
        par(mar = c(5, 4, 2, 2))
        if (logRich) {
            plot(midtimes, div, type = "l", log = "y", xlim = c(max(midtimes), 
                min(midtimes)), xlab = "Time (MY)", ylab = "Lineage Richness (Log Scale)")
        }
        else {
            plot(midtimes, div, type = "l", xlim = c(max(midtimes), 
                min(midtimes)), ylim = c(0, max(div) + 2), xlab = "Time (MY)", 
                ylab = "Lineage Richness")
        }
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
