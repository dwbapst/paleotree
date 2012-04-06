\name{multiDiv}
\alias{multiDiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
multiDiv(data, tblen = 1, calc.median = F, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{tblen}{
%%     ~~Describe \code{tblen} here~~
}
  \item{calc.median}{
%%     ~~Describe \code{calc.median} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (data, tblen = 1, calc.median = F, plot = T) 
{
    require(ape)
    dclass <- sapply(data, class)
    dclass1 <- numeric()
    dclass1[dclass == "matrix"] <- 1
    dclass1[dclass == "list"] <- 2
    dclass1[dclass == "phylo"] <- 3
    if (any(dclass1 == 1)) {
        lims1 <- sapply(data[dclass1 == 1], function(x) c(min(x), 
            max(x)))
    }
    else {
        lims1 <- NA
    }
    if (any(dclass1 == 2)) {
        lims2 <- sapply(data[dclass1 == 2], function(x) c(min(x[[1]][max(x[[2]]), 
            ]), max(x[[1]][min(x[[2]]), ])))
    }
    else {
        lims2 <- NA
    }
    if (any(dclass1 == 3)) {
        lims3 <- numeric()
        for (i in which(dclass1 == 3)) {
            ttree <- data[[i]]
            if (is.null(ttree$root.time)) {
                ntime <- dist.nodes(ttree)[, Ntip(ttree) + 1]
                ntime <- max(ntime) - ntime
            }
            else {
                ntime <- dist.nodes(ttree)[, Ntip(ttree) + 1]
                ntime <- ttree$root.time - ntime
            }
            lims3 <- c(lims3, c(min(ntime), max(ntime)))
        }
    }
    else {
        lims3 <- NA
    }
    end <- min(c(lims1, lims2, lims3), na.rm = T)
    start <- max(c(lims1, lims2, lims3), na.rm = T)
    midtimes <- seq(start + 2 * tblen, end - 2 * tblen, by = -tblen)
    int.start <- midtimes + (tblen/2)
    int.end <- midtimes - (tblen/2)
    int.times <- cbind(int.start, int.end)
    div <- matrix(, length(midtimes), 1)
    if (any(dclass1 == 1)) {
        divs1 <- sapply(data[dclass1 == 1], function(x) taxicDivCont(timeData = x, 
            int.times = int.times, output = T, plot = T)[, 3])
        div <- cbind(div, divs2)
    }
    if (any(dclass1 == 2)) {
        divs2 <- sapply(data[dclass1 == 2], function(x) taxicDivDisc(timeList = x, 
            int.times = int.times, output = T, plot = T)[, 3])
        div <- cbind(div, divs2)
    }
    if (any(dclass1 == 3)) {
        divs3 <- sapply(data[dclass1 == 3], function(x) phyloDiv(ttree = x, 
            int.times = int.times, output = T, plot = T)[, 3])
        div <- cbind(div, divs3)
    }
    div <- div[, -1]
    median <- apply(div, 1, median)
    q1 <- apply(div, 1, quantile, probs = 0.025)
    q2 <- apply(div, 1, quantile, probs = 0.975)
    median.curve <- cbind(median = median, low.95quantile = q1, 
        high.95quantile = q2)
    res <- list(int.times, div, median.curve)
    if (plot) {
        plot_multiDiv(res)
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
