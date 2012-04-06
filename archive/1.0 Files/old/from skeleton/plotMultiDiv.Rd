\name{plotMultiDiv}
\alias{plotMultiDiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotMultiDiv(result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
%%     ~~Describe \code{result} here~~
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
function (result) 
{
    times <- apply(results[[1]], 1, mean)
    mdiv <- results[[3]]
    y_lim <- c(min(mdiv), max(mdiv))
    plot(times, mdiv[, 3], type = "n", ylim = y_lim, xlim = c(max(times), 
        min(times)), xlab = "Time (MYA)", ylab = "Lineage Richness", 
        main = paste("Median Diversity Curve from", ntr, "Trees"))
    polygon(c(times, rev(times)), c(mdiv[, 3], rev(mdiv[, 2])), 
        col = "gray", border = NA)
    lines(times, mdiv[, 1], lwd = 3)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
