\name{bin_srcTimePaleoPhy}
\alias{bin_srcTimePaleoPhy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bin_srcTimePaleoPhy(tree, timeList, sampRate, ntrees = 1, sites = NULL, anc.wt = 1, node.mins = NULL, rand.obs = F, root.max = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{timeList}{
%%     ~~Describe \code{timeList} here~~
}
  \item{sampRate}{
%%     ~~Describe \code{sampRate} here~~
}
  \item{ntrees}{
%%     ~~Describe \code{ntrees} here~~
}
  \item{sites}{
%%     ~~Describe \code{sites} here~~
}
  \item{anc.wt}{
%%     ~~Describe \code{anc.wt} here~~
}
  \item{node.mins}{
%%     ~~Describe \code{node.mins} here~~
}
  \item{rand.obs}{
%%     ~~Describe \code{rand.obs} here~~
}
  \item{root.max}{
%%     ~~Describe \code{root.max} here~~
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
function (tree, timeList, sampRate, ntrees = 1, sites = NULL, 
    anc.wt = 1, node.mins = NULL, rand.obs = F, root.max = 200) 
{
    require(ape)
    if (is.null(sites)) {
        sites <- matrix(1:(Ntip(tree) * 2), , 2)
    }
    else {
        sites[, 1] <- sapply(sites[, 1], function(x) which(x == 
            sort(unique(as.vector(sites)))))
        sites[, 2] <- sapply(sites[, 2], function(x) which(x == 
            sort(unique(as.vector(sites)))))
    }
    rownames(sites) <- rownames(timeList[[2]])
    ttrees <- rmtree(ntrees, 3)
    for (ntrb in 1:ntrees) {
        no_bad_dates <- F
        while (!no_bad_dates) {
            siteTime <- matrix(, max(sites), 2)
            for (i in unique(as.vector(sites))) {
                go <- timeList[[2]][which(sites == i)[1]]
                siteTime[i, ] <- timeList[[1]][go, ]
            }
            siteDates <- apply(siteTime, 1, function(x) runif(1, 
                x[2], x[1]))
            timeData <- cbind(siteDates[sites[, 1]], siteDates[sites[, 
                2]])
            rownames(timeData) <- rownames(timeList[[2]])
            no_bad_dates <- all((timeData[, 1] - timeData[, 2]) >= 
                0)
        }
        if (rand.obs) {
            timeData[, 2] <- apply(timeData, 1, function(x) runif(1, 
                x[2], x[1]))
        }
        ttrees[[ntrb]] <- srcTimePaleoPhy(tree, timeData, sampRate, 
            ntrees = 1, anc.wt = anc.wt, node.mins = node.mins, 
            root.max = root.max, rand.obs = F)
    }
    if (ntrees == 1) {
        ttrees <- ttrees[[1]]
    }
    return(ttrees)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
