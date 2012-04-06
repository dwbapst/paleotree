\name{taxa2cladogram}
\alias{taxa2cladogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
taxa2cladogram(taxad, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxad}{
%%     ~~Describe \code{taxad} here~~
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
function (taxad, plot = F) 
{
    require(ape)
    tlabs <- paste("t", taxad[, 1], sep = "")
    desc <- lapply(taxad[, 1], function(x) (taxad[taxad[, 2] == 
        x, 1])[!is.na(taxad[taxad[, 2] == x, 1])])
    ndesc <- sapply(desc, length)
    rank <- numeric(length(ndesc))
    rank[ndesc == 0] <- 1
    rank[rank == 0] <- NA
    while (any(is.na(rank))) {
        rank <- sapply(1:length(rank), function(x) ifelse(!is.na(rank[x]), 
            rank[x], 1 + max(rank[sapply(desc[[x]], function(y) which(y == 
                taxad[, 1]))])))
    }
    comp <- numeric(length(ndesc))
    lab <- list()
    lab[rank == 1] <- tlabs[rank == 1]
    comp[rank == 1] <- 1
    while (any(comp == 0)) {
        tpot <- comp == 0
        tpot2 <- rank == min(rank[tpot])
        tpick <- which(tpot & tpot2)[1]
        dlab <- paste(unlist(lab[desc[[tpick]]]), ",", sep = "", 
            collapse = "")
        lab[[tpick]] <- paste("(", dlab, tlabs[tpick], ")", sep = "")
        comp[tpick] <- 1
    }
    tree1 <- paste(lab[[1]], ";", sep = "")
    tree2 <- read.tree(text = tree1)
    if (plot) {
        plot(ladderize(tree2), show.tip.label = F)
    }
    return(tree2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
