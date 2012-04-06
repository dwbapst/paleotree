\name{simFossilTaxa}
\alias{simFossilTaxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simFossilTaxa(p, q, w = 0, u = 0, nruns = 1, mintaxa = 1, maxtaxa = 1000, maxtime = 100, nExtant = NULL, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
}
  \item{nruns}{
%%     ~~Describe \code{nruns} here~~
}
  \item{mintaxa}{
%%     ~~Describe \code{mintaxa} here~~
}
  \item{maxtaxa}{
%%     ~~Describe \code{maxtaxa} here~~
}
  \item{maxtime}{
%%     ~~Describe \code{maxtime} here~~
}
  \item{nExtant}{
%%     ~~Describe \code{nExtant} here~~
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
function (p, q, w = 0, u = 0, nruns = 1, mintaxa = 1, maxtaxa = 1000, 
    maxtime = 100, nExtant = NULL, plot = F) 
{
    results <- list()
    if (maxtaxa > 10000 & maxtime > 10000) {
        warning("Please provide a realistic maxtaxa or maxtime limit...")
    }
    if (!is.null(nExtant)) {
        if (nExtant > maxtaxa) {
            warning("Error: nExtant is set higher than maxtaxa")
        }
    }
    for (i in 1:nruns) {
        taxad <- matrix(c(1, NA, 0, NA), 1, )
        pqw <- p + q + w
        maxtime1 <- maxtime
        continue <- T
        while (any(is.na(taxad[, 4])) & continue) {
            tpot <- is.na(taxad[, 4])
            tpot2 <- min(taxad[tpot, 3]) == taxad[, 3]
            tpick <- which(tpot & tpot2)[1]
            tpick_FO <- taxad[tpick, 3]
            wait <- 0
            while (is.na(taxad[tpick, 4]) & continue) {
                wait <- rexp(1, rate = pqw) + wait
                type <- sample(1:3, 1, prob = c(p/pqw, q/pqw, 
                  w/pqw))
                if (type == 1) {
                  type1 <- sample(1:2, 1, prob = c(1 - u, u))
                  if (type1 == 1) {
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA))
                  }
                  if (type1 == 2) {
                    taxad[tpick, 4] <- wait + tpick_FO
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA))
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA))
                  }
                }
                if (type == 2) {
                  taxad[tpick, 4] <- wait + tpick_FO
                }
                if (type == 3) {
                  taxad[tpick, 4] <- wait + tpick_FO
                  taxad <- rbind(taxad, c(max(taxad[, 1]) + 1, 
                    tpick, wait + tpick_FO, NA))
                }
                if (nrow(taxad) >= (maxtaxa + 1)) {
                  maxtime1 <- taxad[maxtaxa + 1, 3]
                }
                if (!is.null(nExtant)) {
                  if (nExtant > 0) {
                    if (sum(is.na(taxad[, 4])) >= (nExtant + 
                      1) & ifelse(!is.null(mintaxa), nrow(taxad) > 
                      mintaxa, T)) {
                      maxtime1 <- max(taxad[is.na(taxad[, 4]), 
                        3])
                    }
                  }
                }
                continue <- ifelse(any(is.na(taxad[, 4])), any(taxad[is.na(taxad[, 
                  4]), 3] <= maxtime1), F)
            }
            mint <- ifelse(!is.null(mintaxa), nrow(taxad) >= 
                mintaxa, T)
            rext <- ifelse(!is.null(nExtant), nExtant == (sum(is.na(taxad[, 
                4])) + sum(taxad[!is.na(taxad[, 4]), 4] > maxtime1)), 
                T)
            if ((!rext | !mint) & !continue) {
                taxad <- matrix(c(1, NA, 0, NA), 1, )
                continue <- T
                maxtime1 <- maxtime
            }
        }
        taxad3 <- matrix(taxad[taxad[, 3] < maxtime1, ], sum(taxad[, 
            3] < maxtime1), )
        taxad1 <- cbind(taxad3, (taxad3[, 4] >= maxtime1))
        taxad1[, 3:4] <- maxtime1 - taxad1[, 3:4]
        taxad1[, 3] <- round(taxad1[, 3], digits = 4)
        taxad1[, 4] <- round(taxad1[, 4], digits = 4)
        taxad1[taxad1[, 4] < 0, 4] <- 0
        taxad2 <- cbind(matrix(1:nrow(taxad1), , 1), matrix(match(taxad1[, 
            2], taxad1[, 1]), , 1), matrix(taxad1[, 3:5], , 3))
        results[[i]] <- taxad2
        if (plot) {
            taxicDivCont(results[[i]][, 3:4, drop = F])
        }
    }
    if (nruns == 1) {
        results <- results[[1]]
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
