\name{sampleRanges}
\alias{sampleRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleRanges(taxad, r, min.taxa = 2, ranges.only = T, modern.samp = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxad}{
%%     ~~Describe \code{taxad} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{min.taxa}{
%%     ~~Describe \code{min.taxa} here~~
}
  \item{ranges.only}{
%%     ~~Describe \code{ranges.only} here~~
}
  \item{modern.samp}{
%%     ~~Describe \code{modern.samp} here~~
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
function (taxad, r, min.taxa = 2, ranges.only = T, modern.samp = T) 
{
    redo <- T
    while (redo) {
        samp_occ <- list()
        for (i in 1:nrow(taxad)) {
            samps <- taxad[i, 3]
            while (min(samps) > taxad[i, 4]) {
                samps <- c(samps, min(samps) - rexp(1, rate = r))
            }
            if (taxad[i, 5] == 1) {
                samps <- c(samps, 0)
            }
            samps <- samps[-c(1, length(samps))]
            if (length(samps) > 0) {
                samp_occ[[i]] <- samps
            }
            else {
                samp_occ[[i]] <- NA
            }
        }
        redo <- min.taxa > sum(sapply(samp_occ, function(x) !any(is.na(x))))
    }
    names(samp_occ) <- paste("t", taxad[, 1], sep = "")
    if (ranges.only) {
        ranges <- cbind(sapply(samp_occ, max), sapply(samp_occ, 
            min))
        rownames(ranges) <- names(samp_occ)
        colnames(ranges) <- c("FAD", "LAD")
        res <- ranges
    }
    else {
        res <- samp_occ
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
