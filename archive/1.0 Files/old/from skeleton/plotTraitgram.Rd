\name{plotTraitgram}
\alias{plotTraitgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTraitgram(trait, tree, trait.name = "'trait'", conf.int = T, lwd = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trait}{
%%     ~~Describe \code{trait} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{trait.name}{
%%     ~~Describe \code{trait.name} here~~
}
  \item{conf.int}{
%%     ~~Describe \code{conf.int} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
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
function (trait, tree, trait.name = "'trait'", conf.int = T, 
    lwd = 1.5) 
{
    require(geiger)
    if (is.null(tree$root.time)) {
        tree$root.time <- max(dist.nodes(tree)[Ntip(tree) + 1, 
            1:Ntip(tree)])
    }
    times <- tree$root.time - dist.nodes(tree)[Ntip(tree) + 1, 
        ]
    if (conf.int) {
        asr <- ace(trait, tree, method = "pic")
        tr1 <- c(trait, asr$ace)
        edges <- tree$edge
        ci <- asr$CI95
        tr2 <- c(tr1, ci)
        plot(1, 1, type = "n", xlim = c(min(tr2) - 0.1, max(tr2) + 
            0.1), ylim = c(max(times), min(times)), xlab = "Trait Values", 
            ylab = "Time from Root", main = paste("Traitgram of", 
                trait.name))
        for (i in 1:nrow(edges)) {
            anc <- edges[i, 1]
            desc <- edges[i, 2]
            lines(c(tr1[anc], tr1[desc]), c(times[anc], times[desc]), 
                lwd = lwd)
        }
        for (i in 1:Nnode(tree)) {
            lines(c(ci[i, 1], ci[i, 2]), c(times[i + Ntip(tree)], 
                times[i + Ntip(tree)]), lwd = lwd)
        }
    }
    else {
        tr1 <- c(trait, getAncStates(trait, tree))
        edges <- tree$edge
        plot(1, 1, type = "n", xlim = c(min(tr1) - 0.1, max(tr1) + 
            0.1), ylim = c(max(times), min(times)), xlab = "Trait Values", 
            ylab = "Time from Root", main = paste("Traitgram of", 
                trait.name))
        for (i in 1:nrow(edges)) {
            anc <- edges[i, 1]
            desc <- edges[i, 2]
            lines(c(tr1[anc], tr1[desc]), c(times[anc], times[desc]), 
                lwd = lwd)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
