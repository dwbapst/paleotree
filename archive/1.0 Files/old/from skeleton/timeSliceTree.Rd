\name{timeSliceTree}
\alias{timeSliceTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
timeSliceTree(ttree, sliceTime, drop.extinct = F, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttree}{
%%     ~~Describe \code{ttree} here~~
}
  \item{sliceTime}{
%%     ~~Describe \code{sliceTime} here~~
}
  \item{drop.extinct}{
%%     ~~Describe \code{drop.extinct} here~~
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
function (ttree, sliceTime, drop.extinct = F, plot = T) 
{
    require(ape)
    if (is.null(ttree$root.time)) {
        ttree$root.time <- max(dist.nodes(ttree)[1:Ntip(ttree), 
            Ntip(ttree) + 1])
        print("Warning: no ttree$root.time! Assuming latest tip is at present (time=0)")
    }
    tslice <- ttree$root.time - sliceTime
    dnode <- dist.nodes(ttree)[, Ntip(ttree) + 1]
    cedge <- which(sapply(1:Nedge(ttree), function(x) any(ttree$edge[x, 
        1] == which(dnode < tslice)) & any(ttree$edge[x, 2] == 
        which(dnode >= tslice))))
    droppers <- numeric()
    for (i in 1:length(cedge)) {
        desc <- ttree$edge[cedge[i], 2]
        if (desc > Ntip(ttree)) {
            desctip <- unlist(prop.part(ttree)[desc - Ntip(ttree)])
            droppers <- c(droppers, desctip[-1])
        }
    }
    stree <- drop.tip(ttree, droppers)
    dnode <- dist.nodes(stree)[, Ntip(stree) + 1]
    cedge <- sapply(1:Nedge(stree), function(x) any(stree$edge[x, 
        2] == which(dnode >= tslice)))
    cnode_depth <- dnode[stree$edge[cedge, 1]]
    stree$edge.length[cedge] <- tslice - cnode_depth
    stree$root.time <- ttree$root.time
    if (drop.extinct) {
        dnode <- dist.nodes(stree)[1:Ntip(stree), Ntip(stree) + 
            1]
        droppers <- which((dnode + 0.01) < tslice)
        stree1 <- drop.tip(stree, droppers)
        stree1$root.time <- sliceTime + max(dist.nodes(stree1)[1:Ntip(stree1), 
            Ntip(stree1) + 1])
    }
    else {
        stree1 <- stree
    }
    if (plot) {
        layout(matrix(1:2, , 2))
        plot(ladderize(ttree), show.tip.label = F)
        axisPhylo()
        plot(stree1, show.tip.label = F)
    }
    return(stree1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
