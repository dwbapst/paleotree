\name{timePaleoPhy}
\alias{timePaleoPhy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
timePaleoPhy(tree, timeData, type = "basic", vartime = NULL, add.term = F, node.mins = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{timeData}{
%%     ~~Describe \code{timeData} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{vartime}{
%%     ~~Describe \code{vartime} here~~
}
  \item{add.term}{
%%     ~~Describe \code{add.term} here~~
}
  \item{node.mins}{
%%     ~~Describe \code{node.mins} here~~
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
function (tree, timeData, type = "basic", vartime = NULL, add.term = F, 
    node.mins = NULL) 
{
    require(ape)
    ntime <- sapply(1:Nnode(tree), function(x) max(timeData[tree$tip.label[unlist(prop.part(tree)[x])], 
        1]))
    ntime <- c(timeData[tree$tip.label, 1], ntime)
    if (length(node.mins) > 0) {
        require(phangorn)
        for (i in (Ntip(tree) + 1):length(ntime)) {
            desc_all <- unlist(Descendants(tree, i, type = "all"))
            desc_nodes <- c(desc_all[desc_all > Ntip(tree)], 
                i) - Ntip(tree)
            node_times <- node.mins[desc_nodes]
            ntime[i] <- max(ntime[i], node_times[!is.na(node_times)])
        }
    }
    if (type == "equal" & length(vartime) > 0) {
        ntime[Ntip(tree) + 1] <- vartime + ntime[Ntip(tree) + 
            1]
    }
    ttree <- tree
    ttree$edge.length <- sapply(1:Nedge(ttree), function(x) ntime[ttree$edge[x, 
        1]] - ntime[ttree$edge[x, 2]])
    if (type == "aba") {
        if (is.na(vartime)) {
            stop("No All Branch Additive Value Given!")
        }
        ttree$edge.length <- ttree$edge.length + vartime
    }
    if (type == "zlba") {
        if (is.na(vartime)) {
            stop("No Branch Additive Value Given!")
        }
        ttree$edge.length[ttree$edge.length < 1e-04] <- ttree$edge.length[ttree$edge.length < 
            1e-04] + vartime
    }
    if (type == "mbl") {
        require(phangorn)
        if (is.na(vartime)) {
            stop("No Minimum Branch Length Value Given!")
        }
        root_node <- Ntip(ttree) + 1
        mbl <- vartime
        while (any(ttree$edge.length < mbl)) {
            mom <- ttree$edge[((1:Nedge(ttree))[ttree$edge.length < 
                mbl])[sample(length((1:Nedge(ttree))[ttree$edge.length < 
                mbl]), 1)], 1]
            mom <- c(mom, Ancestors(ttree, mom))
            debt <- mbl - min(ttree$edge.length[ttree$edge[, 
                1] == mom[1]])
            ttree$edge.length[mom[1] == ttree$edge[, 1]] <- ttree$edge.length[mom[1] == 
                ttree$edge[, 1]] + debt[1]
            if (length(mom) > 1) {
                for (i in 2:length(mom)) {
                  small <- min(ttree$edge.length[ttree$edge[, 
                    1] == mom[i]])
                  mom_blen <- ttree$edge.length[ttree$edge[, 
                    1] == mom[i] & ttree$edge[, 2] == mom[i - 
                    1]]
                  debt[i] <- max(debt[i - 1] - max(mom_blen - 
                    mbl, 0), 0) + max(mbl - small, 0)
                  ttree$edge.length[ttree$edge[, 1] == mom[i] & 
                    ttree$edge[, 2] == mom[i - 1]] <- mom_blen - 
                    max(min(max(mom_blen - mbl, 0), debt[i - 
                      1]), 0) + max(mbl - small, 0)
                  ttree$edge.length[ttree$edge[, 1] == mom[i] & 
                    ttree$edge[, 2] != mom[i - 1]] <- ttree$edge.length[ttree$edge[, 
                    1] == mom[i] & ttree$edge[, 2] != mom[i - 
                    1]] + debt[i]
                }
            }
        }
    }
    if (type == "equal") {
        zbr <- cbind(1:Nedge(ttree), node.depth(ttree)[ttree$edge[, 
            2]])
        zbr <- zbr[ttree$edge.length == 0, ]
        zbr <- zbr[order(zbr[, 2]), 1]
        for (i in zbr) {
            if (ttree$edge.length[i] == 0) {
                brs <- ttree$edge[i, 2]
                mom <- which(ttree$edge[i, 1] == ttree$edge[, 
                  2])
                while (ttree$edge[mom, 1] != (Ntip(ttree) + 1) & 
                  ttree$edge.length[mom] == 0) {
                  brs[length(brs) + 1] <- ttree$edge[mom, 2]
                  mom <- which(ttree$edge[mom, 1] == ttree$edge[, 
                    2])
                }
                brs[length(brs) + 1] <- ttree$edge[mom, 2]
                totbl <- sum(ttree$edge.length[match(brs, ttree$edge[, 
                  2])])
                ntime[brs[-1]] <- ntime[brs[-1]] + cumsum(rep(totbl/length(brs), 
                  length(brs) - 1))
                ttree$edge.length <- sapply(1:Nedge(ttree), function(x) ntime[ttree$edge[x, 
                  1]] - ntime[ttree$edge[x, 2]])
            }
        }
    }
    tree$root.time <- max(timeData[ttree$tip.label, 1]) + min(dist.nodes(tree)[1:Ntip(tree), 
        Ntip(tree) + 1])
    if (add.term) {
        obs_ranges <- timeData[, 1] - timeData[, 2]
        term_id <- ttree$tip.label[ttree$edge[ttree$edge[, 2] <= 
            Ntip(ttree), 2]]
        term_add <- sapply(term_id, function(x) obs_ranges[x])
        ttree$edge.length[ttree$edge[, 2] <= Ntip(ttree)] <- ttree$edge.length[ttree$edge[, 
            2] <= Ntip(ttree)] + term_add
    }
    return(ttree)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
