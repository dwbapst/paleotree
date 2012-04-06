\name{taxa2phylo}
\alias{taxa2phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
taxa2phylo(taxad2, obs_time = NULL, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxad2}{
%%     ~~Describe \code{taxad2} here~~
}
  \item{obs_time}{
%%     ~~Describe \code{obs_time} here~~
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
function (taxad2, obs_time = NULL, plot = F) 
{
    require(ape)
    taxad1 <- taxad2[, 1:4]
    if (any((taxad1[, 4] - taxad1[, 3]) < 0)) {
        taxad1[, 3:4] <- max(taxad1[, 3:4]) - taxad1[, 3:4]
    }
    if (any((taxad1[, 4] - taxad1[, 3]) < 0)) {
        stop("Time Error!")
    }
    if (is.null(obs_time)) {
        obs <- taxad1[, 4]
    }
    else {
        obs <- max(taxad2[, 3:4]) - obs_time
    }
    if (nrow(taxad1) != length(obs)) {
        stop("#obs != #taxa !!")
    }
    fake_taxa <- matrix(sapply((1:nrow(taxad1))[!is.na(obs)], 
        function(x) c(nrow(taxad1) + x, taxad1[x, 1], obs[x], 
            obs[x])), , 4, byrow = T)
    fake_taxa[, 1] <- (1:nrow(fake_taxa)) + nrow(taxad1)
    taxad <- rbind(taxad1, fake_taxa)
    ntaxa <- nrow(taxad)
    desc <- lapply(taxad[, 1], function(x) taxad[c(F, taxad[-1, 
        2] == x), 1])
    births2 <- lapply(desc, function(x) if (length(x) > 0) {
        sapply(x, function(y) taxad[y == taxad[, 1], 3])
    })
    desc <- lapply(1:length(desc), function(x) if (length(desc[[x]]) > 
        0) {
        desc[[x]][match(1:length(births2[[x]]), rank(births2[[x]], 
            ties.method = "random"))]
    }
    else {
        numeric()
    })
    births <- lapply(desc, function(x) if (length(x) > 0) {
        sapply(x, function(y) taxad[y == taxad[, 1], 3])
    })
    events <- lapply(1:ntaxa, function(x) if (length(desc[[x]]) > 
        0) {
        c(taxad[x, 1], desc[[x]])
    }
    else {
        c(taxad[x, 1])
    })
    times <- lapply(1:ntaxa, function(x) if (length(births[[x]]) > 
        0) {
        c(taxad[x, 3], births[[x]], taxad[x, 4])
    }
    else {
        c(taxad[x, 3], taxad[x, 4])
    })
    nseg <- sapply(times, length) - 1
    seg_labs <- lapply(1:ntaxa, function(x) taxad[x, 1] + (1:nseg[x]/(nseg[x] + 
        1)))
    seg_labs <- lapply(seg_labs, function(x) c(floor(x[1]), x[-1]))
    taxa_anc <- c(0, sapply(2:ntaxa, function(x) unlist(seg_labs[taxad[x, 
        2] == taxad[, 1]])[which(unlist(events[taxad[x, 2] == 
        taxad[, 1]]) == taxad[x, 1]) - 1]))
    moms2 <- lapply(seg_labs, function(x) x[-length(x)])
    moms <- lapply(1:ntaxa, function(x) c(taxa_anc[x], unlist(moms2[[x]])))
    lengths <- lapply(times, diff)
    term <- sapply(unlist(seg_labs), function(x) !any(unlist(moms) == 
        x))
    edgeD <- data.frame(id = unlist(seg_labs), anc = unlist(moms), 
        brlen = unlist(lengths), term = term)
    MRCA <- min(edgeD[sapply(edgeD$id, function(x) 1 < sum(x == 
        edgeD$anc)), 1])
    edgeD <- edgeD[-which(edgeD[, 1] <= MRCA), ]
    droppers <- which(edgeD[, 4] & edgeD[, 1] < (nrow(taxad1) + 
        1))
    while (length(droppers) > 0) {
        edgeD <- edgeD[-droppers, ]
        edgeD[, 4] <- sapply(edgeD[, 1], function(x) !any(edgeD[, 
            2] == x))
        droppers <- which(as.logical(edgeD[, 4]) & edgeD[, 1] < 
            (nrow(taxad1) + 1))
    }
    edgeD[, 4] <- sapply(edgeD[, 1], function(x) !any(edgeD[, 
        2] == x))
    ndesc <- sapply(edgeD[, 1], function(x) sum(x == edgeD[, 
        2]))
    while (any(ndesc == 1)) {
        epick <- edgeD[ndesc == 1, ]
        if (is.data.frame(epick)) {
            epick <- epick[1, ]
        }
        edesc <- edgeD[edgeD$anc == epick$id, ]
        edgeD <- edgeD[-which(edgeD$id == edesc$id), ]
        newe <- data.frame(id = edesc$id, anc = epick$anc, brlen = (epick$brlen + 
            edesc$brlen), term = edesc$term)
        edgeD[edgeD$id == epick$id, ] <- newe
        ndesc <- sapply(edgeD$id, function(x) sum(x == edgeD$anc))
    }
    ndesc <- sapply(edgeD$id, function(x) sum(x == edgeD$anc))
    e_fix <- numeric()
    e_fix[edgeD$term] <- sapply(edgeD[edgeD$term, 1], function(x) which(fake_taxa[, 
        1] == x))
    e_fix[!edgeD$term] <- sum(edgeD$term) + 1 + (1:sum(!edgeD$term))
    ea_fix <- sapply(edgeD$anc, function(x) ifelse(x != MRCA, 
        e_fix[edgeD$id == x], sum(edgeD$term) + 1))
    tlabs <- paste("t", taxad1[!is.na(obs), 1], sep = "")
    edgf <- cbind(ea_fix, e_fix)
    colnames(edgf) <- NULL
    tree1 <- list(edge = edgf, tip.label = tlabs, edge.length = edgeD[, 
        3], Nnode = length(unique(edgf[, 1])))
    class(tree1) <- "phylo"
    tree <- reorder(collapse.singles(tree1), "cladewise")
    if (plot) {
        plot(ladderize(tree), show.tip.label = F)
        axisPhylo()
    }
    first_obs_time <- max(taxad2[, 3:4]) - min(obs, na.rm = T)
    tree$root.time <- first_obs_time + min(dist.nodes(tree)[1:Ntip(tree), 
        Ntip(tree) + 1])
    return(tree)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
