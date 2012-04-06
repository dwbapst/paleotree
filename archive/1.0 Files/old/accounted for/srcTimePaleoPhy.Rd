\name{srcTimePaleoPhy}
\alias{srcTimePaleoPhy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
srcTimePaleoPhy(tree, timeData, sampRate, ntrees = 1, anc.wt = 1, rand.obs = F, node.mins = NULL, root.max = 200, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{timeData}{
%%     ~~Describe \code{timeData} here~~
}
  \item{sampRate}{
%%     ~~Describe \code{sampRate} here~~
}
  \item{ntrees}{
%%     ~~Describe \code{ntrees} here~~
}
  \item{anc.wt}{
%%     ~~Describe \code{anc.wt} here~~
}
  \item{rand.obs}{
%%     ~~Describe \code{rand.obs} here~~
}
  \item{node.mins}{
%%     ~~Describe \code{node.mins} here~~
}
  \item{root.max}{
%%     ~~Describe \code{root.max} here~~
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
function (tree, timeData, sampRate, ntrees = 1, anc.wt = 1, rand.obs = F, 
    node.mins = NULL, root.max = 200, plot = F) 
{
    require(ape)
    require(phangorn)
    if (length(sampRate) == 1) {
        sampRate <- rep(sampRate, Ntip(tree))
        names(sampRate) <- tree$tip.label
    }
    else {
        if (length(sampRate) != Ntip(tree)) {
            stop("SR Length != Ntip!")
        }
    }
    if (length(node.mins) != Nnode(tree) & !is.null(node.mins)) {
        stop("node.mins length != Nnode!")
    }
    ttree1 <- timePaleoPhy(tree, timeData, type = "basic", node.mins = node.mins, 
        add.term = F)
    if (length(node.mins) > 0) {
        locked_nodes <- which(!is.na(node.mins)) + +Ntip(tree)
    }
    else {
        locked_nodes <- NA
    }
    ttree1 <- collapse.singles(ttree1)
    ttrees <- rmtree(ntrees, 3)
    for (ntr in 1:ntrees) {
        if (rand.obs) {
            timeData1 <- cbind(timeData[, 1], apply(timeData, 
                1, function(x) runif(1, x[2], x[1])))
        }
        else {
            timeData1 <- timeData
        }
        ktree <- ttree1
        nodes <- (1:Nnode(ktree)) + Ntip(ktree)
        nodes <- nodes[order(-node.depth(ktree)[-(1:Ntip(ktree))])]
        anags <- character()
        budds <- character()
        while (length(nodes) > 0) {
            node <- nodes[1]
            tipl <- ktree$tip.label
            tipd <- cbind(ID = (1:Ntip(ttree1)), FAD = (timeData1[tipl, 
                1]), LAD = (timeData1[tipl, 2]), SR = (sampRate[tipl]))
            if (node == (Ntip(ktree) + 1)) {
                min_zip <- (-root.max)
                root_push <- -seq(min_zip, 0, by = 0.1)
            }
            else {
                min_zip <- (-ktree$edge.length[ktree$edge[, 2] == 
                  node])
                stem_len <- ktree$edge.length[ktree$edge[, 2] == 
                  node]
            }
            dnodes <- ktree$edge[ktree$edge[, 1] == node, 2]
            dlen <- ktree$edge.length[match(dnodes, ktree$edge[, 
                2])]
            minlocked <- ifelse(!all(is.na(locked_nodes)), any(node == 
                locked_nodes), F)
            if (length(dnodes) > 2) {
                dSR <- numeric()
                drng <- numeric()
                for (i in dnodes) {
                  dtips <- match(unlist(Descendants(ktree, i)), 
                    tipd[, 1])
                  dearly <- which(tipd[dtips, 2] == max(tipd[dtips, 
                    2]))
                  dSR[length(dSR) + 1] <- tipd[dearly, 4]
                  drng[length(drng) + 1] <- ifelse(length(dtips) > 
                    1, NA, diff(unlist(tipd[dtips, 3:2])))
                }
                zip_wt <- (dSR * exp(-dSR * dlen))/sum(dSR * 
                  exp(-dSR * dlen))
                dnode1 <- sample(dnodes, 1, prob = zip_wt)
                if (node == (Ntip(ktree) + 1)) {
                  root_prob <- dSR[dnode1 == dnodes] * exp(-dSR[dnode1 == 
                    dnodes] * root_push)
                  stem_len <- sample(root_push, 1, prob = root_prob)
                }
                plin <- c(dnode1, (dlen[dnode1 == dnodes] + stem_len), 
                  drng[dnode1 == dnodes], NA, 0, dlen[dnode1 == 
                    dnodes] + stem_len, dlen[dnode1 == dnodes] + 
                    drng[dnode1 == dnodes] + stem_len)
                plin <- matrix(plin, 1, )
                colnames(plin) <- c("node", "brl", "rng", "anc", 
                  "tSpec", "tFO", "tLO")
                add_nodes <- cbind(dnodes, dlen + stem_len, drng, 
                  dSR)[-which(dnodes == dnode1), ]
                add_nodes <- add_nodes[order(dlen[-which(dnodes == 
                  dnode1)]), ]
                colnames(add_nodes) <- c("node", "dstem", "rng", 
                  "SR")
                for (i in 1:nrow(add_nodes)) {
                  zips <- matrix(, 1, 3)
                  for (j in 1:nrow(plin)) {
                    min_zip <- plin[j, 5]
                    max_zip <- ifelse(anc.wt > 0 & !is.na(plin[j, 
                      7]), min(add_nodes[i, 2], plin[j, 7]), 
                      min(add_nodes[i, 2], plin[j, 6]))
                    if (minlocked) {
                      max_zip <- 0
                    }
                    poss_zip <- seq(min_zip, max_zip, by = 0.1)
                    gap <- add_nodes[i, 2] - poss_zip
                    like <- ifelse(poss_zip > plin[j, 6], anc.wt * 
                      add_nodes[i, 4] * exp(-add_nodes[i, 4] * 
                      gap), add_nodes[i, 4] * exp(-add_nodes[i, 
                      4] * gap))
                    new_zip <- cbind(like, plin[j, 1], poss_zip)
                    zips <- rbind(zips, new_zip)
                  }
                  if (nrow(zips) < 3) {
                    zips <- matrix(zips[-1, ], 1, 3)
                  }
                  else {
                    zips <- zips[-1, ]
                  }
                  colnames(zips) <- c("like", "anc", "tzip")
                  zip_prob <- zips[, 1]/sum(zips[, 1])
                  ch_zip <- sample(1:nrow(zips), 1, prob = zip_prob)
                  ch_anc <- zips[ch_zip, 2]
                  ch_tzip <- zips[ch_zip, 3]
                  if (!is.na(plin[ch_anc == plin[, 1], 7])) {
                    if (plin[ch_anc == plin[, 1], 7] == ch_tzip) {
                      anags <- c(anags, ktree$tip.label[ch_anc])
                    }
                    if (plin[ch_anc == plin[, 1], 6] < ch_tzip) {
                      budds <- c(budds, ktree$tip.label[ch_anc])
                    }
                  }
                  new_lin <- c(add_nodes[i, 1], add_nodes[i, 
                    2] - ch_tzip, add_nodes[i, 3], ch_anc, ch_tzip, 
                    add_nodes[i, 2], add_nodes[i, 2] + add_nodes[i, 
                      3])
                  plin <- rbind(plin, new_lin)
                }
                taxad_o <- t(apply(plin, 1, function(x) c(x[1], 
                  x[4], x[5], ifelse(is.na(x[7]), x[6], x[7]))))
                new_anc <- sapply(taxad_o[, 2], function(x) ifelse(is.na(x), 
                  NA, which(taxad_o[, 1] == x)))
                taxad_n <- cbind(1:nrow(taxad_o), new_anc, taxad_o[, 
                  3:4])
                taxad_n[, 3:4] <- max(taxad_n[, 3:4]) - taxad_n[, 
                  3:4]
                subtree <- taxa2phylo(taxad_n)
                subtree$tip.label <- taxad_o[match(subtree$tip.label, 
                  paste("t", 1:nrow(taxad_o), sep = "")), 1]
                new_stem <- diff(sort(plin[, 5]))[1]
                for (i in dnodes) {
                  dtip <- which(subtree$tip.label == i)
                  if (i > Ntip(ktree)) {
                    subclade <- extract.clade(ktree, i)
                    subtree <- bind.tree(subtree, subclade, where = dtip)
                    subtree <- collapse.singles(subtree)
                  }
                  else {
                    subtree$tip.label[dtip] <- ktree$tip.label[i]
                  }
                }
                if (node != (Ntip(ktree) + 1)) {
                  drtips <- unlist(Descendants(ktree, node))
                  tip_lab <- ktree$tip.label[drtips[1]]
                  droptree <- collapse.singles(drop.tip(ktree, 
                    drtips[-1]))
                  droptree$edge.length[droptree$edge[, 2] == 
                    which(droptree$tip.label == tip_lab)] <- new_stem
                  droptree <- bind.tree(droptree, subtree, where = which(droptree$tip.label == 
                    tip_lab))
                  ktree1 <- droptree
                }
                else {
                  ktree1 <- subtree
                }
                if (length(nodes) > 1) {
                  d_o <- lapply(Descendants(ktree, nodes[-1]), 
                    function(x) ktree$tip.label[x])
                  d_n <- lapply(Descendants(ktree1)[-(1:Ntip(ktree1))], 
                    function(x) ktree1$tip.label[x])
                  nodes1 <- sapply(d_o, function(x) which(sapply(d_n, 
                    function(y) ifelse(length(y) == length(x), 
                      all(sort(y) == sort(x)), F))))
                  nodes1 <- Ntip(ktree1) + nodes1
                  nodes1 <- nodes1[order(-node.depth(ktree1)[nodes1])]
                  if (!all(is.na(locked_nodes))) {
                    d_ol <- lapply(Descendants(ktree, locked_nodes), 
                      function(x) ktree$tip.label[x])
                    locked_nodes <- sapply(d_ol, function(x) which(sapply(d_n, 
                      function(y) ifelse(length(y) == length(x), 
                        all(sort(y) == sort(x)), F))))
                    locked_nodes <- Ntip(ktree1) + locked_nodes
                  }
                }
                else {
                  nodes1 <- numeric()
                }
                ktree1 <- collapse.singles(ktree1)
                ktree <- ktree1
                nodes <- nodes1
            }
            else {
                dlen1 <- min(dlen)
                dlen2 <- max(dlen)
                dnode1 <- dnodes[which(dlen == dlen1)[1]]
                dnode2 <- dnodes[dnodes != dnode1]
                d2FADs <- tipd[match(unlist(Descendants(ktree, 
                  dnode2)), tipd[, 1]), 2]
                d2early <- which(d2FADs == max(d2FADs))
                d2early <- ifelse(length(d2early) > 1, sample(d2early, 
                  1), d2early)
                d2SR <- (tipd[match(unlist(Descendants(ktree, 
                  dnode2)), tipd[, 1]), 4])[d2early]
                d1rng <- ifelse(dnode1 <= Ntip(ktree), diff(unlist(tipd[match(dnode1, 
                  tipd[, 1]), 3:2])), NA)
                d2rng <- ifelse(dnode2 <= Ntip(ktree), diff(unlist(tipd[match(dnode2, 
                  tipd[, 1]), 3:2])), NA)
                max_zip <- ifelse(anc.wt > 0 & !is.na(d1rng), 
                  min(dlen1 + d1rng, dlen2), dlen1)
                if (minlocked) {
                  max_zip <- 0
                }
                poss_zip <- seq(min_zip, max_zip, by = 0.1)
                gap2 <- dlen2 - poss_zip
                prob_zip <- d2SR * exp(-d2SR * gap2)
                prob_zip <- ifelse(poss_zip > dlen1, prob_zip * 
                  anc.wt, prob_zip)
                prob_zip <- prob_zip/sum(prob_zip)
                ch_zip <- sample(poss_zip, 1, prob = prob_zip)
                new_dlen1 <- ifelse(ch_zip > dlen1, NA, dlen1 - 
                  ch_zip)
                new_dlen2 <- dlen2 - ch_zip
                if (is.na(new_dlen1)) {
                  if (ch_zip == max_zip & dlen1 + d1rng < dlen2) {
                    anags <- c(anags, ktree$tip.label[dnode1])
                    new_dlen1 <- 0
                  }
                  else {
                    budds <- c(budds, ktree$tip.label[dnode1])
                    new_dlen1 <- d1rng + dlen1 - ch_zip
                  }
                }
                else {
                  new_dlen1 <- ifelse(!is.na(d1rng), new_dlen1 + 
                    d1rng, new_dlen1)
                }
                new_dlen2 <- ifelse(!is.na(d2rng), new_dlen2 + 
                  d2rng, new_dlen2)
                if (node != (Ntip(ktree) + 1)) {
                  ktree$edge.length[ktree$edge[, 2] == node] <- ch_zip - 
                    min_zip
                }
                ktree$edge.length[match(dnode1, ktree$edge[, 
                  2])] <- new_dlen1
                ktree$edge.length[match(dnode2, ktree$edge[, 
                  2])] <- new_dlen2
                nodes <- nodes[-1]
            }
        }
        ktree <- reorder(collapse.singles(ktree), "cladewise")
        ktree$anag.tips <- anags
        ktree$budd.tips <- budds
        ktree$root.time <- max(timeData1[ktree$tip.label, 2]) + 
            min(dist.nodes(ktree)[1:Ntip(ktree), Ntip(ktree) + 
                1])
        names(ktree$edge.length) <- NULL
        names(ktree$tip.label) <- NULL
        names(ktree$budd.tips) <- NULL
        names(ktree$anag.tips) <- NULL
        tipdiffs <- cbind(diff(sort(-timeData1[, 2])), diff(sort(dist.nodes(ktree)[1:Ntip(ktree), 
            Ntip(ktree) + 1])), diff(sort(-timeData1[, 2])) - 
            diff(sort(dist.nodes(ktree)[1:Ntip(ktree), Ntip(ktree) + 
                1])))
        test1 <- all(tipdiffs[, 3] < (10^-10))
        test2 <- identical(names(sort(-timeData1[, 2])), ktree$tip.label[order(dist.nodes(ktree)[1:Ntip(ktree), 
            Ntip(ktree) + 1])])
        if (length(unique(timeData1[, 2])) < Ntip(tree)) {
            test2 <- T
        }
        if (all(c(test1, test2))) {
            ktree$test <- "passed"
        }
        else {
            warning("ERROR: TERMINAL TIPS IMPROPERLY ALIGNED, CAUSE UNKNOWN")
        }
        if (plot) {
            par(mar = c(2.5, 2.5, 1, 2.5))
            layout(matrix(1:3, 3, ))
            plot(ladderize(tree), show.tip.label = T, use.edge.length = F)
            plot(ladderize(ttree1), show.tip.label = T)
            axisPhylo()
            plot(ladderize(ktree), show.tip.label = T)
            axisPhylo()
        }
        ttrees[[ntr]] <- ktree
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
