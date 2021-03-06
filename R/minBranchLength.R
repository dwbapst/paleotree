#' Scales Edge Lengths of a Phylogeny to a Minimum Branch Length
#' 
#' Rescales a tree with edge lengths so that all edge lengths
#' are at least some minimum branch length 
#' (sometimes abbreviated as "\code{MBL}" or "\code{mbl}").
#' Edge lengths are transformed so they are
#' greater than or equal to the input minimum branch length, by
#' subtracting edge length from more root-ward edges
#' and added to later branches. 
#' This may or may not change the age of the root divergence, depending on the
#' distribution of short branch lengths close to the root.

#  OLD: 'without changing the relative distance of the tips from the root node.'
# 03-13-19 THAT IS AN EVIL LIE DAVID, THIS TOTALLY CHANGES THE ROOT AGE

#' @details
#' This function was formally an internal segment in
#' \code{\link{timePaleoPhy}}, and now is called by \code{timePaleoPhy}
#' instead, allowing users to apply \code{minBranchLength}
#' to trees that already have edge lengths.

#' @param tree A phylogeny with edge lengths of class \code{phylo}.

#' @param mbl The minimum branch length 

#' @param modifyRootAge If \code{TRUE} (the default), the input tree is checked for
#' a root age given as \code{$root.time} and if present it is checked
#' and fixed for any possible movement backwards due to short
#' branches close to the root node.


#' @return
#' A phylogeny with edge lengths of class \code{phylo}.

#' @seealso
#' This function was originally an internal
#' piece of \code{\link{timePaleoPhy}},
#' which implements the minimum branch
#' length time-scaling method along with others,
#' which may be what you're looking for
#' (instead of this miscellaneous function).

#' @author 
#' David W. Bapst

#' @examples
#' 
#' #simulation with an example non-ultrametric tree
#' 
#' tree <- rtree(20)
#' # randomly replace edges with ZLBs
#'    # similar to multi2di output
#' tree <- degradeTree(tree,0.3,
#'    leave.zlb = TRUE) 	
#' 
#' tree2 <- minBranchLength(tree,0.1)
#' 
#' layout(1:2)
#' 
#' plot(tree)
#' axisPhylo()
#' plot(tree2)
#' axisPhylo()
#' 
#' layout(1)
#' 
#' 
#' #now let's try it with an ultrametric case
#' 
#' # get a random tree
#' tree <- rtree(30)
#' # randomly replace edges with ZLBs
#'    # similar to multi2di output
#' tree <- degradeTree(tree,0.5,leave.zlb = TRUE) 
#' # now randomly resolve	
#' tree <- di2multi(tree)
#' # give branch lengths so its ultrametric
#' tree <- compute.brlen(tree)
#' 
#' # and we have an ultrametric tree with polytomies, yay!
#' plot(tree) 
#' 
#' # now randomly resolve
#' tree2 <- multi2di(tree)
#' # get new branch lengths as would with real data
#' tree2 <- minBranchLength(tree2,0.1)
#' 
#' layout(1:2)
#' plot(tree,show.tip.label = FALSE)
#' axisPhylo()
#' plot(tree2,show.tip.label = FALSE)
#' axisPhylo()
#' 
#' layout(1)
#' 
#' # check that root ages aren't being left unmodified
#'    # create a tree with lots of ZBLs at the root
#' x <- stree(10)
#' x$edge.length <- runif(Nedge(x))
#' x <- multi2di(x)
#' # give it a root age
#' x$root.time <- max(node.depth.edgelength(x))
#' 
#' z <- minBranchLength(tree = x, mbl = 1)
#' plot(z)
#' 

#' @name minBranchLength 
#' @aliases minBranchLen minimumBranchLen minimumBranchLength 
#' @rdname minBranchLength
#' @export
minBranchLength <- function(tree, mbl, modifyRootAge = TRUE){	
	#require(phangorn)
	#test arguments
	#tree - a tree with edge lengths
	if(!inherits(tree, "phylo")){
		stop("tree is not of class phylo")
		}
	if(is.null(tree$edge.length)){stop("Tree has no edge lengths")}
	timetree <- tree
	#mbl - a single numeric value
	if(!is.numeric(mbl) | length(mbl) != 1){
		stop("mbl is a not a single numeric value")}
	#
	#root_node <- Ntip(timetree)+1
	while(any(timetree$edge.length<mbl)){
		# get every edge that is short
		shortEdge <- (1:Nedge(timetree))[timetree$edge.length<mbl]
		# pick the shortest one, if multiple of that length, pick first one
		shortLength <- timetree$edge.length[shortEdge]
		shortestLength <- shortEdge[shortLength == min(shortLength)]
		mom <- timetree$edge[shortestLength[1],1]
		#make vector of every mom node that is ancestral
		mom <- c(mom,Ancestors(timetree,mom))
		selNodes <- (mom[1] == timetree$edge[,1])
		debt <- mbl - min(timetree$edge.length[selNodes])
		timetree$edge.length[selNodes] <- timetree$edge.length[selNodes] + debt[1]
		#make vector of smallest brlen with each mom node as anc
		# calculate, simultaneously, the changes in debt
			# and branch lengthening required as go down tree
		# change branch lengths; hypothetically, debt should then equal zero...
		if(length(mom)>1){
			for(i in 2:length(mom)){
				selNodes_mom <- timetree$edge[,1] == mom[i]
				selNodes_child <- timetree$edge[,2] == mom[i-1]
				selNodes_MC <- selNodes_mom & selNodes_child
				selNodes_MnotC <- selNodes_mom & !selNodes_child
				#
				small <- min(timetree$edge.length[selNodes_mom])
				mom_blen <- timetree$edge.length[selNodes_MC]
				#
				maxMBLdiff <- max(mom_blen-mbl, 0)
				maxMBLsmall <- max(mbl-small, 0) 
				#
				debt[i] <- max(debt[i-1] - maxMBLdiff, 0) + maxMBLsmall
				#
				smallDiff <- max(0, min(maxMBLdiff, debt[i-1]))
				#
				timetree$edge.length[selNodes_MC] <- mom_blen + maxMBLsmall - smallDiff
				# was i on drugs when i wrote that??
				#
				timetree$edge.length[selNodes_MnotC] <-  timetree$edge.length[selNodes_MnotC] + debt[i]
				}
			}
		}
	############################################
	# Fix the root age, if present... 
		# if that has been pushed back further
	if(!is.null(timetree$root.time) & modifyRootAge){
		timetree <- fixRootTime(
			treeNew = timetree,
			treeOrig = tree,
			fixingMethod = "rescaleUsingTipToRootDist",
			testConsistentDepth = TRUE)
		}
	#############################
	return(timetree)
	}

	