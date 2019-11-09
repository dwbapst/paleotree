#' Time-Slicing a Phylogeny
#' 
#' Removes the portion of a tree after a set point in time, as if the tree
#' after that moment had been sliced away.

#' @details 
#' The function assumes that the input \code{ttree} will generally have an element called
#' \code{$root.time}, which is the time before present that the root divergence
#' occurred. If \code{$root.time} is not present as an element of \code{ttree}, then it is
#' assumed the tip furthest from the root is at time zero (present-day) and a new
#' \code{$root.time} is calculated (a warning will be issued in this case).
#' 
#' The \code{sliceTime} is always calculated as on the same scale as \code{ttree$root.time}.
#' In other words, if \code{root.time = 100}, then \code{timeSlice = 80} will slice the tree 20
#' time units after the root.
#' 
#' If \code{drop.extinct = TRUE}, then extinct tips are dropped and (if present) the
#' $root.time of ttree is adjusted. This is done using the \code{paleotree} function
#' \code{dropExtinct}.
#' 

#' @param ttree A time-scaled phylogeny of class \code{phylo}.

#' @param sliceTime Time to 'slice' the tree at. See details.

#' @param drop.extinct If \code{TRUE}, drops tips that go extinct before timeSlice.

#' @param tipLabels What sort of tip labels should be placed on cropped branches
#' which had multiple descendants? The default option, \code{"earliestDesc"} labels
#' a clipped branch with the earliest appearing tip descendant of that branch. 
#' Alternatively, if \code{tipLabels = "allDesc"},
#' these tips can instead be labeled with a compound label consisting of
#' all descendants that were on the cropped branch, seperated by semi-colons.

#' @param plot If \code{TRUE}, plots input and output trees for comparison.

#' @return Returns the modified phylogeny as an object of class \code{phylo}. 
#' See argument \code{tipLabels} for how the labeling of the tips for 
#' cut branches is controlled.

#' @author David W. Bapst, with modification of code by Klaus Schliep to avoid use of
#' function \code{dist.nodes}, which has difficulty with large trees, and greatly
#' benefiting the run time of this function.

#' @seealso \code{\link{phyloDiv}}, \code{\link{dropExtinct}},
#' \code{\link{dropExtant}}
#' 
#' Also see the function \code{treeSlice} in the library \code{phytools}, which will slice a
#' tree at some point in and return all the subtrees which remain after the
#' slicing time. (Effectively the reversed \emph{opposite} of \code{timeSliceTree}.)

#' @note
#' Note that the default behavior of \code{tiplabels = "earliestDesc"} labels
#' cut branches with the tip label for the earliest tip descendant.
#' This is somewhat arbitrary; the actual morphotaxon present at that time might have
#' been a different taxon that the earliest appearing tip. For simulated datasets where
#' morphotaxon identity is known throughout and not limited to tip observations,
#' slice the taxon data in that more detailed form, and then transform that morphotaxon
#' data to a tree, perhaps with \code{taxa2phylo}.

#' @examples
#' 
#' # a neat example of using phyloDiv with timeSliceTree 
#'     # to simulate doing extant-only phylogeny studies 
#'     # of diversification...in the past!
#' set.seed(444)
#' record <- simFossilRecord(
#'     p = 0.1, q = 0.1, nruns = 1,
#'     nTotalTaxa = c(30,40), 
#'     nExtant = 0)
#' taxa <- fossilRecord2fossilTaxa(record)
#' taxicDivCont(taxa)
#' 
#' # that's the whole diversity curve
#'    # now let's do it for a particular time-slide
#' tree <- taxa2phylo(taxa)
#' # use timeSliceTree to make tree of relationships
#'     # up until time = 950 
#' tree950 <- timeSliceTree(
#'     tree,
#'     sliceTime = 950,
#'     plot = TRUE,
#'     drop.extinct = FALSE
#'     )
#'
#' # compare tip labels when we use tipLabels = "allDesc"
#' tree950_AD <- timeSliceTree(
#'     tree,
#'     sliceTime = 950,
#'     plot = TRUE,
#'     tipLabel = "allDesc",
#'     drop.extinct = FALSE
#'     )
#'     
#' cbind(tree950$tip.label, tree950_AD$tip.label)
#' 
#' # with timeSliceTree we could
#'     # look at the lineage accumulation curve 
#'     # we would recover from the species extant
#'     # at that point in time
#' 
#' # use drop.extinct = T to only get the
#'     # tree of lineages extant at time = 950
#' tree950 <- timeSliceTree(
#'     tree,
#'     sliceTime = 950,
#'     plot = FALSE,
#'     drop.extinct = TRUE
#'     )
#' # now its an ultrametric tree with many fewer tips...
#'     # lets plot the lineage accumulation plot on a log scale
#' phyloDiv(tree950,
#'     plotLogRich = TRUE
#'     )
#' 

#' @export timeSliceTree
timeSliceTree <- function(ttree,
                          sliceTime,
                          drop.extinct = FALSE,
                          tiplabels = "earliestDesc",
                          plot = TRUE
                          ){
	#take a phylogeny and produce a phylogenetic 'slice' at time X (w/respect to root.time)
		#lineages extant at sliceTime sliced to end at that point
		#if no root.time, then it is assumed the tip furthest from the root is at 0 (present-day)
			#a warning will be issued if this is assumed
		#extinct tips will be dropped if drop.extinct = TRUE
	#require(ape)
	if(!inherits(ttree, "phylo")){
		stop("ttree is not of class phylo")
		}
	if(is.null(ttree$root.time)){
		ttree$root.time <- max(node.depth.edgelength(ttree)[1:Ntip(ttree)])
		message(
		  "Warning: no ttree$root.time! Assuming latest tip is at present (time = 0)"
		  )
		}
  #time from root to slice time
  tslice <- ttree$root.time-sliceTime	
	#first let's drop all edges that branch later than the slice
	#make them all single lineages by dropping all but one taxon
	dnode <- node.depth.edgelength(ttree)
	#identify the ancestor nodes of edges which cross the tslice
	cedge <- which(
	  (dnode[ ttree$edge[, 1] ] < tslice) & (dnode[ttree$edge[, 2] ]   >=  tslice)
	  )
	# test that slice time isn't past the latest tip
	if(length(cedge) == 0){
		stop(
		  "sliceTime is later than latest tip"
		  )
	  }
	droppers <- numeric()
	propPartTree <- prop.part(ttree)
	for(i in 1:length(cedge)){
		desc <- ttree$edge[cedge[i],2]
		if(desc>Ntip(ttree)){	
		  # if an internal edge that goes past the tslice
		  # drop all but one tip
			desctip <- propPartTree[[desc-Ntip(ttree)]]	
			droppers <- c(droppers,desctip[-1])
			#if(tiplabels == "earliestDesc"){
      #  # I don't think anything needs to be done!
			#  }
			if(tiplabels == "allDesc"){
			  ttree$tip.label[desctip[1]] <- paste0(
			    ttree$tip.label[desctip],
			    collapse = ";"
			    )
			  }
		}}
	stree <- drop.tip(ttree,droppers)
	# which edges cross over tslice?
	dnode <- node.depth.edgelength(stree)
	cedge <- (dnode[stree$edge[, 2] ]   >=  tslice)
	cnode_depth <- dnode[stree$edge[cedge,1]]
	stree$edge.length[cedge] <- tslice-cnode_depth
	# root.time shouldn't (can't?) change!
	stree$root.time <- ttree$root.time		
	if(drop.extinct){
		stree1 <- dropExtinct(stree,ignore.root.time = TRUE)
	}else{
		stree1 <- stree
		}
	if(plot){
		layout(1:2)
		plot(ladderize(ttree),show.tip.label = FALSE)
			axisPhylo()
		plot(ladderize(stree1),show.tip.label = FALSE)
			axisPhylo()
		layout(1)
		}
	return(stree1)
	}
