#' Modify, Drop or Bind Terminal Branches of Various Types (Mainly for Paleontological Phylogenies)
#' 
#' Modifying a dated tree with \code{$root.time} elements often
#' changes the actual timing of the root relative to the 
#' tips, such as when dropping tips, extending branches, or
#' shift node ages backwards. When such modifications occur, 
#' the function \code{fixRootTime} can be used to find
#' the correct root age. 

#' This function is mainly used as a
#' utility function called by other tree-modifying functions discussed
#' in the manual page for \code{\link{modifyTerminalBranches}}.

#' This is typically perfomed via the function \code{\link{fixRootTime}}.

#' @param treeOrig A \code{phylo} object of a time-scaled
#' phylogeny with a \code{$root.time} element.

#' @param treeNew A \code{phylo} object containing a
#' modified form of \code{treeOrig} (with
#' no extra tip taxa added, but possibly with
#' some tip taxa removed).

#' @param consistentDepth A logical, either \code{TRUE}
#' or \code{FALSE}. If \code{TRUE} (the default)
#' the tree's root-to-furthest-tip depth is tested
#' to make sure this depth is
#' not greater than the new \code{$root.time}
#' appended to the output tree.

#' @param nodeAgeTransfer A logical. If \code{TRUE}, the
#' \code{$root.time} of the new tree is determined by
#' comparing the clades of taxa between the two input
#' trees. The new root age assigned is the age of
#' (\emph{1}) the \code{treeOrig} clade that contains \emph{all}
#' taxa present in \code{treeNew} and, if the set of (1)
#' contains multiple clades, (\emph{2}) the clade in the
#' first set that contains the fewest taxa not in
#' \code{treeNew}. If \code{FALSE}, the \code{root.time} assigned to
#' \code{treeNew} is the \code{$root.time} of \code{treeOrig}, adjusted
#' based on the change in total tree depth between \code{treeOrig} and
#' \code{treeNew}, as measured between the root and
#' the first matching taxon in both trees. The later is 
#' how \code{fixRootTime} functioned by default
#' prior to paleotree v2.3.


#' @return Gives back a modified phylogeny as a \code{phylo} object, with a
#' modified \code{$root.time} element.

#' @author David W. Bapst

#' @seealso 
#' \code{\link{modifyTerminalBranches}},
#' \code{\link{minBranchLength}} 

#' @examples
#' 
#' \donttest{
#' 
#' #testing dropPaleoTip... and fixRootTime by extension
#'
#' #simple example
#' tree <- read.tree(text = "(A:3,(B:2,(C:5,D:3):2):3);")
#' tree$root.time <- 10
#' plot(tree,no.margin = FALSE)
#' axisPhylo()
#'
#' # now a series of tests, dropping various tips
#' (test <- dropPaleoTip(tree,"A")$root.time) #  = 7
#' (test[2] <- dropPaleoTip(tree,"B")$root.time) #  = 10
#' (test[3] <- dropPaleoTip(tree,"C")$root.time) #  = 10
#' (test[4] <- dropPaleoTip(tree,"D")$root.time) #  = 10
#' (test[5] <- dropPaleoTip(tree,c("A","B"))$root.time) #  = 5
#' (test[6] <- dropPaleoTip(tree,c("B","C"))$root.time) #  = 10
#' (test[7] <- dropPaleoTip(tree,c("A","C"))$root.time) #  = 7
#' (test[8] <- dropPaleoTip(tree,c("A","D"))$root.time) #  = 7
#' 
#' # is it all good? if not, fail so paleotree fails...
#' if(!identical(test,c(7,10,10,10,5,10,7,7))){stop("fixRootTime fails!")}
#' 
#' }

#' @rdname fixRootTime
#' @export
fixRootTime <- function(
		treeOrig, treeNew,
		consistentDepth = TRUE,
		nodeAgeTransfer = TRUE){
	#####################################
	#require(ape)
	if(!inherits(treeOrig, "phylo")){
		stop("treeOrig is not of class phylo")}
	if(!inherits(treeNew, "phylo")){
		stop("treeNew is not of class phylo")}
	if(is.null(treeOrig$root.time)){
		stop("ERROR: treeOrig passed to fixRootTime with no $root.time??")}
	#also need a warning message if taxa present in treeNew that aren't in treeOrig
	taxaNewNM <- treeNew$tip.label[sapply(treeNew$tip.label,function(x) !any(x == treeOrig$tip.label))]
	if(length(taxaNewNM)>0){
		stop(paste("taxa:",taxaNewNM,
			"are present in treeNew but not treeOrig"))
			}
	#two different ways to fix the root time
	if(nodeAgeTransfer){
		##NEW WAY 11-28-14
			#If TRUE, the root.time of the new tree is determined by
			#comparing the clades of taxa between the two input trees.
			#The new root age assigned is the age of
			#(\emph{1}) the treeOrig clade that contains *all* taxa
			#present in treeNew and, if the set of (1)
			#contains multiple clades, (\emph{2}) the clade in the
			#(1) set that contains the fewest taxa not in
			#treeNew.
		dates <- dateNodes(treeOrig,labelDates = FALSE,tolerance = 0.001)
		treeDesc <- lapply(Descendants(treeOrig),function(x) sort(treeOrig$tip.label[x]))
		#treeRootNew <- sort(treeNew$tip.label[Descendants(treeNew)[[Ntip(treeNew)+1]]]) #no
		#the descendants of treeNew's root are ALL the taxa in treeNEW
		#So which treeOrig clades contain ALL taxa in treeNew?
		allNewTaxa <- sapply(treeDesc,function(x) all(sapply(treeNew$tip.label,function(y) any(y == x)))) #logical vector
		#now, if more than one contains all-new-taxa, which of these treeOrig clades minimizes not-shared taxa?
		if(sum(allNewTaxa)>1){
			nUnshared <- sapply(treeDesc,function(x) sum(sapply(x,function(y) all(y != treeNew$tip.label)))) #numeric
			matchRootNew <- which(allNewTaxa & nUnshared == min(nUnshared[allNewTaxa]))
		}else{
			matchRootNew <- which(allNewTaxa)
			}
		if(length(matchRootNew)>1){stop("More than one node contains these taxa")} #maybe sort by age
		if(length(matchRootNew)<1){stop("No nodes match the new tree's root, a root age can not be obtained")}
		treeNew$root.time <- unname(dates[matchRootNew])
	}else{
		##OLD WAY
			#If FALSE, the root.time assigned to treeNew is the root.time of treeOrig, adjusted
			# based on the change in total tree depth between treeOrig and treeNew,
				# as measured between the root and the first matching taxon in both trees.
			# The later is how fixRootTime functioned by default prior to paleotree v2.3.
		orig_dist <- node.depth.edgelength(treeOrig)[
			which(treeNew$tip.label[1] == treeOrig$tip.label)
			]
		new_dist <- node.depth.edgelength(treeNew)[1]
		treeNew$root.time <- treeOrig$root.time-(orig_dist-new_dist)
		}
	if(consistentDepth){
		#if(round(max(node.depth.edgelength(treeNew)) - treeNew$root.time)>0){
		#	stop("fixRootTime isn't fixing correctly, root.time less than max tip-to-root length!")}
		#
		checkRootTime(treeNew, stopIfFail)
		}
	return(treeNew)
	}
	
	

# not exported
checkRootTime <- function(tree, stopIfFail = FALSE){	
	# check that the tree and its root age makes sense
	if(!is.null(tree$root.time) & !is.null(tree$edge.length)){
		if(round(max(node.depth.edgelength(tree)) - tree$root.time)>0){
			#tree$root.time < max(node.depth.edgelength(tree))
			if(stopIfFail){
				stop(paste0("Max tip-to-root tree depth (",
					round(max(node.depth.edgelength(tree))),
					") is greater than root age(",
					round(tree$root.time),
					"), so tips are in the future",
					"\nSomething has probably gone very wrong"))
				}else{
				warning(paste0("Max tip-to-root tree depth (",
					round(max(node.depth.edgelength(tree))),
					") is greater than root age(",
					round(tree$root.time),
					"), so tips are in the future",
					"\nSomething has probably gone very wrong"))
				}
			res <- FALSE
		}else{
			res <- TRUE
			}
		}
	}			