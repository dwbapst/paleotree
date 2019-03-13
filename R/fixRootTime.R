#' @rdname modifyTerminalBranches
#' @export
fixRootTime <- function(treeOrig,treeNew,consistentDepth = TRUE,nodeAgeTransfer = TRUE){
	# UNUSED FUNCTION?
	#treeDepth <- function(tree){
	#	#require(ape)
	#	max(dist.nodes(tree)[,Ntip(tree)+1])
	#	}
	#
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
		stop(paste("taxa:",taxaNewNM,"are present in treeNew but not treeOrig"))}
	#two different ways to fix the root time
	if(nodeAgeTransfer){
		##NEW WAY 11-28-14
			#If TRUE, the root.time of the new tree is determined by
			#comparing the clades of taxa between the two input trees. The new root age assigned is the age of
			#(\emph{1}) the treeOrig clade that contains *all* taxa present in treeNew and, if the set of (1)
			#contains multiple clades, (\emph{2}) the clade in the (1) set that contains the fewest taxa not in
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
		if(round(max(node.depth.edgelength(treeNew)) - treeNew$root.time)>0){
			stop("fixRootTime isn't fixing correctly, root.time less than max tip-to-root length!")}
		}
	return(treeNew)
	}