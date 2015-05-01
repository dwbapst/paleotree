#' Test the Edge Matrix of a 'phylo' Phylogeny Object for Inconsistencies
#'
#' \code{testEdgeMat} is a small simple function which tests the $edge matrix of 'phylo' objects for
#' inconsistencies that can cause downstream analytical problems.
#' The associated function, \code{cleanNewPhylo} puts an input
#' phylo object, presumably freshly created or reconstituted by some function, through a series
#" of post-processing, This includes having singles collapsed,
#' nodes reordered and being written out as a Newick string and read back in,
#' to ensure functionality with ape functions
#' and ape-derived functions. 

#' @aliases cleanNewPhylo cleanTree

#' @details
#' Useful when doing complex manipulations and reconstitutions of 'phylo' objects (or their
#' de novo construction), and thus is used by a number of paleotree functions.

#' @param tree A phylogeny object of type phylo

#' @return
#' For \code{testEdgeMat}, if all the checks in the function pass correctly, the logical TRUE is returned.
#'
#' For \code{cleanNewPhylo}, an object of class 'phylo' is returned.

#' @author David W. Bapst

#' @examples
#'
#' set.seed(444)
#' tree<-rtree(10)
#' # should return TRUE
#' testEdgeMat(tree)
#'
#' # should also work on trees with two taxa
#' testEdgeMat(rtree(2))
#'
#' #running cleanNewPhylo on this tree should have little effect
#' 		#beyond ladderizing it...
#' tree<-cleanNewPhylo(tree)
#'
#' #compare outputs
#' layout(1:2)
#' plot(tree)
#' plot(tree1)
#' layout(1)


#' @name testEdgeMat
#' @rdname testEdgeMat
#' @export testEdgeMat
testEdgeMat<-function(tree){
	if(!is(tree,"phylo")){stop("tree is not of type 'phylo'")}
	#test edge matrix 
	if(length(tree$edge[,2])!=length(unique(tree$edge[,2]))){
		stop(paste("Some nodes are listed as a descendant twice in the edge matrix",
		paste0(tree$edge[duplicated(tree$edge[,2]),2],collapse=", ")))}
	if(Nnode(tree)!=(max(tree$edge[,1])-Ntip(tree))){stop("Number of nodes is incorrect based on edge[,1]?")}
	if(Ntip(tree)>2){
		if(Nnode(tree)!=(max(tree$edge)-Ntip(tree))){
			stop("Number of nodes is incorrect based on edge numbering?")}
		#is every internal node listed as a descendant and ancestor, in edge[,2] and edge[,1]?
		if(!all(sapply((2:Nnode(tree))+Ntip(tree),function(x) any(x==tree$edge[,1])))){
			stop("Not all internal nodes (except root) listed in edge[,1]?")}
		if(!all(sapply((2:Nnode(tree))+Ntip(tree),function(x) any(x==tree$edge[,2])))){
			stop("Not all internal nodes (except root) listed in edge[,2]?")}
	#}else{
	#	if(identical(sort(unique(tree$edge[,2])),c(1L,2L))){stop("Number of nodes is incorrect based on edge[,2]?")}
		}
	return(TRUE)
	}

#' @rdname testEdgeMat
#' @export 
cleanNewPhylo<-function(tree){ 
		#CHECKS
		if(class(tree)!="phylo"){stop("Must be class 'phylo'")}
		if(any(is.na(match(c("edge","tip.label","Nnode","edge.length"),names(tree))))){
			stop("Missing key required elements of a 'phylo' object")}
		oldNtip<-length(tree$tip.label)
		#make it a good tree
		#collapse singles
		tree<-collapse.singles(tree)
		#check it
		if(!testEdgeMat(tree)){stop("Edge matrix has inconsistencies")}
		tree<-reorder(tree,"cladewise") 	#REORDER IT
		tree<-read.tree(text=write.tree(tree))
		tree<-ladderize(tree)
		if(oldNtip!=Ntip(tree)){stop("Final tip taxon number different from original number of tip taxon names")}
		return(tree)
		}