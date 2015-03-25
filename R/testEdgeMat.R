#' Test the edge matrix of a phylo object for inconsistencies
#'
#' This is a small simple function which tests the $edge matrix of 'phylo' objects for
#' inconsistencies that can cause downstream analytical problems.

#' @details
#' Useful when doing complex manipulations and reconstitutions of 'phylo' objects (or their
#' de novo construction), and thus is used by a number of paleotree functions.

#' @param tree A phylogeny object of type phylo

#' @return
#' If all the checks in the function pass correctly, the logical TRUE is returned.

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


#' @name testEdgeMat
#' @rdname testEdgeMat
#' @export testEdgeMat
testEdgeMat<-function(tree){
	if(!is(tree,"phylo")){stop("tree is not of type 'phylo'")}
	#test edge matrix
	if(length(tree$edge[,2])!=length(unique(tree$edge[,2]))){
		stop("Some nodes are listed as a descendant twice in the edge matrix")}
	if(Nnode(tree)!=(max(tree$edge[,1])-Ntip(tree))){stop("Number of nodes is incorrect based on edge[,1]?")}
	if(Ntip(tree)>2){
		if(Nnode(tree)!=(max(tree$edge[,2])-Ntip(tree))){stop("Number of nodes is incorrect based on edge[,2]?")}
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