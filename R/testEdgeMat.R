#' Test the Edge Matrix of a 'phylo' Phylogeny Object for Inconsistencies
#'
#' \code{testEdgeMat} is a small simple function which tests the $edge matrix of 'phylo' objects for
#' inconsistencies that can cause downstream analytical problems.
#' The associated function, \code{cleanNewPhylo} puts an input
#' phylo object, presumably freshly created or reconstituted by some function, through a series
#' of post-processing, This includes having singles collapsed,
#' nodes reordered and being written out as a Newick string and read back in,
#' to ensure functionality with ape functions
#' and ape-derived functions. 

#' @aliases cleanNewPhylo cleanTree

#' @details
#' Useful when doing complex manipulations and reconstitutions of 'phylo' objects (or their
#' de novo construction), and thus is used by a number of paleotree functions.

#' @param tree A phylogeny object of type phylo

#' @param reorderTree A logical indicating whether a step of \code{reorder.phylo()} will be applied.
#' Reordering may cause more problems than it is worth.

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
#' # should also work on star trees
#' testEdgeMat(stree(10))
#'
#' # should also work on trees with two taxa
#' testEdgeMat(rtree(2))
#'
#' # should also work on trees with one taxon
#' testEdgeMat(stree(1))
#'
#' #running cleanNewPhylo on this tree should have little effect
#' 		#beyond ladderizing it...
#' tree1<-cleanNewPhylo(tree)
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
	#test that all but one node has an ancestor
	parentMatch<-match(unique(tree$edge[,1]),tree$edge[,2])
	if(sum(is.na(parentMatch))>1){
		stop(paste("More than one apparent root; \n",
			"more than one internal node without an ancestor listed"))}
	#trace all tips to a single ancestor
	ultimateAnc<-sapply(unique(c(tree$edge[,1],tree$edge[,2])),function(taxa){
		while(any(tree$edge[,2]==taxa)){
			taxa<-tree$edge[tree$edge[,2]==taxa,1]
			if(length(taxa)>1){
				stop("Some nodes are listed as a descendant twice in the edge matrix")}
			}
		return(taxa)
		})
	if(length(unique(ultimateAnc))!=1){
		stop("Tip and internal node IDs in $edge trace back to more than one unique common ancestor")}
	if(Nnode(tree)!=(max(tree$edge[,1])-Ntip(tree))){
		stop("Nnode is lower than number implied by edge[,1]?")}
	if(Ntip(tree)>2){
		if(Nnode(tree)!=(max(tree$edge)-Ntip(tree))){
			stop("Number of nodes is incorrect based on edge numbering?")}
		if(Nnode(tree)>1){
			#is every internal node listed as a descendant and ancestor, in edge[,2] and edge[,1]?
			if(!all(sapply((1:Nnode(tree))+Ntip(tree),function(x) any(x==tree$edge[,1])))){
				stop("Not all internal nodes (including root) listed in edge[,1]?")}
			if(sum(!sapply((1:Nnode(tree))+Ntip(tree),function(x) any(x==tree$edge[,2])))>1){
				stop("Not all internal nodes (except root) listed in edge[,2]?")}
			}
		}
	#if(identical(sort(unique(tree$edge[,2])),c(1L,2L))){stop("Number of nodes is incorrect based on edge[,2]?")}
	return(TRUE)
	}

#' @rdname testEdgeMat
#' @export 
cleanNewPhylo<-function(tree,reorderTree=TRUE){ 
		#CHECKS
		if(class(tree)!="phylo"){stop("Must be class 'phylo'")}
		if(any(is.na(match(c("edge","tip.label","Nnode"),names(tree))))){
			stop("Missing key required elements of a 'phylo' object")}
		oldNtip<-length(tree$tip.label)
		#make it a good tree
		#check it
		if(!testEdgeMat(tree)){stop("Edge matrix has inconsistencies")}
		#collapse singles
			#count number of single nodes
		Nsingle<-sum(sapply(unique(tree$edge[,1]),function(x) sum(x==tree$edge[,1])==1))
		if(Nsingle>0){
			treePrev<-tree
			while(Nsingle>0){
				tree1<-collapse.singles(treePrev)
				if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
				if((Nnode(treePrev)-Nnode(tree1))>Nsingle){
					stop("collapse.singles dropped too many nodes")}
				Nsingle<-sum(sapply(unique(tree1$edge[,1]),function(x) sum(x==tree1$edge[,1])==1))
				treePrev<-tree1
				#print("for counting how many times singles need to be dropped")
				}
		}else{
			tree1<-tree
			}
		#
		if(reorderTree){
			#reorder
			tree1<-reorder(tree1,"cladewise") 	#REORDER IT
			if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
			}
		tree1<-read.tree(text=write.tree(tree1))
		if(!testEdgeMat(tree1)){stop("Edge matrix has inconsistencies")}
		tree1<-ladderize(tree1)
		if(oldNtip!=Ntip(tree1)){stop("Final tip taxon number different from original number of tip taxon names")}
		return(tree1)
		}