addTermBranchLength<-function(tree,addtime=0.001){
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]<-tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]+addtime
	if(any(tree$edge.length<0)){stop("Error: tree has negative branch lengths!")}
	if(!is.null(tree$root.time)){tree$root.time<-tree$root.time+addtime}
	return(tree)
	}