dropZLB<-function(tree){
	#drops terminal branches that are zero length
		#adjusts tree$root.time if necessary
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	drop_e<-(tree$edge[,2]<(Ntip(tree)+1)) & (tree$edge.length==0)
	drop_t<-(tree$edge[,2])[drop_e]
	if((Ntip(tree)-length(drop_t))>1){
		tree1<-drop.tip(tree,drop_t)
		if(!is.null(tree$root.time)){tree1<-fixRootTime(tree,tree1)}
		res<-tree1
	}else{res<-NA}
	return(res)
	}