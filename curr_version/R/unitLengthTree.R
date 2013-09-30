unitLengthTree<-function(tree){
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	tree$edge.length<-rep(1,Nedge(tree))
	tree$root.time<-NULL
	return(tree)
	}