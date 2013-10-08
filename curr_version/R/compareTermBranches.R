compareTermBranches<-function(tree1,tree2){
	#output vector of shifts in terminal branch lengths
	#require(ape)
	if(!is(tree1, "phylo")){stop("Error: tree1 is not of class phylo")}
	if(!is(tree2, "phylo")){stop("Error: tree2 is not of class phylo")}
	tree1<-drop.tip(tree1,tree1$tip.label[is.na(match(tree1$tip.label,tree2$tip.label))])
	tree2<-drop.tip(tree2,tree2$tip.label[is.na(match(tree2$tip.label,tree1$tip.label))])
	term1<-tree1$edge.length[tree1$edge[,2]<=Ntip(tree1)]
	term1<-term1[order(tree1$edge[tree1$edge[,2]<=Ntip(tree1),2])]
	term2<-tree2$edge.length[tree2$edge[,2]<=Ntip(tree2)]
	term2<-term2[order(tree2$edge[tree2$edge[,2]<=Ntip(tree2),2])]
	term2<-term2[match(tree1$tip.label,tree2$tip.label)]
	term_diff<-term2-term1
	names(term_diff)<-tree1$tip.label
	return(term_diff)
	}