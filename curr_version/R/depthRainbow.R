depthRainbow<-function(tree){
	#plots a tree with edges color-coded to depth
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	tree<-ladderize(tree)
	ndepth<-dist.nodes(tree)[,Ntip(tree)+1]
	#nodelabels(ceiling(ndepth[(Ntip(tree):Nedge(tree))+1]),node=(Ntip(tree):Nedge(tree))+1)
	edepth<-ceiling((ndepth[(Ntip(tree):Nedge(tree))+1])[tree$edge[,1]-Ntip(tree)])+1
	col_edge<-rainbow(max(edepth))[edepth]
	plot(ladderize(tree),show.tip.label=FALSE,edge.color=col_edge);axisPhylo()
	}