fixRootTime<-function(treeOrig,treeNew){
	treeDepth<-function(tree){
		#require(ape)
		max(dist.nodes(tree)[,Ntip(tree)+1])
		}
	#require(ape)
	if(class(treeOrig)!="phylo"){stop("Error: treeOrig is not of class phylo")}
	if(class(treeNew)!="phylo"){stop("Error: treeNew is not of class phylo")}
	if(is.null(treeOrig$root.time)){stop("ERROR: treeOrig passed to fixRootTime with no $root.time??")}
	orig_dist<-dist.nodes(treeOrig)[
		which(treeNew$tip.label[1]==treeOrig$tip.label),Ntip(treeOrig)+1
		]
	new_dist<-dist.nodes(treeNew)[1,Ntip(treeNew)+1]
	treeNew$root.time<-treeOrig$root.time-(orig_dist-new_dist)
	if(round(max(dist.nodes(treeNew)[, Ntip(treeNew) + 1]) - treeNew$root.time)>0){
		stop("Error: fixRootTime isn't fixing correctly, root.time less than max tip-to-root length!")}
	return(treeNew)
	}