trueTermTaxaTree<-function(TermTaxaRes,time.obs){
		#for a term-taxa (candle) tree datasets
	#given a per-taxon vector of observation times, returns the true tree
	#time.obs must have taxon names
	#require(ape)
	#first, check if the times of observations are outside of original taxon ranges
	taxR<-TermTaxaRes$taxonRanges
	nameMatch<-match(names(time.obs),rownames(taxR))
	if(any(is.na(nameMatch))){stop("ERROR: names on time.obs and in TermTaxaRes don't match")}
	if(any(sapply(1:length(time.obs),function(x) if(is.na(time.obs[x])){FALSE}else{
			(time.obs[x]>taxR[nameMatch[x],1])|(time.obs[x]<taxR[nameMatch[x],2])}))){
		stop("ERROR: Given time.obs are outside of the original taxon ranges")}
	#now onwards with the actual function
	tree1<-TermTaxaRes$tree
	termEdge<-sapply(tree1$edge[,2],function(x) any(x==(1:Ntip(tree1))))	
	newDurations<-taxR[nameMatch,1]-time.obs
	if(is.null(names(time.obs))){stop("ERROR: No taxon names on observation vector?")}
	tipMatch<-sapply(1:Ntip(tree1),function(x) which(tree1$tip.label[x]==names(time.obs)))
	dropTaxa<-character()
	for(i in 1:Ntip(tree1)){
		newDur<-newDurations[tipMatch[i]]
		if(!is.na(newDur)){
			if(tree1$edge.length[tree1$edge[,2]==i]<newDur){stop("Error: New duration longer than original taxon ranges?")}
			tree1$edge.length[tree1$edge[,2]==i]<-newDur
		}else{
			dropTaxa<-c(dropTaxa,tree1$tip.label[i])
			}
		}
	treeNoDrop<-tree1		#the root.time can't change cause the term branches got shifted
	if(length(dropTaxa)>0){tree1<-drop.tip(tree1,dropTaxa)}
	#need to correct root.time if basal outgroups were removed
	tree1<-fixRootTime(treeNoDrop,tree1)	#use the tree after adjusting the term branch lengths
	return(tree1)
	}