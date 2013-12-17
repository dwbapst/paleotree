simTermTaxa<-function(ntaxa,sumRate=0.2){
		#previously known as 'candle'
	#This function will make ideal cladist datasets
	#all taxa will be descendants
		#all evolution prior to branching or differentiation is unsampled
		#taxa do not differentiatiate over those times
	#require(ape)
	tree<-deadTree(ntaxa=ntaxa,sumRate=sumRate)
	#taxonNames<-tree$tip.label
	termEdge<-sapply(tree$edge[,2],function(x) any(x==(1:ntaxa)))
	#termAnc<-tree$edge[termEdge,1]
	taxonDurations<-tree$edge.length[termEdge]
	nodeDist<-dist.nodes(tree)
	taxonLADs<-tree$root.time-nodeDist[ntaxa+1,1:ntaxa]
	taxonFADs<-taxonLADs+taxonDurations
	taxonRanges<-cbind(taxonFADs,taxonLADs)
	rownames(taxonRanges)<-tree$tip.label[tree$edge[termEdge,2]]
	res<-list(taxonRanges=taxonRanges,tree=tree)
	return(res)
	}