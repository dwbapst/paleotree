simTermTaxaAdvanced<-function(p=0.1,q=0.1,mintaxa=1,maxtaxa=1000,mintime=1,maxtime=1000,
		minExtant=0,maxExtant=NULL,min.cond=TRUE){
		#previously known as 'candle'
	#This function will make ideal cladist datasets
	#all taxa will be monophyletic descendants
		#all evolution prior to branching or differentiation is unsampled
		#taxa do not differentiate over those times
	#extant example
		#p=0.1;q=0.1;mintaxa=50;maxtaxa=100;mintime=1;maxtime=1000;minExtant=10;maxExtant=20;min.cond=FALSE
	#require(ape)
	taxa<-simFossilTaxa(p=p,q=q,mintaxa=mintaxa,maxtaxa=maxtaxa,mintime=mintime,maxtime=maxtime,
		minExtant=minExtant,maxExtant=maxExtant,min.cond=min.cond,nruns=1,
		anag.rate=0,prop.bifurc=0,prop.cryptic=0,count.cryptic=FALSE,print.runs=FALSE,plot=FALSE)
	tree<-dropZLB(taxa2phylo(taxa))
	ntaxa<-Ntip(tree)
	#taxonNames<-tree$tip.label
	termEdge<-sapply(tree$edge[,2],function(x) any(x==(1:ntaxa)))
	termNodes<-tree$edge[termEdge,2]
	#termAnc<-tree$edge[termEdge,1]
	taxonDurations<-tree$edge.length[termEdge]
	nodeDist<-dist.nodes(tree)
	taxonLADs<-tree$root.time-nodeDist[ntaxa+1,termNodes]
	taxonFADs<-taxonLADs+taxonDurations
	taxonRanges<-cbind(taxonFADs,taxonLADs)
	rownames(taxonRanges)<-tree$tip.label[termNodes]
	res<-list(taxonRanges=taxonRanges,tree=tree)
	return(res)
	}