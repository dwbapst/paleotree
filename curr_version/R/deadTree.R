deadTree<-function(ntaxa,sumRate=0.2){
	#all taxa will always be extinct
		#assuming the clade is extinct, p and q are meaningless
		#all trees would have us assume that p=q if analyzed...
		#we'll just have a 'rate' of events (branching or ext; p+q)
		#the probability of a branch ending in branching or ext is meaningless!
			#by definition equal proportion must branch and go extinct
			#anything else should be indep in a b-d process
	#lets make the phylogeny 
	#require(ape)
	tree<-rtree(ntaxa)
	tree$edge.length<-rexp(ntaxa+ntaxa-2,sumRate)
	tree$root.time<-max(dist.nodes(tree)[ntaxa+1,])+200
	return(tree)
	}