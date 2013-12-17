cladogeneticTraitCont<-function(taxa,rate=1,meanChange=0,rootTrait=0){
	#simulate speciational trait evolution for datasets from simFossilTaxa
	#idiot proofing
	taxa<-taxa[order(taxa[,1]),]
	if(any(taxa[-1,2]>=taxa[-1,1])){stop("Ancestors have higher IDs than Descendants?")}
	anctaxa<-sapply(taxa[-1,2],function(x) which(x==taxa[,1]))
	traits<-rootTrait
	for(i in 2:nrow(taxa)){
		if(taxa[i,1]==taxa[i,6]){
			traits[i]<-traits[anctaxa[i-1]]+rnorm(1,mean=meanChange,sd=sqrt(rate))
		}else{
			traits[i]<-traits[anctaxa[i-1]]
			}
		}
	names(traits)<-rownames(taxa)
	return(traits)
	}