taxa2cladogram<-function(taxad,drop.cryptic=FALSE,plot=FALSE){
	#take a taxad and turn it into an unscaled cladogram
		#do this by forming the tree as newick format first
		#essentially, this algorithm works by defining clades as only those sets of taxa which start with the FAD of the first taxon
		#This could be a 'clade' of a single OTU
		#or a clade of a bunch of taxa where one is a budding ancestor and the other are budding descendants or whatever
		#the important thing is that there's only so many nodes as there are instances where morphotaxa end/start allowing for homologies
	#require(ape)
	if(any(taxad[,6]!=taxad[,1])){
		for(i in which(taxad[,1]!=taxad[,6])){
			#reset descendants of cryptic taxa so all bud off of first cryptic species		
			taxad[taxad[,2]==taxad[i,1],2]<-taxad[i,6]
			}
		}
	tlabs<-rownames(taxad)
	desc<-lapply(taxad[,1],function(x) (taxad[taxad[,2]==x,1])[!is.na(taxad[taxad[,2]==x,1])])
	ndesc<-sapply(desc,length)
	rank<-numeric(length(ndesc))
	rank[ndesc==0]<-1
	rank[rank==0]<-NA
	while(any(is.na(rank))){
		rank<-sapply(1:length(rank),function(x) ifelse(!is.na(rank[x]),rank[x],
				1+max(rank[sapply(desc[[x]],function(y) which(y==taxad[,1]))])))}
	#okay, now all taxa are ranked by their depth from the tips
	comp<-numeric(length(ndesc))
	lab<-list()
	lab[rank==1]<-tlabs[rank==1]
	comp[rank==1]<-1
	while(any(comp==0)){
		tpot<-comp==0
		tpot2<-rank==min(rank[tpot])
		tpick<-which(tpot & tpot2)[1]
		dlab<-paste(unlist(lab[desc[[tpick]]]),",",sep="",collapse="")
		lab[[tpick]]<-paste("(",dlab,tlabs[tpick],")",sep="")
		comp[tpick]<-1
		}
	tree1<-paste(lab[[1]],";",sep="")
	tree2<-read.tree(text=tree1)
	if(drop.cryptic & any(taxad[,6]!=taxad[,1])){
		tree2<-drop.tip(tree2,tlabs[taxad[,6]!=taxad[,1]])
		tree2<-collapse.singles(tree2)
		}
	if(plot){plot(ladderize(tree2),show.tip.label=FALSE)}
	return(tree2)
	}