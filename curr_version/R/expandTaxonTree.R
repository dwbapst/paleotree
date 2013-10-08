expandTaxonTree<-function(taxonTree,taxaData,collapse=NULL,keepBrLen=FALSE,plot=FALSE){
	#this function takes a higher-level taxon tree and
		#expands it to a lower level species-level tree
		#using a species list
	#"taxa" here represents the groups to be replaced on the taxonTree
	#taxonTree = tree with taxon IDs as tips
	#taxaData = character vector of higher taxon ids for each new tip, tip labels as vector names
	#collapse = if present, vector of taxa names to be collapsed
	#should be possible to take a tree of mixed species/genera
		#and just replace the genera
	#taxonTree<-rtree(10);taxonTree$tip.label<-as.character(1:10);collapse<-sample(taxonTree$tip.label,5)
	#taxaData<-as.character(sample(1:10,100,replace=TRUE));names(taxaData)<-paste("t",1:100,sep="")
	#require(ape)
	if(!is(taxonTree, "phylo")){stop("Error: taxonTree is not of class phylo")}
	if(any(is.na(taxaData))){stop("Error: some values of taxonData missing!!")}
	if(!is.null(collapse) & keepBrLen){
		message("Warning: collapsed branch lengths turned to zero when keepBrLen is TRUE!")}
	tree<-taxonTree
	if(!keepBrLen){tree$edge.length<-rep(1,Nedge(tree))}		#get rid of all branch lengths
	#first, expand all higher taxa to lower taxon polytomies
	for(i in unique(taxaData)){				#loop through all 	
		tip<-which(tree$tip.label==i)
		if(length(collapse)>0){if(any(collapse==i)){
			tree$edge.length[which(tree$edge[,2]==tip)]<-0
			}}
		cotaxa<-names(taxaData)[taxaData==i]	#which species do I want? These...
		repTree<-stree(length(cotaxa))		#replacement polytomy
		if(keepBrLen){repTree$edge.length<-rep(0,length(cotaxa))
			}else{repTree$edge.length<-rep(1,length(cotaxa))}
		repTree$tip.label<-cotaxa			#replace names,edge.lengths
		tree<-bind.tree(tree,repTree,tip)	#replace the right tip	
		}
		#now collapse non-monophyletic groupings
	if(!keepBrLen){
		tree1<-di2multi(tree);tree1$edge.length<-NULL;tree1<-collapse.singles(tree1)
		tree1<-read.tree(text=write.tree(tree1))
		}else{tree1<-tree}
	if(plot==TRUE){layout(1:2);plot(taxonTree);plot(tree1);layout(1)}
	return(tree1)
	}