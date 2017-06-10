



# function for scaling posterior trees to absolute time (get root.time)
#


setRootAge<-function(tree,fixedAges=NULL{
	#For all trees to be comparable, we will use the $root.time convention from paleotree (Bapst, 2012)
	if(!is.phylo()){
		stop("tree must be of type 'phylo'")
	if(is.null(fixedAges)){
		if(is.null(attr(tree,"fixedTable")){
			stop("fixedAges must be supplied")
		}else{
			fixedAges<-attr(tree,"fixedTable")
			}
		}
	fixedTaxa<-as.character(fixedAges[,1])
	fixedAges<-as.numeric(fixedAges[,2])
	taxaTree<-tree$tip.labels
	# drop unshared taxa
	missingAge<-sapply(fixedTaxa,function(x) all(x!=taxaTree))
	if(length(missingAge)==length(fixedTaxa)){
		stop("None of the taxa in fixedAges found as OTU tip labels on tree")
		}	
	fixedAges<-fixedAges[!missingAge]
	fixedTaxa<-fixedTaxa[!missingAge]
	# get first taxon at youngest age
	youngest<-fixedAges==min(fixedAges)[1]
	youngDate<-fixedAges[youngest]
	youngTipDepth<-node.depth.edgelength(tree)[taxaTree==treefixedTaxa[youngest]]
	tree$root.time<-youngTipDepth+youngDate	
	return(tree)
	}

'