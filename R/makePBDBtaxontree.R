#' Creating a Taxon-Tree from Taxonomic Data Downloaded from the Paleobiology Database
#'
#' This function creates phylogeny-like object of type phylo from the taxonomic hierarchy
#' recorded in a taxonomy download from the PBDB for a given group.

#' @details
#' This function should not be taken too seriously. Many groups in the Paleobiology Database have
#' out-of-date or very incomplete taxonomic information. This function is meant to help visualize
#' what information is present, and by use of time-scaling functions, allow us to visualize the intersection
#' of temporal and phylogenetic, mainly to look for incongruence due to either incorrect taxonomic placements,
#' erroneous occurrence data or both.
#'
#' Note however that, contrary to common opinion among some paleontologists, taxon-trees may be just as useful for 
#' macroevolutionary studies as reconstructed phylogenies (Soul and Friedman, in press.)

#' @param data A table of taxonomic data collected from the Paleobiology Database, using the taxa list option
#' with show=phylo. 

#' @param rank The selected taxon rank; must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'.

#' @return
#' A phylogeny of class 'phylo', where each tip is a taxon of the given 'rank'. Edges are scaled so that
#' the distance from one taxon rank to another 1, then merged to remove singleton nodes. As not all
#' taxa have parents at the immediate taxon level above, this leads to some odd cases. For example,
#' two genera emanating from a node representing a class but with a very short (length=1) branch
#' and a long branch (length=3) means one genus is simply placed in the class, with no family or order listed
#' wile the one on the long branch is within an order and family that is otherwise monogeneric.

#' @seealso
#' Other functions for manipulating PBDB data can be found at \code{\link{taxonSortPBDBocc}},
#' \code{\link{occData2timeList}}, and the example data at \code{\link{graptPBDB}}.

#' @author David W. Bapst

#' @references
#' Soul, L. C., and M. Friedman. In Press. Taxonomy and Phylogeny Can Yield
#' Comparable Results in Comparative Palaeontological Analyses.\emph{Systematic Biology}
#' (http://sysbio.oxfordjournals.org/content/early/2015/03/23/sysbio.syv015.abstract)

#' @examples
#' \dontrun{
#' 
#' getTaxaPBDB<-function(taxon){
#' 	#let's get some taxonomic data
#' 	taxaData<-read.csv(paste0("http://paleobiodb.org/",
#' 		"data1.1/taxa/list.txt?base_name=",taxon,
#' 		"&rel=all_children&show=phylo,img,nav&status=senior"))
#' 	return(taxaData)
#' 	}
#' 
#' #graptolites
#' data<-getTaxaPBDB("Graptolithina")
#' graptTree<-makePBDBtaxontree(data,"genus")
#' plot(graptTree)
#' 
#' #conodonts
#' data<-getTaxaPBDB("Conodonta")
#' conoTree<-makePBDBtaxontree(data,"genus")
#' plot(conoTree)
#' 
#' #asaphid trilobites
#' data<-getTaxaPBDB("Asaphida")
#' asaTree<-makePBDBtaxontree(data,"genus")
#' plot(asaTree)
#' 
#' }
#' 
#' ###################################
#' 
#' #let's try time-scaling the graptolite tree
#' 
#' #get some example occurrence and taxonomic data
#' data(graptPBDB)
#' 
#' #get taxontree
#' graptTree<-makePBDBtaxontree(graptTaxaPBDB,"genus")
#' plot(graptTree)
#' 
#' #get time data from occurrences
#' graptOccGenus<-taxonSortPBDBocc(graptOccPBDB,rank="genus",onlyFormal=FALSE)
#' graptTimeGenus<-occData2timeList(occList=graptOccGenus)
#' 
#' #let's time-scale this tree with paleotree
#' timeTree<-bin_timePaleoPhy(graptTree,timeList=graptTimeGenus,
#' 	nonstoch.bin=TRUE,type="mbl",vartime=3)
#' 
#' #drops a lot of taxa; some of this is due to mispellings, etc
#' 
#' \dontrun{
#' 
#' #make pretty plot with library strap
#' library(strap)
#' geoscalePhylo(timeTree, ages=timeTree$ranges.used)
#' 
#' }
#' 


#' @name makePBDBtaxontree
#' @rdname makePBDBtaxontree
#' @export
makePBDBtaxontree<-function(data,rank){	
	if(!is(data,"data.frame")){stop("data isn't a data.frame")}
	if(length(rank)!=1){stop("rank must be a single value")}
	if(!any(sapply(c("species","genus","family","order","class","phylum"),function(x) x==rank))){
		stop("rank must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")}
	#if(!any(colnames(data)=="taxon_name")){stop("Data must be a taxonomic download under vocab='pbdb'")}
	if(!any(colnames(data)=="family")){stop("Data must be a taxonomic download with show=phylo")}
	# Do some translation
	#if com vocab
	if(any("rnk"==colnames(data))){	
		colnames(data)[colnames(data)=="rnk"]<-"taxon_rank"
		colnames(data)[colnames(data)=="fml"]<-"family"
		colnames(data)[colnames(data)=="odl"]<-"order"
		colnames(data)[colnames(data)=="cll"]<-"class"	
		colnames(data)[colnames(data)=="phl"]<-"phylum"	
		# taxon rank translation vectors for compact vocab
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")
		taxRankCOM<-2:26
		#change contents of "identified_rank" and "accepted_rank"
		data$taxon_rank<-sapply(data$taxon_rank,function(x) taxRankPBDB[x==taxRankCOM])
		message("compact vocab detected, relevant fields will be translated")
		}
	#if 1.1
	if(any(colnames(data)=="rank")){
		colnames(data)[colnames(data)=="rank"]<-"taxon_rank"
		}
	#
	#filter on rank
	data<-data[data[,"taxon_rank"]==rank,]
	#
	#if 1.2 and there is an accepted_name column
	if(any(colnames(data)=="accepted_name")){
		#fill empty accepted_name values with taxon_name
		nameFormal<-data[,"accepted_name"]
		nameFormal[is.na(nameFormal)]<-as.character(
			data[is.na(nameFormal),"taxon_name"])
		if(length(nameFormal)!=length(unique(nameFormal))){
			stop("Duplicated taxon names??")}
		#replace taxon_name
		data[,"taxon_name"]<-nameFormal	
		}
	#
	#get the fields you want
	taxonFields<-c("kingdom","phylum","class","order","family",
		"taxon_name")
	taxonData<-data[,taxonFields]
	taxonData<-apply(taxonData,2,as.character)
	#remove constant columns
	constantCol<-apply(taxonData,2,function(x) all(x==x[[1]]))
	constantCol[rev(which(constantCol))[1]]<-FALSE
	taxonData<-taxonData[,!constantCol]
	taxonData[taxonData==""]<-NA
	#are there any missing accepted names
	if(any(is.na(taxonData[,ncol(taxonData)]))){stop("Missing taxon names in dataset?")}
	#head(taxonData)
	labels<-taxonData[,ncol(taxonData)]
	edge<-matrix(NA,,2)
	#need to define columns BACKWARDS
	levels<-rev(1:ncol(taxonData))[-1]
	for(level in levels){
		newNodes<-unique(taxonData[,level])[!is.na(unique(taxonData[,level]))]
		for(node in newNodes){
			labels<-c(labels,node)
			nodeID<-length(labels)
			#find the descendant nodes
			whichDesc<-which(taxonData[,level]==node)
			descTaxRow<-lapply(whichDesc,function(x) (taxonData[x,-(1:level)]))
			descTax<-sapply(descTaxRow,function(x) x[!is.na(x)][1])
			descTax<-unique(descTax)
			descID<-sapply(descTax,function(x) which(x==labels))
			names(descID)<-NULL
			newEdge<-cbind(nodeID,descID)
			edge<-rbind(edge,newEdge)
			}
		}
	edge<-edge[-1,]
	edge.length<-rep(1,nrow(edge))
	Nnode<-length(unique(edge[,1]))
	tip.label<-taxonData[,ncol(taxonData)]
	node.label<-labels[-(1:length(tip.label))]
	#the node labels get lost when I do read.tree(write.tree())
		#I have no idea what to do about it...
	Ntip<-length(tip.label)
	#need to flip node numbers
	nodes<-sort(unique(edge[,1]))
	nodes<-cbind(c(1:Ntip,nodes),c(1:Ntip,rev(nodes)))
	edge[,1]<-sapply(edge[,1],function(x) nodes[x==nodes[,1],2])
	edge[,2]<-sapply(edge[,2],function(x) nodes[x==nodes[,1],2])
	#check root
	nodes<-sort(unique(edge[,1]))
	root<-nodes[sapply(nodes,function(x) all(x!=edge[,2]))]
	if(root!=(Ntip+1)){stop("Root isn't renumbering correctly")}
	#reorder edge
	edge<-edge[order(edge[,1],edge[,2]),]
	#make the tree
	tree<-list(edge=edge,tip.label=tip.label,edge.length=edge.length,
		Nnode=Nnode)	#,node.label=rev(node.label)
	class(tree)<-"phylo"
	#collapse singles
	tree<-collapse.singles(tree)
	#check it
	if(!testEdgeMat(tree)){stop("Edge matrix has inconsistencies")}
	#make it a good tree
	tree1<-reorder(tree,"cladewise") 	#REORDER IT
	tree1<-read.tree(text=write.tree(tree))
	tree1<-ladderize(tree1)
	#plot(tree1)
	return(tree1)
	}