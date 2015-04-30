#' Creating a Taxon-Tree from Taxonomic Data Downloaded from the Paleobiology Database
#'
#' This function creates phylogeny-like object of type \code{phylo} from the taxonomic information
#' recorded in a taxonomy download from the PBDB for a given group. Two different algorithms are provided,
#' the default being based on parent-child taxon relationships, the other based on the nested Linnean hierarchy.

#' @details
#' This function should not be taken too seriously. Many groups in the Paleobiology Database have
#' out-of-date or very incomplete taxonomic information. This function is meant to help visualize
#' what information is present, and by use of time-scaling functions, allow us to visualize the intersection
#' of temporal and phylogenetic, mainly to look for incongruence due to either incorrect taxonomic placements,
#' erroneous occurrence data or both. 
#'
#' Note however that, contrary to common opinion among some paleontologists, taxon-trees may be just as useful for 
#' macroevolutionary studies as reconstructed phylogenies (Soul and Friedman, in press.).

#' @param data A table of taxonomic data collected from the Paleobiology Database, using the taxa list option
#' with show=phylo. Should work with versions 1.1-1.2 of the API, with either the 'pbdb' or 'com' vocab. However,
#' as 'accepted_name' is not available in API v1.1, the resulting tree will have a taxon's *original* name and not
#' any formally updated name.

#' @param rank The selected taxon rank; must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'.

#' @param cleanTree By default, the tree is run through a series of post-processing, including having singles collapsed,
#' nodes reordered and being written out as a Newick string and read back in, to ensure functionality with ape functions
#' and ape-derived functions. If FALSE, none of this post-processing is done and users should beware, as such trees can
#' lead to hard-crashes of R.

#' @param method Controls which algorithm is used for calculating the taxon-tree: either \code{method="parentChild"}
#' (the default option) which converts the listed binary parent-child taxon relationships from the input PBDB data,
#' or \code{method="Linnean"}, which converts a taxon-tree by creating a table of the Linnean
#' taxonomic assignments (family, order, etc), which are provided when
#' option 'show=phylo' is used in PBDB API calls.

#' @param tipSet This argument only impacts analyses where the argument
#' \code{method="parentChild"} is also used. This \code{tipSet} controls
#' which taxa are selected as tip taxa for the
#' output tree. The default \code{tipSet="nonParents"} selects all child taxa which
#' are not listed as parents in \code{parentChild}. Alternatively, \code{tipSet="all"}
#' will add a tip to every internal node with the parent-taxon name encapsulated in
#' parentheses.

# @param cleanDuplicate If TRUE (\emph{not} the default), duplicated taxa of a taxonomic rank *not* selected by argument \code{rank}
# will be removed silently. Only duplicates of the taxonomic rank of interest will actually result in an error message.

#' @return
#' A phylogeny of class 'phylo', where each tip is a taxon of the given 'rank'. See additional details
#' regarding branch lengths can be found in the sub-algorithms used to create the taxon-tree by this function:
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}.
#'
#' Please note that when applied to output from the taxa option of the API version 1.1, the taxon names
#' returned are the \emph{original} taxon names as 'accepted_name' is not available in API v1.1, while
#' under API v1.2, the returned taxon names should be the most up-to-date formal names for those taxa.
#' Similar issues also effect the identification of parent taxa, as the accepted name of the
#' parent ID number is only provided in version 1.2 of the API.

#' @seealso
#' Two sub-algorithms in paleotree are used to create the taxon-tree within this function:
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}.
#'
#' Other functions for manipulating PBDB data can be found at \code{\link{taxonSortPBDBocc}},
#' \code{\link{occData2timeList}}, and the example data at \code{\link{graptPBDB}}.

#' @author David W. Bapst

#' @references
#' Soul, L. C., and M. Friedman. In Press. Taxonomy and Phylogeny Can Yield
#' Comparable Results in Comparative Palaeontological Analyses. \emph{Systematic Biology}
#' (http://sysbio.oxfordjournals.org/content/early/2015/03/23/sysbio.syv015.abstract)

#' @examples
#' \dontrun{
#' 
#' easyGetPBDBtaxa<-function(taxon,show=c("phylo","img","app")){
#' 	#let's get some taxonomic data
#' 	taxaData<-read.csv(paste0("http://paleobiodb.org/",
#' 		"data1.1/taxa/list.txt?base_name=",taxon,
#' 		"&rel=all_children&show=,
#'		paste0(show,collapse=","),"&status=senior"))
#' 	return(taxaData)
#' 	}
#' 
#' #graptolites
#' graptData<-easyGetPBDBtaxa("Graptolithina")
#' graptTree<-makePBDBtaxonTree(graptData,"genus")
#' plot(graptTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(graptTree$node.label,adj=c(0,1/2))
#' 
#' #conodonts
#' conoData<-easyGetPBDBtaxa("Conodonta")
#' conoTree<-makePBDBtaxonTree(conoData,"genus")
#' plot(conoTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(conoTree$node.label,adj=c(0,1/2))
#' 
#' #asaphid trilobites
#' asaData<-easyGetPBDBtaxa("Asaphida")
#' asaTree<-makePBDBtaxonTree(asaData,"genus")
#' plot(asaTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(asaTree$node.label,adj=c(0,1/2))
#' 
#' #Ornithischia
#' ornithData<-easyGetPBDBtaxa("Ornithischia")
#' #need to drop repeated taxon first: Hylaeosaurus
#' ornithData<-ornithData[-(which(ornithData[,"taxon_name"]=="Hylaeosaurus")[1]),]
#' ornithTree<-makePBDBtaxonTree(ornithData,"genus")
#' plot(ornithTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(ornithTree$node.label,adj=c(0,1/2))
#' 
#' #Rhynchonellida
#' rynchData<-easyGetPBDBtaxa("Rhynchonellida")
#' #need to drop repeated taxon first: Rhynchonelloidea
#' rynchData<-rynchData[-(which(rynchData[,"taxon_name"]=="Rhynchonelloidea")[1]),]
#' rynchTree<-makePBDBtaxonTree(rynchData,"genus")
#' plot(rynchTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(rynchTree$node.label,adj=c(0,1/2))
#' 
#' #some of these look pretty messy!
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
#' #get the taxon tree
#' graptTree<-makePBDBtaxonTree(graptTaxaPBDB,"genus")
#' plot(graptTree)
#' nodelabels(graptTree$node.label,cex=0.5)
#' 
#' #get time data from occurrences
#' graptOccGenus<-taxonSortPBDBocc(graptOccPBDB,rank="genus",onlyFormal=FALSE)
#' graptTimeGenus<-occData2timeList(occList=graptOccGenus)
#' 
#' #let's time-scale this tree with paleotree
#'		# use minimum branch length for visualization
#' 		# and nonstoch.bin so we plot maximal ranges
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
#' nodelabels(timeTree$node.label,cex=0.5)
#' }
#' 

#' @name makePBDBtaxonTree
#' @rdname makePBDBtaxonTree
#' @export
makePBDBtaxonTree<-function(data,rank,method="parentChild",tipSet="nonParents",cleanTree=TRUE){		# ,cleanDuplicate=FALSE
	#CHECKS
	if(length(method)!=1 | !is.character(method)){
		stop("method must be a single character value")}
	if(!any(method==c("Linnean","parentChild"))){
		stop("method must be one of either 'Linnean' or 'parentChild'")}
	if(!is(data,"data.frame")){stop("data isn't a data.frame")}
	if(length(rank)!=1 | !is.character(rank)){
		stop("rank must be a single character value")}
	if(!any(sapply(c("species","genus","family","order","class","phylum"),function(x) x==rank))){
		stop("rank must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")}
	#
	#translate to a common vocabulary
	data<-translatePBDBtaxa(data)
	#
	if(method=="parentChild"){
		# translate rank and taxon_rank to a number
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")	#keep informal as high, never drop!
		rank<-which(rank==taxRankPBDB)
		numTaxonRank<-sapply(data[,"taxon_rank"],function(x) which(x==taxRankPBDB))		
		#drop taxa below specified rank
		data<-data[rank>=numTaxonRank,]
		#get parent name
		parentName<-sapply(data[,"parent_no"],function(x){ 
			z<-which(data[,"taxon_no"]==x)
			#check for multiple matches
			if(length(z)>1){
				stop("Multiple taxon matches to taxon IDs? Check if you used '&status=senior' in API call")}
			#if zero matches
			if(length(z)<1){
				z<-as.character(x)
			}else{
				z<-data[,"taxon_name"]}
			return(z)
			})
		#get matrix of all parentChild relationships
		parentChildAll<-cbind(parentName,data[,childName])
		#first pull out all for desired rank
		parentChildMat<-parentChildAll[rank==numTaxonRank,]
		#count parents floating without ancestors of their own
		floatParents<-parentChildMat[sapply(parentChildMat[,1],function(x)
			all(x!=parentChildmat[,2])),1]
		#start tracing back
		while(length(floatParents)>1){	#so only a root can float
			newRelations<-match(floatParents,parentChildAll[,2])
			newRelations<-newRelations[!is.na(newRelations)]
			parentChildMat<-rbind(parentChildMat,parentChildAll[newRelations,])
			#recalculate float
			floatParents2<-parentChildMat[sapply(parentChildMat[,1],function(x)
				all(x!=parentChildmat[,2])),1]
			#put in stopping condition
			if(length(floatParents)==length(floatParents2)){
				stop("Provided PBDB Dataset does not appear to have a monophyletic set of parent-child relationship pairs")
			}else{
				floatParents<-floatParents2}
			}
		tree<-taxonTable2taxonTree(parentChild=parentChildMat,tipSet=tipSet,cleanTree=cleanTree)
		}
	#
	if(method=="Linnean"){
		# TRASH:
		#check to make sure no taxon names are listed twice, first clean then check again
		#nDup<-sapply(data[,"taxon_name"],function(x) sum(data[,"taxon_name"]==x)>1)
		#if(any(nDup) & cleanDuplicate){
		#	#find any taxa of not right rank, remove them
		#	droppers<-which(data[,"taxon_rank"]!=rank & nDup)
		#	message(paste0("Duplicate taxa dropped: ",
		#		paste0("(",droppers,") ",data[droppers,"taxon_name"]," [rank: ",data[droppers,"taxon_rank"],"]",collapse=", ")))
		#	data<-data[-(droppers),]
		#	}
		#
		#Check if show=phylo was used
		if(!any(colnames(data)=="family")){stop("Data must be a taxonomic download with show=phylo for method='Linnean'")}
		#now check and return an error if duplicates remain
		nDup<-sapply(nrow(data),function(x) sum(data[,"taxon_name"]==data[x,"taxon_name"])>1 & data[x,"taxon_rank"]==rank)
		if(any(nDup)){
			stop(paste0("Duplicate taxa of selected rank: ",paste0("(",which(nDup),") ",data[nDup,"taxon_name"],collapse=", ")))}
		#filter on rank
		data<-data[data[,"taxon_rank"]==rank,]
		#
		#get the fields you want
		taxonFields<-c("kingdom","phylum","class","order","family",
			"taxon_name")
		taxonData<-data[,taxonFields]
		taxonData<-apply(taxonData,2,as.character)
		tree<-taxonTable2taxonTree(taxonTable=taxonData,cleanTree=cleanTree)
		}
	return(tree)
	}

#hidden function, don't import
translatePBDBtaxa<-function(data){
	# Do some translation
	#if com vocab
	if(any("rnk"==colnames(data))){	
		#apparently it doesn't matter if these columns *are* present or not
		colnames(data)[colnames(data)=="acn"]<-"accepted_name"
		colnames(data)[colnames(data)=="snp"]<-"senpar_no"
		colnames(data)[colnames(data)=="rnk"]<-"taxon_rank"
		colnames(data)[colnames(data)=="nam"]<-"taxon_name"
		colnames(data)[colnames(data)=="fml"]<-"family"
		colnames(data)[colnames(data)=="odl"]<-"order"
		colnames(data)[colnames(data)=="cll"]<-"class"	
		colnames(data)[colnames(data)=="phl"]<-"phylum"	
		colnames(data)[colnames(data)=="kgl"]<-"kingdom"
		colnames(data)[colnames(data)=="par"]<-"parent_no"
		colnames(data)[colnames(data)=="oid"]<-"taxon_no"
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
	#if 1.1 and vocab is pbdb
	if(any(colnames(data)=="rank")){
		colnames(data)[colnames(data)=="rank"]<-"taxon_rank"
		}
	#
	#if 1.2 and there is an accepted_name column..
	if(any(colnames(data)=="accepted_name")){
		#fill empty accepted_name values with taxon_name
		nameFormal<-data[,"accepted_name"]
		nameFormal[is.na(nameFormal)]<-as.character(data[is.na(nameFormal),"taxon_name"])
		#replace taxon_name
		data[,"taxon_name"]<-nameFormal	
		#also replace parent_no in the same way with senpar_no
		parNum<-data[,"senpar_no"]
		parNum[is.na(parNum)]<-as.character(data[is.na(parNum),"parent_no"])	
		data[,"parent_no"]<-parNum		
		}
	#
	return(data)
	}