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
#' Depending on the \code{method}
#' used, either the element \code{$parentChild} or \code{$taxonTable} is added to the list structure of
#' the output phylogeny object, which was used as input for one of the two algorithms mentioned above.
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
#'		paste0(show,collapse=","),"&status=senior"),
#'		stringsAsFactors=FALSE)
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
#' #get the taxon tree: Linnean method first
#' graptTree<-makePBDBtaxonTree(graptTaxaPBDB,"genus",method="Linnean")
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
makePBDBtaxonTree<-function(data,rank,method="parentChild",queryMissing=FALSE,
					tipSet="nonParents",cleanTree=TRUE){		
	#
	# rank="genus"; method="parentChild"; tipSet="nonParents"; cleanTree=TRUE
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
	data1<-translatePBDBtaxa(data)
	#
	if(method=="parentChild"){
		#need two things: a table of parent-child relationships as IDs
			#and a look-up table of IDs and taxon names
		#
		#filer out lower than selected rank
		# translate rank and taxon_rank to a number
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")	#keep informal as high, never drop!
		rank1<-which(rank==taxRankPBDB)
		numTaxonRank<-sapply(data1[,"taxon_rank"],function(x) which(x==taxRankPBDB))		
		#drop taxa below specified rank
		data1<-data1[rank1<=numTaxonRank,]
		#also recreate numTaxonRank
		numTaxonRank<-sapply(data1[,"taxon_rank"],function(x) which(x==taxRankPBDB))		
		#
		#create lookup table for taxon names
		taxonNameTable<-data1[,c("taxon_no","taxon_name")]
		#add parents not listed
		parentFloat<-unique(data[,"parent_no"])
		parentFloat<-parentFloat[is.na(match(parentFloat,taxonNameTable[,1]))]
		taxonNameTable<-rbind(taxonNameTable,cbind(parentFloat,paste("ID:",as.character(parentFloat)))
		#DONE
		#
		#now need to put together parentChild table
		#first, get table of all parentChild relationships in data
		pcAll<-data[,c("parent_no","taxon_no")]
		#then start creating final matrix: first, those of desired rank
		pcMat<-pcAll[rank1==numTaxonRank,]
		




		#identify IDs of parents floating without ancestors of their own
		getFloat<-function(pcMat){unique(pcMat[sapply(pcMat[,1],function(x) all(x!=pcMat[,2])),1])}
		floatParents<-getFloat(pcMat=parentChildMat)

		#TRASH
		#get parent name
		#parentName<-sapply(data1[,"parent_no"],function(x){ 
			z<-which(data1[,"taxon_no"]==x)
			#check for multiple matches
			if(length(z)>1){
				stop("Multiple taxon matches to taxon IDs? Check if you used '&status=senior' in API call")}
			#if zero matches
			if(length(z)<1){
				z<-paste("ID:",as.character(x))
			}else{
				z<-data1[z,"taxon_name"]
				}
			return(z)
			})
		#check that the same name wasn't returned for multiple ID numbers
		parentNameNo<-unique(cbind(parentName,data1[,"parent_no"]))
		if(any(duplicated(parentNameNo[,2]))){
			stop("multiple taxon numbers assigned to same parent taxon name?")}


		#start tracing back
		nCount<-0
		while(length(floatParents)>1){	#so only a root can float
			nCount<-nCount+1
			newRelations<-match(floatParents,parentChildAll[,2])
			newRelations<-newRelations[!is.na(newRelations)]
			parentChildMat<-rbind(parentChildMat,parentChildAll[newRelations,])
			#recalculate float
			floatParents2<-getFloat(pcMat=parentChildMat)
			#put in stopping condition
			if(length(floatParents2)>1 & identical(sort(floatParents),sort(floatParents2))){
				if(queryMissing){
					getPBDBtaxaID<-function(ID){
						#let's get some taxonomic data
						taxaData<-read.csv(paste0("http://paleobiodb.org/",
							"data1.1/taxa/list.txt?id=",paste0(ID,collapse=","),
							"&rel=self&status=senior&vocab=pbdb"),
							stringsAsFactors=FALSE)
							return(taxaData)
							}
					floatID<-parentNameNo[match(floatParents2,parentNameNo[,1]),2]
					#drop Eukarya, as it won't return if status=senior under 1.1
					floatParents2<-floatParents2[floatID!="1"]
					floatID<-floatID[floatID!="1"]	
					floatData<-getPBDBtaxaID(floatID)
					floatData<-translatePBDBtaxa(floatData)
					#update taxon names in parentNameNo,parentChildMat,parentChildAll
					newTaxa<-cbind(floatData$taxon_name,floatData[,c("taxon_no")])
					newTaxa<-cbind(oldname=paste("ID:",as.character(newTaxa[,2])),newTaxa)
					for(i in 1:nrow(newTaxa)){
						parentNameNo[parentNameNo==newTaxa[i,1]]<-newTaxa[i,2]
						parentChildMat[parentChildMat==newTaxa[i,1]]<-newTaxa[i,2]
						parentChildAll[parentChildAll==newTaxa[i,1]]<-newTaxa[i,2]
						}
					newParents<-floatData[,"parent_no"]
					newParents<-cbind(oldname=paste("ID:",as.character(newParents)),newParents)
					#update parentNameNo
					parentNameNo<-rbind(parentNameNo,newParents)
					#update parentChildMat, parentChildAll
					newEntries<-cbind(newParents[,1],floatParents2)
					parentChildMat<-rbind(parentChildMat,newEntries)
					parentChildAll<-rbind(parentChildAll,newEntries)
					floatParents<-getFloat(pcMat=parentChildMat)
				}else{
					stop(paste0("Provided PBDB Dataset does not appear to have a \n",
						" monophyletic set of parent-child relationship pairs. \n",
						"Multiple taxa appear to be listed as parents, but are not \n",
						"listed themselves so have no parents listed: \n",
						paste0(floatParents,collapse=", ")))
					}
			}else{
				floatParents<-floatParents2}
			}
		tree<-parentChild2taxonTree(parentChild=parentChildMat,tipSet=tipSet,cleanTree=cleanTree)
		tree$parentChild<-parentChildMat
		}
	#
	if(method=="Linnean"){
		#Check if show=phylo was used
		if(!any(colnames(data1)=="family")){stop("Data must be a taxonomic download with show=phylo for method='Linnean'")}
		#message that tipSet and queryMissing are ignored
		message("Linnean taxon-tree option selected, arguments 'tipSet' and 'queryMissing' ignored")
		#now check and return an error if duplicate taxa of selected rank
		nDup<-sapply(nrow(data1),function(x) sum(data1[,"taxon_name"]==data1[x,"taxon_name"])>1 & data1[x,"taxon_rank"]==rank)
		if(any(nDup)){
			stop(paste0("Duplicate taxa of selected rank: ",paste0("(",which(nDup),") ",data1[nDup,"taxon_name"],collapse=", ")))}
		#filter on rank
		data1<-data1[data1[,"taxon_rank"]==rank,]
		#
		#get the fields you want
		taxonFields<-c("kingdom","phylum","class","order","family",
			"taxon_name")
		taxonData<-data1[,taxonFields]
		taxonData<-apply(taxonData,2,as.character)
		tree<-taxonTable2taxonTree(taxonTable=taxonData,cleanTree=cleanTree)
		tree$taxonTable<-taxonData
		}
	return(tree)
	}

#hidden function, don't import
translatePBDBtaxa<-function(data){
	# Do some translation
	#need to replace any empty string values with NAs (due perhaps to use of read.csv with the API)
	data[data==""]<-NA
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
		#
		#replace taxon_no with accepted_no
		taxNum<-data[,"accepted_no"]
		taxNum[is.na(taxNum)]<-as.character(data[is.na(taxNum),"taxon_no"])	
		data[,"taxon_no"]<-taxNum			
		#
		#also replace parent_no in the same way with senpar_no
		parNum<-data[,"senpar_no"]
		parNum[is.na(parNum)]<-as.character(data[is.na(parNum),"parent_no"])	
		data[,"parent_no"]<-parNum		
		}
	#
	return(data)
	}