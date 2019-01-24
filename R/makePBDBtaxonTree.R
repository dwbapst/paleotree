#' Creating a Taxon-Tree from Taxonomic Data Downloaded from the Paleobiology Database
#'
#' The function \code{makePBDBtaxonTree} creates phylogeny-like 
#' object of class \code{phylo} from the taxonomic information
#' recorded in a taxonomy download from the PBDB for
#' a given group. Two different algorithms are provided,
#' the default being based on parent-child taxon relationships,
#' the other based on the nested Linnean hierarchy. The function
#' \code{plotTaxaTreePBDB} is also provided as a minor helper
#' function for optimally plotting the labeled topologies that are
#' output by \code{makePBDBtaxonTree}.
#' 
 

#' @details
#' This function should not be taken too seriously.
#' Many groups in the Paleobiology Database have
#' out-of-date or very incomplete taxonomic information.
#' This function is meant to help visualize
#' what information is present, and by use of time-scaling
#' functions, allow us to visualize the intersection
#' of temporal and phylogenetic, mainly to look for incongruence
#' due to either incorrect taxonomic placements,
#' erroneous occurrence data or both. 
#'
#' Note however that, contrary to common opinion among some
#' paleontologists, taxon-trees may be just as useful for 
#' macroevolutionary studies as reconstructed phylogenies
#' (Soul and Friedman, 2015.).

#' @param taxaDataPBDB A table of taxonomic data collected from
#' the Paleobiology Database, using the taxa list option
#' with \code{show = class}. Should work with versions 1.1-1.2 of
#' the API, with either the \code{pbdb} or \code{com} vocab. However,
#' as \code{accepted_name} is not available in API v1.1, the resulting
#' tree will have a taxon's *original* name and not
#' any formally updated name.

#' @param rank The selected taxon rank; must be one of 'species',
#' 'genus', 'family', 'order', 'class' or 'phylum'.

#' @param cleanTree When \code{TRUE} (the default), the tree is run through a series of
#' post-processing, including having singles collapsed,
#' nodes reordered and being written out as a Newick string and read
#' back in, to ensure functionality with ape functions
#' and ape-derived functions. If \code{FALSE}, none of this post-processing
#' is done and users should beware, as such trees can
#' lead to hard-crashes of R.

#' @param method Controls which algorithm is used for calculating
#' the taxon-tree. The default option is \code{method  = "parentChild"}
#' which converts the listed binary parent-child taxon relationships in
#' the Paleobiology Database- these parent-child relationships (if missing
#' from the input dataset) are autofilled using API calls to the
#' Paleobiology Database. Alternatively, users may use
#' \code{method = "Linnean"}, which converts the table of Linnean taxonomic
#' assignments (family, order, etc as provided by \code{show = class} in
#' PBDB API calls) into a taxon-tree. Two methods formerly both implemented
#' under \code{method  = "parentChild"} are also available as
#' \code{method = "parentChildOldMergeRoot"} and \code{method = "parentChildOldQueryPBDB"}
#' respectively. Both of these use similar algorithms as the current
#' \code{method  = "parentChild"} but differ in how they treat taxa with
#' parents missing from the input taxonomic dataset.
#' \code{method = "parentChildOldQueryPBDB"} behaves most similar
#' to \code{method = "parentChild"}  in that it queries the Paleobiology
#' Database via the API , but repeatedly does so for information on parent
#' taxa of the 'floating' parents, and continues within a \code{while}
#' loop until only one such unassigned parent taxon remains. This latter
#' option may talk a long time or never finish, depending on the
#' linearity and taxonomic structures encountered in the PBDB taxonomic
#' data; i.e. if someone a taxon was ultimately its own indirect child
#' in some grand loop by mistake, then under this option
#' \code{makePBDBtaxonTree} might never finish. In cases where taxonomy
#' is bad due to weird and erroneous taxonomic assignments reported by
#' the PBDB, this routine may search all the way back to a very ancient
#' and deep taxon, such as the \emph{Eukaryota} taxon.
#' \code{method = "parentChildOldMergeRoot"} will combine these disparate
#' potential roots and link them to an artificially-constructed
#' pseudo-root, which at least allows for visualization of the taxonomic
#' structure in a limited dataset. This latter option will be fully
#' offline, as it does nto do any additional API calls
#' of the Paleobiology Database, unlike other options.

#  @param solveMissing Under \code{method  = "parentChild"}, what should \code{makePBDBtaxonTree} do about
#  multiple 'floating' parent taxa, listed without their own parent taxon information in the input
#   dataset under \code{taxaDataPBDB}? Each of these is essentially a separate root taxon, for a different set
#   of parent-child relationships, and thus poses a problem as far as returning a single phylogeny is
#   concerned. If \code{solveMissing = NULL} (the default), nothing is done and the operation halts with
#   an error, reporting the identity of these taxa. Two alternative solutions are offered: first,
#   \code{solveMissing  = "mergeRoots"} will combine these disparate potential roots and link them to an
#   artificially-constructed pseudo-root, which at least allows for visualization of the taxonomic
#   structure in a limited dataset. Secondly, \code{solveMissing  = "queryPBDB"} queries the Paleobiology
#   Database repeatedly via the API for information on parent taxa of the 'floating' parents, and continues
#   within a \code{while()} loop until only one such unassigned parent taxon remains. This latter option may
#   talk a long time or never finish, depending on the linearity and taxonomic structures encountered in the
#   PBDB taxonomic data; i.e. if someone a taxon was ultimately its own indirect child in some grand loop by
#   mistake, then under this option \code{makePBDBtaxonTree} might never finish. In cases where taxonomy is
#   bad due to weird and erroneous taxonomic assignments reported by the PBDB, this routine may search all
#   the way back to a very ancient and deep taxon, such as the Eukaryota taxon.
#  Users should thus use \code{solveMissing  = "queryPBDB"} only with caution.

#' @param tipSet This argument only impacts analyses where the argument
#' \code{method  = "parentChild"} is also used. This \code{tipSet} controls
#' which taxa are selected as tip taxa for the
#' output tree. The default \code{tipSet  = "nonParents"} selects all child taxa which
#' are not listed as parents in \code{parentChild}. Alternatively, \code{tipSet = "all"}
#' will add a tip to every internal node with the parent-taxon name encapsulated in
#' parentheses.

#' @param APIversion Version of the Paleobiology Database API used by
#' \code{makePBDBtaxonTree} when \code{method  = "parentChild"} or
#' \code{method  = "parentChildOldQueryPBDB"} is used. The current default
#' is \code{APIversion = "1.2"}, the most recent API version as of 12/11/2018.

# @param cleanDuplicate If \code{TRUE} (\emph{not} the default), duplicated taxa of a
# taxonomic rank \emph{not} selected by argument \code{rank}
# will be removed silently. Only duplicates of the taxonomic rank of interest
# will actually result in an error message.

#' @param annotatedDuplicateNames A logical determining whether duplicate taxon names,
#' when found in the Paleobiology Database for taxa (presumably reflecting an issue with
#' taxa being obsolete but with incomplete seniority data), should be annotated to include
#' sequential numbers so to modify them, via function\code{base}'s
#' \code{\link[base]{make.unique}}. This only applies to
#' \code{method = "parentChild"}, with the default option being
#' \code{annotatedDuplicateNames = TRUE}. If more than 26 duplicates are found, an error
#' is issued. If this argument is \code{FALSE}, an error is issued if duplicate taxon
#' names are found.

#' @param taxaTree A phylogeny of class \code{phylo}, presumably a taxon tree as output from
#' \code{makePBDBtaxonTree} with higher-taxon names as node labels.

#' @param edgeLength The edge length that the plotted tree should be plotted
#' with (\code{plotTaxaTreePBDB} plots phylogenies as non-ultrametric,
#' not as a cladogram with aligned tips).

#' @return
#' A phylogeny of class \code{phylo}, where each tip is a taxon of the given 'rank'. See additional details
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
#' Two other functions in paleotree are used as sub-algorithms by \code{makePBDBtaxonTree}
#' to create the taxon-tree within this function,
#' and users should consult their manual pages for additional details:
#'
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}
#' 
#' Closely related functions for 
#'
#' Other functions for manipulating PBDB data can be found at \code{\link{taxonSortPBDBocc}},
#' \code{\link{occData2timeList}}, and the example data at \code{\link{graptPBDB}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.
#' 
#' Soul, L. C., and M. Friedman. 2015. Taxonomy and Phylogeny Can Yield
#' Comparable Results in Comparative Palaeontological Analyses. \emph{Systematic Biology} 
#' (\href{http://sysbio.oxfordjournals.org/content/early/2015/03/23/sysbio.syv015.abstract}{Link})


#' @examples
#' \dontrun{
#' 
#' 
#' #graptolites
#' graptData <- easyGetPBDBtaxa("Graptolithina")
#' graptTree <- makePBDBtaxonTree(graptData,"genus",
#' 	   method = "parentChild")
#' #try Linnean
#' graptTree <- makePBDBtaxonTree(graptData,"genus",
#' 	   method = "Linnean")
#' 
#' # plot it!
#' plot(graptTree,show.tip.label = FALSE,
#'     no.margin = TRUE,edge.width = 0.35)
#' nodelabels(graptTree$node.label,adj = c(0,1/2))
#' 
#' #################
#' #conodonts
#' conoData <- easyGetPBDBtaxa("Conodonta")
#' conoTree <- makePBDBtaxonTree(conoData,"genus",
#' 	   method = "parentChild")
#' 
#' # plot it!
#' plot(conoTree,show.tip.label = FALSE,
#'    no.margin = TRUE,edge.width = 0.35)
#' nodelabels(conoTree$node.label,adj = c(0,1/2))
#' 
#' ############################
#' #asaphid trilobites
#' asaData <- easyGetPBDBtaxa("Asaphida")
#' asaTree <- makePBDBtaxonTree(asaData,"genus",
#' 	   method = "parentChild")
#' 
#' # plot it!
#' plot(asaTree,show.tip.label = FALSE,
#'    no.margin = TRUE,edge.width = 0.35)
#' nodelabels(asaTree$node.label,adj = c(0,1/2))
#' 
#' #Ornithischia
#' ornithData <- easyGetPBDBtaxa("Ornithischia")
#' ornithTree <- makePBDBtaxonTree(ornithData,"genus",
#' 	  method = "parentChild")
#' #try Linnean
#' #need to drop repeated taxon first: Hylaeosaurus
#' ornithData <- ornithData[
#'    -(which(ornithData[,"taxon_name"] == "Hylaeosaurus")[1]),]
#' ornithTree <- makePBDBtaxonTree(ornithData,"genus",
#' 	  method = "Linnean")
#' plot(ornithTree,show.tip.label = FALSE,
#'    no.margin = TRUE,edge.width = 0.35)
#' nodelabels(ornithTree$node.label,adj = c(0,1/2))
#' 
#' #Rhynchonellida
#' rynchData <- easyGetPBDBtaxa("Rhynchonellida")
#' rynchTree <- makePBDBtaxonTree(rynchData,"genus",
#' 	  method = "parentChild")
#' plot(rynchTree, show.tip.label = FALSE, 
#'    no.margin = TRUE, edge.width = 0.35)
#' nodelabels(rynchTree$node.label,adj = c(0,1/2))
#' 
#' #some of these look pretty messy!
#' 
#' }
#' 
#' ###################################
#' \donttest{
#' 
#' #let's try time-scaling the graptolite tree
#' 
#' #get some example occurrence and taxonomic data
#' data(graptPBDB)
#' 
#' #get the taxon tree: Linnean method
#' graptTree <- makePBDBtaxonTree(graptTaxaPBDB,
#'     "genus", method = "Linnean")
#' plot(graptTree, cex = 0.4)
#' nodelabels(graptTree$node.label, cex = 0.5)
#'
#' #get the taxon tree: parentChild method
#' graptTree <- makePBDBtaxonTree(graptTaxaPBDB,
#'     "genus", method = "parentChild")
#' plot(graptTree,cex = 0.4)
#' nodelabels(graptTree$node.label,cex = 0.5)
#' 
#' #get time data from occurrences
#' graptOccGenus <- taxonSortPBDBocc(graptOccPBDB,
#'     rank = "genus", onlyFormal = FALSE)
#' graptTimeGenus <- occData2timeList(occList = graptOccGenus)
#' 
#' #let's time-scale the parentChild tree with paleotree
#' 		# use minimum branch length for visualization
#' 		# and nonstoch.bin so we plot maximal ranges
#' timeTree <- bin_timePaleoPhy(graptTree,
#'     timeList = graptTimeGenus,
#'     nonstoch.bin = TRUE,
#'     type = "mbl", vartime = 3)
#' 
#' #drops a lot of taxa; some of this is due to mispellings, etc
#' 
#' }
#' \dontrun{
#' 
#' #make pretty plot with library strap
#' library(strap)
#' geoscalePhylo(timeTree,
#'     ages = timeTree$ranges.used)
#' nodelabels(timeTree$node.label,cex = 0.5)
#'
#' }
#' 



#' @name makePBDBtaxonTree
#' @rdname makePBDBtaxonTree
#' @export
makePBDBtaxonTree <- function(taxaDataPBDB, rank,
					method = "parentChild", #solveMissing = NULL,
					tipSet = "nonParents", cleanTree = TRUE,
					annotatedDuplicateNames = TRUE,
					APIversion = "1.2"){		
	############################################################
	############################################################
	# library(paleotree);data(graptPBDB);
	# taxaDataPBDB <- graptTaxaPBDB; rank = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
	# taxaDataPBDB <- graptTaxaPBDB; rank = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
	# taxaDataPBDB <- graptTaxaPBDB; rank = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
	# taxaDataPBDB <- graptTaxaPBDB; rank = "genus"; method = "Linnean"; 
	#
	#CHECKS
	if(length(method) != 1 | !is.character(method)){
		stop("method must be a single character value")
		}
	if(!any(method == c("Linnean", "parentChild"))){
		stop("method must be one of either 'Linnean' or 'parentChild'")
		}
	#if(!is.null(solveMissing)){
	#	if(length(solveMissing)>1 | !is.character(solveMissing)){
	#		stop("solveMissing must be either NULL or a single character value")
	#		}
	#	if(is.na(match(solveMissing,c("queryPBDB","mergeRoots")))){
	#		stop('solveMissing but be either NULL or "queryPBDB" or "mergeRoots"')
	#		}
	#	}
	if(!is(taxaDataPBDB,"data.frame")){
		stop("taxaDataPBDB isn't a data.frame")
		}
	if(length(rank) != 1 | !is.character(rank)){
		stop("rank must be a single character value")
		}
	if(!any(sapply(c("species","genus","family","order","class","phylum"),
			function(x) x == rank))){
		stop("rank must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")
		}
	#
	#########################################
	# CLEAN DATA
	#
	#translate to a common vocabulary
	dataTransform <- translatePBDBtaxa(taxaDataPBDB)
	dataTransform <- unique(dataTransform)
	#
	if(method == "parentChild"){
		# need two things: a table of parent-child relationships as IDs
			#and a look-up table of IDs and taxon names
		# 
		##############################
		# FIND ALL PARENTS FIRST
			# three column matrix with taxon name, taxon ID, parent ID
			# (in that order)
		parData<- getAllParents(dataTransform,
			status="all", annotatedDuplicateNames = annotatedDuplicateNames)
		#
		#######################################
		# NOW FILTER OUT TIP TAXA WE WANT
		#tipIDs <- getTaxaIDsDesiredRank(data=dataTransform, rank=rank)
		tipIDs <- dataTransform$taxon_no[dataTransform$taxon_rank==rank]
		# figure out which taxon numbers match tip IDs
		whichTip <- match(tipIDs, parData$taxon_no)
		#
		###############################
		# BUILD PARENT-CHILD matrix
		# get parent-child matrix for just desired OTUs 
			# start matrix with those parent-child relationships
			# subset these from parData using the taxon ID numbers
		pcMat <- subsetParDataPBDB(subsetNum = tipIDs,
			parData = parData)			
		# starting from desired tip OTUs, work backwards to a common ancestor from the full parData
		pcMat<-constructParentChildMatrixPBDB(initPCmat=pcMat,
			parData = parData)	
		#################################
		# convert parent-child matrix to accepted taxon names
		#print(pcMat)
		pcMat <- apply(pcMat, 1:2, 
			convertParentChildMatNames,
			parData = parData)
		# remove "NODE" root-stem 
		pcMat <- pcMat[pcMat[,1] != "NODE", ]
		############################
		# Calculate the taxon tree!
		tree <- parentChild2taxonTree(
			parentChild = pcMat,
			tipSet = tipSet,
			cleanTree = cleanTree)
		#
		tree$parentChild <- pcMat
		}
	##############
	#
	if(method == "Linnean"){
		dataTransform <- apply(dataTransform, 2, as.character)
		#Check if show = class was used
		if(!any(colnames(dataTransform) == "family")){
			stop("taxaDataPBDB must be a taxonomic download with show = class for method = 'Linnean'")
			}
		#message that tipSet (and solveMissing) is ignored
		message("Linnean taxon-tree option selected, argument 'tipSet' is ignored")
		#now check and return an error if duplicate taxa of selected rank
		nDup <- sapply(nrow(dataTransform),function(x)
			sum(dataTransform[,"taxon_name"] == dataTransform[x,"taxon_name"])>1
			 & dataTransform[x,"taxon_rank"] == rank
			)
		if(any(nDup)){
			stop(
				paste0(
					"Duplicate taxa of selected rank: ",
					paste0("(",which(nDup),") ",
					dataTransform[nDup,"taxon_name"],
					collapse = ", ")
					)
				)
			}
		#filter on rank
		dataTransform <- dataTransform[dataTransform[,"taxon_rank"] == rank,]
		#
		#get the taxonomic fields you want
		#from API documentation: "phylum","class","order","family","genus"  
			# apparently 'kingdom' dropped from v1.2 API
		taxonFields <- c(
			#"kingdom",
			"phylum", "class",
			"order", "family",
			"taxon_name")
		#print(colnames(dataTransform))
		taxonData <- dataTransform[,taxonFields]
		taxonData <- apply(taxonData,2,as.character)
		tree <- taxonTable2taxonTree(taxonTable = taxonData,cleanTree = cleanTree)
		tree$taxonTable <- taxonData
		}
	#########
	#
	if(method == "parentChildOldMergeRoot" | method == "parentChildOldQueryPBDB"){
		dataTransform <- apply(dataTransform, 2, as.character)
		# need two things: a table of parent-child relationships as IDs
			#and a look-up table of IDs and taxon names
		# 
		# filer out lower than selected rank
		# 
		# translate rank and taxon_rank to a number
		# taxon rank translation vectors for compact vocab
		taxRankPBDB <- getTaxRankPBDB()
		rankID <- which(rank == taxRankPBDB)
		numTaxonRank <- sapply(dataTransform[,"taxon_rank"],
			function(x) which(x == taxRankPBDB))		
		#drop taxa below specified rank
		dataTransform <- dataTransform[rankID <= numTaxonRank,]
		#also recreate numTaxonRank
		numTaxonRank <- sapply(dataTransform[,"taxon_rank"],
			function(x) which(x == taxRankPBDB))		
		#
		#create lookup table for taxon names
		taxonNameTable <- cbind(
			as.numeric(dataTransform[,"taxon_no"]), 
			as.character(dataTransform[,"accepted_name"])
			)
		#add parents not listed
		parentFloat <- unique(dataTransform[,"parent_no"])
		parentFloat <- parentFloat[is.na(match(parentFloat, taxonNameTable[,1]))]
		taxonNameTable <- rbind(taxonNameTable,
			cbind(parentFloat,
				paste("ID:", as.character(parentFloat))
				)
			)
		#DONE
		#
		#now need to put together parentChild table
		#first, get table of all parentChild relationships in taxaDataPBDB
		pcAll <- cbind(as.numeric(dataTransform[,"parent_no"]), 
			as.numeric(dataTransform[,"taxon_no"]))
		#then start creating final matrix: first, those of desired rank
		pcMat <- pcAll[rankID == numTaxonRank,]
		#identify IDs of parents floating without ancestors of their own
		floaters <- getFloatAncPBDB(pcDat = pcMat)
		#
		nCount <- 0	#start tracing back
		while(length(floaters)>1){	#so only ONE root can float
			nCount <- nCount+1
			#get new relations: will 'anchor' the floaters
			anchors <- match(floaters,pcAll[,2])
			pcMat <- rbind(pcMat, pcAll[anchors[!is.na(anchors)],])	#bind to pcMat
			floatersNew <- getFloatAncPBDB(pcDat = pcMat)	#recalculate float
			#stopping condition, as this is a silly while() loop...
			if(length(floatersNew)>1 & identical(sort(floaters),sort(floatersNew))){
				#if(!is.null(solveMissing)){
					if(method == "parentChildOldQueryPBDB"){
						floatData <- queryMissingParents(
							taxaID = floatersNew, 
							APIversion = APIversion
							)	
						#update taxon names in taxonNameTable
						whichUpdate <- match(floatData[,"taxon_no"],taxonNameTable[,1])
						taxonNameTable[whichUpdate[!is.na(whichUpdate)],2] <- floatData[
							!is.na(whichUpdate),"taxon_name"]
						#add any new parent taxa to taxonNameTable
						parentFloat <- unique(floatData[,"parent_no"])
						parentFloat <- parentFloat[is.na(match(parentFloat,taxonNameTable[,1]))]
						if(length(parentFloat)>0){
							taxonNameTable <- rbind(taxonNameTable,
								cbind(parentFloat,paste("ID:",as.character(parentFloat))))
							}
						#update parentChildMat, parentChildAll
						newEntries <- floatData[,c("parent_no","taxon_no")]
						pcMat <- rbind(pcMat,newEntries)
						pcAll <- rbind(pcAll,newEntries)
						}
					if(method == "parentChildOldMergeRoot"){
						pcMat <- rbind(pcMat,cbind("ArtificialRoot",floaters))
						taxonNameTable <- rbind(taxonNameTable,c("ArtificialRoot","ArtificialRoot"))
						message(paste0(
							"Multiple potential root-taxa artificially merged at a common root: \n",
							paste0(taxonNameTable[match(floaters,taxonNameTable[,1]),2]
								,collapse = ", ")))
						}
					#regardless of method used, recalculate floaters
					floaters <- getFloatAncPBDB(pcDat = pcMat)
				if(length(floatersNew)>1 & identical(sort(floaters),sort(floatersNew))){
					stop(
						paste0("Provided PBDB Dataset does not appear to have a \n",
							"  monophyletic set of parent-child relationship pairs. \n",
							"Multiple taxa appear to be listed as parents, but are not \n",
							"  listed themselves so have no parents listed: \n",
							paste0(
								taxonNameTable[match(floaters,taxonNameTable[,1]),2],
								collapse = ", "
								)
							)
						)
					}
			}else{
				floaters <- floatersNew
				}
			}
		pcMat <- apply(pcMat,2,as.character)
		#
		# 
		tree <- parentChild2taxonTree(parentChild = pcMat,
			tipSet = tipSet,
			cleanTree = cleanTree)
		#convert tip.label and node.label to taxon names from taxonNameTable
		tree$tip.label <- taxonNameTable[match(tree$tip.label,
			taxonNameTable[,1]),2]
		tree$node.label <- taxonNameTable[match(tree$node.label, 
			taxonNameTable[,1]),2]
		tree$parentChild <- pcMat
		}
	####################
	tree$taxaDataPBDB <- taxaDataPBDB
	return(tree)
	}

	

	
	
#' @rdname makePBDBtaxonTree
#' @export	
plotTaxaTreePBDB<-function(taxaTree, edgeLength = 1){
	taxaTree$edge.length <- rep(edgeLength ,Nedge(taxaTree))
	taxaTree$edge.length [taxaTree$edge[,2] <= Ntip(taxaTree)] <- edgeLength/5
	taxaTree$root.edge <- edgeLength*1.5
	plot(taxaTree,
		show.tip.label = FALSE,
    		no.margin = TRUE,
		root.edge=TRUE,
		edge.width = 0.2)
	nodelabels(taxaTree$node.label,
		cex=0.5, 
		adj = c(1.1,0.5))
	}
	
	
	
#hidden function, don't import
translatePBDBtaxa <- function(taxaDataPBDB){
	# Do some translation
	#need to replace any empty string values with NAs (due perhaps to use of read.csv with the API)
	taxaDataPBDB[taxaDataPBDB == ""] <- NA
	#if com vocab
	if(any("rnk" == colnames(taxaDataPBDB))){	
		#apparently it doesn't matter if these columns *are* present or not
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "acn"] <- "accepted_name"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "snp"] <- "senpar_no"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "rnk"] <- "taxon_rank"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "nam"] <- "taxon_name"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "fml"] <- "family"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "odl"] <- "order"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "cll"] <- "class"	
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "phl"] <- "phylum"	
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "kgl"] <- "kingdom"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "par"] <- "parent_no"
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "oid"] <- "taxon_no"
		# taxon rank translation vectors for compact vocab
		taxRankPBDB <- getTaxRankPBDB()
		taxRankCOM <- 2:26
		#change contents of "identified_rank" and "accepted_rank"
		taxaDataPBDB$taxon_rank <- sapply(taxaDataPBDB$taxon_rank,function(x)
			taxRankPBDB[x == taxRankCOM])
		message("compact vocab detected, relevant fields will be translated")
		}
	###########
	# following are closet cases that mostly only apply to OLD API calls
	#
	if(any(colnames(taxaDataPBDB) == "rank")){
		#if 1.1 and vocab is pbdb
		colnames(taxaDataPBDB)[colnames(taxaDataPBDB) == "rank"] <- "taxon_rank"
		}	
	#
	if(any(colnames(taxaDataPBDB) == "accepted_name")){
		#if 1.2 and there is an accepted_name column..
			#fill empty accepted_name values with taxon_name
		nameFormal <- taxaDataPBDB[,"accepted_name"]
		nameFormal[is.na(nameFormal)] <- as.character(taxaDataPBDB[is.na(nameFormal),"taxon_name"])
		#replace taxon_name
		taxaDataPBDB[,"taxon_name"] <- nameFormal
		#
		#replace taxon_no with accepted_no
		taxNum <- taxaDataPBDB[,"accepted_no"]
		taxNum[is.na(taxNum)] <- as.character(taxaDataPBDB[is.na(taxNum),"taxon_no"])	
		taxaDataPBDB[,"taxon_no"] <- taxNum			
		#
		}
	if(any(colnames(taxaDataPBDB)=="senpar_no")){
		#if this is OLD v1.2 taxaDataPBDB, and there is a senpar_no column
			# replace parent_no in the same way with senpar_no
		parNum <- taxaDataPBDB[,"senpar_no"]
		parNum[is.na(parNum)] <- as.character(taxaDataPBDB[is.na(parNum),"parent_no"])	
		taxaDataPBDB[,"parent_no"] <- parNum
		}
	#
	return(taxaDataPBDB)
	}

getTaxRankPBDB<-function(){
	return(c("subspecies","species","subgenus","genus",
		"subtribe","tribe","subfamily",
		"family","superfamily","infraorder",
		"suborder","order","superorder","infraclass",
		"subclass","class","superclass","subphylum",
		"phylum","superphylum","subkingdom",
		"kingdom","unranked clade","informal"
		#keep informal as high, never drop!
		))
	}

convertParentChildMatNames <- function(taxID, parData){
	taxMatch <- match(taxID, parData$taxon_no)
	if(is.na(taxMatch)){
		newName <- "NODE"
	}else{
		newName <- parData$taxon_name[taxMatch]
		newName <- as.character(newName)
		}
	return(newName)
	}
	
queryMissingParents <- function(taxaID,
		APIversion = "1.2",
		status = "all"){
	# find missing parents by access API
	#
	# drop Eukarya, as it won't return if status = senior under 1.1
	# taxaID <- as.numeric(taxaID[taxaID != "1"])
	#
	# let's get some taxonomic data
	floatData <- read.csv(
		paste0("http://paleobiodb.org/data",APIversion,
			"/taxa/list.txt?taxon_id=",paste0(taxaID,collapse=","),
				### should we take all or only ACCEPTED parents?
			 "&rel=exact&status=",status,
			 "&vocab=pbdb"
			),
		stringsAsFactors = FALSE)
	if(nrow(floatData) == 0){
		stop(
			paste("Current PBDB API would not return info for the following taxon IDs: \n",
				paste0(taxaID,collapse = ", ")
				)
			)
		}
	floatData <- translatePBDBtaxa(floatData)
	#parse down to just taxon_name, taxon_no, parent_no
	floatData <- parseParentPBDBData(floatData)
	return(floatData)
	}


parseParentPBDBData <- function(parentData){
	#parse down to just taxon_name, taxon_no, parent_no
	if(any(colnames(parentData)=="accepted_name")){
		taxon_name = parentData$accepted_name
	}else{
		taxon_name = parentData$taxon_name
		}
	#
	result <- data.frame(
		taxon_name = as.character(taxon_name),
		parent_no = as.numeric(parentData$parent_no),
		taxon_no = as.numeric(parentData$taxon_no)
		)
	rownames(result) <- NULL
	return(result)
	}

getAllParents<-function(inputData, status, annotatedDuplicateNames = TRUE){
	parData<-parseParentPBDBData(inputData)
	noParentMatch<-findNoParentMatch(parData)
	while(sum(noParentMatch)>1){
		floatingParentNum <- unique(parData$parent_no[noParentMatch])
		dataNew<-queryMissingParents(floatingParentNum, status=status)
		parData<-rbind(parData,dataNew)
		noParentMatch<-findNoParentMatch(parData)
		}
	# TESTS
	if(sum(noParentMatch)!=1){
		stop("Cannot find a single common ancestor by tracing parents")
		}
	if(nrow(parData) != nrow(unique(parData))){
		print(parData[duplicated(parData),])
		stop("getAllParents added duplicate parent-child relationships, see print-out above")
		}
	if(length(parData$taxon_no) != length(unique(parData$taxon_no))){
		print(parData[duplicated(parData$taxon_no),])
		stop("getAllParents added duplicate child taxa, see print-out above")
		}
	if(length(parData$taxon_name) != length(unique(parData$taxon_name))){
		if(annotatedDuplicateNames){
			parData$taxon_name <- make.unique(as.character(parData$taxon_name))
		}else{
			print(parData[duplicated(parData$taxon_name),])
			stop("getAllParents added duplicate accepted taxon names, see print-out above")
			}
		}
	return(parData)
	}


findNoParentMatch<-function(parData){
	res <- is.na(match(parData$parent_no, parData$taxon_no))
	return(res)
	}

	
getFloatAncPBDB <- function(pcDat){
	#identify IDs of parents floating without ancestors of their own
	# which parents have no children?
	noChildren <- sapply(pcDat$parent_no, function(x) 
			all(x != pcDat$child_no)
			)
	# list all unique no-child parents
	res <- unique(pcDat$parent_no[noChildren])
	return(res)
	}	
	
constructParentChildMatrixPBDB <- function(initPCmat, parData){
	# starting from desired tip OTUs, work backwards to a common ancestor from the full parData
	if(nrow(initPCmat) != nrow(unique(initPCmat))){
		print(initPCmat[duplicated(initPCmat),])
		stop("initial parent-child matrix had duplicate parent-child relationships, see print-out above")
		}
	newPCmat <- initPCmat
	# find floating parents in current newPCmat
	#identify IDs of parents floating without ancestors of their own
	floaters <- getFloatAncPBDB(pcDat = newPCmat)
	# use a while loop to complete the parent-child matrix
	while(length(floaters)>1){	#so only ONE root can float
		# get new relations: will 'anchor' the floaters
		anchorMat <- subsetParDataPBDB(subsetNum = floaters, parData = parData)
		# bind to newPCmat
		newPCmat <- rbind(newPCmat,anchorMat)
		if(nrow(newPCmat) != nrow(unique(newPCmat))){
			print(newPCmat[duplicated(newPCmat),])
			stop("annotated pcMat had duplicate parent-child relationships, see print-out above")
			}
		# recalculate floater taxa
		floatersNew <- getFloatAncPBDB(pcDat = newPCmat)	#recalculate float
		#
		# put in a stopping condition for the love of god
		if( length(floatersNew)>1 & identical(floaters,floatersNew) ){
			stop(paste0(
				"Provided PBDB Dataset does not appear to have a \n",
				"  monophyletic set of parent-child relationship pairs. \n",
				"Multiple taxa appear to be listed as parents, but are not \n",
				"  listed themselves so have no parents listed: \n",
				paste0(floaters,
					collapse = ", "
					),
				paste0(
					parData$taxon_name[match(floaters,
						parData$taxon_no)],
					collapse = ", "
					)
				))	
		}else{
			floaters <- floatersNew
			}
		}
	return(newPCmat)
	}	
	
subsetParDataPBDB <- function(subsetNum,parData){
	# pull parent-child relationships for a set of taxon numbers from parData
	subsetRows <- match(subsetNum, parData$taxon_no)
	# remove the non matching subset (usually the root)
	subsetRows <- subsetRows[!is.na(subsetRows)]
	# pull the subset of parent-child relationships from parData
	#subsetMat <- parData[subsetRows, c("parent_no", "taxon_no")]
	subsetMat <- data.frame(
		parent_no = parData$parent_no[subsetRows],
		child_no = parData$taxon_no[subsetRows]
		)
	rownames(subsetMat)<-as.character(parData$taxon_no[subsetRows])
	return(subsetMat)
	}
	
	
#getTaxaIDsDesiredRank<-function(data, rank){
#	# filter out lower than selected rank (for tip taxa)
#		# so need to know which ranks are lower/higher
#	# get taxon rank translation vectors for compact vocab
#	taxRankPBDB <- getTaxRankPBDB()
#	# translate rank to a number
#	# translate taxon_rank to a number
#	numTaxonRank <- sapply(data[,"taxon_rank"],
#		function(x) which(x == taxRankPBDB))		
#	#now need to put together parentChild table
#	# get taxon ID numbers of just those of desired rank
#	desiredIDs <- data[rankID == numTaxonRank, "parent_no"]
#	desiredIDs <- as.numeric(desiredIDs)
#	return(desiredIDs)
#	}


	