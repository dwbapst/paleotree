#' Creating a Taxon-Tree from Taxonomic Data Downloaded from the Paleobiology Database
#'
#' This function creates phylogeny-like object of type
#' \code{phylo} from the taxonomic information
#' recorded in a taxonomy download from the PBDB for
#' a given group. Two different algorithms are provided,
#' the default being based on parent-child taxon relationships,
#' the other based on the nested Linnean hierarchy.

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

#' @param data A table of taxonomic data collected from
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
#   dataset under \code{data}? Each of these is essentially a separate root taxon, for a different set
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
#' Two other functions in paleotree are used as sub-algorithms by \code{makePBDBtaxonTree}
#' to create the taxon-tree within this function,
#' and users should consult their manual pages for additional details:
#'
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}
#'
#' Other functions for manipulating PBDB data can be found at \code{\link{taxonSortPBDBocc}},
#' \code{\link{occData2timeList}}, and the example data at \code{\link{graptPBDB}}.

#' @author David W. Bapst

#' @references
#' Soul, L. C., and M. Friedman. 2015. Taxonomy and Phylogeny Can Yield
#' Comparable Results in Comparative Palaeontological Analyses. \emph{Systematic Biology} 
#' (\href{http://sysbio.oxfordjournals.org/content/early/2015/03/23/sysbio.syv015.abstract}{Link})
#' 

#' @examples
#' \dontrun{
#' 
#' easyGetPBDBtaxa <- function(taxon,
#' 		show = c("class","img","app","parent"),
#' 		status = "accepted"){
#' 	##########################################
#' 	# 12-30-18: modified for API version 1.2
#' 	#let's get some taxonomic data
#' 	taxaData <- read.csv(
#' 		paste0("http://paleobiodb.org/",
#' 			"data1.2/taxa/list.txt?base_name=", taxon,
#' 			"&show=",paste0(show,collapse = ","),
#' 			# status -> all, accepted, valid
#' 				# accepted -> only senior synonyms
#' 				# valid -> snior synonyms + valid subjective synonyms
#' 				# all -> valid taxa + repressed invalid taxa
#' 			"&taxon_status=",status
#' 			),
#' 	    stringsAsFactors = FALSE)
#' 	#######################################
#' 	return(taxaData)
#' 	}
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
makePBDBtaxonTree <- function(data, rank,
					method = "parentChild", #solveMissing = NULL,
					tipSet = "nonParents", cleanTree = TRUE,
					APIversion = "1.2"){		
	############################################################
	############################################################
	# library(paleotree);data(graptPBDB);
	# data <- graptTaxaPBDB; rank = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
	# data <- graptTaxaPBDB; rank = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
	# data <- graptTaxaPBDB; rank = "genus"; method = "parentChild"; tipSet = "nonParents"; cleanTree = TRUE
	# data <- graptTaxaPBDB; rank = "genus"; method = "Linnean"; 
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
	if(!is(data,"data.frame")){
		stop("data isn't a data.frame")
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
	dataTransform <- translatePBDBtaxa(data)
	dataTransform <- apply(dataTransform, 2, as.character)
	#
	if(method == "parentChild"){
		# need two things: a table of parent-child relationships as IDs
			#and a look-up table of IDs and taxon names
		# 
		##############################
		# FIND ALL PARENTS FIRST
			# three column matrix with taxon name, taxon ID, parent ID
			# (in that order)
		parData<- getAllParents(dataTransform)
		#
		#######################################
		# NOW FILTER OUT TIP TAXA WE WANT
		tipIDs <- getTaxaIDsDesiredRank(data=dataTransform, rank=rank)
		# figure out which taxon numbers match tip IDs
		whichTip <- match(tipIDs, parData$taxon_no)
		#
		###############################
		# BUILD PARENT-CHILD matrix
		# get parent-child matrix for just desired OTUs 
			# start matrix with those parent-child relationships
		pcMat <- parData[whichTip, c("parent_no","taxon_no")]
		# find floating parents in current pcMat
		#identify IDs of parents floating without ancestors of their own
		floaters <- getFloatAncPBDB(pcDat = pcMat)
		# use a while loop to complete the parent-child matrix
		while(length(floaters)>1){	#so only ONE root can float
			# get new relations: will 'anchor' the floaters
			anchors <- match(floaters,parData[,2])
			anchorMat <- parData[anchors, c("parent_no","taxon_no")]
			# bind to pcMat
			pcMat <- rbind(pcMat,anchorMat)
			# recalculate floater taxa
			floatersNew <- getFloatAncPBDB(pcDat = pcMat)	#recalculate float
			if( identical(floaters,floatersNew) ){
				stop(paste0(
					"Provided PBDB Dataset does not appear to have a \n",
					"  monophyletic set of parent-child relationship pairs. \n",
					"Multiple taxa appear to be listed as parents, but are not \n",
					"  listed themselves so have no parents listed: \n",
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
		#################################
		# convert parent-child matrix to accepted taxon names
		pcMat <- apply(pcMat, 2, function(x) 
			as.character(parData$taxon_name[x])
			)
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
		#Check if show = class was used
		if(!any(colnames(dataTransform) == "family")){
			stop("Data must be a taxonomic download with show = class for method = 'Linnean'")
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
		#first, get table of all parentChild relationships in data
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
			pcMat <- rbind(pcMat,pcAll[anchors[!is.na(anchors)],])	#bind to pcMat
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
	return(tree)
	}

#hidden function, don't import
translatePBDBtaxa <- function(data){
	# Do some translation
	#need to replace any empty string values with NAs (due perhaps to use of read.csv with the API)
	data[data == ""] <- NA
	#if com vocab
	if(any("rnk" == colnames(data))){	
		#apparently it doesn't matter if these columns *are* present or not
		colnames(data)[colnames(data) == "acn"] <- "accepted_name"
		colnames(data)[colnames(data) == "snp"] <- "senpar_no"
		colnames(data)[colnames(data) == "rnk"] <- "taxon_rank"
		colnames(data)[colnames(data) == "nam"] <- "taxon_name"
		colnames(data)[colnames(data) == "fml"] <- "family"
		colnames(data)[colnames(data) == "odl"] <- "order"
		colnames(data)[colnames(data) == "cll"] <- "class"	
		colnames(data)[colnames(data) == "phl"] <- "phylum"	
		colnames(data)[colnames(data) == "kgl"] <- "kingdom"
		colnames(data)[colnames(data) == "par"] <- "parent_no"
		colnames(data)[colnames(data) == "oid"] <- "taxon_no"
		# taxon rank translation vectors for compact vocab
		taxRankPBDB <- getTaxRankPBDB()
		taxRankCOM <- 2:26
		#change contents of "identified_rank" and "accepted_rank"
		data$taxon_rank <- sapply(data$taxon_rank,function(x) taxRankPBDB[x == taxRankCOM])
		message("compact vocab detected, relevant fields will be translated")
		}
	###########
	# following are closet cases that mostly only apply to OLD API calls
	#
	if(any(colnames(data) == "rank")){
		#if 1.1 and vocab is pbdb
		colnames(data)[colnames(data) == "rank"] <- "taxon_rank"
		}	
	#
	if(any(colnames(data) == "accepted_name")){
		#if 1.2 and there is an accepted_name column..
			#fill empty accepted_name values with taxon_name
		nameFormal <- data[,"accepted_name"]
		nameFormal[is.na(nameFormal)] <- as.character(data[is.na(nameFormal),"taxon_name"])
		#replace taxon_name
		data[,"taxon_name"] <- nameFormal
		#
		#replace taxon_no with accepted_no
		taxNum <- data[,"accepted_no"]
		taxNum[is.na(taxNum)] <- as.character(data[is.na(taxNum),"taxon_no"])	
		data[,"taxon_no"] <- taxNum			
		#
		}
	if(any(colnames(data)=="senpar_no")){
		#if this is OLD v1.2 data, and there is a senpar_no column
			# replace parent_no in the same way with senpar_no
		parNum <- data[,"senpar_no"]
		parNum[is.na(parNum)] <- as.character(data[is.na(parNum),"parent_no"])	
		data[,"parent_no"] <- parNum
		}
	#
	return(data)
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
		paste0("http://paleobiodb.org/",
			"data",APIversion,
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
	result <- data.frame(
		if(any(colnames(parentData)=="accepted_name")){
			"taxon_name" = as.character(parentData[,"accepted_name"])
		}else{
			"taxon_name" = as.character(parentData[,"taxon_name"])
			}
		,
		"parent_no" = as.numeric(parentData[,"parent_no"]),
		"taxon_no" = as.numeric(parentData[,"taxon_no"]))
	return(result)
	}



getFloatAncPBDB <- function(pcDat){
	#identify IDs of parents floating without ancestors of their own
	res <- unique(pcDat[
		sapply(pcDat[,1], function(x) 
			all(x != pcDat[,2])
			)
		,1])
	return(res)
	}


getAllParents<-function(inputData){
	parData<-parseParentPBDBData(inputData)
	noParentMatch<-findNoParentMatch(parData)
	while(sum(noParentMatch)>1){
		floatingParentNum <- unique(parData$parent_no[noParentMatch])
		dataNew<-queryMissingParents(floatingParentNum)
		parData<-rbind(parData,dataNew)
		noParentMatch<-findNoParentMatch(parData)
		}
	if(sum(noParentMatch)!=1){
		stop("Cannot find a single common ancestor by tracing parents")
		}
	return(parData)
	}


findNoParentMatch<-function(parData){
	res <- is.na(match(parData$parent_no, parData$taxon_no))
	return(res)
	}


getTaxaIDsDesiredRank<-function(data,rank){
	# filter out lower than selected rank (for tip taxa)
		# so need to know which ranks are lower/higher
	# get taxon rank translation vectors for compact vocab
	taxRankPBDB <- getTaxRankPBDB()
	# translate rank to a number
	rankID <- which(rank == taxRankPBDB)
	# translate taxon_rank to a number
	numTaxonRank <- sapply(data[,"taxon_rank"],
		function(x) which(x == taxRankPBDB))		
	#now need to put together parentChild table
	# get taxon ID numbers of just those of desired rank
	desiredIDs <- data[rankID == numTaxonRank, "parent_no"]
	desiredIDs <- as.numeric(desiredIDs)
	return(desiredIDs)
	}


	