#' Get the Sample of Posterior Trees from a Dated Phylogenetic Analysis with MrBayes (Or a Summary Tree, such as the MCCT)
#' 
#' MrBayes is not great for getting samples of dated posterior phylogenies, or for obtaining certain summary trees from
#' the posterior (specifically the MCCT and MAP, which are specific trees in the posterior). This is because the tree
#' samples as returned are scaled relative to rate parameters in a separate file. This function attempts to automate
#' the handling of multiple files (both .t tree files and .p parameter files), as well as multiple files
#' associated with separate runs, to obtain samples of posterior trees, or summary trees such as the MCCT or MAP.
#' These resulting trees are now scaled to units of time, but not be placed correctly on an absolute time-scale
#' if all tips are extinct. See details of output below.

#' @details
#' This function is most useful for dealing with dating analyses in MrBayes, particularly when tip-dating
#' a tree with fossil taxa, as the half-compatibility and all-compatibility summary trees offered by the
#' 'sumt' command in MrBayes can have issues properly portraying summary trees from such datasets.

#' @param runFile A filename in the current directory, or a path to a file that is either a .p 
#' or .t file from a MrBayes analysis. This filename and path will be used for finding additional 
#' .t and .p files, via the \code{nRuns} settings and assuming that files are in the
#' same directory \emph{and} these files are named under
#' typical MrBayes file naming conventions. (In other words, if you have renamed your .p or .t files,
#' this function probably won't be able to find them.)

#' @param nRuns The number of runs in your analysis. This variable is used for figuring out what 
#' filenames will be searched for: if you specify that you have less runs than you
#' actually ran in reality, then some runs won't be examined in thi function. Conversely,
#' specify too many, and this function will throw an error when it cannot find files it expects
#' but do not exist. The default for this argument (\emph{two} runs) is based on the default number of runs in MrBayes.

#' @param burnin The fraction of trees sampled in the posterior discarded  and not returned
#' by this function directly, nor included in calculation of summary trees. Must be a numeric
#' value greater than 0 and less than 1.

#' @param getFixedTimes If \code{TRUE}, this function will also look for, scan, and parse an
#' associated NEXUS file. Ignoring any commented lines (ie. anything between "[   ]" ), commands
#' for fixing taxa will be identified, parsed and returned to the user, either as a message
#' printed to the R console if output is read to a file, or as a attribute named 'fixed ages'
#' if output as an R object (formatted as a two-column table of OTU names and their respective fixed ages).
#' If the output is an R object, these objects with 
#' 
#' Please note: the code for \code{getFixedTimes = TRUE} contains a \code{while()}
#' loop in it for removing nested series of
#' square brackets (i.e. treated as comments in NEXUS files). Thus files with
#' ridiculously nested series of brackets may cause this code to take a while
#' to complete, or may even cause it to hang.

#' @param originalNexusFile Filename (and possibly path too) to the original NEXUS file for this analysis.
#' Only tried if \code{getFixedTimes = TRUE}. If \code{NULL} (the default), then this function will
#' instead try to find a NEXUS file with the same name as implied by the filename used in other inputs. If
#' this file cannot be found, the function will fail. 

#' @param labelPostProb Logical. If \code{TRUE}, then nodes of the
#' output tree will be labeled with their respective posterior 
#' probabilities, as calculated based on the frequency of a clade
#' occurring across the post-burnin posterior tree sample. If \code{FALSE},
#' this is skipped.
	
#' @param outputTrees Determines the output trees produced; for format of output, see section
#' on returned Value below. Must be of length one, and either \code{"all"},
#' which means all trees from the post-burnin posterior will
#' returned, a number greater than zero, which will be the number of trees
#' randomly sampled from across the post-burning posterior and returned,
#' or \code{"MCCT"} and \code{"MAP"}, which stand for 'maximum clade compatibility tree'
#' and 'maximum a posteri tree' respectively. The MAP is the single tree from the post-burnin
#' posterior with the highest marginal likelihood. The MCCT is the single tree from the post-burnin
#' posterior which contains clades with the highest product of posterior probabilities for its
#' component clades. Thus, the MAP is the best overall tree, while the MCCT may be the best
#' tree for summarizing topological support.

#' @param file Filename (possibly with path) as a character string
#' leading to a file which will be overwritten with the output trees (or summary tree),
#' as a NEXUS file. If \code{NULL} (the default), the output will
#' instead be directly returned by this function.

#' @return
#' Depending on argument \code{file}, the output tree or trees is either
#' returned directly, or instead written out in NEXUS format via
#' ape's \code{write.NEXUS} function to an external file. The output
#' will consist either of multiple trees sampled from the post-burnin posterior,
#' or will consist of a single phylogeny (a summary tree, either the MCCT or the MAP).
#'' 
#' Users are warned that the resulting dated trees will not have \code{$root.time} elements necessary
#' for comparison against an absolute time-scale. Wile the trees may be scaled to units of absolute
#' time now, rather than with branch lengths expressed in the rate of character change, the dates
#' estimated by some phylogenetics functions in R may give inaccurate estimates of when events
#' occur on the absolute time-scale if all tips are extinct.
#' This is because most functions for phylogenetics in R (and
#' elsewhere) will instead presume that the latest tip will be at time 0 (the modern), which
#' may be wrong if you are using \code{paleotree} for analyzing \emph{paleontological} datasets
#' consisting of entirely extinct taxa. This can be solved by using argument \code{getFixedTimes = TRUE}
#' to obtain fixed tip ages, and then scaling the resulting output to absolute time using
#' the function \code{\link{setRootAges}} to obtain a \code{$root.time} element for each tree.
#'


#' @seealso
#' When the argument \code{getFixedTimes = TRUE} is used, the resulting output can easily be scaled to absolute time 
#' if you used a taxon with a fixed age, and function \code{\link{setRootAges}}.
#' 
#' Maximum Clade Credibility trees are estimated using the function\code{\link[phangorn]{maxCladeCred}} in package phangorn.

#' @author
#' David Bapst, with rescaling of raw output
#' trees via code originally written by Nicholas Crouch.

#' @examples
#' \dontrun{
#' 
#' MCCT <- obtainDatedPosteriorTreesMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, burnin = 0.5,
#' 		outputTrees = "MCCT", file = NULL)
#'
#' MAP <- obtainDatedPosteriorTreesMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, burnin = 0.5, getFixedTimes = TRUE,
#' 		outputTrees = "MAP", file = NULL)
#'
#' # get a root age from the fixed ages for tips
#' setRootAge(tree = MAP)
#' 
#' #pull a hundred trees randomly from the posterior
#' hundredRandomlySelectedTrees <- obtainDatedPosteriorTreesMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, burnin = 0.5, getFixedTimes = TRUE,
#' 		outputTrees = 100, file = NULL)
#'
# # set the root age using the fixed ages
#' setRootAges(trees = hundredRandomlySelectedTrees)
#' 
#' }
#'



#' @name obtainDatedPosteriorTreesMrB
#' @rdname obtainDatedPosteriorTreesMrB
#' @export
obtainDatedPosteriorTreesMrB <- function(runFile,nRuns = 2,burnin = 0.5,
	outputTrees,labelPostProb = FALSE,
	getFixedTimes = FALSE,originalNexusFile = NULL,
	file = NULL){
	#checks
	if(length(outputTrees) != 1){
		stop("outputTrees must be of length 1")
		}
	if(all(outputTrees != c("all","MCCT","MAP")) & !is.numeric(outputTrees)){
		stop("outputTrees must be one of 'all', 'MCCT', or a numeric value indicating the number of trees to randomly sample from the posterior")
		}
	if(is.numeric(outputTrees)){	
		if(outputTrees<1){
			stop("If numeric, outputTrees must be greater than 0")
			}
		}
	if(burnin>1 | burnin<0){
		stop("burnin must be a value between 0 and 1")
		}
	# Load tree, paramater and mcmc files produced by MrBayes
	#
	# # pick either a .p or .t file; can be either. Used for getting path. 
	# all needed .p and .t files must be in the same directory as this file
	#
	# take indicated file and get the basic run name
	runPath <- strsplit(runFile,split = "\\.run")[[1]]
	if(length(runPath) != 2){
		stop("Unable to parse the runPath correctly")
		}
	runPath <- runPath[[1]]
	# get list of tree files
	treeFiles <- lapply(1:nRuns,function(x){
		read.nexus(file = paste0(runPath,".run",x,".t"))}
		)
	# get list of .p files
	parFiles <- lapply(1:nRuns,function(x)
		read.table(file = paste0(runPath,".run",x,".p"),
			 header = T, skip = 1)
		)
	#############################
	# checks for length 
	for(i in 1:nRuns){
		if(length(treeFiles[[i]]) != nrow(parFiles[[i]])){
			stop("Parameter data and tree data not of the same length")
			}
		}
	####################################
	#	
	# @param getFixedTimes If \code{TRUE}, this function will also look for, scan, and parse an
	# associated NEXUS file. Ignoring any commented lines (ie. anything between "[   ]" ), commands
	# for fixing taxa will be identifiedd, parsed and returned to the user, either as a message
	# printed to the R console if output is read to a file, or as a attribute named 'fixed ages'
	# if output as an R object (formatted as a two-column table of OTU names and their respective fixed ages).
	# The search for the NEXUS file is controlled with argument \code{originalNexusFile}
	#
	# Please note: this has a while() loop in it for removing nested series of
	# square brackets (i.e. treated as comments in NEXUS files) then files with
	# extremely nested series of brackets may cause this code to take a while
	# to complete, or may even cause it to hang.
	#
	if(getFixedTimes){
		if(is.null(originalNexusFile)){
			# file name, if presuming shares run name with 
			#originalNexusFile <- paste0(runPath,".nex")
			originalNexusFile <- runPath
			}
		fixedTable <- getMrBFixedAgesFromNexus(originalNexusFile)
	}else{
		fixedTable <- NULL
		}
	####################
	#
	# Specify sample start based on burnin - here 50%
	sampleStart <- floor(length(treeFiles[[1]])*burnin)
	#
	# remove burnin from both
	treeFilesBurn <- lapply(treeFiles,function(x){
		startSamp <- floor(length(x)*burnin)
		x[startSamp:length(x)]
		}
	)
	parFilesBurn <- lapply(parFiles,function(x){
		startSamp <- floor(nrow(x)*burnin)
		x[startSamp:nrow(x),]
		}
	)
	# checks for length again
	for(i in 1:nRuns){
		if(length(treeFilesBurn[[i]]) != nrow(parFilesBurn[[i]])){
			stop("Parameter data and tree data not of the same length")
			}
		}
	###
	#
	# The branch lengths of the trees need to be scaled
		# by the estimated clock rate in the par files
	#
	# Function to rescale branch lengths of a phylogeny by clockrate
	rescaleMrBTree <- function(phy, rate){
		phy$edge.length <- phy$edge.length / rate$clockrate
		return(phy)
		}
	# now apply it
	rescaledTrees <- lapply(1:nRuns,function(x){ 
		lapply(1:length(treeFilesBurn[[x]]),function(y){
			rescaleMrBTree(treeFilesBurn[[x]][[y]],parFilesBurn[[x]][y,])
			})
		})
	#################################
	#
	# attach marginal likelihoods to each tree
	for(i in 1:nRuns){
		for(j in 1:length(rescaledTrees[[i]])){
			rescaledTrees[[i]][[j]]$LnPr <- parFilesBurn[[i]][j,]$LnPr
			}
		}
	# concatanate trees from each run
	lumpTrees <- unlist(rescaledTrees,recursive = FALSE)
	class(lumpTrees) <- "multiPhylo"
	#
	##########################################
	#
	if(outputTrees == "all"){
		outTree <- lumpTrees
		if(labelPostProb){
			# assign posterior probabilities as node labels
			outTree <- lapply(outTree,getPosteriorProbabiities,postBurninTrees = lumpTrees)
			}
		}
	#
	if(outputTrees == "MAP"){
		# get MAP
		LnPr <- sapply(lumpTrees,function(x) x$LnPr)
		whichMAP <- which(LnPr == max(LnPr))
		outTree <- lumpTrees[[whichMAP]]
		if(labelPostProb){
			# assign posterior probabilities as node labels
			outTree <- getPosteriorProbabiities(tree = outTree,postBurninTrees = lumpTrees)
			}
		}
	#
	if(outputTrees == "MCCT"){
		# turns out the MCCT isn't the tree with the highest likelihood
		outTree <- phangorn::maxCladeCred(lumpTrees)
		if(labelPostProb){
			# assign posterior probabilities as node labels
			outTree <- getPosteriorProbabiities(tree = outTree,postBurninTrees = lumpTrees)
			}
		}
	#
	if(is.numeric(outputTrees)){
		# randomly sample N trees from the posterior
		# first check
		if(outputTrees>length(lumpTrees)){
			stop(paste0("Numeric value for outputTrees (",outputTrees,
				") greater than the total number of post-burning trees (",
				length(lumpTrees),")"))
			}
		whichOutput <- sample(length(lumpTrees),outputTrees,replace = FALSE)
		outTree <- lumpTrees[whichOutput]
		if(labelPostProb){
			# assign posterior probabilities as node labels
			if(outputTrees>1){
				outTree <- lapply(outTree,getPosteriorProbabiities,postBurninTrees = lumpTrees)
			}else{
				outTree <- getPosteriorProbabiities(tree = outTree,postBurninTrees = lumpTrees)
				}
			}
		}
	########################################
	if(!is.null(file)){
		write.nexus(outTree, file = file)
		if(!is.null(fixedTable)){
			print("You may need to know which tip ages were fixed to a precise date")
			print("in order to assign absolute dates to the tree, as follows:")
			print(fixedTable)
			}
	}else{
		if(!is.null(fixedTable)){
			attr(outTree,"fixedTable") <- fixedTable
			}
		return(outTree)
		}
	}

getPosteriorProbabiities <- function(tree,postBurninTrees){
	cladeFreq <- prop.clades(tree,postBurninTrees)
	postProbs <- cladeFreq/length(postBurninTrees)
	tree$node.label <- postProbs
	return(tree)
	}

