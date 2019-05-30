#' Get the Sample of Posterior Trees from a Dated Phylogenetic Analysis
#' with MrBayes (Or a Summary Tree, such as the MCCT)
#' 
#' MrBayes is not great for getting samples of dated posterior
#' phylogenies, or for obtaining certain summary trees from
#' the posterior (specifically the MCCT and MAP, which are specific
#' trees in the posterior). This is because the tree
#' samples as returned are scaled relative to rate parameters in a
#' separate file. This function attempts to automate
#' the handling of multiple files (both \emph{.t} tree files and \emph{.p}
#' parameter files), as well as multiple files
#' associated with separate runs, to obtain samples of posterior
#' trees, or summary trees such as the MCCT or MAP.
#' These resulting trees are now scaled to units of time, but
#' not be placed correctly on an absolute time-scale
#' if all tips are extinct. See details of output below.

#' @details
#' This function is most useful for dealing with dating
#' analyses in MrBayes, particularly when tip-dating
#' a tree with fossil taxa, as the half-compatibility
#' and all-compatibility summary trees offered by the
#' '\code{sumt}' command in MrBayes can have issues properly
#' portraying summary trees from such datasets.

#' @param runFile A filename in the current directory,
#' or a path to a file that is either a \emph{.p} 
#' or \emph{.t} file from a MrBayes analysis. This filename
#' and path will be used for finding additional 
#' \emph{.t} and \emph{.p} files, via the \code{nRuns}
#' settings and assuming that files are in the
#' same directory \emph{and} these files are named under
#' typical MrBayes file naming conventions. (In other words,
#' if you have renamed your \emph{.p} or \emph{.t} files,
#' this function probably won't be able to find them.)

#' @param nRuns The number of runs in your analysis. This variable is used for figuring out what 
#' filenames will be searched for: if you specify that you have less runs than you
#' actually ran in reality, then some runs won't be examined in thi function. Conversely,
#' specify too many, and this function will throw an error when it cannot find files it expects
#' but do not exist. The default for this argument
#' (\emph{two} runs) is based on the default number of runs in MrBayes.

#' @param burnin The fraction of trees sampled in the posterior discarded  and not returned
#' by this function directly, nor included in calculation of summary trees. Must be a numeric
#' value greater than 0 and less than 1.

#' @param getFixedTimes If \code{TRUE}, this function will also look for, scan, and parse an
#' associated NEXUS file. Ignoring any commented lines (ie. anything between "[   ]" ), commands
#' for fixing taxa will be identified, parsed and returned to the user, either as a message
#' printed to the R console if output is read to a file, or as a attribute named 'fixed ages'
#' if output as an R object (formatted as a two-column
#' table of OTU names and their respective fixed ages).
#' If the output is an R object, these objects with 
#' 
#' Please note: the code for \code{getFixedTimes = TRUE} contains a \code{while()}
#' loop in it for removing nested series of
#' square brackets (i.e. treated as comments in NEXUS files). Thus files with
#' ridiculously nested series of brackets may cause this code to take a while
#' to complete, or may even cause it to hang.

#' @param getRootAges \code{FALSE} by default. 
#' If \code{TRUE}, and \code{getFixedTimes = TRUE}
#' as well as \code{file = NULL}
#' (such that trees will be assigned within the R memory
#' rather than saved to an external file), the
#' functions \code{setRootAge} and its wrapper function
#' \code{setRootAges} will be applied to the output
#' so that all output trees have \code{root.time}
#' elements for use with other functions in \code{paleotree}
#' as well as other packages.

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
#' or a label for a type of summary tree selected from the posterior based on various
#' properties. The two most commonly seen such point-estimate-summaries are
#' the \emph{MCCT} tree, which stands for the 'maximum clade compatibility tree',
#' and the \emph{MAP} tree, which stands for the 'maximum a posteri tree'. 
#' The MCCT is the single tree from the post-burnin posterior which
#' contains the set of bifurcations (clades) with the highest product of posterior
#' probabilities (i.e. are found on the most trees in the post-burnin posterior).
#' The MCCT tree is returned if the argument \code{outputTrees = "MCCT"} is used.
#' The MAP is the single tree from the post-burnin posterior with the highest
#' posterior probabilty associated with it. Unfortunately, versions of
#' \code{paleotree} prior to version 3.2.1 did not use the posterior probability
#' to select the supposed 'MAP' tree. MrBayes provides two values
#' for each sampled tree and corresponding parameters: 
#' \code{LnPr}, the log prior probability of the current parameter 
#' proposals under the specified prior distributions, 
#' and \code{LnL}, the log likelihood of the cold chain, i.e. the log-likelihood
#' of the sampled parameter values, given the observed data and specified models.
#' Neither of these are the posterior probability. 
#' The true posterior probability (as given by Bayes Theorem) is 
#' the prooduct of the likelihood and the prior probability, divided by
#' the likelihood of the model, the latter of which is very rarely known.
#' More commonly, the calculatable portion of the posterior probability is
#' the product of the likelihood and the prior probability; or, here, easily
#' calculated as the log posterior probability, as the sum of the
#' log likelihood and log prior probability. Given confusion over application of 'MAP'
#' trees in previous version of paleotree, three options are available:
#' \code{"MAPosteriori"} for the Maximum A Posteriori tree (the MAP tree,
#' or the single tree in the posterior with the highest
#' posterior probability, as given by \code{LnPr + LnL}),
#' \code{"MAPriori"} for the Maximum A Priori tree (the tree in the
#' posterior sample with the highest prior probability, indep of
#' the data), and \code{"MaxLikelihood"}, the
#' tree with the highest model likelihood given the data, ignoring the prior probability.
#' The former option of \code{outputTrees = "MAP"} is deprecated,
#' as its previous implementation only examine \code{LnPr} and thus
#' returned the tree now referred to here as the \code{"MAPriori"} tree.
#' Overall, the Maximum A Posteriori tree is the "best" tree based
#' on the metric most directly considered by Bayesian analysis for
#' proposal acceptance, but the MCCT may be
#' the best tree for summarizing topological support.
#' In either case, point estimates of topology are often
#' problematic summaries for phylogenetic analyses.

#' @param file Filename (possibly with path) as a character string
#' leading to a file which will be overwritten with the output trees (or summary tree),
#' as a NEXUS file. If \code{NULL} (the default), the output will
#' instead be directly returned by this function.



#' @return
#' Depending on argument \code{file}, the output tree or trees is either
#' returned directly, or instead written out in NEXUS format via
#' ape's \code{write.NEXUS} function to an external file. The output
#' will consist either of multiple trees sampled from the post-burnin posterior,
#' or will consist of a single phylogeny (a summary tree, either
#' the MCCT or the MAP - see the details for the argument \code{outputTrees}).
#' 
#' If the argument \code{setRootAges = TRUE} is not used,
#' users are warned that the resulting dated trees will
#' not have \code{$root.time} elements necessary
#' for comparison against an absolute time-scale. Wile the
#' trees may be scaled to units of absolute
#' time now, rather than with branch lengths expressed in
#' the rate of character change, the dates
#' estimated by some phylogenetics functions in 
#' R may give inaccurate estimates of when events
#' occur on the absolute time-scale if all tips are extinct.
#' This is because most functions for phylogenetics in R (and
#' elsewhere) will instead presume that the latest tip
#' will be at time 0 (the modern), which
#' may be wrong if you are using \code{paleotree} for
#' analyzing \emph{paleontological} datasets
#' consisting of entirely extinct taxa. This can be
#' solved by using argument \code{getFixedTimes = TRUE}
#' to obtain fixed tip ages, and then scaling the resulting output to absolute time using
#' the argument \code{setRootAges = TRUE}, which obtains
#' a \code{$root.time} element for each tree
#' using the functions \code{\link{setRootAge}} and 
#'\code{\link{setRootAges}} (for single and multiple phylogenies).
#' 


#' @seealso 
#' When the arguments \code{getFixedTimes = TRUE} and 
#' \code{setRootAges = TRUE} are used, the resulting output will be scaled to absolute time 
#' with the available fixed ages using functions \code{\link{setRootAge}}
#' and \code{\link{setRootAges}} (for single and multiple phylogenies).
#' This is only done if fixed ages are available and if the tree is not
#' being saved to an external file.
#' 
#' Maximum Clade Credibility trees are estimated using the function
#' \code{\link[phangorn]{maxCladeCred}} in package phangorn.
#' 
#' See function \code{link{tipDatingCompatabilitySummaryMrB}} for additional
#' ways of solely evaluating the topoligical information
#' in trees taken from MrBayes posterior samples.

#' @author 
#' David Bapst, with rescaling of raw output
#' trees via code originally written by Nicholas Crouch.

#' @examples
#' \dontrun{
#' 
#' MCCT <- obtainDatedPosteriorTreesMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, 
#'  	burnin = 0.5,
#'  	outputTrees = "MCCT", 
#'  	file = NULL)
#' 
#' MAP <- obtainDatedPosteriorTreesMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, 
#'  	burnin = 0.5, 
#'  	getFixedTimes = TRUE,
#'  	outputTrees = "MAPosteriori", 
#'  	file = NULL)
#' 
#' # get a root age from the fixed ages for tips
#' setRootAge(tree = MAP)
#' 
#' #pull a hundred trees randomly from the posterior
#' hundredRandomlySelectedTrees <- obtainDatedPosteriorTreesMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, 
#'  	burnin = 0.5, 
#'  	getFixedTimes = TRUE,
#'  	getRootAges = TRUE,
#'  	outputTrees = 100, 
#'  	file = NULL)
#' 
#' 
#' }
#' 



#' @name obtainDatedPosteriorTreesMrB
#' @rdname obtainDatedPosteriorTreesMrB
#' @export
obtainDatedPosteriorTreesMrB <- function(
	runFile,nRuns = 2,burnin = 0.5,
	outputTrees,labelPostProb = FALSE,
	getFixedTimes = FALSE,
	getRootAges = FALSE,
	originalNexusFile = NULL,
	file = NULL){
	#
	#########################################################
	#checks
	if(length(outputTrees) != 1){
		stop("outputTrees must be of length 1")
		}
	if(all(
			outputTrees != c("all","MCCT","MAPosteriori","MAPriori","MaxLikelihood")
			)
		& !is.numeric(outputTrees)){
		########
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
	if(getRootAges & !getFixedTimes){
		stop("Root ages cannot be set for output trees if fixed tip dates are not obtained")
		}
	if(getRootAges & !is.null(file)){
		stop("Root ages cannot be set when files are read to an external file format, please use file = NULL (the default)")
		}		
	#
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
			# does runPath end with .nex already?
			endNex <- grepl("\\.nex$", runPath , ignore.case=TRUE)
			#
			# if the nexus file does not end with .nex, add it
			if(endNex){
				originalNexusFile <- runPath
			}else{
				originalNexusFile <- paste0(runPath,".nex")
				}			
			}
		fixedTable <- getMrBFixedAgesFromNexus(origNexusFile = originalNexusFile)
		if(nrow(fixedTable)==0 & getRootAges){
			stop("No fixed ages found for obtaining $root.time, cannot use argument getRootAges = TRUE")
			}
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
			rescaledTrees[[i]][[j]]$LnL <- parFilesBurn[[i]][j,]$LnL
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
	if(outputTrees == "MAPriori"){
		# get Maximum A Priori tree
		LnPr <- sapply(lumpTrees,function(x) x$LnPr)
		whichMaxPrior <- which(LnPr == max(LnPr))
		outTree <- lumpTrees[[whichMaxPrior]]
		if(labelPostProb){
			# assign posterior probabilities as node labels
			outTree <- getPosteriorProbabiities(tree = outTree,postBurninTrees = lumpTrees)
			}
		}
	#
	if(outputTrees == "MaxLikelihood"){
		# get Maximum Model Likelihood tree
		LnL <- sapply(lumpTrees,function(x) x$LnL)
		whichMaxLikelihood <- which(LnL == max(LnL))
		outTree <- lumpTrees[[whichMaxLikelihood]]
		if(labelPostProb){
			# assign posterior probabilities as node labels
			outTree <- getPosteriorProbabiities(tree = outTree,postBurninTrees = lumpTrees)
			}
		}		
	#
	if(outputTrees == "MAPosteriori"){
		# get Maximum A Posteriori tree
		#
		# get posterior probability
		LnL <- sapply(lumpTrees,function(x) x$LnL)
		LnPr <- sapply(lumpTrees,function(x) x$LnPr)
		# Posterior Prob = Prior Prob * Likelihood(data) / likelihood (model)
			# likelihood (model) is unknown (but is constant!)
		# the log posterior probability would be the sum of LnL and LnPr
		#print(LnPr)
		#print(LnL)
		LnPost <- LnPr + LnL
		#
		whichMaxPosterior <- which(LnPost == max(LnPost))
		outTree <- lumpTrees[[whichMaxPosterior]]
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
			#message("For use in R with paleotree and other packages, you may want to set the root age")
			if(is(outTree,"multiPhylo")){
				#message("Simply use function setRootAges() next")
				outTree<-setRootAges(outTree)
			}else{
				#message("Simply use function setRootAge() next")
				outTree<-setRootAge(outTree)
				}
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

