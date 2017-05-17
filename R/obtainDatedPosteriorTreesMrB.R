#' Get the Sample of Posterior Trees from a Dated Phylogenetic Analysis with MrBayes (Or a Summary Tree, such as the MCCT)
#' 
#' 

#' @details



# pick either a .p or .t file; can be either. Used for getting path. 
# all needed .p and .t files must be in the same directory as this file


# number of runs in MrBayes, default is 2

#' @param runfile A filename in the current directory, or a path to a file that is either a .p 
#' or .t file from a MrBayes analysis. This filename and path will be used for finding additional 
#' .t and .p files, via the \code{nRuns} settings and assuming that files are in the
#' same directory \emph{and} these files are named under
#' typical MrBayes file naming conventions. (In other words, if you've renamed your .p or .t files,
#' this function probably won't be able to find them,.)

#' @param nRuns The number of runs 

#' @param burnin 

#' @param outputTrees 

#' @param file 


#' @return
#' Returns either a single phylo object (the maximum clade 

#' @aliases

#' @seealso
#' \code{\link{maxCladeCred}}

#' @author
#' David Bapst, with rescaling of raw output
#' trees via code originally written by Nicholas Crouch.

#' @examples
#' \dontrun{
#' 
#' MCCT<-obtainDatedPosteriorTreesMrB(
#'  	runFile="C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns=2,burnin=0.5,
#' 	outputTrees="MCCT",file=NULL)
#' 
#' }
#'



#' @name obtainDatedPosteriorTreesMrB
#' @rdname obtainDatedPosteriorTreesMrB
#' @export
obtainDatedPosteriorTreesMrB<-function(runFile,nRuns=2,burnin=0.5,outputTrees,file=NULL){
	#checks
	if(length(outputTrees)!=1){
		stop("outputTrees must be of length 1")
		}
	if(all(outputTrees!=c("all","MCCT","MAP")) & !is.numeric(outputTrees)){
		stop("outputTrees must be one of 'all', 'MCCT', or a numeric value indicating the number of trees to randomly sample from the posterior")
		}
	if(is.numeric(outputTrees)){	
		if(outputTrees<1){
			stop("If numeric, outputTrees must be greater than 1")
			}
		}
	if(burnin>1 | burnin<0){
		stop("burnin must be a value between 0 and 1")
		}
	# Load tree, paramater and mcmc files produced by MrBayes
	#
	# take indicated file and get the basic run name
	runPath<-strsplit(runFile,split="\\.run")[[1]]
	if(length(runPath)!=2){
		stop("Unable to parse the runPath correctly")
		}
	runPath<-runPath[[1]]
	# get list of tree files
	treeFiles<-lapply(1:nRuns,function(x){
		read.nexus(file=paste0(runPath,".run",x,".t"))}
		)
	# get list of .p files
	parFiles<-lapply(1:nRuns,function(x)
		read.table(file=paste0(runPath,".run",x,".p"),
			 header = T, skip=1)
		)
	#############################
	# checks for length 
	for(i in 1:nRuns){
		if(length(treeFiles[[i]])!=nrow(parFiles[[i]])){
			stop("Parameter data and tree data not of the same length")
			}
		}
	####################
	#
	# Specify sample start based on burnin - here 50%
	sampleStart <- floor(length(treeFiles[[1]])*burnin)
	#
	# remove burnin from both
	treeFilesBurn<-lapply(treeFiles,function(x){
		startSamp<-floor(length(x)*burnin)
		x[startSamp:length(x)]
		}
	)
	parFilesBurn<-lapply(parFiles,function(x){
		startSamp<-floor(nrow(x)*burnin)
		x[startSamp:nrow(x),]
		}
	)
	# checks for length again
	for(i in 1:nRuns){
		if(length(treeFilesBurn[[i]])!=nrow(parFilesBurn[[i]])){
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
	rescaledTrees<-lapply(1:nRuns,function(x){ 
		lapply(1:length(treeFilesBurn[[x]]),function(y){
			rescaleMrBTree(treeFilesBurn[[x]][[y]],parFilesBurn[[x]][y,])
			})
		})
	#################################
	#
	# attach marginal likelihoods to each tree
	for(i in 1:nRuns){
		for(j in 1:length(rescaledTrees[[i]])){
			rescaledTrees[[i]][[j]]$LnPr<-parFilesBurn[[i]][j,]$LnPr
			}
		}
	# concatanate trees from each run
	lumpTrees<-unlist(rescaledTrees,recursive=FALSE)
	class(lumpTrees)<-"multiPhylo"
	#
	##########################################
	#
	if(outputTrees=="all"){
		outTree<-lumpTrees
		}
	#
	if(outputTrees=="MAP"){
		# get MAP
		LnPr<-sapply(lumpTrees,function(x) x$LnPr)
		whichMAP<-which(LnPr==max(LnPr))
		outTree<-lumpTrees[[whichMAP]]
		}
	#
	if(outputTrees=="MCCT"){
		# turns out the MCCT isn't the tree with the highest likelihood
		outTree<-phangorn::maxCladeCred(lumpTrees)
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
		whichOutput<-sample(length(lumpTrees),outputTrees,replace=FALSE)
		outTree<-lumpTrees[whichOutput]
		}
	########################################
	if(!is.null(file)){
		write.nexus(outTree, file=file)
	}else{
		return(outTree)
		}
	}



