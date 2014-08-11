#' Simulating Un-Conditioned Trees of Fossil Taxa
#' 
#' Creates sets of paleo-trees with minimal conditioning and with sampling such
#' that lineages may be missing and terminal branches are truncated, but
#' node-times are known perfectly
#' 
#' @details
#' This function is essentially a complex wrapper for simFossilTaxa, sampleRanges 
#' and taxa2phylo, and users should refer to the respective documentation for those
#' functions in addition to this manual page: \code{\link{simFossilTaxa}}, 
#' \code{\link{taxa2phylo}} and \code{\link{sampleRanges}}.
#'
#' This function will output simulated phylogenies of fossil taxa
#' where the divergence times are perfectly known, only sampled taxa are
#' included on the tree and tips are located at the last observed time for the
#' species (the apparent time of extinction, except for living taxa).
#' 
#' simPaleoTrees essentially uses simFossilTaxa with no.cond as TRUE and other
#' minimal conditioning, so as to get as unbiased a sample of simulations as
#' possible (without exceeding the maximum constraints). This is useful for
#' birth-death analyses, although a number of options available in
#' simFossilTaxa are thus unavailable in simPaleoTrees. By default, there is no
#' conditioning on the number of extant taxa, living taxa are sampled perfectly
#' at time 0 and zero-length branches are dropped. Unlike simFossilTaxa, you
#' cannot condition on a certain number of extant taxa, only whether they are
#' allowed or not (via all.extinct). As of version 1.6, there are now options
#' relating to speciation modes. By default, taxa are only simulated under
#' budding cladogenesis but this can be changed with the arguments anag.rate,
#' prop.bifurc and prop.cryptic. Trees with cryptic taxa are always returned
#' with cryptic taxa unmerged (see sampleRanges).
#' 
#' Because the divergence times are known perfectly, yet tips are at the
#' apparent time of extinction and unsampled taxa are dropped, one should not
#' use the output of this analysis except for very specialized simulation
#' analyses. The results are probably not anything like real datasets of
#' paleontological phylogenies, at least in most aspects.
#' 
#' The print.runs argument does not work precisely as in simFossilTaxa; it only
#' counts how many accepted datasets from simFossilTaxa are acceptable for
#' output after the simulation of sampling.

#' @inheritParams simFossilTaxa
#' @inheritParams simFossilTaxa_SRCond

#' @param ntrees Number of trees to simulate

#' @param all.extinct Condition on all taxa being extinct by modern? Default is
#' false

#' @param modern.samp.prob Probability of sampling living taxa at the present
#' day (time=0), see documentation for sampleRanges

#' @param prop.cryptic Proportion of branching events with no morphological
#' differentiation (i.e. cryptic speciation) relative to branching events
#' associated with morphological differentiation (budding and/or bifurcating 
#' cladogenesis). Trees with cryptic taxa are always returned with cryptic 
#' taxa unmerged (see sampleRanges).

#' @param drop.zlb Should zero-length branches be dropped?

#' @param ranges.only If TRUE (the default), the ranges returned in $ranges
#' are given as taxon first and last occurrences only. If
#' FALSE, $ranges returns the times of all sampling events for each taxon
#' as vectors within a list.

#' @param plot Should data be plotted as it is simulated?

#' @return Output is an object of class multiphylo containing the simulated
#' phylogenies, unless ntrees is one in which case the output is a phylogeny of
#' class 'phylo'.
#' 
#' Additionally, each of these simulated phylogenies will have the original
#' simulated taxa data (from simFossilTaxa) and sampled ranges (from
#' sampleRanges) attached as the elements $taxa and $ranges to each phylo
#' object.

#' @seealso \code{\link{simFossilTaxa}}, \code{\link{taxa2phylo}},
#' \code{\link{sampleRanges}}

#' @examples
#' 
#' set.seed(444)
#' #simulate trees conditioned to have no living descendants
#' trees <- simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=10,all.extinct=TRUE,maxtime=100,
#'     print.runs=TRUE,plot=TRUE)
#' #number of tips
#' sapply(trees,Ntip)
#' 
#' #simulate trees conditioned to have possible living taxa and perfect sampling at modern
#' trees <- simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=10,all.extinct=FALSE,maxtime=100,
#'     modern.samp.prob=TRUE,print.runs=TRUE,plot=TRUE)
#' #number of tips
#' sapply(trees,Ntip)
#' 
#' @export simPaleoTrees
simPaleoTrees<-function(p,q,r,ntrees=1,all.extinct=FALSE,modern.samp.prob=1.0,mintime=1,maxtime=100,
	mintaxa=2,maxtaxa=500,anag.rate=0,prop.bifurc=0,prop.cryptic=0,drop.zlb=TRUE,print.runs=FALSE,ranges.only=TRUE,
	plot=FALSE){
	#this is a wrapper which will create many paleo trees with at least two observed tips
		#uses simFossilTaxa, sampRanges,taxa2phylo, etc
		#good if you want to simulate many many trees with extinct taxa
		#divergence times for nodes will be perfectly known
	#by default: (minimal conditioning)
		#no conditioning on the number of extant taxa
		#living taxa are sampled perfectly at the present
		#zero-length branches are dropped
	#simPaleoTrees(p=0.1,q=0.1,r=0.1,ntrees=10)
	#require(ape)
	if(mintaxa<2){stop("Error: Need at least two taxa per tree; increase mintaxa")}
	if(ntrees<1){stop("Error: number of trees to simulate is <1")}
	res<-rmtree(ntrees,2)
	ntries<-0
	for(i in 1:ntrees){
		rerun<-TRUE
		while(rerun){
			ntries<-ntries+1
			taxa<-suppressMessages(simFossilTaxa(p=p,q=q,anag.rate=anag.rate,prop.bifurc=prop.bifurc,prop.cryptic=prop.cryptic,nruns=1,mintaxa=mintaxa,
				maxtaxa=maxtaxa,maxtime=maxtime,maxExtant=ifelse(all.extinct,0,maxtaxa),min.cond=FALSE,plot=plot))
			ranges<-sampleRanges(taxa,r,min.taxa=0,modern.samp.prob=modern.samp.prob,merge.cryptic=FALSE,ranges.only=ranges.only)
			if(ranges.only){
				ranges1<-ranges
			}else{
				ranges1<-cbind(sapply(ranges,max),sapply(ranges,min))
				rownames(ranges1)<-names
				colnames(ranges1)<-c("FAD","LAD")
				}
			if(sum(!is.na(ranges1[,1]))>1){
				tree<-taxa2phylo(taxa,obs_time=ranges1[,2],plot=plot)
				if(drop.zlb){tree<-dropZLB(tree)}
				if(all(!is.na(tree))){
					numext<-sum(ranges1[,2]==0)
					minta<-Ntip(tree)>mintaxa
					minti<-(max(ranges1,na.rm=TRUE)-min(ranges1,na.rm=TRUE))>mintime
					if(minta & minti){
						rerun<-FALSE
						}
					}
				}
			}
		tree$taxa<-taxa		#attach original taxon matrix
		tree$ranges<-ranges	#attach sampled ranges
		res[[i]]<-tree
		if(ntrees>10){if((i%%(ntrees/10))==0){message(cat(round(i*100/ntrees),"% ",sep=""))}}
		}
	if(print.runs){message(paste(ntrees," trees accepted from ",ntries," total runs (",signif(ntrees/ntries,2)," Acceptance Probability)",sep=""))}
	if(ntrees==1){res<-res[[1]]}
	return(res)
	}
