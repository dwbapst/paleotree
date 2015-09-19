#' Plot a Traitgram for Continuous Traits
#' 
#' This function uses maximum-likelihood ancestral trait estimation to plot a
#' 'traitgram' (Ackerly, 2009) given a tree and a set of continuous trait
#' values.
#' 
#' @details By default, this function uses \code{\link{ace}} from the library ape to
#' reconstruct ancestral traits and confidence intervals using the PIC method.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' @param trait a vector of continuous trait values
#' @param tree a phylo object
#' @param trait.name The name of the trait plotted, used for the figure's main
#' title
#' @param conf.int if true, confidence intervals are plotted
#' @param lwd The line width used for the figure
#' @return Returns no value, just plots the traitgram.
#' @note One should probably never do ancestral trait estimation without
#' looking at the confidence intervals, as these reconstructed estimates tend
#' to be very uncertain.
#' @author David W. Bapst
#' @seealso \code{\link{ace}}
#' 
#' Also see the functions \code{traitgram} in the library picante and
#' \code{phenogram} in the library phytools.
#' @references Ackerly, D. 2009 Conservatism and diversification of plant
#' functional traits: Evolutionary rates versus phylogenetic signal.
#' \emph{Proceedings of the National Academy of Sciences} \bold{106}(Supplement
#' 2):19699--19706.
#' @examples
#' 
#' set.seed(444)
#' tree <- rtree(10)
#' trait <- rTraitCont(tree)
#' 
#' #first, traitgram without conf intervals
#' plotTraitgram(trait,tree,conf.int=FALSE)
#' 
#' #now, with
#' plotTraitgram(trait,tree)
#' #not much confidence, eh?
#' 
#' @export plotTraitgram
plotTraitgram<-function(trait,tree,trait.name="'trait'",conf.int=TRUE,lwd=1.5){
	#traitgram plotted using ML ASR from geiger (or ace() from ape if ci=TRUE)
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	#sort trait, if not sorted already
	if(is.null(names(trait))){
		message("No names for trait data, assuming in same order as tree$tip.label")
	}else{
		trait<-trait[tree$tip.label]
		}
	#get root time
	if(is.null(tree$root.time)){tree$root.time<-max(dist.nodes(tree)[Ntip(tree)+1,1:Ntip(tree)])}
	times<-tree$root.time-dist.nodes(tree)[Ntip(tree)+1,]
	if(conf.int){
		asr<-ace(trait,tree,method="pic")
		tr1<-c(trait,asr$ace);edges<-tree$edge
		ci<-asr$CI95;tr2<-c(tr1,ci)
		plot(1,1,type="n",
			xlim=c(min(tr2)-0.1,max(tr2)+0.1),
			ylim=c(max(times),min(times)),
			xlab="Trait Values",ylab="Time (Before Present)",
			main=paste("Traitgram of",trait.name))
		for(i in 1:nrow(edges)){
			anc<-edges[i,1];desc<-edges[i,2]
			lines(c(tr1[anc],tr1[desc]),c(times[anc],times[desc]),lwd=lwd)
			}
		for(i in 1:Nnode(tree)){
			lines(c(ci[i,1],ci[i,2]),c(times[i+Ntip(tree)],times[i+Ntip(tree)]),lwd=lwd)
			}
	}else{
		#require(geiger)	#not anymore, not as of 04-03-13
		tr1<-c(trait,ace(trait,tree,method="pic")$ace);edges<-tree$edge
		plot(1,1,type="n",xlim=c(min(tr1)-0.1,max(tr1)+0.1),ylim=c(max(times),min(times)),
			xlab="Trait Values",ylab="Time  (Before Present)",
			main=paste("Traitgram of",trait.name))
		for(i in 1:nrow(edges)){anc<-edges[i,1];desc<-edges[i,2]
			lines(c(tr1[anc],tr1[desc]),c(times[anc],times[desc]),lwd=lwd)}
		}
	}
