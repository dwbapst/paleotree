

#' Modify or Drop Terminal Branches of Various Types
#' 
#' These functions modify terminal branches or drop certain terminal branches
#' based on various criteria.
#' 
#' DropZLB drops tip-taxa that are attached to the tree via zero-length
#' terminal branches ("ZLBs"). This is sometimes useful for paleo-trees, as
#' various time-scaling methods often produce these ZLBs, taxa whose early
#' appearance causes them to be functionally interpreted as ancestors in some
#' time-scaling methods. Removing ZLBs is advised for analyses of
#' diversification/diversity, as these will appear as simultaneous
#' speciation/extinction events. Note this function only drops tips attached to
#' a terminal zero-length branch; if you want to collapse internal zero-length
#' branches, see the ape function \code{\link{di2multi}}.
#' 
#' DropExtinct drops all terminal branches which end before the modern (i.e.
#' extinct taxa). DropExtant drops all terminal branches which end at the
#' modern (i.e. extant/still-living taxa). In both cases, the modern is defined
#' based on tree$root.time if available, or the modern is inferred to be the
#' point in time when the tip furthest from the root (the latest tip)
#' terminates.
#' 
#' If the input tree has a $root.time element, as expected for most paleo-tree
#' objects handled by this library, that root.time is adjusted if the relative
#' time of the root divergence changes when terminal branches are dropped.
#' Adjusted root.times are only given if the input tree has root.times.
#' 
#' addTermBranchLength adds an amount equal to the argument 'addtime' to the
#' terminal branch lengths of the tree. If there is a $root.time element, this
#' is increased by an amount equal to addtime. A negative amount can be input
#' to reduce the length of terminal branches. However, if negative branch
#' lengths are produced, the function fails and a warning is produced.
#' 
#' When a tree is modified, such as having tips dropped or branches extended,
#' fixRootTime can be used to find the new $root.time. It is mainly used as a
#' utility function called by the other functions discussed in this help file.
#' 
#' @aliases dropZLB dropExtinct dropExtant addTermBranchLength fixRootTime
#' @param tree A phylogeny as a phylo object
#' @param tol Tolerance for determining modern age; used for distinguishing
#' extinct from extant taxa. Tips which end within 'tol' of the furthest
#' distance from the root will be treated as 'extant' taxa for the purpose of
#' keeping or dropping.
#' @param ignore.root.time Ignore root.time in calculating which tips are
#' extinct? root.time will still be adjusted
#' @param addtime Extra amount of time to add to all terminal branch lengths.
#' @param treeOrig A phylo object of a time-scaled phylogeny with a $root.time
#' element
#' @param treeNew A phylo object containing a modified form of treeOrig (with
#' no extra tips added)
#' @return Gives back a modified phylogeny as a phylo object, generally with a
#' modified $root.time element.
#' @author David W. Bapst
#' @seealso \code{\link{phyloDiv}}, \code{\link{drop.tip}},
#' \code{\link{compareTermBranches}}
#' @examples
#' 
#' set.seed(444)
#' #Simulate some fossil ranges with simFossilTaxa
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont<-sampleRanges(taxa,r=0.5)
#' #Now let's make a tree using taxa2phylo
#' tree <- taxa2phylo(taxa,obs_time=rangesCont[,2])
#' #compare the two trees
#' layout(1:2)
#' plot(ladderize(tree))
#' plot(ladderize(dropZLB(tree)))
#' layout(1)
#' 
#' #example using dropExtinct and dropExtant
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=40,maxtime=1000,maxExtant=20)
#' tree<-taxa2phylo(taxa)
#' phyloDiv(tree)
#' tree1 <- dropExtinct(tree)
#' phyloDiv(tree1)
#' tree2 <- dropExtant(tree)
#' phyloDiv(tree2)
#' 
#' #example using addTermBranchLength
#' set.seed(444)
#' treeA <- rtree(10)
#' treeB <- addTermBranchLength(treeA,1)
#' compareTermBranches(treeA,treeB)
#' 

dropZLB<-function(tree){
	#drops terminal branches that are zero length
		#adjusts tree$root.time if necessary
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	drop_e<-(tree$edge[,2]<(Ntip(tree)+1)) & (tree$edge.length==0)
	drop_t<-(tree$edge[,2])[drop_e]
	if((Ntip(tree)-length(drop_t))>1){
		tree1<-drop.tip(tree,drop_t)
		if(!is.null(tree$root.time)){tree1<-fixRootTime(tree,tree1)}
		res<-tree1
	}else{res<-NA}
	return(res)
	}

dropExtinct<-function(tree,tol=0.01,ignore.root.time=FALSE){
	#drop all terminal taxa that are less than 0.001 from the modern
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(is.null(tree$root.time)){
		message("No tree$root.time: Assuming latest tip is at present (time=0)")
		}
	dnode<-dist.nodes(tree)[1:Ntip(tree),Ntip(tree)+1]
	dnode<-round(dnode,6)
	if(!is.null(tree$root.time) & !ignore.root.time){if(round(tree$root.time,6)>max(dnode)){stop("Error: all tips are extinct based on tree$root.time!")}}
	droppers<-which((dnode+tol)<max(dnode))
	if((Ntip(tree)-length(droppers))<2){stop("Error: Less than 2 tips are extant on the tree!")}
	stree<-drop.tip(tree,droppers)
	if(!is.null(tree$root.time)){
		#now need to add $root.time given the droppers
		#should be root.time MINUS distance from furthest tip in tree PLUS distance from latest tip to root of stree
		#stree$root.time<-tree$root.time-max(dnode)+max(dist.nodes(stree)[1:Ntip(stree),Ntip(stree)+1])
		stree<-fixRootTime(tree,stree)
		}
	return(stree)
	}

dropExtant<-function(tree,tol=0.01){
	#drop all terminal taxa that are more than 0.001 from the modern
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(is.null(tree$root.time)){
		message("Warning: no tree$root.time! Assuming latest tip is at present (time=0)")
		}
	dnode<-dist.nodes(tree)[1:Ntip(tree),Ntip(tree)+1]
	dnode<-round(dnode,6)
	if(!is.null(tree$root.time)){if(round(tree$root.time,6)>max(dnode)){stop("Error: all tips are extinct based on tree$root.time!")}}
	droppers<-which((dnode+tol)>max(dnode))
	if((Ntip(tree)-length(droppers))<2){stop("Error: Less than 2 tips extinct on the tree!")}
	stree<-drop.tip(tree,droppers)
	if(!is.null(tree$root.time)){
		#now need to add $root.time given the droppers
		#should be root.time MINUS distance from earilest tip in tree PLUS distance from earliest tip to root of stree
		#stree$root.time<-tree$root.time-min(dnode)+min(dist.nodes(stree)[1:Ntip(stree),Ntip(stree)+1])
		stree<-fixRootTime(tree,stree)
		}
	return(stree)
	}

addTermBranchLength<-function(tree,addtime=0.001){
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]<-tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]+addtime
	if(any(tree$edge.length<0)){stop("Error: tree has negative branch lengths!")}
	if(!is.null(tree$root.time)){tree$root.time<-tree$root.time+addtime}
	return(tree)
	}

fixRootTime<-function(treeOrig,treeNew){
	treeDepth<-function(tree){
		#require(ape)
		max(dist.nodes(tree)[,Ntip(tree)+1])
		}
	#require(ape)
	if(!is(treeOrig, "phylo")){stop("Error: treeOrig is not of class phylo")}
	if(!is(treeNew, "phylo")){stop("Error: treeNew is not of class phylo")}
	if(is.null(treeOrig$root.time)){stop("ERROR: treeOrig passed to fixRootTime with no $root.time??")}
	orig_dist<-dist.nodes(treeOrig)[
		which(treeNew$tip.label[1]==treeOrig$tip.label),Ntip(treeOrig)+1
		]
	new_dist<-dist.nodes(treeNew)[1,Ntip(treeNew)+1]
	treeNew$root.time<-treeOrig$root.time-(orig_dist-new_dist)
	if(round(max(dist.nodes(treeNew)[, Ntip(treeNew) + 1]) - treeNew$root.time)>0){
		stop("Error: fixRootTime isn't fixing correctly, root.time less than max tip-to-root length!")}
	return(treeNew)
	}