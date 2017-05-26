


#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso
#' See \code{phangorn}'s function for calculating the Robinson-Foulds distance: \code{\link{phangorn::treedist}}

#' @author
#' David W. Bapst. This code was produced as part of a project 
#' funded by National Science Foundation grant EAR-1147537 to S. J. Carlson.

#' @references
#' This measure was introduced in our recent Systematic Biology paper:
#' 
#' Bapst, D. W., H. A. Schreiber, and S. J. Carlson. In press. Combined analysis of extant Rhynchonellida
#' (Brachiopoda) using morphological and molecular data. \emph{Systematic Biology} doi: 10.1093/sysbio/syx049
#'

#' @examples
#' 
#' # let's simulate two trees
#' 
#' set.seed(1)
#' treeA<-rtree(30,br=NULL)
#' treeB<-rtree(30,br=NULL)
#' 
#' \dontrun{
#' 
#' # visualize the difference between these two trees
#' library(phytools)
#' plot(cophylo(treeA,treeB))
#' 
#' # what is the Robinson-Foulds (RF) distance between these trees?
#' library(phangorn)
#' treedist(treeA,treeB)
#' 
#' }
#' 
#' # The RF distance is less intuitive when 
#'     # we consider a tree that isn't well-resolved
#' 
#' # let's simulate the worst resolved tree possible: a star tree
#' treeC<-stree(30)
#' 
#' # plot the tanglegram between A and C
#' plot(cophylo(treeA,treeC))
#' 
#' # however the RF distance is *not* zero
#' # even though the only difference is a difference in resolution
#' \dontrun{
#' treedist(treeA,treeC)
#' }
#' 
#' # the contradiction difference (CD) ignores differences in resolution
#' 
#' # Tree C (the star tree) has zero CD between it and trees A and B
#' identical(treeContradiction(treeA,treeC),0)  # should be zero distance
#' identical(treeContradiction(treeB,treeC),0)  # should be zero distance
#' 
#' # two identical trees also have zero CD between them (as you'd hope) 
#' identical(treeContradiction(treeA,treeA),0)  # should be zero distance
#' 
#' #' and here's the CD between A and B
#' treeContradiction(treeA,treeB)  # should be non-zero distance
#' 
#' # a less ideal property of the CD is that two taxon on opposite ends of the 
#' # moving from side of the topology to the other of an otherwise identical tree
#' # will return the maximum contradiction difference possible (i.e., `= 1`)
#' 
#' # an example
#' treeAA<-read.tree(text="(A,(B,(C,(D,(E,F)))));")
#' treeBB<-read.tree(text="(E,(B,(C,(D,(A,F)))));")
#' plot(cophylo(treeAA,treeBB))
#' treeContradiction(treeAA,treeBB)
#' 
#' \dontrun{
#' # Note however also a property of RF distance too:
#' treedist(treeAA,treeBB)
#' }
#'


#' @name treeContradiction
#' @rdname treeContradiction
#' @export
treeContradiction<-function(tree1,tree2,rescale=TRUE){
  # checks
  if(!inherits(tree1, "phylo")){
		stop("tree1 is not of class phylo")
    }
  if(!inherits(tree2, "phylo")){
		stop("tree2 is not of class phylo")
    }
  #
  tree1<-drop.tip(tree1,setdiff(tree1$tip.label,tree2$tip.label))
  tree2<-drop.tip(tree2,setdiff(tree2$tip.label,tree1$tip.label))
  #
  # more checks
  if(Ntip(tree1)!=Ntip(tree2)){
      stop("Trees do not contain same number of tips after pruning to tips with identical labels (?!)")}
  if(Ntip(tree1)<2 | Ntip(tree2)<2){
      stop("Trees contain less than one tip after pruning")}
  #
  # now measure number of contraditions
  nUnshared<-sum(1==attr(prop.part(tree1,tree2),"number"))
  part1<-lapply(prop.part(tree1),function(x) tree1$tip.label[x])
  part2<-lapply(prop.part(tree2),function(x) tree2$tip.label[x])
  nContra1<-nContradiction(part1,part2)
  nContra2<-nContradiction(part2,part1)
  res<-nContra1+nContra2
  #
  # rescale to 0-1 scale?
  if(rescale){
    #number of possible nodes that could contradict on one unrooted tree
    nPossNodes<-Ntip(tree1)-2   # per tree   
    res<-res/(2*nPossNodes)     # per two trees
    }
  return(res)
  }
  
  
testContradiction<-function(namesA,namesB){
	matchA<-namesA %in% namesB
	matchB<-namesB %in% namesA
	if(any(matchB)){
		res<-!(all(matchA) | all(matchB))
	}else{
		res<-FALSE
		}
	return(res)
	}
  
  
nContradiction<-function(partA,partB){
  partContra<-sapply(partA,function(x) 
      any(sapply(partB,function(y) 
        testContradiction(x,y))))  
  res<-sum(partContra)
  return(res)
  }
