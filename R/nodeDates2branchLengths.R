#' Obtaining Edge Lengths for Undated Phylogenies Using Known Branching Node and Tip Ages
#' 
#' This function takes some undated phylogenetic topology, a set of ages in
#' absolute time, for the internal nodes and (by default) the terminal tips
#' of that phylogeny,  and returns a dated phylogeny consistent with those input ages.

#' @details
#' The function \code{\link{compute.brtime}} in package \code{ape} does
#' a very similar functionality, but is limited in its application for
#' only ultrametric trees, as it does not allow for tips to have
#' incongruent ages. It also only accepts node ages as on the relative
#' scale where the latest tips are at zero, as assumed in general
#' elsewhere in package \code{ape}.

#' @param nodeDates Under default \code{allTipsModern = FALSE} conditions,
#' \code{nodeDates} should be a vector of length \code{Ntip(tree) + Nnode(tree)}
#' which contains the dates for all terminal tip nodes and internal nodes for
#' the tree, in that order, as numbered in the \code{tree$edge} matrix. Such
#' a vector is produced as output by \code{\link{dateNodes}}. If \code{allTipsModern = TRUE},
#' then the vector should only be as long as the number of nodes,
#' and contain the dates only for those same internal nodes in \code{tree}. 
#' These dates should always on a descending scale (i.e. time before present), with respect to an
#' absolute time-scale. It is possible for the time 0 date to represent
#' a date far in the future from the latest tip.

#' @param tree An undated phylogeny object, of class \code{phylo}, lacking edge lengths.
#' If the tree appears to be dated (i.e. has edge lengths), the function will issue a warning.

#' @param  allTipsModern A logical, default is \code{FALSE}. If \code{FALSE},
#' then the function expects \code{nodeDates} to contain ages for all
#' 'nodes' - both internal branching nodes and terminal tips. If \code{TRUE}, then
#' the function will expect \code{nodeDates} to contain ages only for internal
#' branching nodes, and all tips will be assumed to be at time 0. (Thus, if your
#' tree is ultrametric but tips aren't all at the modern, do \emph{not}
#' use {allTipsModern = TRUE}).

#' @return
#' A dated tree as a list of class \code{phylo}, with a \code{$root.time}
#' element for referencing the tree against absolute time.

#' @seealso
#' This function will likely often be used in conjunction with
#' \code{\link{dateNodes}}, such as for summarizing node and tip age
#' estimates from a sample of trees, to produce a single dated tree
#' to act as a point estimate. Beware however that point estimates of
#' tree samples may have little resemblance to any individual tree in that sample.
#' 
#' This function should perform identically for ultrametric trees as package \code{ape}'s function \code{\link{compute.brtime}}.

#' @author David W. Bapst

#' @references

#' @examples
#' set.seed(444)
#' 
#' # we'll do a number of tests, let's check at the end that all are TRUE
#' tests<-logical()
#' 
#' # with a non-ultrametric tree
#' chrono<-rtree(10)
#' # make an undated tree
#' notChrono<-chrono
#' notChrono$edge.length<-NULL
#' 
#' # now lets try with dateNodes in paleotree
#' nodeTimes<-dateNodes(chrono)
#' # need to use allTipsModern = FALSE because tip ages are included
#' chronoRedux <-  nodeDates2branchLengths(tree=notChrono,
#'     nodeDates=nodeTimes, allTipsModern=FALSE)
#' # test that its the same
#' (tests<-c(tests,all.equal.numeric(chrono$edge.length,chronoRedux$edge.length))
#' 
#' ######################################
#' # modern ultrametric tree
#' chrono<-rcoal(10)
#' # make an undated tree
#' notChrono<-chrono
#' notChrono$edge.length<-NULL
#' 
#' # with ultrametric trees, you could just use ape's compute.brtime 
#' 
#' # getting branching times with ape
#' branchingTimes <- branching.times(chrono)	
#' # setting those branching times with ape
#' chronoRedux <-  compute.brtime(notChrono, branchingTimes)
#' # test that its the same
#' (tests<-c(tests,all.equal.numeric(chrono$edge.length,chronoRedux$edge.length))
#' 
#' # lets do the same thing but with nodeDates2branchLengths
#' 
#' # can use branching.times from ape 
#' 	# (but only for ultrametric trees!)
#' chronoRedux <-  nodeDates2branchLengths(tree=notChrono,
#'     nodeDates=branchingTimes, allTipsModern=TRUE)
#' # test that its the same
#' (tests<-c(tests,all.equal.numeric(chrono$edge.length,chronoRedux$edge.length))
#' 
#' # now lets try with dateNodes in paleotree
#' nodeTimes<-dateNodes(chrono)
#' # need to use allTipsModern = FALSE because tip ages are included
#' chronoRedux <-  nodeDates2branchLengths(tree=notChrono,
#'     nodeDates=nodeTimes, allTipsModern=FALSE)
#' # test that its the same
#' (tests<-c(tests,all.equal.numeric(chrono$edge.length,chronoRedux$edge.length))
#' 
#' # get just the node times (remove tip dates)
#' nodeOnlyTimes<-nodeTimes[-(1:Ntip(chrono))]
#' # let's use the allTipsModern = TRUE setting
#' chronoRedux <-  nodeDates2branchLengths(tree=notChrono,
#'     nodeDates=nodeOnlyTimes, allTipsModern=TRUE)
#' # test that its the same
#' (tests<-c(tests,all.equal.numeric(chrono$edge.length,chronoRedux$edge.length))
#' 
#' # did all tests come out as TRUE?
#' if(!all(tests)){stop("nodeDates2branchLengths isn't functioning correctly")}
#' 



#' @name nodeDates2branchLengths
#' @rdname nodeDates2branchLengths
#' @export
nodeDates2branchLengths<-function(nodeDates, tree, allTipsModern=FALSE){
    # checks
	if(!inherits(tree,"phylo")){
		stop("tree is not of class phylo")
		}
	if(!is.null(tree$edge.lengths)){
        warning("input tree has $edge.lengths present, these will be replaced")
		}
    #
	if(allTipsModern){
	    if(length(nodeDates)!=Nnode(tree)){
			stop("nodeDates must be same length as number of nodes on input tree if allTipsModern=TRUE")
			}
		#add zero ages for tips
		allAges<-c(rep(0,Ntip(tree)),nodeDates)
		stop("nodeDates2branchLengths doesn't handle non-ultrametric trees or trees where not all tips are at the modern day... yet")
	}else{
	    if(length(nodeDates)!=(Nnode(tree)+Ntip(tree))){
			stop("nodeDates must be same length as number of nodes AND tips on input tree if allTipsModern=FALSE")
			}
		allAges<-nodeDates
		}
	# check that all ages are provided
	if(any(is.na(allAges)) | any(is.null(allAges))){
		stop("Some input ages appear to be NA or NULL - ages for all nodes MUST be provided")
		}	
    ######################################
    # get mother node age for each edge
    momAges<-allAges[tree$edge[,1]]
    # get node ages for child nodes of each edge
    childAges<-allAges[tree$edge[,2]]
    #edge lengths = mom - child
    edgeLengths<-momAges-childAges
	# check edgeLengths
	if(any(edgeLengths<0)){
		stop("Check ages - some edges are calculated as having negative lengths!")
		}
	if(any(edgeLengths==0)){
		message("Caution: some edges are calculated as being zero-length - is this expected?")
		}
	#
    tree$edge.length<-edgeLengths
	# set root age
	tree$root.time<-allAges[(Ntip(tree)+1]
    return(tree)
    }
