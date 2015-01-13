#' Randomly Collapse a Portion of Nodes on a Phylogeny
#' 
#' \code{degradeTree} removes a proportion of the total nodes in a tree, chosen
#' randomly, collapsing the nodes to produce a less-resolved tree. The related function \code{collapseNodes} 
#'given a tree and a vector of nodes to collapse, removes those nodes from a tree, creating a polytomy.

#' @details In the function \code{degradeTree}, the nodes are removed at random
#' using the basic R function \code{sample}. \code{degradeTree}
#' can be conditioned to remove nodes of a particular depth with greater
#' probability/frequency by setting node.depth to a value between zero
#' (favoring the removal of deep nodes close to the root) or one (shallow nodes
#' far from the root). Depth is evaluated based on the number of descendant
#' tips. If node.depth is not NA, the relative proportion of descendants from
#' each node is calculated, summed to 1 and the node.depth value subtracted
#' from this proportion. These values are then squared, normalized again to
#' equal to 1 and then used as the probabilities for sampling nodes for
#' removal.
#' 
#' By default, branch lengths are removed from the input tree prior to
#' degradation and entirely absent from the output tree. This is changed if
#' argument leave.zlb is TRUE.
#' 
#' @param tree A phylogeny of class 'phylo'
#' @param prop_collapse Proportion of nodes to collapse
#' @param nCollapse Number of nodes to collapse, can be supplied as an
#' alternative to prop_collapse
#' @param node.depth A number between 0 to 1, which conditions the depth of
#' nodes removed. If NA, no conditioning (this is the default).
#' @param leave.zlb If FALSE, the default option, the original branch length
#' distribution is destroyed and branches set to zero by this function will
#' return polytomies. If TRUE, then the original edge lengths are kept for
#' unmodified edges, and modified edges are changed to zero length, and are not
#' collapsed into polytomies. The removed branch length is not shifted to other
#' edges.
#' @param nodeID The node ID number(s) to be collapsed into a polytomy, as identified in
#' the \code{$edge} matrix of the'phylo' object. Must be a vector of one or more ID numbers.

#' @param collapseDir Whether to collapse the edge leading the listed node (if "forward"), or to
#' collapse the child edges leading away from the node (if "backward"). Collapsing a node into a
#' polytomy conceptually could be either and users should heed this option carefully. 

#' @return Returns the modified tree as an object of class phylo, with no edge
#' lengths by default.
#' @seealso \code{\link{di2multi}},\code{\link{timeLadderTree}},
#' @examples
#' 
#' set.seed(444)
#' tree <- rtree(100)
#' tree1 <- degradeTree(tree,0.5) 
#' 
#' #let's compare the input and output
#' layout(matrix(1:2,,2))
#' plot(tree,show.tip.label=FALSE,use.edge.length=FALSE)
#' plot(tree1,show.tip.label=FALSE,use.edge.length=FALSE)
#' 
#' #now with collapseNodes
#' tree <- rtree(10)
#' #collapse nodes backwards
#'    #let's collapse lucky node number 13!
#' tree1 <- collapseNodes(nodeID=13,tree=tree,collapseDir="backward")  
#' #collapse nodes forwards 
#' tree2 <- collapseNodes(nodeID=13,tree=tree,collapseDir="forward")
#' 
#' #let's compare
#' layout(1:3)
#' plot(tree,use.edge.length=FALSE,main="original")
#' plot(tree1,use.edge.length=FALSE,main="backward collapse")
#' plot(tree2,use.edge.length=FALSE,main="forward collapse")
#' 
#' layout(1)
#' 
#' @name degradeTree
#' @export degradeTree
degradeTree<-function(tree,prop_collapse,nCollapse=NULL,node.depth=NA,leave.zlb=FALSE){
	#collapses a given proportion of internal edges, creating polytomies
		#node.depth conditions on depth of edge in tree
			# 1 removes more shallow nodes, 0 removes deeper nodes
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	edge<-(1:length(tree$edge))[which(tree$edge[,2]>Ntip(tree))]	#internal edges
	if(is.null(nCollapse)){nCollapse<-round(prop_collapse*length(edge))}
	if(is.na(node.depth)){
		cedge<-sample(edge,nCollapse)	#edges chosen to collapse
	}else{
		node_pdesc<-sapply(prop.part(tree),length)/Ntip(tree)	#prop desc per int node
		edge_pdesc<-node_pdesc[tree$edge[edge,2]-Ntip(tree)]
		edge_prob<-(edge_pdesc-node.depth)^2;edge_prob<-edge_prob/sum(edge_prob)
		cedge<-sample(edge,nCollapse,prob=edge_prob)	#chosen edges	
		}
	if(leave.zlb){
		tree$edge.length[cedge]<-0
	}else{
		tree$edge.length<-NULL
		tree<-collapse.singles(tree)
		tree$edge.length<-rep(1,Nedge(tree))
		tree$edge.length[cedge]<-0
		tree<-di2multi(tree)
		tree<-collapse.singles(tree)
		tree$edge.length<-NULL		
		}
	return(tree)
	}

#' @rdname degradeTree
#' @export 
collapseNodes<-function(tree,nodeID,collapseDir,leave.zlb=FALSE){
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(!all(nodeID>0) | !all(nodeID<(Ntip(tree)+Nnode(tree)+1))){
		stop("Error: some nodeID values outside the range of tip and node IDs")}
	if(!all(nodeID>Ntip(tree))){
		message("Warning: Some nodeID values indicate terminal tips; collapsing these generally doesn't do anything")}
	if(!any(collapseDir==c("forward","backward"))){stop("collapseDir must be either 'forward' or 'backward'")}
	if(collapseDir=="backward" | collapseDir=="forward"){
		if(collapseDir=="backward"){
			# 01-15-14 if the edge is collapsed backward, then its the edges who have that node as a descendant!!!
			cedge<-which(sapply(tree$edge[,2],function(x) any(x==nodeID)))
			}
		if(collapseDir=="forward"){
			cedge<-which(sapply(tree$edge[,1],function(x) any(x==nodeID)))
			}
	}else{
			stop("collapseDir is not either 'backward' or 'forward' ??")
			}
	if(leave.zlb){
		tree$edge.length[cedge]<-0
	}else{
		tree$edge.length<-NULL
		tree<-collapse.singles(tree)
		tree$edge.length<-rep(1,Nedge(tree))
		tree$edge.length[cedge]<-0
		tree<-di2multi(tree)
		tree<-collapse.singles(tree)
		tree$edge.length<-NULL		
		}
	return(tree)
	}