#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @author 

#' @examples
#' 
#' #simulation with an example non-ultrametric tree
#' 
#' tree<-rtree(20)
#' tree<-degradeTree(tree,0.3,leave.zlb=TRUE) 	#randomly replaces edges with ZLBs, similar to multi2di output
#' 
#' tree2<-minBranchLength(tree,0.1)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree2);axisPhylo()
#' 
#' #now let's try it with an ultrametric case
#' 
#' tree<-rtree(30)
#' tree<-degradeTree(tree,0.5,leave.zlb=TRUE) 	#randomly replaces edges with ZLBs, similar to multi2di output
#' tree<-di2multi(tree)
#' tree<-compute.brlen(tree)
#' 
#' plot(tree) #ultrametric tree with polytomies, yay
#' 
#' #now randomly resolve, get new branch lengths as would with real data
#' tree2<-multi2di(tree)
#' tree2<-minBranchLength(tree2,0.1)
#' 
#' layout(1:2)
#' plot(tree);axisPhylo()
#' plot(tree2);axisPhylo()
#' 
#' layout(1)

#' @name minBranchLength
#' @rdname minBranchLength
#' @export
minBranchLength<-function(tree, mbl){	
	#require(phangorn)
	#test arguments
	#tree - a tree with edge lengths
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(is.null(tree$edge.length)){stop("Error: Tree has no edge lengths")}
	timetree<-tree
	#mbl - a single numeric value
	if(!is.numeric(mbl) | length(mbl)!=1){
		stop("Error: mbl is a not a single numeric value")}
	#
	root_node<-Ntip(timetree)+1
	while(any(timetree$edge.length<mbl)){
		#pick one at random, make vector of every mom node that is ancestral
		mom<-timetree$edge[((1:Nedge(timetree))[timetree$edge.length<mbl])
			[sample(length((1:Nedge(timetree))[timetree$edge.length<mbl]),1)],1]
		mom<-c(mom,Ancestors(timetree,mom))
		debt<-mbl-min(timetree$edge.length[timetree$edge[,1]==mom[1]])
		timetree$edge.length[mom[1]==timetree$edge[,1]]<-timetree$edge.length[mom[1]==timetree$edge[,1]] + debt[1]
		#make vector of smallest brlen with each mom node as anc
		#calculate, simulatenously, the changes in debt and branch lengthening required as go down tree
		#change branch lengths; hypothetically, debt should then equal zero...
		if(length(mom)>1){for(i in 2:length(mom)){
			small<-min(timetree$edge.length[timetree$edge[,1]==mom[i]])
			mom_blen<-timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]==mom[i-1]]
			debt[i]<-max(debt[i-1] - max(mom_blen-mbl,0),0) + max(mbl-small,0) 
			timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]==mom[i-1]] <- 
			mom_blen - max(min(max(mom_blen-mbl,0),debt[i-1]),0) + max(mbl-small,0)
			timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]!=mom[i-1]] <-  
			timetree$edge.length[timetree$edge[,1]==mom[i] & timetree$edge[,2]!=mom[i-1]] + debt[i]
			}}
		}
	return(timetree)
	}
