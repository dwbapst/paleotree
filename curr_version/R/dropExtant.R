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