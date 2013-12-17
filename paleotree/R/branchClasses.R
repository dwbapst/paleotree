branchClasses<-function(tree,whichExtant=NULL,tol=0.01){
	#names will be time of origin for each branch
	#require(ape)
	if(!is(tree,"phylo")){stop("Error: tree is not of class phylo")}
	dists<-dist.nodes(tree)[Ntip(tree)+1,]
	if(is.null(whichExtant)){
		dists1<-dists[1:Ntip(tree)]
		if(is.null(tree$root.time)){modern<-max(dists1)}else{modern<-tree$root.time}
		modern.tips<-which(dists1>(modern-tol))
		dead.tips<-which(dists1<=(modern-tol))
	}else{
		modern.tips<-which(sapply(tree$tip.label,function(x) any(x==names(which(whichExtant==1)))))
		dead.tips<-which(sapply(tree$tip.label,function(x) any(x==names(which(whichExtant==0)))))
		}
	int.edges<-tree$edge[,2]>Ntip(tree)
	live.edges<-sapply(tree$edge[,2],function(x) any(x==modern.tips))
	dead.edges<-sapply(tree$edge[,2],function(x) any(x==dead.tips))
	if(is.null(tree$root.time)){
		depths<-max(dists)-dists[tree$edge[,1]]
		}else{
		depths<-tree$root.time-dists[tree$edge[,1]]
		}
	brlen.all<-tree$edge.length
		names(brlen.all)<-depths
	brlen.int<-tree$edge.length[int.edges]
		names(brlen.int)<-depths[int.edges]
	brlen.live<-tree$edge.length[live.edges]
		names(brlen.live)<-depths[live.edges]		
	brlen.dead<-tree$edge.length[dead.edges]
		names(brlen.dead)<-depths[dead.edges]
	list(brlen.all=brlen.all,brlen.int=brlen.int,brlen.live=brlen.live,brlen.dead=brlen.dead)
	}