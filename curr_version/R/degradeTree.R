degradeTree<-function(tree,prop_collapse,nCollapse=NULL,node.depth=NA,leave.zlb=FALSE){
	#collapses a given proportion of internal edges, creating polytomies
		#node.depth conditions on depth of edge in tree
			# 1 removes more shallow nodes, 0 removes deeper nodes
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
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