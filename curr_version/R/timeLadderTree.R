timeLadderTree<-function(tree,timeData){
	#resolves all polytomies in a tree as ladders to match FADs in timeData
	#only applicable to continuous time data
	#require(ape)	
	#first sanitize data
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	if(ncol(timeData)==6){timeData<-timeData[,3:4,drop=FALSE]}	#also allow it to accept taxad objects
	#first clean out all taxa which are NA or missing in timeData
	#remove taxa that are NA or missing in timeData
	tree<-drop.tip(tree,tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))])
	if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data")}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Error: Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	#node IDs of polytomies
	nodes<-Ntip(tree)+1:Nnode(tree)
	polys<-nodes[sapply(nodes,function(x) sum(tree$edge[,1]==x)>2)]
	#error if no polytomies
	if(length(polys)==0){stop("Error: No Polytomies in Cropped Tree?")}
	while(length(polys)>0){
		node<-polys[1]
		#make a list of the first FAD of each descendant
		desc<-tree$edge[tree$edge[,1]==node,2]
		dTips<-lapply(desc,function(x) if(x>Ntip(tree)){
			tree$tip.label[prop.part(tree)[[x-Ntip(tree)]]]
			}else{tree$tip.label[x]})
		#get max FAD
		dFADs<-sapply(dTips,function(x) max(timeData[x,1]))
		#build a ladderized subtree of the node
		subtree<-stree(length(dFADs),type="left")
		#figure out the order of tips, use a second vector of random numbers for ties
		subtree$tip.label<-desc[order(-dFADs,sample(1:length(dFADs)))]
		#stick on tip labels
		for(i in desc){
			dtip<-which(subtree$tip.label==i)
			if(i>Ntip(tree)){		#if its a clade
				subclade<-extract.clade(tree,i)
				subtree<-bind.tree(subtree,subclade,where=dtip)
				subtree<-collapse.singles(subtree)
				}else{subtree$tip.label[dtip]<-tree$tip.label[i]}	#if its a tip
			}
		#replace original node with new, resolved, scaled node
		if(node!=(Ntip(tree)+1)){	#if it isn't the root node
			drtips<-prop.part(tree)[[node-Ntip(tree)]]
			tip_lab<-tree$tip.label[drtips[1]]	#cut out all but one tip, just to put it back together later
			droptree<-collapse.singles(drop.tip(tree,drtips[-1]))
			droptree<-bind.tree(droptree,subtree,where=which(droptree$tip.label==tip_lab))	#put in subtree at tip
			tree<-droptree
		}else{				#if it is the root node
			tree<-subtree
			}	
		nodes<-Ntip(tree)+1:Nnode(tree)
		polys<-nodes[sapply(nodes,function(x) sum(tree$edge[,1]==x)>2)]
		}
	tree$edge.length<-NULL
	return(tree)
	}