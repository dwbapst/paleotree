compareNodeAges<-function(tree1,tree2,dropUnshared=FALSE){
	#output vector of shifts in node dates
		#08-02-12: Allows multiple trees2 to be multiple trees
			#will produce a matrix, each row is a tree in tree2, each column a different but commonly shared clade
	#require(ape)
	if(class(tree1)!="phylo"){stop("Error: tree1 is not of class phylo")}
	tree1orig<-tree1
	if(class(tree2)!="phylo"){
		if(class(tree2)!="multiPhylo"){stop("Error: tree2 is not of class phylo or multiphylo")}
		trees2<-tree2
		#if it isn't multiphylo, make it into one!
		}else{trees2<-list(tree2);class(trees2)<-"multiPhylo"}
	#okay, need to find all matches common to tree1 and tree2
		#we'll make a MATRIX of all clades held in common between each tree in trees2 and to tree1
		#each row will be a different tree, each column a different clade
		#node dates for each tree will get added as a new column or to an old column
	matchMat<-NULL
	for(i in 1:length(trees2)){
		tree2<-trees2[[i]]
		tree1<-tree1orig
		#incredibly, all of the following is necessary to properly adjust node dates (??!!)
		matches1<-which(!is.na(match(tree1$tip.label,tree2$tip.label)))[1]
		if(length(matches1)<1){stop(paste("Error: No shared taxa between tree2 and tree1[[",i,"]]!"))}
		tipmatch<-tree1$tip.label[matches1]
		mtimeA<-dist.nodes(tree1)[matches1,Ntip(tree1)+1]
		mtimeB<-dist.nodes(tree2)[match(tipmatch,tree2$tip.label),Ntip(tree2)+1]
		tree1<-drop.tip(tree1,tree1$tip.label[is.na(match(tree1$tip.label,tree2$tip.label))])
		tree2<-drop.tip(tree2,tree2$tip.label[is.na(match(tree2$tip.label,tree1$tip.label))])
		ntime1<-dist.nodes(tree1)[,Ntip(tree1)+1]
		ntime2<-dist.nodes(tree2)[,Ntip(tree2)+1]
		mtime1<-ntime1[match(tipmatch,tree1$tip.label)]
		mtime2<-ntime2[match(tipmatch,tree2$tip.label)]
		if(!is.null(tree1$root.time)){
			tree1$root.time<-tree1$root.time-(mtimeA-mtime1)
			ntime1<-tree1$root.time-ntime1
			ntime1<-round(ntime1,6)
			if(min(ntime1)<0){stop(paste("Error: tree1$root.time is less than total depth of tree1!"))}
		}else{
			ntime1<-max(ntime1)-ntime1
			}
		if(!is.null(tree2$root.time)){
			tree2$root.time<-tree2$root.time-(mtimeB-mtime2)
			ntime2<-tree2$root.time-ntime2
			ntime2<-round(ntime2,6)
			if(min(ntime2)<0){stop("Error: tree2[",i,"]$root.time is less than total depth of that tree!")}
		}else{
			ntime2<-max(ntime2)-ntime2
			}
		clades1<-lapply(prop.part(tree1),function(x) sort(tree1$tip.label[x]))
		clades2<-lapply(prop.part(tree2),function(x) sort(tree2$tip.label[x]))
		matches<-match(clades1,clades2)
		cladesMatches<-clades2[matches[!is.na(matches)]]
		if(length(matches[!is.na(matches)])==1){cladesMatches<-list(cladesMatches)}
		ages1<-ntime1[Ntip(tree1)+which(!is.na(matches))]
		ages2<-ntime2[Ntip(tree2)+matches[!is.na(matches)]]
		age_diff<-ages1-ages2
		names(age_diff)<-NULL
		#okay, need to find all matches common to tree1 and tree2
			#we'll make a MATRIX of all clades held in common between each tree in trees2 and to tree1
			#each row will be a different tree, each column a different clade
			#node dates for each tree will get added as a new column or to an old column
		cladesMatches<-sapply(cladesMatches,function(x) paste(x,collapse=","))
		if(is.null(matchMat)){	#if the first tree examined...
			matchMat<-matrix(age_diff,1,)
			colnames(matchMat)<-cladesMatches
		}else{
			currMatches<-match(cladesMatches,colnames(matchMat))
			matchMat<-rbind(matchMat,rep(NA,ncol(matchMat)))
			for(j in 1:length(currMatches)){
				if(!is.na(currMatches[j])){
					matchMat[i,currMatches[j]]<-age_diff[j]
				}else{
					matchMat<-cbind(matchMat,c(rep(NA,i-1),age_diff[j]))
					colnames(matchMat)[ncol(matchMat)]<-cladesMatches[j]
					}
				}
			}
		}
	if(dropUnshared){
		matchMat<-matchMat[,apply(matchMat,2,function(x) all(!is.na(x))),drop=FALSE]
		}
	rownames(matchMat)<-names(trees2)
	if(length(trees2)==1){
		matchMat<-matchMat[1,]
		names(matchMat)<-NULL
		}
	return(matchMat)
	}