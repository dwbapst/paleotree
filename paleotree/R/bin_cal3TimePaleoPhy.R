bin_cal3TimePaleoPhy<-function(tree,timeList,brRate,extRate,sampRate,ntrees=1,nonstoch.bin=FALSE,
		sites=NULL,point.occur=FALSE,anc.wt=1,node.mins=NULL,rand.obs=FALSE,FAD.only=FALSE,
		adj.obs.wt=TRUE,root.max=200,step.size=0.1,randres=FALSE,plot=FALSE){
	#see the bin_cal3 function for more notation...
	#require(ape)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(class(timeList[[1]])!="matrix"){if(class(timeList[[1]])=="data.frame"){timeList[[1]]<-as.matrix(timeList[[1]])
		}else{stop("Error: timeList[[1]] not of matrix or data.frame format")}}
	if(class(timeList[[2]])!="matrix"){if(class(timeList[[2]])=="data.frame"){timeList[[2]]<-as.matrix(timeList[[2]])
		}else{stop("Error: timeList[[2]] not of matrix or data.frame format")}}
	if(ntrees<1){stop("Error: ntrees<1")}
	if(ntrees==1){message("Warning: Do not interpret a single cal3 time-scaled tree")}
	if(ntrees==1 & !nonstoch.bin){
		message("Warning: Do not interpret a single tree; dates are stochastically pulled from uniform distributions")}
	if(rand.obs & FAD.only){stop("Error: rand.obs and FAD.only cannot both be true")}
	if(!is.null(sites) & point.occur){stop("Error: Inconsistent arguments, point.occur=TRUE will replace input 'sites' matrix")}
	#clean out all taxa which are NA or missing for timeData
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeList[[2]][,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		tree<-drop.tip(tree,droppers)
		if(is.null(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data!")}
		timeList[[2]][which(!sapply(rownames(timeList[[2]]),function(x) any(x==tree$tip.label))),1]<-NA
		}
	timeList[[2]]<-timeList[[2]][!is.na(timeList[[2]][,1]),]
	if(any(is.na(timeList[[2]]))){stop("Weird NAs in Data??")}
	if(any(apply(timeList[[1]],1,diff)>0)){stop("Error: timeList[[1]] not in intervals in time relative to modern")}
	if(any(timeList[[1]][,2]<0)){stop("Error: Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeList[[2]],1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeList[[2]][,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	if(is.null(sites)){
		if(point.occur){
			if(any(timeList[[2]][,1]!=timeList[[2]][,2])){
				stop("Error: point.occur=TRUE but some taxa have FADs and LADs listed in different intervals?!")}
			sites<-matrix(c(1:Ntip(tree),1:Ntip(tree)),Ntip(tree),2)
		}else{
			sites<-matrix(1:(Ntip(tree)*2),,2)
			}
	}else{	#make sites a bunch of nicely behaved sorted integers
		sites[,1]<-sapply(sites[,1],function(x) which(x==sort(unique(as.vector(sites)))))
		sites[,2]<-sapply(sites[,2],function(x) which(x==sort(unique(as.vector(sites)))))
		}
	ttrees<-rmtree(ntrees,3)
	siteTime<-matrix(,max(sites),2)
	for (i in unique(as.vector(sites))){		#build two-col matrix of site's FADs and LADs
		go<-timeList[[2]][which(sites==i)[1]]	#find an interval for this site
		siteTime[i,]<-timeList[[1]][go,]
		}
	for(ntrb in 1:ntrees){
		if(!nonstoch.bin){
			bad_sites<-unique(as.vector(sites))
			siteDates<-apply(siteTime,1,function(x) runif(1,x[2],x[1]))
			while(length(bad_sites)>0){
				siteDates[bad_sites]<-apply(siteTime[bad_sites,],1,function(x) runif(1,x[2],x[1]))
				bad_sites<-unique(as.vector(sites[(siteDates[sites[,1]]-siteDates[sites[,2]])<0,]))
				}
			timeData<-cbind(siteDates[sites[,1]],siteDates[sites[,2]])
		}else{
			timeData<-cbind(siteTime[sites[,1],1],siteTime[sites[,2],2])
			}
		rownames(timeData)<-rownames(timeList[[2]])
		#if(rand.obs){timeData[,2]<-apply(timeData,1,function(x) runif(1,x[2],x[1]))}
		#if(FAD.only){timeData[,2]<-timeData[,1]}
		tree2<-suppressMessages(cal3TimePaleoPhy(tree,timeData,brRate=brRate,extRate=extRate,sampRate=sampRate,
			ntrees=1,anc.wt=anc.wt,node.mins=node.mins,adj.obs.wt=adj.obs.wt,root.max=root.max,step.size=step.size,
			FAD.only=FAD.only,rand.obs=rand.obs,randres=randres,plot=plot))
		tree2$ranges.used<-timeData
		names(tree2$edge.length)<-NULL
		ttrees[[ntrb]]<-tree2
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}