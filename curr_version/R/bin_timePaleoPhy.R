bin_timePaleoPhy<-function(tree,timeList,type="basic",vartime=NULL,ntrees=1,nonstoch.bin=FALSE,randres=FALSE,timeres=FALSE,
	sites=NULL,point.occur=FALSE,add.term=FALSE,inc.term.adj=FALSE,rand.obs=FALSE,node.mins=NULL,plot=FALSE){
	#wrapper for applying non-SRC time-scaling to timeData where FADs and LADs are given as bins 
		#see timePaleoPhy function for more details
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
	#sites is a matrix, used to indicate if binned FADs or LADs of multiple species were obtained from the locality / time point
			#i.e. the first appearance of species A, B and last appearance of C are all from the same lagerstatten
			#this will fix these to always have the same date relative to each other across many trees
			#this will assume that species listed for a site all are listed as being from the same interval...
				#this function also assumes that the sites matrix is ordered exactly as the timeList data is
	#if rand.obs=TRUE, the the function assumes that the LADs in timeList aren't where you actually want the tips
		#instead, tips will be randomly placed anywhere in that taxon's range with uniform probability
		#thus, tip locations will differ slightly for each tree in the sample
		#this is useful when you have a specimen or measurement but you don't know its placement in the species' range
	#default options
	#type="basic";vartime=NULL;ntrees=1;nonstoch.bin=FALSE;randres=FALSE;timeres=FALSE
	#sites=NULL;point.occur=FALSE;add.term=FALSE;inc.term.adj=FALSE;rand.obs=FALSE;node.mins=NULL;plot=FALSE
	#require(ape)
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	if(class(timeList[[1]])!="matrix"){if(class(timeList[[1]])=="data.frame"){timeList[[1]]<-as.matrix(timeList[[1]])
		}else{stop("Error: timeList[[1]] not of matrix or data.frame format")}}
	if(class(timeList[[2]])!="matrix"){if(class(timeList[[2]])=="data.frame"){timeList[[2]]<-as.matrix(timeList[[2]])
		}else{stop("Error: timeList[[2]] not of matrix or data.frame format")}}
	if(ntrees==1 & !nonstoch.bin){
		message("Warning: Do not interpret a single tree; dates are stochastically pulled from uniform distributions")}
	if(ntrees<1){stop("Error: ntrees<1")}
	#clean out all taxa which are NA or missing for timeData
	if(ntrees==1 & randres){message("Warning: Do not interpret a single randomly-resolved tree")}
	if(randres & timeres){stop(
		"Error: Inconsistent arguments: You cannot randomly resolve polytomies and resolve with respect to time simultaneously!")}
	if(!is.null(sites) & point.occur){stop("Error: Inconsistent arguments, point.occur=TRUE will replace input 'sites' matrix")}
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeList[[2]][,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		tree<-drop.tip(tree,droppers)
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
				#print(length(bad_sites))
				}
			timeData<-cbind(siteDates[sites[,1]],siteDates[sites[,2]])
		}else{
			timeData<-cbind(siteTime[sites[,1],1],siteTime[sites[,2],2])
			}
		rownames(timeData)<-rownames(timeList[[2]])
		#if(rand.obs){timeData[,2]<-apply(timeData,1,function(x) runif(1,x[2],x[1]))}
		tree1<-tree
		if(!is.binary.tree(tree)){
			if(randres){tree1<-multi2di(tree)}
			if(timeres){tree1<-timeLadderTree(tree,timeData)}	
			}
		tree2<-suppressMessages(timePaleoPhy(tree1,timeData,type=type,vartime=vartime,ntrees=1,
			randres=FALSE,add.term=add.term,inc.term.adj=inc.term.adj,rand.obs=rand.obs,
			node.mins=node.mins,plot=plot))
		tree2$ranges.used<-timeData
		names(tree2$edge.length)<-NULL
		ttrees[[ntrb]]<-tree2
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}