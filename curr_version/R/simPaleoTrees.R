simPaleoTrees<-function(p,q,r,ntrees=1,all.extinct=FALSE,modern.samp.prob=1.0,mintime=1,maxtime=100,
	mintaxa=2,maxtaxa=500,anag.rate=0,prop.bifurc=0,prop.cryptic=0,drop.zlb=TRUE,print.runs=FALSE,
	plot=FALSE){
	#this is a wrapper which will create many paleo trees with at least two observed tips
		#uses simFossilTaxa, sampRanges,taxa2phylo, etc
		#good if you want to simulate many many trees with extinct taxa
		#divergence times for nodes will be perfectly known
	#by default: (minimal conditioning)
		#no conditioning on the number of extant taxa
		#living taxa are sampled perfectly at the present
		#zero-length branches are dropped
	#simPaleoTrees(p=0.1,q=0.1,r=0.1,ntrees=10)
	#require(ape)
	if(mintaxa<2){stop("Error: Need at least two taxa per tree; increase mintaxa")}
	if(ntrees<1){stop("Error: number of trees to simulate is <1")}
	res<-rmtree(ntrees,2)
	ntries<-0
	for(i in 1:ntrees){
		rerun<-TRUE
		while(rerun){
			ntries<-ntries+1
			taxa<-suppressMessages(simFossilTaxa(p=p,q=q,anag.rate=anag.rate,prop.bifurc=prop.bifurc,prop.cryptic=prop.cryptic,nruns=1,mintaxa=mintaxa,
				maxtaxa=maxtaxa,maxtime=maxtime,maxExtant=ifelse(all.extinct,0,maxtaxa),min.cond=FALSE,plot=plot))
			ranges<-sampleRanges(taxa,r,min.taxa=0,modern.samp.prob=modern.samp.prob,merge.cryptic=FALSE)
			if(sum(!is.na(ranges[,1]))>1){
				tree<-taxa2phylo(taxa,obs_time=ranges[,2],plot=plot)
				if(drop.zlb){tree<-dropZLB(tree)}
				if(all(!is.na(tree))){
					numext<-sum(ranges[,2]==0)
					minta<-Ntip(tree)>mintaxa
					minti<-(max(ranges,na.rm=TRUE)-min(ranges,na.rm=TRUE))>mintime
					if(minta & minti){
						rerun<-FALSE
						}
					}
				}
			}
		tree$taxa<-taxa		#attach original taxon matrix
		tree$ranges<-ranges	#attach sampled ranges
		res[[i]]<-tree
		if(ntrees>10){if((i%%(ntrees/10))==0){message(cat(round(i*100/ntrees),"% ",sep=""))}}
		}
	if(print.runs){message(paste(ntrees," trees accepted from ",ntries," total runs (",signif(ntrees/ntries,2)," Acceptance Probability)",sep=""))}
	if(ntrees==1){res<-res[[1]]}
	return(res)
	}