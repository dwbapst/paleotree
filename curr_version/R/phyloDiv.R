phyloDiv<-function(tree,int.length=1,int.times=NULL,plot=TRUE,plotLogRich=FALSE,
		drop.ZLB=TRUE,timelims=NULL){
	#function that computes a diversity curve from a tree file
		#aka lineage-through-time plot
	#root.time
		#ttree$root.time is used to place the diversity curve in time
		#if no root.time, then it is assumed latest tip is at 0 time (present day)
	#time interval starts and ends can be pre-input as a 2 column matrix
		#int.length is ignored in this case
	#this function will randomly resolve any tree it is given using multi2di()
		#this shouldn't affect anything to my knowledge
	#this function also automatically drops zero-length branches from the tree
		#this is advised for paleo-tree analyses of diversification
	#output (if TRUE) is 3 col matrix of bin-start, bit-end, div
	#plotLogRich just decides if the div plot if log-scale or not on the y axis
	#require(ape)
	ttree<-tree
	if(class(ttree)!="phylo"){stop("Error: ttree is not of class phylo")}
	tblen<-int.length
	if(drop.ZLB){ttree<-dropZLB(ttree)}
	savetree<-ttree
	if(!is.binary.tree(ttree)){ttree<-multi2di(ttree)}
	if(is.null(ttree$root.time)){
		ntime<-dist.nodes(ttree)[,Ntip(ttree)+1]
		ntime<-max(ntime)-ntime
	}else{
		ntime<-dist.nodes(ttree)[,Ntip(ttree)+1]
		ntime<-ttree$root.time-ntime
		ntime<-round(ntime,6)
		if(min(ntime)<0){stop("Error: tree$root.time is less than total depth of tree!")}
		}
	if(is.null(int.times)){
		midtimes<-seq(max(ntime)+3*tblen,min(ntime)-2*tblen,by=-tblen)
		midtimes<-midtimes[midtimes>0]
		int.start<-midtimes+(tblen/2)
		int.end<-midtimes-(tblen/2)
	}else{
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	LAD<-ntime[1:Ntip(ttree)]				#death
	FAD<-ntime[(Ntip(ttree)+1):length(ntime)]		#birth
	div<-sapply(1:length(midtimes),function(x) 1+sum(FAD>=int.end[x])-sum(LAD>int.start[x]))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		layout(matrix(1:2,2,1))
		parOrig<-par(mar=c(1,4,1,1))
		plot(ladderize(savetree),show.tip.label=FALSE)
		par(mar=c(5,4,2,2))
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Lineage Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				ylim=c(0,max(div1)+1),
				xlab="Time (Before Present)",ylab="Lineage Richness")
			}
		par(parOrig);layout(1)
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}