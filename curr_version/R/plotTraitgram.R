plotTraitgram<-function(trait,tree,trait.name="'trait'",conf.int=TRUE,lwd=1.5){
	#traitgram plotted using ML ASR from geiger (or ace() from ape if ci=TRUE)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(is.null(tree$root.time)){tree$root.time<-max(dist.nodes(tree)[Ntip(tree)+1,1:Ntip(tree)])}
	times<-tree$root.time-dist.nodes(tree)[Ntip(tree)+1,]
	if(conf.int){
		asr<-ace(trait,tree,method="pic")
		tr1<-c(trait,asr$ace);edges<-tree$edge
		ci<-asr$CI95;tr2<-c(tr1,ci)
		plot(1,1,type="n",
			xlim=c(min(tr2)-0.1,max(tr2)+0.1),
			ylim=c(max(times),min(times)),
			xlab="Trait Values",ylab="Time (Before Present)",
			main=paste("Traitgram of",trait.name))
		for(i in 1:nrow(edges)){
			anc<-edges[i,1];desc<-edges[i,2]
			lines(c(tr1[anc],tr1[desc]),c(times[anc],times[desc]),lwd=lwd)
			}
		for(i in 1:Nnode(tree)){
			lines(c(ci[i,1],ci[i,2]),c(times[i+Ntip(tree)],times[i+Ntip(tree)]),lwd=lwd)
			}
	}else{
		#require(geiger)	#not anymore, not as of 04-03-13
		tr1<-c(trait,ace(trait,tree,method="pic")$ace);edges<-tree$edge
		plot(1,1,type="n",xlim=c(min(tr1)-0.1,max(tr1)+0.1),ylim=c(max(times),min(times)),
			xlab="Trait Values",ylab="Time  (Before Present)",
			main=paste("Traitgram of",trait.name))
		for(i in 1:nrow(edges)){anc<-edges[i,1];desc<-edges[i,2]
			lines(c(tr1[anc],tr1[desc]),c(times[anc],times[desc]),lwd=lwd)}
		}
	}