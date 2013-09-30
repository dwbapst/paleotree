timeSliceTree<-function(ttree,sliceTime,drop.extinct=FALSE,plot=TRUE){
	#take a phylogeny and produce a phylogenetic 'slice' at time X (w/respect to root.time)
		#lineages extant at sliceTime sliced to end at that point
		#if no root.time, then it is assumed the tip furthest from the root is at 0 (present-day)
			#a warning will be issued if this is assumed
		#extinct tips will be dropped if drop.extinct=TRUE
	#require(ape)
	if(class(ttree)!="phylo"){stop("Error: ttree is not of class phylo")}
	if(is.null(ttree$root.time)){
		ttree$root.time<-max(dist.nodes(ttree)[1:Ntip(ttree),Ntip(ttree)+1])
		message("Warning: no ttree$root.time! Assuming latest tip is at present (time=0)")
		}
	tslice<-ttree$root.time-sliceTime	#time from root to slice time
	#first let's drop all edges that branch later than the slice
	#make them all single lineages by dropping all but one taxon
	dnode<-dist.nodes(ttree)[,Ntip(ttree)+1]
	#identify the ancestor nodes of edges which cross the tslice
	cedge<-which(sapply(1:Nedge(ttree),function(x) any(ttree$edge[x,1]==which(dnode<tslice))
			& any(ttree$edge[x,2]==which(dnode>=tslice))))
	droppers<-numeric()
	for(i in 1:length(cedge)){
		desc<-ttree$edge[cedge[i],2]
		if(desc>Ntip(ttree)){	#if an internal edge that goes past the tslice
			desctip<-unlist(prop.part(ttree)[desc-Ntip(ttree)])	#drop all but one tip
			droppers<-c(droppers,desctip[-1])
		}}
	stree<-drop.tip(ttree,droppers)
	#which edges cross over tslice?
	dnode<-dist.nodes(stree)[,Ntip(stree)+1]
	cedge<-sapply(1:Nedge(stree),function(x) any(stree$edge[x,2]==which(dnode>=tslice)))
	cnode_depth<-dnode[stree$edge[cedge,1]]
	stree$edge.length[cedge]<-tslice-cnode_depth
	stree$root.time<-ttree$root.time		#root.time shouldn't (can't?) change!
	if(drop.extinct){
		stree1<-dropExtinct(stree,ignore.root.time=TRUE)
	}else{stree1<-stree}
	if(plot){layout(1:2);plot(ladderize(ttree),show.tip.label=FALSE);axisPhylo();
		plot(ladderize(stree1),show.tip.label=FALSE);layout(1)}
	return(stree1)
	}