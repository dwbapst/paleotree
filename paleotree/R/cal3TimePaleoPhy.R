cal3TimePaleoPhy<-function(tree,timeData,brRate,extRate,sampRate,ntrees=1,anc.wt=1,node.mins=NULL,
	rand.obs=FALSE,FAD.only=FALSE,adj.obs.wt=TRUE,root.max=200,step.size=0.1,randres=FALSE,plot=FALSE){
	#see SRC function for more notation...
	#function for Ps - use pqr2Ps
	#example data
	#tree<-rtree(10);tree$edge.length<-sample(0:1,Nedge(tree),replace=TRUE);tree<-di2multi(tree)
	#ntrees=2;anc.wt=1;node.mins=NULL;sampRate=rep(0.1,Ntip(tree));names(sampRate)<-tree$tip.label
	#brRate<-extRate<-sampRate
	#timeData<-runif(Ntip(tree),200,400);timeData<-cbind(timeData,timeData-runif(Ntip(tree),1,80))
	#rownames(timeData)<-tree$tip.label;root.max=200;plot=TRUE;rand.obs=FALSE;FAD.only=TRUE
	#step.size=0.1;adj.obs.wt=TRUE;randres=FALSE
	#
	#taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=50,maxtaxa=100,maxtime=1000,maxExtant=0)
	#cladogram<-taxa2cladogram(taxa);timeData<-sampleRanges(taxa,r=0.1)
	#
	###trying to see if adj.wts works
	#tree<-rtree(2);anc.wt=1;node.mins=NULL;brRate<-extRate<-sampRate<-0.1;timeData<-cbind(c(100,95),c(90,85));adj.obs.wt=TRUE
	#rownames(timeData)<-tree$tip.label;root.max=200;plot=TRUE;rand.obs=FALSE;FAD.only=TRUE;ntrees=1;randres=FALSE;step.size=0.1
	#
	#add.zombie=FALSE;node.mins<-c(-sort(-runif(1,600,900)),rep(NA,Nnode(tree)-1))	#assume two very deep divergences
	#
	#require(ape)#;require(phangorn)
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	if(rand.obs & FAD.only){stop("Error: rand.obs and FAD.only cannot both be true")}
	#first clean out all taxa which are NA or missing in timeData
	if(ntrees==1){message("Warning: Do not interpret a single cal3 time-scaled tree")}
	if(ntrees<1){stop("Error: ntrees<1")}
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		tree<-drop.tip(tree,droppers)
		if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data!")}
		timeData[which(!sapply(rownames(timeData),function(x) any(x==tree$tip.label))),1]<-NA
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	#make sure all taxa have a sampRate, brRate and exRate
	if(length(sampRate)==1){sampRate<-rep(sampRate,Ntip(tree));names(sampRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(sampRate))))){
			stop("Sampling Rates Not Given For All Taxa on Tree!")}}
	if(length(brRate)==1){brRate<-rep(brRate,Ntip(tree));names(brRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(brRate))))){
			stop("Branching Rates Not Given For All Taxa on Tree!")}}
	if(length(extRate)==1){extRate<-rep(extRate,Ntip(tree));names(extRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(extRate))))){
			stop("Extinction Rates Not Given For All Taxa on Tree!")}}
	#allow per-taxon anc.wt values
	if(length(anc.wt)==1){anc.wt<-rep(anc.wt,Ntip(tree));names(anc.wt)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(anc.wt))))){
			stop("Ancestral Weights Not Given For All Taxa on Tree!")}}
	Ps<-sapply(tree$tip.label,function(x) pqr2Ps(brRate[x],extRate[x],sampRate[x]))
	names(Ps)<-tree$tip.label
	if(length(node.mins)!=Nnode(tree) & !is.null(node.mins)){stop("node.mins length != Nnode!")}
	ttree1<-timePaleoPhy(tree,timeData,type="basic",node.mins=node.mins,add.term=FALSE,inc.term.adj=FALSE)
	#identify which nodes are min-locked; make sure to update when resolving polytomies
	if(length(node.mins)>0){locked_nodes<-which(!is.na(node.mins))++Ntip(tree)}else{locked_nodes<-NA}
	ttree1<-collapse.singles(ttree1)
	ttrees<-rmtree(ntrees,3)
	for(ntr in 1:ntrees){
		#10/30/12: get FAD, new LAD (time of observation), and then calculate difference between t.obs and LAD
		if(rand.obs | FAD.only){
			if(FAD.only){
				timeData1<-cbind(timeData[,1],timeData[,1],timeData[,1]-timeData[,2])
				}
			if(rand.obs){
				timeData1<-cbind(timeData[,1],apply(timeData,1,function(x) runif(1,x[2],x[1])))
				timeData1<-cbind(timeData1,timeData1[,2]-timeData[,2])
				}
		}else{
			timeData1<-cbind(timeData,0)
			}
		ktree<-ttree1
		if(randres){ktree<-multi2di(ktree)}
		nodes<-(1:Nnode(ktree))+Ntip(ktree)		#get a vector of all internal nodes	
		nodes<-nodes[order(-node.depth(ktree)[-(1:Ntip(ktree))])]	#order by depth
		anags<-character();budds<-character();nAdjZip<-0
		while(length(nodes)>0){		#can't use a for() because # of nodes may change
			#save_tree<-ktree;dev.new();plot(ktree)
			node<-nodes[1]
			tipl<-ktree$tip.label
			#put together tip data: diffLAD is the difference between the time of observation and the true LAD
			tipd<-cbind(ID=(1:Ntip(ttree1)),FAD=(timeData1[tipl,1]),time.obs=(timeData1[tipl,2]),SR=sampRate[tipl],
				BR=brRate[tipl],ER=extRate[tipl],Ps=Ps[tipl],ancWt=anc.wt[tipl],diffLAD=(timeData1[tipl,3]))
			if(node==(Ntip(ktree)+1)){
				min_zip<-(-root.max)	#if root, allow to be push back up to root.max
				stem_len<-root.max
				root_push<--seq(min_zip,0,by=step.size)
			}else{									#if not root, push down to lower node
				min_zip<-(-ktree$edge.length[ktree$edge[,2]==node])
				stem_len<-ktree$edge.length[ktree$edge[,2]==node]
				}		
			dnodes<-ktree$edge[ktree$edge[,1]==node,2]	#find the daughter nodes
			dlen<-ktree$edge.length[match(dnodes,ktree$edge[,2])]	#find the dedges lengths
			minlocked<-ifelse(!all(is.na(locked_nodes)),any(node==locked_nodes),FALSE)	#is this node min-locked?
			#if(any(dnodes==which(ktree$tip.label=="t10"))){break()}
			if(length(dnodes)>2){		#if node is a polytomy, use PARALLEL ZIPPER
				#pick a starting lineage
				dSR<-drng<-dBR<-dER<-dPs<-danc.wt<-ddfLAD<-numeric()
				for(i in dnodes){		#for each desc, get vector of SR for earliest and range if desc is a tip
					dtips<-match(unlist(Descendants(ktree,i)),tipd[,1])
					dearly<-which(tipd[dtips,2]==max(tipd[dtips,2]))[1]
					dSR[length(dSR)+1]<-tipd[dearly,4]
					dBR[length(dBR)+1]<-tipd[dearly,5]
					dER[length(dER)+1]<-tipd[dearly,6]
					dPs[length(dPs)+1]<-tipd[dearly,7]
					danc.wt[length(danc.wt)+1]<-tipd[dearly,8]
					ddfLAD[length(ddfLAD)+1]<-tipd[dearly,9]		#10-30-12: diff between t.obs and LAD
					drng[length(drng)+1]<-ifelse(length(dtips)>1,NA,diff(unlist(tipd[dtips,3:2])))
					}
				#08-01-12: choice of starting lineage doesn't matter (see SRC method)
				dnode1<-dnodes[which(dlen==min(dlen))[1]]		#just pick first appearing
				#make sure to include stem length in calculations!
				#08-03-12: as with new SRC above, the stem length will JUST be the -min_zip: max root.push!
				#if(node==(Ntip(ktree)+1)){
				#	#07-29-12: this treats the stem length as a single gap
				#		#given that this should be a single gap and not gamma distributed
				#	root_density<-dexp(root_push,rate=dSR[dnode1==dnodes]+(dBR[dnode1==dnodes]*dPs[dnode1==dnodes]))
				#	stem_len<-sample(root_push,1,prob=root_density)
				#	}
				#make data structure for placed lineages; anc= row of anc lineage, events in time-from-stem 
				plin<-c(dnode1,(dlen[dnode1==dnodes]+stem_len),drng[dnode1==dnodes],NA,
					0,dlen[dnode1==dnodes]+stem_len,dlen[dnode1==dnodes]+drng[dnode1==dnodes]+stem_len,
					danc.wt[dnode1==dnodes],ddfLAD[dnode1==dnodes])
				plin<-matrix(plin,1,)
				colnames(plin)<-c("node","brl","rng","anc","tSpec","tFO","tLO","ancWt","diffLAD")
				#place additional lineages, in order of max zip, using parallel zippers along placed lineages
				#add_nodes will hold necessary information on dlen, rng, SR for unplaced nodes
				add_nodes<-cbind(dnodes,dlen+stem_len,drng,dSR,dBR,dER,dPs,danc.wt,ddfLAD)[-which(dnodes==dnode1),]
				add_nodes<-add_nodes[order(dlen[-which(dnodes==dnode1)]),]
				colnames(add_nodes)<-c("node","dstem","rng","SR","BR","ER","Ps","ancWt","diffLAD")	#dstem is distance from stem to FAD
				for(i in 1:nrow(add_nodes)){
					#identify each currently placed lineage, produce zipper for each relative to stem point
					zips<-matrix(,1,3)
					for(j in 1:nrow(plin)){
						min_zip<-plin[j,5]	#time of speciation (stem time for each placed lineage)
						max_zip<-ifelse(plin[j,8]>0 & !is.na(plin[j,7]),	#max zip is complicated, dependent on anc.wt
							min(add_nodes[i,2],plin[j,7]),	#min of LAD of the placed lineage or FAD of the lineage to be placed
							min(add_nodes[i,2],plin[j,6]))	#min of FAD of the placed lineage or FAD of the lineage to be placed
						if(minlocked){max_zip<-stem_len}		#if minlocked, node must be not earlier than original node time
						#and on a stemTime=0 scale as used for the parallel zipper, stem_len is the original node time, unlike single zipper below
						poss_zip<-seq(min_zip,max_zip,by=step.size)
						#08-01-12: getting the density via the Cal3 algorithm
						gap1<-plin[j,6]-poss_zip		#waiting time from FAD1 to zip
						gap1<-ifelse(gap1>0,gap1,0)		#waiting time from branching point to FAD1
						gap2<-add_nodes[i,2]-poss_zip			#waiting time from branching point to FAD2
						gapStem2zip<-ifelse(plin[j,6]>poss_zip,poss_zip,plin[j,6])	#waiting time from stem to FAD1 OR br node
						totalgap<-gap1+gap2+gapStem2zip
						#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
							#under different combinations with p=q=r,p=q>r and p=q<r
							#admittedly not a perfect fit, though!
							#use rates from the lineages being added
						linDense<-dgamma(totalgap,shape=2,rate=add_nodes[i,4]+(add_nodes[i,5]*add_nodes[i,7]))
						linDense<-ifelse(poss_zip>plin[j,6],
							plin[j,8]*linDense,	#with anc.wt (of the PLACED LINEAGE)
								    linDense)	#without anc.wt, but with gap1 probability
						#10-30-12 need to correct linDense value for time of observation if difference from LAD =/= 0
							#use the diffLAD of the PLACED LINEAGE (CAUSE ITS THE ANCESTOR)
							#if adj.obs.wt, if anc.wt>0 & diffLAD>0 & plin's LAD is before the add_nodes's FAD...
						if(ifelse(is.na(plin[j,7]),FALSE,add_nodes[i,2]>(plin[j,7]+step.size)) & 
								adj.obs.wt & plin[j,8]>0 & plin[j,9]>0){	
							adj_max<-max(add_nodes[i,2],plin[j,9]+plin[j,7])
							adj_zips<-seq(plin[j,7]+step.size,adj_max,by=step.size) 	#the zips we ain't looking at
							#08-01-12: getting the density via the Cal3 algorithm
							gap1<-plin[j,6]-adj_zips		#waiting time from FAD1 to zip
							gap1<-ifelse(gap1>0,gap1,0)		#waiting time from branching point to FAD1
							gap2<-add_nodes[i,2]-adj_zips			#waiting time from branching point to FAD2
							gapStem2zip<-ifelse(plin[j,6]>adj_zips,adj_zips,plin[j,6])	#waiting time from stem to FAD1 OR br node
							adj_totalgap<-gap1+gap2+gapStem2zip
							#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
								#under different combinations with p=q=r,p=q>r and p=q<r
								#admittedly not a perfect fit, though!
							adj_linDense<-dgamma(adj_totalgap,shape=2,rate=add_nodes[i,4]+(add_nodes[i,5]*add_nodes[i,7]))
							linDense[length(linDense)]<-linDense[length(linDense)]+sum(adj_linDense*plin[j,8])
							nAdjZip<-nAdjZip+1
							}
						new_zip<-cbind(linDense,plin[j,1],poss_zip)
						zips<-rbind(zips,new_zip)
						}
					if(nrow(zips)<3){zips<-matrix(zips[-1,],1,3)}else{zips<-zips[-1,]}
					colnames(zips)<-c("linDensity","anc","tzip")
					#new as of 08-21-12
					zip_prob<-zips[,1]
					zip_prob[is.na(zip_prob)]<-0
					if(sum(zip_prob)==0){zip_prob<-rep(1,length(zip_prob))}
					zip_prob<-zip_prob/sum(zip_prob)
					ch_zip<-sample(1:nrow(zips),1,prob=zip_prob)	#sample zips
					ch_anc<-zips[ch_zip,2]
					ch_tzip<-zips[ch_zip,3]
					#if anagenesis, add to anags; if budding, add to budds
					if(!is.na(plin[ch_anc==plin[,1],7])){	#if the anc is terminal
						if(plin[ch_anc==plin[,1],7]==ch_tzip){anags<-c(anags,ktree$tip.label[ch_anc])}	#if anagenetic
						if(plin[ch_anc==plin[,1],6]<ch_tzip){budds<-c(budds,ktree$tip.label[ch_anc])}		#if budding
						}
					new_lin<-c(add_nodes[i,1],add_nodes[i,2]-ch_tzip,add_nodes[i,3],
						ch_anc,ch_tzip,add_nodes[i,2],add_nodes[i,2]+add_nodes[i,3],add_nodes[i,8],add_nodes[i,9])
					plin<-rbind(plin,new_lin)	#put in plin
					}
				#turn into a subtree using taxa2phylo()
				taxad_o<-t(apply(plin,1,function(x) c(x[1],x[4],x[5],ifelse(is.na(x[7]),x[6],x[7]))))
				new_anc<-sapply(taxad_o[,2],function(x) ifelse(is.na(x),NA,which(taxad_o[,1]==x)))
				taxad_n<-cbind(1:nrow(taxad_o),new_anc,taxad_o[,3:4])
				taxad_n[,3:4]<-max(taxad_n[,3:4])-taxad_n[,3:4]
				rownames(taxad_n)<-paste("t",1:nrow(taxad_o),sep="")
				subtree<-taxa2phylo(taxad_n)
				subtree$tip.label<-taxad_o[match(subtree$tip.label,paste("t",1:nrow(taxad_o),sep="")),1]
				new_stem<-diff(sort(plin[,5]))[1]	#time to stem (branch length for stem)
				#stick desc tips onto the subtree
				for(i in dnodes){
					dtip<-which(subtree$tip.label==i)
					if(i>Ntip(ktree)){		#if its a clade
						subclade<-extract.clade(ktree,i)
						subtree<-bind.tree(subtree,subclade,where=dtip)
						subtree<-collapse.singles(subtree)
					}else{				#if its a tip
						subtree$tip.label[dtip]<-ktree$tip.label[i]
					}}
				#replace original node with new, resolved, scaled node
				if(node!=(Ntip(ktree)+1)){	#if it isn't the node
					drtips<-prop.part(ktree)[[node-Ntip(ktree)]]
					tip_lab<-ktree$tip.label[drtips[1]]	#I need to cut out all but one tip, for the sole purpose of putting it all back later)
					droptree<-collapse.singles(drop.tip(ktree,drtips[-1]))
					droptree$edge.length[droptree$edge[,2]==which(droptree$tip.label==tip_lab)]<-new_stem	#reset edge length leading to remaining tip to new_stem
					droptree<-bind.tree(droptree,subtree,where=which(droptree$tip.label==tip_lab))	#put in subtree at tip
					ktree1<-droptree
				}else{				#if it is the node
					ktree1<-subtree
					}				
				#once you've changed the structure of the tree find original nodes in new tree (surprisingly frustating to code!)
				if(length(nodes)>1){
					d_o<-lapply(Descendants(ktree,nodes[-1]),function(x) ktree$tip.label[x])
					d_n<-lapply(Descendants(ktree1)[-(1:Ntip(ktree1))],function(x) ktree1$tip.label[x])
					nodes1<-sapply(d_o,function(x) which(sapply(d_n,function(y) 
						ifelse(length(y)==length(x),all(sort(y)==sort(x)),FALSE))))
					nodes1<-Ntip(ktree1)+nodes1
					nodes1<-nodes1[order(-node.depth(ktree1)[nodes1])]	#order by depth
					if(!all(is.na(locked_nodes))){	#update locked_nodes, can re-use d_n
						d_ol<-lapply(Descendants(ktree,locked_nodes),function(x) ktree$tip.label[x])
						locked_nodes<-sapply(d_ol,function(x) which(sapply(d_n,function(y) 
							ifelse(length(y)==length(x),all(sort(y)==sort(x)),FALSE))))
						locked_nodes<-Ntip(ktree1)+locked_nodes
						}					
				}else{nodes1<-numeric()}				#don't bother if no more nodes left...
				#layout(matrix(1:2,2,));plot(save_tree);plot(ktree1);layout(1)
				#update tipd and nodes (tree str will have changed)
				ktree1<-collapse.singles(ktree1)
				ktree<-ktree1
				nodes<-nodes1
			}else{	#if node is NOT a polytomy, then use regular zipper
				dlen1<-min(dlen)				#the shortest branch (nothing can be done about this one)
				dlen2<-max(dlen)				#the longest branch
				dnode1<-dnodes[which(dlen==dlen1)[1]]
				dnode2<-dnodes[dnodes!=dnode1]
				#get the SR of earliest tip of d2
				d2FADs<-tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),2]	#need the NODE for prop.part, idiot!
				d2early<-which(d2FADs==max(d2FADs))	
				d2early<-ifelse(length(d2early)>1,sample(d2early,1),d2early)	#if more than one of same FAD, just randomly choose one
				d2SR<-(tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),4])[d2early]
				d2BR<-(tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),5])[d2early]
				d2Ps<-(tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),7])[d2early]
				d1rng<-ifelse(dnode1<=Ntip(ktree),diff(unlist(tipd[match(dnode1,tipd[,1]),3:2])),NA)	#if clade, range = NA	
				d2rng<-ifelse(dnode2<=Ntip(ktree),diff(unlist(tipd[match(dnode2,tipd[,1]),3:2])),NA)
				d1ancWt<-ifelse(dnode1<=Ntip(ktree),unlist(tipd[match(dnode1,tipd[,1]),8]),0)
				d1diffLAD<-ifelse(dnode1<=Ntip(ktree),unlist(tipd[match(dnode1,tipd[,1]),9]),0)
				#ZIPPPER: first create a list of scenarios, treat the position of the node like a zipper which can be moved up or down
				max_zip<-ifelse(d1ancWt>0 & !is.na(d1rng),min(dlen1+d1rng,dlen2),dlen1)	#if d1 isn't a clade and there can be ancestors	
				if(minlocked){max_zip<-0}		#in single zipper, unlike parallel zipper, 0 is original node time
				poss_zip<-seq(min_zip,max_zip,by=step.size)
				#waiting time from zip (branching point) to FAD1
				gap1<-ifelse(poss_zip>dlen1,0,dlen1-poss_zip)		#restrict to be dlen1 if longer than dlen1 (FAD1-time, prob 0 anyway!)
				gap2<-dlen2-poss_zip			#waiting time from branching point to FAD2
				#waiting time from stem to FAD1 OR branch node (zip)
				gapStem2zip<-ifelse((stem_len+dlen1)>(stem_len+poss_zip),
					stem_len+poss_zip,stem_len+dlen1)	
				totalgap<-gap1+gap2+gapStem2zip
				#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
					#under different combinations with p=q=r,p=q>r and p=q<r
						#admittedly not a perfect fit, though!
				linDense<-dgamma(totalgap,shape=2,rate=d2SR+(d2BR*d2Ps))
				linDensity<-ifelse((stem_len+dlen1)<(stem_len+poss_zip),
					d1ancWt*linDense,	#with anc.wt (of first taxon)
					 	  linDense)	#without anc.wt, but with gap1 probability
				#10-30-12 need to correct linDense value for time of observation if difference from LAD =/= 0
					#use the diffLAD of the potential ancestor
					#if adj.obs.wt, if anc.wt>0 & diffLAD>0 & d1's LAD is before d2's FAD...
				if(adj.obs.wt & d1ancWt>0 & d1diffLAD>0 & (dlen1+d1rng+step.size)<dlen2){
					adj_max<-max(dlen2,d1diffLAD+(dlen1+d1rng))
					adj_zips<-seq((dlen1+d1rng)+step.size,adj_max,by=step.size) 	#the zips we ain't looking at
					#08-01-12: getting the density via the Cal3 algorithm
					#waiting time from zip (branching point) to FAD1
					gap1<-ifelse(adj_zips>dlen1,0,dlen1-adj_zips)		#restrict to be dlen1 if longer than dlen1 (FAD1-time, prob 0 anyway!)
					gap2<-dlen2-adj_zips			#waiting time from branching point to FAD2
					#waiting time from stem to FAD1 OR branch node (zip)
					gapStem2zip<-ifelse((stem_len+dlen1)>(stem_len+adj_zips),
						stem_len+adj_zips,stem_len+dlen1)
					adj_totalgap<-gap1+gap2+gapStem2zip
					#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
						#under different combinations with p=q=r,p=q>r and p=q<r
						#admittedly not a perfect fit, though!
					adj_linDense<-dgamma(totalgap,shape=2,rate=d2SR+(d2BR*d2Ps))
					linDensity[length(linDensity)]<-linDensity[length(linDensity)]+sum(adj_linDense*d1ancWt)
					nAdjZip<-nAdjZip+1
					}
				#new as of 08-21-12
				linDensity[is.na(linDensity)]<-0
				if(sum(linDensity)==0){linDensity[length(linDensity)]<-1}
				linDensity1<-linDensity/sum(linDensity)
				ch_zip<-sample(poss_zip,1,prob=linDensity1)	#pick zipper location
				#calculate new branch lengths, adding terminal ranges to tips
				new_dlen1<-ifelse(ch_zip>dlen1,NA,dlen1-ch_zip)			#If not budding or anagenesis
				new_dlen2<-dlen2-ch_zip		
				if(is.na(new_dlen1)){							#if budding or anagensis
					if(ch_zip==max_zip & dlen1+d1rng<dlen2){ 			#if anagensis
						anags<-c(anags,ktree$tip.label[dnode1])
						new_dlen1<-0
					}else{								#if budding
						budds<-c(budds,ktree$tip.label[dnode1])
						new_dlen1<-d1rng+dlen1-ch_zip
				}}else{new_dlen1<-ifelse(!is.na(d1rng),new_dlen1+d1rng,new_dlen1)}	#if tip, add term range
				new_dlen2<-ifelse(!is.na(d2rng),new_dlen2+d2rng,new_dlen2)		#if tip, add rng to new dlen
				#rescale branches according to their new lengths
				if(node!=(Ntip(ktree)+1)){	#if not root
					ktree$edge.length[ktree$edge[,2]==node]<-ch_zip-min_zip	#change stem length
					}
				ktree$edge.length[match(dnode1,ktree$edge[,2])]<-new_dlen1
				ktree$edge.length[match(dnode2,ktree$edge[,2])]<-new_dlen2
				#print(c(node,ch_zip-min_zip,new_dlen1,new_dlen2))
				nodes<-nodes[-1]			#update nodes
			}}
		ktree<-reorder(collapse.singles(ktree),"cladewise")
		ktree$anag.tips<-anags	#record the number of anagenetic ancestors
		ktree$budd.tips<-budds	#record the number of budding ancestors	
		ktree$nAdjZip<-nAdjZip
		#now add root.time: because NO TIPS ARE DROPPED (due to anagenesis) can calculate this now
			#must be calculated on LADs because the terminal ranges are added to the TREE!!!
			#should be time of earliest LAD + distance of root from earliest tip
		ktree$root.time<-max(timeData1[ktree$tip.label,2])+min(dist.nodes(ktree)[1:Ntip(ktree),Ntip(ktree)+1])	
		names(ktree$edge.length)<-NULL;names(ktree$tip.label)<-NULL;names(ktree$budd.tips)<-NULL;names(ktree$anag.tips)<-NULL
		#stuff for checking if things are correct
		tipdiffs<-cbind(diff(sort(-timeData1[,2])),diff(sort(dist.nodes(ktree)[1:Ntip(ktree),Ntip(ktree)+1]))
			,diff(sort(-timeData1[,2]))-diff(sort(dist.nodes(ktree)[1:Ntip(ktree),Ntip(ktree)+1])))	
		test1<-all(tipdiffs[,3]<(10^-10))
		test2<-identical(names(sort(-timeData1[,2])),ktree$tip.label[order(dist.nodes(ktree)[1:Ntip(ktree),Ntip(ktree)+1])])
		if(length(unique(timeData1[,2]))<Ntip(tree)){test2<-TRUE}	#test 2 does not work if any LADS are same
		if(all(c(test1,test2))){ktree$test<-"passed"}else{warning("Warning: Terminal tips improperly aligned, cause unknown. Use ouput with care.")}
		if(plot){
			parOrig<-par(mar=c(2.5,2.5,1,2.5));layout(matrix(1:3,3,));plot(ladderize(tree),show.tip.label=TRUE,use.edge.length=FALSE)
			plot(ladderize(ttree1),show.tip.label=TRUE);axisPhylo()			
			plot(ladderize(ktree),show.tip.label=TRUE);axisPhylo();
			layout(1);par(parOrig)		
			}
		names(ktree$edge.length)<-NULL
		ttrees[[ntr]]<-ktree
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}