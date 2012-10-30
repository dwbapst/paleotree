#srcTimeScaling_deprecated_10-29-12.R

srcTimePaleoPhy<-function(tree,timeData,sampRate,ntrees=1,anc.wt=1,node.mins=NULL,
	rand.obs=FALSE,FAD.only=FALSE,root.max=200,randres=FALSE,old.src=FALSE,plot=FALSE){
	#Samp Rate Conditioned Time Scaling via the stochastic ZIPPER method!
		#NOW with polytomy-resolving power via the PARALLEL ZIPPER method!
	#resolves relationships and time-scales trees using sampling-rate calibration
		#PRODUCES A STOCHASTIC SAMPLE, so run further analyses over a sample of MULTIPLE trees
		#initially time-scales a tree by making node time = earliest FAD of tip taxa (i.e. type='basic' in time_paleophy
		#will resolve polytomies and infer ancestors according to a budding model of speciation
	#NOTE: this function will return ALL tips put into it (it does not drop anagenetic tips)
		#for diversification analyses, use drop.zbl() to get rid of anagenetic tips (term zbls)
	#ARGUMENTS: 
		#'tree' is a phylogeny of taxa without branch lengths (should not be randomly resolved prior to time-scaling)
		#timeData is a matrix of FADs and 'observation datums' with rownames = species IDs
			#time is expected to be in standard paleo reference, such as MYA (i.e. 'larger' date is older)
		#observed taxon ranges (LAD-FAD) are ALWAYS added to the tree, with the second column LADs as the location of the tips
			#to allow for tips to be at range midpoints (recc. for trait evol analyses), replace LADs in timeData with mid-range dates
			#to place tips at FADs, simply set second column equal to the first column prior to running this function
		#sampling rate can be (a) a vector for all species with the instantaneous per-time-unit sampling rate (called 'r' by Foote,'97)
			#or (b) a single value of r for all species
		#ntrees is number of trees in output sample
		#if rand.obs=TRUE, the the function assumes that the LADs in timeData aren't where you actually want the tips
			#instead, tips will be randomly placed anywhere in that taxon's range with uniform probability
			#thus, tip locations will differ slightly for each tree in the sample
			#this will when you have a specimen or measurement but you don't know its placement in the species' range
		#Allows some or all mininum node times to be set pre-analysis
			#node.mins = vector of minimum time estimates for ind nodes, numbered as in edges, minus Ntip(ptree)
			#nodes with node.mins will be 'locked' so they can only be pushed back further than the given value
		#anc.wt will modify the likelihood weights of inferring direct ancestor-descendant relationships
			#anc.wt=0 will essentially shut this process off
			#forcing assumption that all taxa are terminal taxa resulting from bifurcating speciation
		#root.max is the maximum amount of time that the root could possibly be pushed to, from either:
				#(a) the first FAD or (b) the earliest node.mins
			#default is set to an arbitrary high number
		#plot will create a plot of the input tree, the basic timescaled tree 
			#and the new SRC time-scaled trees as it makes each of the later
	#additional information added to tree structure
		#$root.time
			#ALL TREES ARE OUTPUT WITH ELEMENTs "$root.time"
			#this is the time of the root on the tree, which is important for comparing across trees
			#this must be calculated prior to adding anything to terminal branches
		#$budd.tips - vector of tip.labels for all taxa placed as ancestors via budding
		#$anag.tips - vector of tip.labels for all taxa placed as ancestors via anagensis
	#example data
	#tree<-rtree(10);tree$edge.length<-sample(0:1,Nedge(tree),replace=TRUE);tree<-di2multi(tree)
	#ntrees=2;anc.wt=1;add.zombie=FALSE;node.mins=NULL;sampRate=rep(0.1,Ntip(tree));names(sampRate)<-tree$tip.label
	#timeData<-runif(Ntip(tree),200,400);timeData<-cbind(timeData,timeData-runif(Ntip(tree),1,80))
	#rownames(timeData)<-tree$tip.label;root.max=200;plot=TRUE;rand.obs=TRUE;FAD.only=FALSE;old.src=FALSE;randres=FALSE
	#node.mins<-c(-sort(-runif(1,600,900)),rep(NA,Nnode(tree)-1))	#assume two very deep divergences
	#
	#better example data?
	#set.seed(444);taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
	#timeData <- sampleRanges(taxa,r=0.5);tree <- taxa2cladogram(taxa,plot=TRUE);SRres<-getSampRateCont(rangesCont)
	#sampRate<-SRres[[2]][2];ntrees=1;anc.wt=1;node.mins=NULL;rand.obs=FALSE;FAD.only=FALSE;root.max=200
	#randres=FALSE;old.src=FALSE;plot=FALSE
	#
	require(ape)#;require(phangorn)
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	if(rand.obs & FAD.only){stop("Error: rand.obs and FAD.only cannot both be true")}
	#first clean out all taxa which are NA or missing in timeData
	if(ntrees==1){message("Warning: Do not interpret a single SRC time-scaled tree")}
	if(ntrees<1){stop("Error: ntrees<1")}
	tree<-drop.tip(tree,tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))])
	if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data")}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(length(sampRate)==1){sampRate<-rep(sampRate,Ntip(tree));names(sampRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(sampRate))))){
			stop("Sampling Rates Not Given For All Taxa!")}}
	#allow per-taxon anc.wt values
	if(length(anc.wt)==1){anc.wt<-rep(anc.wt,Ntip(tree));names(anc.wt)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(anc.wt))))){
			stop("Ancestral Weights Not Given For All Taxa on Tree!")}}
	if(length(node.mins)!=Nnode(tree) & !is.null(node.mins)){stop("node.mins length != Nnode!")}
	ttree1<-timePaleoPhy(tree,timeData,type="basic",node.mins=node.mins,add.term=FALSE)
	#identify which nodes are min-locked; make sure to update when resolving polytomies
	if(length(node.mins)>0){locked_nodes<-which(!is.na(node.mins))++Ntip(tree)}else{locked_nodes<-NA}
	ttree1<-collapse.singles(ttree1)
	ttrees<-rmtree(ntrees,3)
	for(ntr in 1:ntrees){
		if(rand.obs){
			timeData1<-cbind(timeData[,1],apply(timeData,1,function(x) runif(1,x[2],x[1])))
			}else{timeData1<-timeData}
		if(FAD.only){timeData1<-cbind(timeData[,1],timeData[,1])}else{timeData1<-timeData}
		ktree<-ttree1
		if(randres){ktree<-multi2di(ktree)}
		nodes<-(1:Nnode(ktree))+Ntip(ktree)		#get a vector of all internal nodes	
		nodes<-nodes[order(-node.depth(ktree)[-(1:Ntip(ktree))])]	#order by depth
		anags<-character();budds<-character()
		while(length(nodes)>0){		#can't use a for() because # of nodes may change
			#save_tree<-ktree;dev.new();plot(ktree)
			node<-nodes[1]
			tipl<-ktree$tip.label
			tipd<-cbind(ID=(1:Ntip(ttree1)),FAD=(timeData1[tipl,1]),LAD=(timeData1[tipl,2])
				,SR=(sampRate[tipl]),ancWt=(anc.wt[tipl]))
			if(node==(Ntip(ktree)+1)){
				min_zip<-(-root.max)	#if root, allow to be push back up to root.max
				stem_len<-root.max
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
				dSR<-drng<-danc.wt<-numeric()		
				for(i in dnodes){		#for each desc, get vector of SR for earliest and range if desc is a tip
					dtips<-match(unlist(Descendants(ktree,i)),tipd[,1])
					dearly<-which(tipd[dtips,2]==max(tipd[dtips,2]))[1]
					dSR[length(dSR)+1]<-tipd[dearly,4]
					danc.wt[length(danc.wt)+1]<-tipd[dearly,5]
					drng[length(drng)+1]<-ifelse(length(dtips)>1,NA,diff(unlist(tipd[dtips,3:2])))
					}
				if(old.src){	#old SRC has first lineage picked by probability of gap
					#first, randomly pick one desc lineage, weighted by implied unobs evol history of max zip
					zip_wt<-(dSR*exp(-dSR*dlen))/sum(dSR*exp(-dSR*dlen))	#get weights from density distribution, no call to anc.wt necc
					dnode1<-sample(dnodes,1,prob=zip_wt)
				}else{	#07-29-12: in new SRC method with 3-gap Y, choice of starting lineage doesn't matter
					dnode1<-dnodes[which(dlen==min(dlen))[1]]		#just pick first appearing
					}
				#make sure to include stem length in calculations!
				if(old.src & node==(Ntip(ktree)+1)){
					#07-29-12: this treats the stem length as a single gap
						#given that this should be a single gap and not gamma distributed, no old.src distinction
					#08-02-12: WRONG, why do this at all? just let root.max be stem.len!!!
					root_push<--seq(min_zip,0,by=0.1)
					root_prob<-dSR[dnode1==dnodes]*exp(-dSR[dnode1==dnodes]*root_push)
					stem_len<-sample(root_push,1,prob=root_prob)
					#in new SRC, the stem_len will remain the -min_zip
					}
				#make data structure for placed lineages; anc= row of anc lineage, events in time-from-stem 
				plin<-c(dnode1,(dlen[dnode1==dnodes]+stem_len),drng[dnode1==dnodes],NA,
					0,dlen[dnode1==dnodes]+stem_len,dlen[dnode1==dnodes]+drng[dnode1==dnodes]+stem_len,
					danc.wt[dnode1==dnodes])
				plin<-matrix(plin,1,)
				colnames(plin)<-c("node","brl","rng","anc","tSpec","tFO","tLO","ancWt")
				#place additional lineages, in order of max zip, using parallel zippers along placed lineages
				#add_nodes will hold necessary information on dlen, rng, SR for unplaced nodes
				add_nodes<-cbind(dnodes,dlen+stem_len,drng,dSR,danc.wt)[-which(dnodes==dnode1),]
				add_nodes<-add_nodes[order(dlen[-which(dnodes==dnode1)]),]
				colnames(add_nodes)<-c("node","dstem","rng","SR","ancWt")	#dstem is distance from stem to FAD
				for(i in 1:nrow(add_nodes)){
					#identify each currently placed lineage, produce zipper for each relative to stem point
					zips<-matrix(,1,3)
					for(j in 1:nrow(plin)){
						min_zip<-plin[j,5]	#time of speciation (stem time for each placed lineage)
						max_zip<-ifelse(plin[j,8]>0 & !is.na(plin[j,7]),	#max zip is complicated, dependent on anc.wt
							min(add_nodes[i,2],plin[j,7]),	#min of LAD of the placed lineage or FAD of the lineage to be placed
							min(add_nodes[i,2],plin[j,6]))	#min of LAD of the placed lineage or FAD of the lineage to be placed
						if(minlocked){max_zip<-stem_len}			#if minlocked, node must be not earlier than original node time
						#and on a stemTime=0 scale as used for the parallel zipper, stem_len is the original node time, unlike single zipper below
						poss_zip<-seq(min_zip,max_zip,by=0.1)
						if(old.src){
							gap<-add_nodes[i,2]-poss_zip	#inferred gap for lineage to be placed
							linDense<-ifelse(poss_zip>plin[j,6],
								plin[j,8]*add_nodes[i,4]*exp(-add_nodes[i,4]*gap),	#with anc.wt (of the PLACED LINEAGE)
									    add_nodes[i,4]*exp(-add_nodes[i,4]*gap))	#without anc.wt
						}else{
							gap1<-plin[j,6]-poss_zip		#waiting time from FAD1 to zip
							gap1<-ifelse(gap1>0,gap1,0)		#waiting time from branching point to FAD1
							gap2<-add_nodes[i,2]-poss_zip			#waiting time from branching point to FAD2
							gapStem2zip<-ifelse(plin[j,6]>poss_zip,poss_zip,plin[j,6])	#waiting time from stem to FAD1 OR br node
							totalgap<-gap1+gap2+gapStem2zip
							#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
								#under different combinations with p=q=r,p=q>r and p=q<r
								#admittedly not a perfect fit, though!
							linDense<-dgamma(totalgap,shape=2,rate=add_nodes[i,4])
							linDense<-ifelse(poss_zip>plin[j,6],
								plin[j,8]*linDense,	#with anc.wt (of the PLACED LINEAGE)
									    linDense)	#without anc.wt, but with gap1 probability
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
					ch_zip<-sample(1:nrow(zips),1,prob=zip_prob)	#sample zips
					ch_anc<-zips[ch_zip,2]
					ch_tzip<-zips[ch_zip,3]
					#if anagenesis, add to anags; if budding, add to budds
					if(!is.na(plin[ch_anc==plin[,1],7])){	#if the anc is terminal
						if(plin[ch_anc==plin[,1],7]==ch_tzip){anags<-c(anags,ktree$tip.label[ch_anc])}	#if anagenetic
						if(plin[ch_anc==plin[,1],6]<ch_tzip){budds<-c(budds,ktree$tip.label[ch_anc])}		#if budding
						}
					new_lin<-c(add_nodes[i,1],add_nodes[i,2]-ch_tzip,add_nodes[i,3],
						ch_anc,ch_tzip,add_nodes[i,2],add_nodes[i,2]+add_nodes[i,3],add_nodes[i,5])
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
				d1rng<-ifelse(dnode1<=Ntip(ktree),diff(unlist(tipd[match(dnode1,tipd[,1]),3:2])),NA)	#if clade, range = NA	
				d2rng<-ifelse(dnode2<=Ntip(ktree),diff(unlist(tipd[match(dnode2,tipd[,1]),3:2])),NA)
				d1ancWt<-ifelse(dnode1<=Ntip(ktree),unlist(tipd[match(dnode1,tipd[,1]),5]),0)				
				#ZIPPPER: first create a list of scenarios, treat the position of the node like a zipper which can be moved up or down
				max_zip<-ifelse(d1ancWt>0 & !is.na(d1rng),min(dlen1+d1rng,dlen2),dlen1)	#if d1 isn't a clade and there can be ancestors	
				if(minlocked){max_zip<-0}		#in single zipper, unlike parallel zipper, 0 is original node time
				poss_zip<-seq(min_zip,max_zip,by=0.1)
				if(old.src){
					#any position of the node posits two gaps: one for d edge 1 and for 2
					#edge1 doesn't matter; it's always same length, from the base to dnode1 so no effect on the probability
					gap2<-dlen2-poss_zip							#calculate gap 2
					prob_zip<-d2SR*exp(-d2SR*gap2)					#get likelihood weights
					prob_zip<-ifelse(poss_zip>dlen1,prob_zip*d1ancWt,prob_zip)
					linDensity<-prob_zip/sum(prob_zip)
				}else{
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
					linDense<-dgamma(totalgap,shape=2,rate=d2SR)
					linDensity<-ifelse((stem_len+dlen1)<(stem_len+poss_zip),
						d1ancWt*linDense,	#with anc.wt (of the first taxon)
						 	  linDense)	#without anc.wt, but with gap1 probability
					}
				#new as of 08-21-12
				linDensity[is.na(linDensity)]<-0
				if(sum(linDensity)==0){linDensity[1]<-1}
				ch_zip<-sample(poss_zip,1,prob=linDensity)				#pick zipper location
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
		ttrees[[ntr]]<-ktree
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}

bin_srcTimePaleoPhy<-function(tree,timeList,sampRate,ntrees=1,nonstoch.bin=FALSE,sites=NULL,anc.wt=1,node.mins=NULL,
	rand.obs=FALSE,FAD.only=FALSE,root.max=200,randres=FALSE,old.src=FALSE,plot=FALSE){
	#wrapper for applying SRC time-scaling to timeData where FADs and LADs are given as bins 
		#see SRC function for more details; SR MUST be instantaneous rate (if R, convert to r using functions in this library)
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
		#assumes time is in millions of years
	#sites is a matrix, used to indicate if binned FADs or LADs of multiple species were obtained from the locality / time point
			#i.e. the first appearance of species A, B and last appearance of C are all from the same lagerstatten
			#this will fix these to always have the same date relative to each other across many trees
			#this will assume that species listed for a site all are listed as being from the same interval...
				#this function also assumes that the sites matrix is ordered exactly as the timeList data is
	#if rand.obs=TRUE, the the function assumes that the LADs in timeList aren't where you actually want the tips
		#instead, tips will be randomly placed anywhere in that taxon's range with uniform probability
		#thus, tip locations will differ slightly for each tree in the sample
		#this is useful when you have a specimen or measurement but you don't know its placement in the species' range
	require(ape)
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	if(class(timeList[[1]])!="matrix"){if(class(timeList[[1]])=="data.frame"){timeList[[1]]<-as.matrix(timeList[[1]])
		}else{stop("Error: timeList[[1]] not of matrix or data.frame format")}}
	if(class(timeList[[2]])!="matrix"){if(class(timeList[[2]])=="data.frame"){timeList[[2]]<-as.matrix(timeList[[2]])
		}else{stop("Error: timeList[[2]] not of matrix or data.frame format")}}
	if(ntrees<1){stop("Error: ntrees<1")}
	if(ntrees==1){message("Warning: Do not interpret a single SRC time-scaled tree")}
	if(ntrees==1 & !nonstoch.bin){
		message("Warning: Do not interpret a single tree; dates are stochastically pulled from uniform distributions")}
	if(rand.obs & FAD.only){stop("Error: rand.obs and FAD.only cannot both be true")}
	#clean out all taxa which are NA or missing for timeData
	tree<-drop.tip(tree,tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeList[[2]][,1])))))])
	if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data")}
	timeList[[2]]<-timeList[[2]][!is.na(timeList[[2]][,1]),]
	if(any(is.na(timeList[[2]]))){stop("Weird NAs in Data??")}
	if(any(apply(timeList[[1]],1,diff)>0)){stop("Error: timeList[[1]] not in intervals in time relative to modern")}
	if(any(timeList[[1]][,2]<0)){stop("Error: Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeList[[2]],1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeList[[2]][,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	if(is.null(sites)){
		sites<-matrix(1:(Ntip(tree)*2),,2)
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
		if(rand.obs){timeData[,2]<-apply(timeData,1,function(x) runif(1,x[2],x[1]))}
		if(FAD.only){timeData[,2]<-timeData[,1]}
		tree2<-suppressMessages(srcTimePaleoPhy(tree,timeData,sampRate,ntrees=1,anc.wt=anc.wt,
			node.mins=node.mins,root.max=root.max,rand.obs=FALSE,randres=randres,old.src=old.src,plot=plot))
		tree2$ranges.used<-timeData
		ttrees[[ntrb]]<-tree2
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}
