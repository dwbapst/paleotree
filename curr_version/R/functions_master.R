#functions_master.R

timePaleoPhy<-function(tree,timeData,type="basic",vartime=NULL,ntrees=1,randres=FALSE,add.term=FALSE,
	rand.obs=FALSE,node.mins=NULL,plot=FALSE){
	#fast time calibration for phylogenies of fossil taxa; basic methods
		#this code inspired by similar code from G. Lloyd and G. Hunt
	#INITIAL: 
		#time-scales a tree by making node time = earliest FAD of tip taxa
		#tree is a phylogeny of taxa without branch lengths
		#timeData is a matrix of FADs and LADs with rownames = species IDs
			#time is expected to be in standard paleo reference, such as MYA (i.e. 'larger' date is older)
		#vartime is a time variable used for time-scaling methods that are not "basic", ignored if "basic"
		#Allows some or all node times to be set pre-analysis
			#node.mins = vector of minimum time estimates for ind nodes, numbered as in edges, minus Ntip(ptree)
		#will make multiple randomly resolved trees if ntrees>1 and randres=T; polytomies resolved with multi2di() from ape
			#not any reason to do this unless you have polytomies
			#do !not! !ever! trust a single tree like that!! ever!!
	#TYPES
		#if (type="basic") just gives initial raw time-scaled tree (vartime is ignored), many zero-length branches
		#if (type="aba") then adds vartime to all branches (All Branch Additive)
		#if (type="zlba") then adds vartime to zero length branches (Zero Length Branch Additive) 
		#if (type="mbl") scales up all branches greater than vartime and subtracts from lower (Min Branch Length)
		#if (type="equal") "equal" method of G. Lloyd, recreated here; vartime is used as time added to root
	#ADDING TERMINAL BRANCHES TO PHYLOGENY 
		#if addterm!=FALSE, then observed taxon ranges (LAD-FAD) are added to the tree, with LADs as the location of the tips
		#to allow for tips to be at range midpoints (recc. for trait evol analyses), replace LADs in timeData with mid-range dates
	#root.time
		#ALL TREES ARE OUTPUT WITH ELEMENTs "$root.time"
		#this is the time of the root on the tree, which is important for comparing across trees
		#this must be calculated prior to adding anything to terminal branches
	#tree<-rtree(10);tree$edge.length<-NULL;type="basic";vartime=NULL;add.term="none";node.mins=NULL
	#timeData<-runif(10,30,200);timeData<-cbind(timeData,timeData-runif(10,1,20));rownames(timeData)<-tree$tip.label
	#node.mins<-runif(9,50,300)
	require(ape)	
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	#remove taxa that are NA or missing in timeData
	if(ntrees<1){stop("Error: ntrees<1")}
	if(!add.term & rand.obs){stop("Error: Inconsistent arguments: add.term must be true for rand.obs to have any effect on output!")}
	if(ntrees>1 & !randres & !rand.obs){stop("Error: Time-scale more trees without randomly resolving or random obs?!")}
	if(ntrees==1 & randres){message("Warning: Do not interpret a single randomly-resolved tree")}
	if(ntrees==1 & rand.obs){message("Warning: Do not interpret a single tree with randomly-placed obs times")}
	tree<-drop.tip(tree,tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))])
	if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data")}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	ttrees<-rmtree(ntrees,2)
	savetree<-tree			#save tree now so can replace with each loop for multi2di()
	saveTD<-timeData
	for(ntr in 1:ntrees){
		#resolve nodes, if tree is not binary
		if(!is.binary.tree(savetree) & randres){tree<-multi2di(savetree)}else{tree<-savetree}
		if(rand.obs){timeData[,2]<-apply(timeData,1,function(x) runif(1,x[2],x[1]))}else{timeData<-saveTD}
		ntime<-sapply(1:Nnode(tree),function(x) 
			max(timeData[tree$tip.label[unlist(prop.part(tree)[x])],1]))	#first, get node times
		ntime<-c(timeData[tree$tip.label,1],ntime)
		if(length(node.mins)>0){	#if there are node.mins, alter ntime as necessary
			#of course, node.mins is referring to nodes in unresolved savetree
			#need to figure out which nodes are which now if randres; remake node.mins
			if(!is.binary.tree(savetree) & randres){
				node_changes<-match(prop.part(savetree),prop.part(tree))
				node.mins1<-rep(NA,Nnode(tree))
				node.mins1[node_changes]<-node.mins
			}else{
				node.mins1<-node.mins
				}
			#require(phangorn)
			for(i in (Ntip(tree)+1):length(ntime)){	#all internal nodes
				desc_all<-unlist(Descendants(tree,i,type="all"))
				desc_nodes<-c(desc_all[desc_all>Ntip(tree)],i)-Ntip(tree)	#INCLUDING ITSELF			
				node_times<-node.mins1[desc_nodes]
				ntime[i]<-max(ntime[i],node_times[!is.na(node_times)])
				}
			}
		if(type=="equal" & !is.null(vartime)){				#add to root, if method="equal"
			ntime[Ntip(tree)+1]<-vartime+ntime[Ntip(tree)+1]
			#anchor_adjust<-vartime+anchor_adjust
			}	
		ttree<-tree
		ttree$edge.length<-sapply(1:Nedge(ttree),function(x) 
			ntime[ttree$edge[x,1]]-ntime[ttree$edge[x,2]])	#finds each edge length easy peasy, based on G. Lloyd's code
		#ttree_basic<-ttree
		##if type=basic, I don't have to do anything but set root.time
		if(type=="aba"){	#if (type="aba") then adds vartime to all branches (All Branch Additive) 
			if(is.na(vartime)){stop("No All Branch Additive Value Given!")}
			ttree$edge.length<-ttree$edge.length+vartime
			}
		if(type=="zlba"){	#if (type="zlba") then adds vartime to zero length branches (Zero Length Branch Additive) 
			if(is.na(vartime)){stop("No Branch Additive Value Given!")}
			ttree$edge.length[ttree$edge.length<0.0001]<-ttree$edge.length[ttree$edge.length<0.0001]+vartime
			}
		if(type=="mbl"){
			#if (type="mbl") scales up all branches greater than vartime and subtracts from lower
				#as long as there are branches smaller than vartime
			#require(phangorn)
			if(is.na(vartime)){stop("No Minimum Branch Length Value Given!")}
			root_node<-Ntip(ttree)+1;mbl<-vartime
			while(any(ttree$edge.length<mbl)){
				#pick one at random, make vector of every mom node that is ancestral
				mom<-ttree$edge[((1:Nedge(ttree))[ttree$edge.length<mbl])
					[sample(length((1:Nedge(ttree))[ttree$edge.length<mbl]),1)],1]
				mom<-c(mom,Ancestors(ttree,mom))
				debt<-mbl-min(ttree$edge.length[ttree$edge[,1]==mom[1]])
				ttree$edge.length[mom[1]==ttree$edge[,1]]<-ttree$edge.length[mom[1]==ttree$edge[,1]] + debt[1]
				#make vector of smallest brlen with each mom node as anc
				#calculate, simulatenously, the changes in debt and branch lengthening required as go down tree
				#change branch lengths; hypothetically, debt should then equal zero...
				if(length(mom)>1){for(i in 2:length(mom)){
					small<-min(ttree$edge.length[ttree$edge[,1]==mom[i]])
					mom_blen<-ttree$edge.length[ttree$edge[,1]==mom[i] & ttree$edge[,2]==mom[i-1]]
					debt[i]<-max(debt[i-1] - max(mom_blen-mbl,0),0) + max(mbl-small,0) 
					ttree$edge.length[ttree$edge[,1]==mom[i] & ttree$edge[,2]==mom[i-1]] <- 
						mom_blen - max(min(max(mom_blen-mbl,0),debt[i-1]),0) + max(mbl-small,0)
					ttree$edge.length[ttree$edge[,1]==mom[i] & ttree$edge[,2]!=mom[i-1]] <-  
						ttree$edge.length[ttree$edge[,1]==mom[i] & ttree$edge[,2]!=mom[i-1]] + debt[i]
					}}
				}
			}
		if(type=="equal"){	#G. Lloyd's "equal" method
			#get a depth-ordered vector that identifies zero-length branches
			zbr<-cbind(1:Nedge(ttree),node.depth(ttree)[ttree$edge[,2]]) 	#Get branch list; 1st col = end-node, 2nd = depth
			zbr<-zbr[ttree$edge.length==0,]						#Parses zbr to just zero-length branches
			zbr<-zbr[order(zbr[,2]),1]							#order zbr by depth
			for(i in zbr){if (ttree$edge.length[i] == 0) {			#starting with most shallow zlb, is this branch a zlb?
				#if zlb, make a vector of mom-zlbs, going down the tree
				brs<-ttree$edge[i,2] 						#branches to rescale, starting with picked branch
				mom<-which(ttree$edge[i,1]==ttree$edge[,2])
				while(ttree$edge[mom,1]!=(Ntip(ttree)+1) & ttree$edge.length[mom]==0){ #keep going while preceding edge is zero len and isn't the root
					brs[length(brs)+1]<-ttree$edge[mom,2]  		#keep adding these branches to brs
					mom<-which(ttree$edge[mom,1]==ttree$edge[,2])	#reset mom
					}
				brs[length(brs)+1]<-ttree$edge[mom,2] 				#Add final branch (which isn't zlb)
				totbl<-sum(ttree$edge.length[match(brs,ttree$edge[,2])]) 	#Amount of time to be shared
				ntime[brs[-1]]<-ntime[brs[-1]]+cumsum(rep(totbl/length(brs),length(brs)-1))
				ttree$edge.length<-sapply(1:Nedge(ttree),function(x) 
					ntime[ttree$edge[x,1]]-ntime[ttree$edge[x,2]])	#update branch lengths using ntime
				}}
			}
		#now add root.time: should be time of earliest FAD + distance of root from earliest tip
		ttree$root.time<-max(timeData[ttree$tip.label,1])+min(dist.nodes(ttree)[1:Ntip(ttree),Ntip(ttree)+1])	
		#if add.term!=FALSE, then taxon observed ranges are added to the tree, with the LADs becoming the location of the tips
		if(add.term){
			obs_ranges<-timeData[,1]-timeData[,2]
			term_id<-ttree$tip.label[ttree$edge[ttree$edge[,2]<=Ntip(ttree),2]]
			term_add<-sapply(term_id,function(x) obs_ranges[x])
			ttree$edge.length[ttree$edge[,2]<=Ntip(ttree)]<-ttree$edge.length[ttree$edge[,2]<=Ntip(ttree)]+term_add
			}
		if(plot){
			parOrig<-par(mar=c(2.5,1,1,0.5));layout(1:2)
			plot(ladderize(tree),show.tip.label=TRUE,use.edge.length=FALSE)
			plot(ladderize(ttree),show.tip.label=TRUE);axisPhylo()
			layout(1);par(parOrig)
			}
		ttrees[[ntr]]<-ttree
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}

bin_timePaleoPhy<-function(tree,timeList,type="basic",vartime=NULL,ntrees=1,nonstoch.bin=FALSE,randres=FALSE,sites=NULL,
	add.term=FALSE,rand.obs=FALSE,node.mins=NULL,plot=FALSE){
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
	require(ape)
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
		if(!is.binary.tree(tree) & randres){tree1<-multi2di(tree)}else{tree1<-tree}
		ttrees[[ntrb]]<-suppressMessages(timePaleoPhy(tree1,timeData,type=type,vartime=vartime,ntrees=1,
			randres=FALSE,add.term=add.term,rand.obs=FALSE,node.mins=node.mins,plot=plot))
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}

srcTimePaleoPhy<-function(tree,timeData,sampRate,ntrees=1,anc.wt=1,node.mins=NULL,rand.obs=FALSE,FAD.only=FALSE,root.max=200,plot=FALSE){
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
	#rownames(timeData)<-tree$tip.label;root.max=200;plot=TRUE;rand.obs=TRUE;FAD.only=FALSE
	#node.mins<-c(-sort(-runif(1,600,900)),rep(NA,Nnode(tree)-1))	#assume two very deep divergences
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
		}else{if(length(sampRate)!=Ntip(tree)){stop("SR Length != Ntip!")}}
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
		nodes<-(1:Nnode(ktree))+Ntip(ktree)		#get a vector of all internal nodes	
		nodes<-nodes[order(-node.depth(ktree)[-(1:Ntip(ktree))])]	#order by depth
		anags<-character();budds<-character()
		while(length(nodes)>0){		#can't use a for() because # of nodes may change
			#save_tree<-ktree;dev.new();plot(ktree)
			node<-nodes[1]
			tipl<-ktree$tip.label
			tipd<-cbind(ID=(1:Ntip(ttree1)),FAD=(timeData1[tipl,1]),LAD=(timeData1[tipl,2]),SR=(sampRate[tipl]))
			if(node==(Ntip(ktree)+1)){
				min_zip<-(-root.max)	#if root, allow to be push back up to root.max
				root_push<--seq(min_zip,0,by=0.1)
			}else{									#if not root, push down to lower node
				min_zip<-(-ktree$edge.length[ktree$edge[,2]==node])
				stem_len<-ktree$edge.length[ktree$edge[,2]==node]
				}		
			dnodes<-ktree$edge[ktree$edge[,1]==node,2]	#find the daughter nodes
			dlen<-ktree$edge.length[match(dnodes,ktree$edge[,2])]	#find the dedges lengths
			minlocked<-ifelse(!all(is.na(locked_nodes)),any(node==locked_nodes),FALSE)#is this node min-locked?
			#if(any(dnodes==which(ktree$tip.label=="t10"))){break()}
			if(length(dnodes)>2){		#if node is a polytomy, use PARALLEL ZIPPER
				#first, randomly pick one desc lineage, weighted by implied unobs evol history of max zip
				dSR<-numeric();drng<-numeric()		
				for(i in dnodes){		#for each desc, get vector of SR for earliest and range if desc is a tip
					dtips<-match(unlist(Descendants(ktree,i)),tipd[,1])
					dearly<-which(tipd[dtips,2]==max(tipd[dtips,2]))[1]
					dSR[length(dSR)+1]<-tipd[dearly,4]
					drng[length(drng)+1]<-ifelse(length(dtips)>1,NA,diff(unlist(tipd[dtips,3:2])))
					}
				zip_wt<-(dSR*exp(-dSR*dlen))/sum(dSR*exp(-dSR*dlen))	#get likelihood weights, no call to anc.wt necc
				dnode1<-sample(dnodes,1,prob=zip_wt)
				#make sure to include stem length in calculations!
				if(node==(Ntip(ktree)+1)){
					root_prob<-dSR[dnode1==dnodes]*exp(-dSR[dnode1==dnodes]*root_push)
					stem_len<-sample(root_push,1,prob=root_prob)
					}
				#make data structure for placed lineages; anc= row of anc lineage, events in time-from-stem 
				plin<-c(dnode1,(dlen[dnode1==dnodes]+stem_len),drng[dnode1==dnodes],NA,
					0,dlen[dnode1==dnodes]+stem_len,dlen[dnode1==dnodes]+drng[dnode1==dnodes]+stem_len)
				plin<-matrix(plin,1,)
				colnames(plin)<-c("node","brl","rng","anc","tSpec","tFO","tLO")
				#place additional lineages, in order of max zip, using parallel zippers along placed lineages
				#add_nodes will hold necessary information on dlen, rng, SR for unplaced nodes
				add_nodes<-cbind(dnodes,dlen+stem_len,drng,dSR)[-which(dnodes==dnode1),]
				add_nodes<-add_nodes[order(dlen[-which(dnodes==dnode1)]),]
				colnames(add_nodes)<-c("node","dstem","rng","SR")	#dstem is distance from stem to FAD
				for(i in 1:nrow(add_nodes)){
					#identify each currently placed lineage, produce zipper for each relative to stem point
					zips<-matrix(,1,3)
					for(j in 1:nrow(plin)){
						min_zip<-plin[j,5]
						max_zip<-ifelse(anc.wt>0 & !is.na(plin[j,7]),
							min(add_nodes[i,2],plin[j,7]),min(add_nodes[i,2],plin[j,6]))
						if(minlocked){max_zip<-0}
						poss_zip<-seq(min_zip,max_zip,by=0.1)
						gap<-add_nodes[i,2]-poss_zip	#inferred gap for lineage to be placed
						like<-ifelse(poss_zip>plin[j,6],
							anc.wt*add_nodes[i,4]*exp(-add_nodes[i,4]*gap),add_nodes[i,4]*exp(-add_nodes[i,4]*gap))
						new_zip<-cbind(like,plin[j,1],poss_zip)
						zips<-rbind(zips,new_zip)
						}
					if(nrow(zips)<3){zips<-matrix(zips[-1,],1,3)}else{zips<-zips[-1,]}
					colnames(zips)<-c("like","anc","tzip")
					zip_prob<-zips[,1]/sum(zips[,1])
					ch_zip<-sample(1:nrow(zips),1,prob=zip_prob)	#sample zips
					ch_anc<-zips[ch_zip,2]
					ch_tzip<-zips[ch_zip,3]
					#if anagenesis, add to anags; if budding, add to budds
					if(!is.na(plin[ch_anc==plin[,1],7])){	#if the anc is terminal
						if(plin[ch_anc==plin[,1],7]==ch_tzip){anags<-c(anags,ktree$tip.label[ch_anc])}	#if anagenetic
						if(plin[ch_anc==plin[,1],6]<ch_tzip){budds<-c(budds,ktree$tip.label[ch_anc])}		#if budding
						}
					new_lin<-c(add_nodes[i,1],add_nodes[i,2]-ch_tzip,add_nodes[i,3],
						ch_anc,ch_tzip,add_nodes[i,2],add_nodes[i,2]+add_nodes[i,3])
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
					drtips<-unlist(Descendants(ktree,node))	#ape needs better tree editting functions
					tip_lab<-ktree$tip.label[drtips[1]]	#I need to cut out all but one tip, for the sole purpose of putting it all back later)
					droptree<-collapse.singles(drop.tip(ktree,drtips[-1]))
					droptree$edge.length[droptree$edge[,2]==which(droptree$tip.label==tip_lab)]<-new_stem	#reset edge length leading to remaining tip to new_stem
					droptree<-bind.tree(droptree,subtree,where=which(droptree$tip.label==tip_lab))	#put in subtree at tip
					ktree1<-droptree
				}else{				#if it is the node
					ktree1<-subtree
					}				
				#once you've changed the structure of the tree find original nodes in new tree (surprisingly frustating to code)
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
			}else{	#if node is not a polytomy, then use regular zipper
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
				#ZIPPPER: first create a list of scenarios, treat the position of the node like a zipper which can be moved up or down
				max_zip<-ifelse(anc.wt>0 & !is.na(d1rng),min(dlen1+d1rng,dlen2),dlen1)	#if d1 isn't a clade and there can be ancestors	
				if(minlocked){max_zip<-0}
				poss_zip<-seq(min_zip,max_zip,by=0.1)
				#any position of the node posits two gaps: one for d edge 1 and for 2
				#edge1 doesn't matter; it's always same length, from the base to dnode1 so no effect on the probability
				gap2<-dlen2-poss_zip							#calculate gap 2
				prob_zip<-d2SR*exp(-d2SR*gap2)					#get likelihood weights
				prob_zip<-ifelse(poss_zip>dlen1,prob_zip*anc.wt,prob_zip)
				prob_zip<-prob_zip/sum(prob_zip)
				ch_zip<-sample(poss_zip,1,prob=prob_zip)				#pick zipper location
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
	rand.obs=FALSE,FAD.only=FALSE,root.max=200,plot=FALSE){
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
		ttrees[[ntrb]]<-suppressMessages(srcTimePaleoPhy(tree,timeData,sampRate,ntrees=1,
			anc.wt=anc.wt,node.mins=node.mins,root.max=root.max,rand.obs=FALSE,plot=plot))
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}

expandTaxonTree<-function(taxonTree,taxaData,collapse=NULL,plot=FALSE){
	#this function takes a higher-level taxon tree and
		#expands it to a lower level species-level tree
		#using a species list
	#"taxa" here represents the groups to be replaced on the taxonTree
	#taxonTree = tree with taxon IDs as tips
	#taxaData = character vector of higher taxon ids for each new tip, tip labels as vector names
	#collapse = if present, vector of taxa names to be collapsed
	#should be possible to take a tree of mixed species/genera
		#and just replace the genera
	#taxonTree<-rtree(10);taxonTree$tip.label<-as.character(1:10);collapse<-sample(taxonTree$tip.label,5)
	#taxaData<-as.character(sample(1:10,100,replace=TRUE));names(taxaData)<-paste("t",1:100,sep="")
	require(ape)
	if(class(taxonTree)!="phylo"){stop("Error: taxonTree is not of class phylo")}
	tree<-taxonTree;tree$edge.length<-rep(1,Nedge(tree))		#get rid of all branch lengths
	#first, expand all higher taxa to lower taxon polytomies
	for(i in unique(taxaData)){				#loop through all 	
		tip<-which(tree$tip.label==i)
		if(length(collapse)>0){if(any(collapse==i)){
			tree$edge.length[which(tree$edge[,2]==tip)]<-0
			}}
		cotaxa<-names(taxaData)[taxaData==i]	#which species do I want? These...
		repTree<-stree(length(cotaxa))		#replacement polytomy
		repTree$edge.length<-rep(1,length(cotaxa))
		repTree$tip.label<-cotaxa			#replace names,edge.lengths
		tree<-bind.tree(tree,repTree,tip)	#replace the right tip	
		}
		#now collapse non-monophyletic groupings
	tree1<-di2multi(tree);tree1$edge.length<-NULL;tree1<-collapse.singles(tree1)
	tree1<-read.tree(text=write.tree(tree1))
	if(plot==TRUE){layout(1:2);plot(taxonTree);plot(tree1);layout(1)}
	return(tree1)
	}

depthRainbow<-function(tree){
	#plots a tree with edges color-coded to depth
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	tree<-ladderize(tree)
	ndepth<-dist.nodes(tree)[,Ntip(tree)+1]
	#nodelabels(ceiling(ndepth[(Ntip(tree):Nedge(tree))+1]),node=(Ntip(tree):Nedge(tree))+1)
	edepth<-ceiling((ndepth[(Ntip(tree):Nedge(tree))+1])[tree$edge[,1]-Ntip(tree)])+1
	col_edge<-rainbow(max(edepth))[edepth]
	plot(ladderize(tree),show.tip.label=FALSE,edge.color=col_edge);axisPhylo()
	}

degradeTree<-function(tree,prop_collapse,node.depth=NA){
	#collapses a given proportion of internal edges, creating polytomies
		#node.depth conditions on depth of edge in tree
			# 1 removes more shallow nodes, 0 removes deeper nodes
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	tree$edge.length<-NULL
	tree<-collapse.singles(tree)
	edge<-(1:length(tree$edge))[which(tree$edge[,2]>Ntip(tree))]	#internal edges
	if(is.na(node.depth)){
		cedge<-sample(edge,round(prop_collapse*length(edge)))	#edges chosen to collapse
	}else{
		node_pdesc<-sapply(prop.part(tree),length)/Ntip(tree)	#prop desc per int node
		edge_pdesc<-node_pdesc[tree$edge[edge,2]-Ntip(tree)]
		edge_prob<-(edge_pdesc-node.depth)^2;edge_prob<-edge_prob/sum(edge_prob)
		cedge<-sample(edge,round(prop_collapse*length(edge)),prob=edge_prob)	#chosen edges	
		}
	tree$edge.length<-rep(1,Nedge(tree))
	tree$edge.length[cedge]<-0
	tree<-di2multi(tree)
	tree<-collapse.singles(tree)
	tree$edge.length<-NULL
	return(tree)
	}

unitLengthTree<-function(tree){
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	tree$edge.length<-rep(1,Nedge(tree))
	tree$root.time<-NULL
	return(tree)
	}

qsProb2Comp<-function(R,q){
	#calculate completeness given R and mu
	res<-numeric()
	for(t in 1:10000){
		res[t]<-(1-(1-R)^t)*(exp(-q*(t-1))-exp(-q*t))
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}

qsRate2Comp<-function(r,q){
	#calculate completeness given r and mu
	res<-r/(r+q)
	names(res)<-NULL
	return(res)
	}

pqsRate2sProb<-function(r,p,q,int.length=1){
	#A more accurate estimat of R given r, p and q
	#assuming p,q,r are constant and the timespan is infinte
		#dt is interval length for R
	#USES equations 26-29 from appendix to Foote (2000)
	#prob of samp for lineages that cross both boundaries
		#note typo in Foote (2000), eq 26, corrected version below
	dt<-int.length
	PDbt<-function(r,dt){1-exp(-r*dt)}
	#prob of samp for lineages that only cross bottom boundary
	PDbL<-function(q,r,dt){
		(((r+(q*exp(-(q+r)*dt)))/(q+r))-exp(-q*dt))/(1-exp(-q*dt))
		}
	#prob of samp for lineages that only cross upper boundary
	PDFt<-function(p,r,dt){
		(((r+(p*exp(-(p+r)*dt)))/(p+r))-exp(-p*dt))/(1-exp(-p*dt))
		}
	#prob of samp for lineages that cross neither boundary
		#29b corrected with addition sign!
	PDFL<-function(p,q,r,dt){
		if(p==q){
			NbNFL<-1/(exp(-q*dt)+(p*dt)-1)		#N(b)/N(FL) based on eq 1b and 6b
			term1<-(r*dt)/(p+r)				#first term in square brackets in eq 29b
			term2<-(1-exp(-p*dt))/p				#second term
			term3<-(p*(1-exp(-(p+r)*dt)))/((p+r)^2)	#third term
			terms<-term1-term2+term3			#full terms in square brackets
			res<-(NbNFL)*p*terms					#P(D|FL)
		}else{
			NbNFL<-1/(((q*exp((p-q)*dt))+((p-q)*exp(-q*dt))-p)/(p-q))
			term1<-(p*r*(exp((p-q)*dt)-1))/((q+r)*(p-q))
			term2<-(p*q*exp(-(q+r)*dt)*(exp((p+r)*dt)-1))/((p+r)*(q+r))
			term3<-exp(-q*dt)*(exp(p*dt)-1)
			terms<-term1+term2-term3
			res<-(NbNFL)*terms
			}
		res
		}
	#need to weight the PDs by the P of those taxon classes
		#use N equations from Foote (2000), relative to Nb to be probs
	Pbt<-exp(-q*dt)	
	PbL<-(1-exp(-q*dt))
	PFt<-exp((p-q)*dt)*(1-exp(-p*dt))
	if(p==q){PFL<-exp(-q*dt)+(p*dt)-1
		}else{PFL<-((q*exp((p-q)*dt))+((p-q)*exp(-q*dt))-p)/(p-q)}
	res<-sum(PDbt(r,dt)*Pbt,PDbL(q,r,dt)*PbL,
		PDFt(p,r,dt)*PFt,PDFL(p,q,r,dt)*PFL)
	names(res)<-NULL
	return(res)
	}

sProb2sRate<-function(R,int.length=1){
	res<-(-log(1-R)/int.length)	#rough estimate
	names(res)<-NULL
	return(res)
	}

sRate2sProb<-function(r,int.length=1){
	res<-1-exp(-r*int.length)	#rough estimate
	names(res)<-NULL
	return(res)
	}

probAnc<-function(p,q,R){	
	#calculates prob of taxa with indirect desc under budding speciation
		#under infinite time, with p=q or p<q
	Pd<-function(q,Ti){exp(-q*(Ti-1))-exp(-q*Ti)}
	PN<-function(p,Ti,N){(exp(-p*Ti)*(p*Ti)^N)/factorial(N)}
	Qm<-function(p,q,M){
		x<-(4*p*q)/((p+q)^2)
		((p+q)/(2*p))*(factorial(2*M)/((2^(2*M))*factorial(M)^2))*((x^M)/((2*M)-1))
		}
	Pinf<-function(p,q,Pp){
		res<-numeric()
		for(M in 1:80){
		res[M]<-Qm(p,q,M)*(1-(1-Pp)^M)
			}
		sum(res)
		}
	Pp<-qsProb2Comp(R,q)
	res<-numeric()
	for(t in 1:2000){
		Nres<-numeric()
		for(N in 1:100){
			Nres[N]<-PN(p,t,N)*(1-(Pp)^N)
			}
		res[t]<-(((1-(1-R)^t)*Pd(q,t))/Pp)*sum(Nres)
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}

getSampRateCont<-function(timeData,n_tbins=1,grp1=NA,grp2=NA,threshold=0.1,est_only=FALSE,iterations=1000000,initial=0.5){
	#this is the multi-parameter maximum likelihood analysis of continuous-time fossil ranges
		#uses a set of timeData (FADs and LADs) to fit models of different samp rates and ext rates
		#can allow for free-moving time windows and different groups
			#these models can then be compared with AIC
		#if est_only=TRUE, then the q and r estimates will be given back per-species
	#x<-runif(100);x<-cbind(x+rexp(100),x);getSampRateCont(x)
	#getSampRateCont(x,grp1=(sample(1:2,100,replace=TRUE)))
	#REQUIRED FUNCTIONS BELOW
	######################
	make_ft<-function(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s){
	qr_predict<-qr_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
	function(par){
		qr<-qr_predict(par)
		q<-qr[,1]; r<-qr[,2]
		ft<-ifelse(dur==0,q/(r+q),q*r*exp(-q*dur)/(r+q))
		-sum(log(ft))
		}
	}
	####################
	qr_predict_multpar<-function(FO,LO,n_tbins,grp1,grp2,g1s,g2s){
	#par consists of a vector, with first n_tbins-1 elements represent bin LENGTHS
		#unless n_tbins==1
		#after that, its a matrix with k rows, first col is q, second is r
			#first comes rows associated with each time bin, 
			#than each group in grp1, than grp2
	#THIS FUNCTION SHOULD OUTPUT A Nx2 MATRIX OF q and r VALUES to the support function
		#MWAHAHAH!
	if(n_tbins>1){	#IF THERE ARE TIME BINS
		if(g1s>0){	#IF THERE ARE GROUPS + TIMEBINS
			if(g2s>0){	#IF THERE IS TWO GROUPS +TIMEBINS
				function(par){
					n_tb<-n_tbins-1
					t_raw<-c(0.5,par[1:n_tb])
					t_prop<-t_raw/sum(t_raw)
					t_bl<-t_prop*(max(FO)-min(LO))
					tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
					mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
					mqrt<-matrix(mqr[1:n_tbins,],,2,byrow=TRUE)
					mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]
					mqrg2<-mqr[(n_tbins+g1s+1):(n_tbins+g1s+g2s),]
					tcat<-sapply(FO,function(x) sum(tbin>=x))
					qrt<-sapply(tcat,function(x) mqrt[x,])
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					qrg2<-sapply(grp2,function(x) mqrg2[x,])
					t((qrt+qrg1+qrg2)/3)
				}
			}else{	#IF THERE IS ONE GROUP +TIMEBINS
				function(par){
					n_tb<-n_tbins-1
					t_raw<-c(0.5,par[1:n_tb])
					t_prop<-t_raw/sum(t_raw)
					t_bl<-t_prop*(max(FO)-min(LO))
					tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
					mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
					mqrt<-matrix(mqr[1:n_tbins,],,2,byrow=TRUE)
					mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]
					tcat<-sapply(FO,function(x) sum(tbin>=x))
					qrt<-sapply(tcat,function(x) mqrt[x,])
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					t((qrt+qrg1)/2)
				}
			}
		}else{	#IF THERE ARE NO GROUPS AND TIMEBINS
			function(par){
				n_tb<-n_tbins-1
				t_raw<-c(0.5,par[1:n_tb])
				t_prop<-t_raw/sum(t_raw)
				t_bl<-t_prop*(max(FO)-min(LO))
				tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
				mqrt<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
				tcat<-sapply(FO,function(x) sum(tbin>=x))
				tcount<-sapply(sort(unique(tcat)),function(x) sum(tcat==x))
				#if(all(tcount>100)){t(sapply(tcat,function(x) mqrt[x,]))
				#}else{matrix(0.99,length(FO),2,byrow=TRUE)}
				t(sapply(tcat,function(x) mqrt[x,]))
			}
		}
	}else{	#IF THERE ARE NO TIMEBINS
		if(g1s>0){	#IF THERE ARE GROUPS AND NO TIMEBINS
			if(g2s>0){	#IF THERE IS TWO GROUPS AND NO TIMEBINS
				function(par){
					mqr<-matrix(par,,2,byrow=TRUE)
					mqrg1<-mqr[1:g1s,]
					mqrg2<-mqr[(g1s+1):(g1s+g2s),]
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					qrg2<-sapply(grp2,function(x) mqrg2[x,])
					t((qrg1+qrg2)/2)
				}
			}else{	#IF THERE IS ONE GROUP AND NO TIMEBINS
				function(par){
					mqr<-matrix(par,,2,byrow=TRUE)
					t(sapply(grp1,function(x) mqr[x,]))
				}
			}
		}else{	#IF THERE ARE NO GROUPS AND NO TIMEBINS (2 param model)
			function(par){
				matrix(par,length(FO),2,byrow=TRUE)
				}	
			}
		}			
	}
	#########################
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	#get rid of any NAs
	if(length(grp1)>1){
		if(length(grp1)!=nrow(timeData)){stop("Error: grp1 is not same length as timeData")}
		grp1<-grp1[!is.na(timeData[,1])]
		}
	if(length(grp2)>1){
		if(length(grp2)!=nrow(timeData)){stop("Error: grp2 is not same length as timeData")}		
		grp2<-grp2[!is.na(timeData[,1])]
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	FO<-timeData[,1];LO<-timeData[,2]
	dur<-(-LO)-(-FO)
	#THRESHOLD DETERMINES RANGES TOO SMALL TO BE CONSIDERED NOT ONE-TIMERS
	dur[dur<threshold]<-0
	#NOW THAT DUR IS GOOD, CONSTRUCT PROPER SUPPORT FUNCTION
	#NOTE THAT THE FIRST BIN MUST ALWAYS BE THE FIRST INTERVAL
	if(length(grp1)>1){g1s<-length(unique(grp1))}else{g1s<-0}
	if(length(grp2)>1){g2s<-length(unique(grp2))}else{g2s<-0}
	support_ft<-make_ft(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s)
	npar<-(n_tbins-1)+(ifelse(n_tbins>1,2*n_tbins,0)+ifelse(g1s>0,g1s*2,0)
		+ifelse(g2s>0,g2s*2,0))
	if(npar==0){npar<-2}
	par_lim<-c(0.0001,10)
	par_init<-rep(initial,npar)
	par_min<-rep(par_lim[1],npar)
	par_max<-rep(par_lim[2],npar)
	#TIME PARAMS WILL BE MADE PROPORTIONAL TO EACH OTHER TO BE BIN LENS OF TOTAL INTERVAL
	answer<-optim(par_init,support_ft,method="L-BFGS-B",lower=par_min,
		upper=par_max,control=list(maxit=iterations))
	#answer
	if(answer$convergence!=0){message("Warning: optimizer did not converge; see help page.")}
	par<-answer$par
	names(par)<-NULL
	if(est_only){
		qr_predict<-qr_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
		qr<-qr_predict(par)
		colnames(qr)<-c("qMax","rMax")
		res<-qr
	}else{
	mes<-answer$message
	conv<-answer$convergence
	SMax<-(-answer$value)
	n<-length(dur)
	aicc<-(2*npar)-(2*SMax)+((2*npar*(npar+1))/(n-npar-1))
	if((n-npar-1)<1){aicc<-"calc failed, npar>(N-1) !"}
	title<-paste("Analysis with",n_tbins,"time bins and",
		sum(c(g1s>0,g2s>0)),"groupings (",g1s,"and",g2s,
		"States),with",npar,"parameters and",n,"taxa")
	title<-c(title,"Note that parameter output is rate components, NOT the typical rate for a member of that group. See help file.")
	if(n_tbins>1){
		n_tb<-n_tbins-1
		t_raw<-c(0.5,par[1:n_tb])
		t_prop<-t_raw/sum(t_raw)
		t_bl<-t_prop*(max(FO)-min(LO))
		t_ends<-c(max(FO),max(FO)-cumsum(t_bl))
		mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
		Pp<-mqr[,2]/(mqr[,1]+mqr[,2])
		mqr<-cbind(mqr,Pp)
		colnames(mqr)<-c("extRate","sampRate","Completeness")
		mqrt<-mqr[1:n_tbins,]
		if(g1s>0){
			mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]	
			if(g2s>0){
				mqrg2<-mqr[(n_tbins+1+g1s):(n_tbins+g1s+g2s),]
				res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqrt[,1:2],pars.grp1=mqrg1[,1:2],
					pars.grp2=mqrg2[,1:2],int.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
			}else{
				res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqrt[,1:2],
					pars.grp1=mqrg1[,1:2],int.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
				}
		}else{
			res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.time=mqrt,
				int.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
			}	
	}else{
		mqr<-matrix(par,,2,byrow=TRUE)
		Pp<-mqr[,2]/(mqr[,1]+mqr[,2])
		mqr<-cbind(mqr,Pp)
		colnames(mqr)<-c("extRate","sampRate","Completeness")
		if(g1s>0){
			mqrg1<-mqr[1:g1s,]
			if(g2s>0){
				mqrg2<-mqr[(g1s+1):(g1s+g2s),]
				res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.grp1=mqrg1[,1:2],pars.grp2=mqrg2[,1:2],
					convergence=conv,message=mes)
			}else{
				res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.grp1=mqrg1,message=mes)
				}
		}else{
			colnames(mqr)<-NULL
			rMax<-mqr[,2]
			qMax<-mqr[,1]
			Pp<-rMax/(rMax+qMax) #From Solow and Smith, 1997
			res<-list(Title=title[1],parameters=c("extRate"=qMax,"sampRate"=rMax,"Completeness"=Pp),
				log.likelihood=-SMax,AICc=aicc,convergence=conv,message=mes)
			}
		}}
	return(res)
	}


getSampProbDisc<-function(timeData,n_tbins=1,grp1=NA,grp2=NA,est_only=FALSE,iterations=10000,initial=0.5){
	#this is the multi-parameter maximum likelihood analysis of binned timeData
		#uses a set of binned-interval timeData (just the by-species first and last intervals matrix) 
			#to fit models of different samp probs and ext rates
			#output from binTimeData() can be input directly (only looks at second matrix)
		#can allow for free-moving time windows and different groups
			#these models can then be compared with AIC
		#if est_only=TRUE, then the q and R estimates will be given back per-species
		#this actually runs very slowly; the current settings are optimized for speed (throttle=1). 
			#Increase throttle to 2-4 for inc accuracy
	#x<-runif(100);x<-cbind(x+rexp(100),x);y<-binTimeData(x);getSampProbDisc(y[[2]],est_only=FALSE)
	#x<-runif(100);x<-cbind(x+rexp(100),x);timeData<-binTimeData(x)[[2]];n_tbins=1;grp1=NA;grp2=NA;est_only=FALSE;throttle=1
	#######################
	#REQUIRED FUNCTIONS BELOW
	#######################
	getPp<-function(R,q){
		#calculate completeness given R and q
		res<-numeric()
		for(Ti in 1:10000){
			res[Ti]<-(1-(1-R)^Ti)*(exp(-q*(Ti-1))-exp(-q*Ti))
			}
		sum(res)
		}
	###########################
	f_dur<-function(q,R,dur){
		TMax<-max(dur)
		Ti<-1:TMax	#unique durations
		N<-sapply(Ti,function(t) sum(t==dur))		#number with those durations
		PDT<-exp(-q*(Ti-1))-exp(-q*Ti)
		Rex<-c(1,rep(2,TMax-1))
		ft<- sapply(Ti,function(t) sum(PDT[t:TMax]*((R^Rex[t])*(Ti[t:TMax]-t+1)*((1-R)^(Ti[t:TMax]-t)))))
		ft<-log(ft/sum(ft))	#divide by sum; same as normalizing by Pp (really? COOL!)
		res<-sum((N*ft)[!is.infinite(ft)])
		return(res)
		}
	####################
	make_ft<-function(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s){
		qR_predict<-qR_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
		function(par){
			qR<-qR_predict(par)
			q<-qR[,1]; R<-qR[,2]/10
			data<-cbind(q,R,dur)
			class<-apply(qR,1,function(x) which(apply(unique(qR),1,function(y) all(x==y))))	#get unique classes of q,R
			f<-sapply(unique(class),function(x) f_dur(q[x==class][1],R[x==class][1],dur[x==class]))
			res<-(-sum(f))
			return(res)
			}
		}
	####################
	qR_predict_multpar<-function(FO,LO,n_tbins,grp1,grp2,g1s,g2s){
	#par consists of a vector, with first n_tbins-1 elements represent bin LENGTHS
		#unless n_tbins==1
		#after that, its a matrix with k rows, first col is q, second is R
			#first comes rows associated with each time bin, 
			#than each group in grp1, than grp2
	#THIS FUNCTION SHOULD OUTPUT A Nx2 MATRIX OF q and R VALUES to the support function
		#MWAHAHAH!
	#note that little r is used here as a shorthand for big R (samp prob, not samp rate)
	if(n_tbins>1){	#IF THERE ARE TIME BINS
		if(g1s>0){	#IF THERE ARE GROUPS + TIMEBINS
			if(g2s>0){	#IF THERE IS TWO GROUPS +TIMEBINS
				function(par){
					n_tb<-n_tbins-1
					t_raw<-c(0.5,par[1:n_tb])
					t_prop<-t_raw/sum(t_raw)
					t_bl<-t_prop*(max(FO)-min(LO))
					tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
					mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
					mqrt<-matrix(mqr[1:n_tbins,],,2,byrow=TRUE)
					mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]
					mqrg2<-mqr[(n_tbins+g1s+1):(n_tbins+g1s+g2s),]
					tcat<-sapply(FO,function(x) sum(tbin>=x))
					qrt<-sapply(tcat,function(x) mqrt[x,])
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					qrg2<-sapply(grp2,function(x) mqrg2[x,])
					t((qrt+qrg1+qrg2)/3)
				}
			}else{	#IF THERE IS ONE GROUP +TIMEBINS
				function(par){
					n_tb<-n_tbins-1
					t_raw<-c(0.5,par[1:n_tb])
					t_prop<-t_raw/sum(t_raw)
					t_bl<-t_prop*(max(FO)-min(LO))
					tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
					mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
					mqrt<-matrix(mqr[1:n_tbins,],,2,byrow=TRUE)
					mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]
					tcat<-sapply(FO,function(x) sum(tbin>=x))
					qrt<-sapply(tcat,function(x) mqrt[x,])
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					t((qrt+qrg1)/2)
				}
			}
		}else{	#IF THERE ARE NO GROUPS AND TIMEBINS
			function(par){
				n_tb<-n_tbins-1
				t_raw<-c(0.5,par[1:n_tb])
				t_prop<-t_raw/sum(t_raw)
				t_bl<-t_prop*(max(FO)-min(LO))
				tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
				mqrt<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
				tcat<-sapply(FO,function(x) sum(tbin>=x))
				tcount<-sapply(sort(unique(tcat)),function(x) sum(tcat==x))
				#if(all(tcount>100)){t(sapply(tcat,function(x) mqrt[x,]))
				#}else{matrix(0.99,length(FO),2,byrow=TRUE)}
				t(sapply(tcat,function(x) mqrt[x,]))
			}
		}
	}else{	#IF THERE ARE NO TIMEBINS
		if(g1s>0){	#IF THERE ARE GROUPS AND NO TIMEBINS
			if(g2s>0){	#IF THERE IS TWO GROUPS AND NO TIMEBINS
				function(par){
					mqr<-matrix(par,,2,byrow=TRUE)
					mqrg1<-mqr[1:g1s,]
					mqrg2<-mqr[(g1s+1):(g1s+g2s),]
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					qrg2<-sapply(grp2,function(x) mqrg2[x,])
					t((qrg1+qrg2)/2)
				}
			}else{	#IF THERE IS ONE GROUP AND NO TIMEBINS
				function(par){
					mqr<-matrix(par,,2,byrow=TRUE)
					t(sapply(grp1,function(x) mqr[x,]))
				}
			}
		}else{	#IF THERE ARE NO GROUPS AND NO TIMEBINS (2 param model)
			function(par){
				matrix(par,length(FO),2,byrow=TRUE)
				}	
			}
		}			
	}
	#############################
	if(length(timeData)==2){timeData<-timeData[[2]]}	#if a timeList matrix...
	#get rid of any NAs
	if(length(grp1)>1){
		if(length(grp1)!=nrow(timeData)){stop("Error: grp1 is not same length as timeData")}
		grp1<-grp1[!is.na(timeData[,1])]
		}
	if(length(grp2)>1){
		if(length(grp2)!=nrow(timeData)){stop("Error: grp2 is not same length as timeData")}		
		grp2<-grp2[!is.na(timeData[,1])]
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(timeData,1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	dur<-apply(timeData,1,diff)+1
	timeData1<-max(timeData)-timeData+1
	FO<-timeData1[,1];LO<-timeData1[,2]
	#NOW CONSTRUCT PROPER SUPPORT FUNCTION
	#NOTE THAT THE FIRST BIN MUST ALWAYS BE THE FIRST INTERVAL
	if(length(grp1)>1){g1s<-length(unique(grp1))}else{g1s<-0}
	if(length(grp2)>1){g2s<-length(unique(grp2))}else{g2s<-0}
	support_ft<-make_ft(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s)
	npar<-(n_tbins-1)+(ifelse(n_tbins>1,2*n_tbins,0)+ifelse(g1s>0,g1s*2,0)
		+ifelse(g2s>0,g2s*2,0))
	if(npar==0){npar<-2}
	par_lim<-c(0.0001,10)
	par_init<-rep(initial,npar)
	par_min<-rep(par_lim[1],npar)
	par_max<-rep(par_lim[2],npar)
	#TIME PARAMS WILL BE MADE PROPORTIONAL TO EACH OTHER TO BE BIN LENS OF TOTAL INTERVAL
	answer<-optim(par_init,support_ft,method="L-BFGS-B",lower=par_min,
		upper=par_max,control=list(maxit=iterations))		#,trace=1
	#answer
	par<-answer$par
	if(answer$convergence!=0){message("Warning: optimizer did not converge; see help page.")}
	if(est_only){
		qR_predict<-qR_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
		qR<-qR_predict(par)
		colnames(qR)<-c("qMax","RMax")
		qR[,2]<-qR[,2]/10
		res<-qR
	}else{
		mes<-answer$message
		conv<-answer$convergence
		SMax<-(-answer$value)
		n<-length(dur)
		aicc<-(2*npar)-(2*SMax)+((2*npar*(npar+1))/(n-npar-1))
		if((n-npar-1)<1){aicc<-"calc failed, npar>(N-1) !"}
		title<-paste("Analysis with",n_tbins,"time bins and",
			sum(c(g1s>0,g2s>0)),"groupings (",g1s,"and",g2s,
			"States),with",npar,"parameters and",n,"taxa")
		title<-c(title,"Note that parameter output is rate components, NOT the typical rate for a member of that group. See help file.")
		if(n_tbins>1){
			n_tb<-n_tbins-1
			t_raw<-c(0.5,par[1:n_tb])
			t_prop<-t_raw/sum(t_raw)
			t_bl<-t_prop*(max(FO)-min(LO))
			t_ends<-c(max(FO),max(FO)-cumsum(t_bl))
			mqR<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
			mqR[,2]<-mqR[,2]/10
			Pp<-sapply(1:nrow(mqR),function(x) getPp(mqR[x,2],mqR[x,1]))
			mqR<-cbind(mqR,Pp)		
			colnames(mqR)<-c("extRate","sampProb","Completeness")
			mqRt<-mqR[1:n_tbins,]
			if(g1s>0){
			mqRg1<-mqR[(n_tbins+1):(n_tbins+g1s),]	
				if(g2s>0){
					mqRg2<-mqR[(n_tbins+1+g1s):(n_tbins+g1s+g2s),]
					res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqRt[,1:2],pars.grp1=mqRg1[,1:2],
						parts.grp2=mqRg2[,1:2],interval.boundaries=t_ends,interval.lengths=t_bl,
						convergence=conv,message=mes)
				}else{
					res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqRt[,1:2],
						pars.grp1=mqRg1[,1:2],interval.boundaries=t_ends,interval.lengths=t_bl,
						convergence=conv,message=mes)
					}
			}else{
				res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.time=mqRt,
					interval.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
				}	
		}else{
			mqR<-matrix(par,,2,byrow=TRUE)
			mqR[,2]<-mqR[,2]/10
			Pp<-sapply(1:nrow(mqR),function(x) getPp(mqR[x,2],mqR[x,1]))
			mqR<-cbind(mqR,Pp)
			colnames(mqR)<-c("extRate","sampProb","Completeness")
			if(g1s>0){
				mqRg1<-mqR[1:g1s,]
				if(g2s>0){
					mqRg2<-mqR[(g1s+1):(g1s+g2s),]
					res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.grp1=mqRg1[,1:2],pars.grp2=mqRg2[,1:2],
						convergence=conv,message=mes)
				}else{
					res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.grp1=mqRg1,convergence=conv,message=mes)
					}
			}else{
				colnames(mqR)<-NULL
				RMax<-mqR[,2]
				qMax<-mqR[,1]
				Pp<-sapply(1:nrow(mqR),function(x) getPp(mqR[x,2],mqR[x,1])) #From Foote, 1997
				res<-list(Title=title[1],parameters=c("extRate"=qMax,"sampProb"=RMax,"Completeness"=Pp),
					log.likelihood=SMax,AICc=aicc,convergence=conv,message=mes)
				}
		}}
	return(res)
	}

binTimeData<-function(timeData,int.length=1,start=NA,int.times=NULL){
	#bin temporal data
	#input: continuous time data (two column of FADs and LADs)
	#output: a list with two 2-col matrices as elements, bin-times and taxon occurences
			#intervals, UNLIKE TIME, always go up (earliest is 1 and increase...)
		#arbitrarily starts bin at the first fad; this can be changed by setting 'start'
			#start must be greater than max(timeData)
			#the last bin is cut off at zero (present day)
	#x<-runif(100);x<-cbind(x+rexp(100),x)
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data?")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	if(is.null(int.times)){
		if(is.na(start)){start<-max(timeData)+int.length}else{if(start<max(timeData)){stop("Error:Start<max(timeData)?")}}
		end<-start-(ceiling((start-min(timeData))/int.length)+1)*int.length
		bins<-seq(start,end,by=-int.length)
		bins<-unique(ifelse(bins<0,0,bins))	#get rid of any extra zeroes or negative numbers
		fads<-sapply(timeData[,1],function(x) which(bins<x)[1]-1)
		lads<-sapply(timeData[,2],function(x) which(bins<x)[1]-1)
		res<-list(int.times=cbind(int.start=bins[1:(length(bins)-1)],int.end=bins[2:length(bins)]),
			taxon.times=cbind(first.int=fads,last.int=lads))
	}else{
		int.durs<-int.times[,1]-int.times[,2]
		if(any(int.durs<=0)){stop("Error: Some input time intervals have zero or negative durations?")}
		int.times<-int.times[order(int.durs),]
		Fint<-sapply(timeData[,1],function(x) which(apply(int.times,1,function(y) y[1]>=x & y[2]<x))[1])
		Lint<-sapply(timeData[,2],function(x) which(apply(int.times,1,function(y) y[1]>=x & y[2]<x))[1])
		taxon.times<-cbind(first.int=Fint,last.int=Lint)
		rownames(taxon.times)<-rownames(timeData)
		taxon.times<-taxon.times[!apply(taxon.times,1,function(x) any(is.na(x))),]
		new.order<-rank(-int.times[,1])	
		taxon.times[,1]<-new.order[taxon.times[,1]]
		taxon.times[,2]<-new.order[taxon.times[,2]]
		int.times<-int.times[order(-int.times[,1]),]
		res<-list(int.times=int.times,taxon.times=taxon.times)
		}
	return(res)
	}

dropZLB<-function(tree){
	#drops terminal branches that are zero length
		#adjusts tree$root.time if necessary
	require(ape)
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	drop_e<-(tree$edge[,2]<(Ntip(tree)+1)) & (tree$edge.length==0)
	drop_t<-(tree$edge[,2])[drop_e]
	if((Ntip(tree)-length(drop_t))>1){
		tree1<-drop.tip(tree,drop_t)
		orig_dist<-dist.nodes(tree)[which(tree1$tip.label[1]==tree$tip.label),Ntip(tree)+1]
		new_dist<-dist.nodes(tree1)[1,Ntip(tree1)+1]
		if(!is.null(tree$root.time)){tree1$root.time<-tree$root.time-(orig_dist-new_dist)}
		res<-tree1
	}else{res<-NA}
	return(res)
	}

sampleRanges<-function(taxad,r,alpha=1,beta=1,rTimeRatio=1,modern.samp.prob=1,min.taxa=2,
	ranges.only=TRUE,minInt=0.01,merge.cryptic=TRUE,alt.method=FALSE,plot=FALSE){
	#sample ranges using a taxad matrix as input 
	#if (ranges.only=TRUE): outputs matrix of FADs/LADs, with NAs for unsampled taxa
	#if (ranges.only=FALSE): outputs per-species list with vectors of dates where that species was sampled
		#ranges and occurance are output on a BACKWORD-moving timescale as expected for paleo data
	#if modern.samp.prob=1, then all still-living taxa (taxa at 0 for LAD) are ALWAYS last observed at zero
		#this approximates the fact that we think the present-day living biota is almost perfectly sampled
			#(well, relative to the modern)
	#if r is a vector, then it is considered to be a vector of per-species sampling rates
		#r is the mean sampling rate for a species
	#rTimeRatio is the ratio of the latest taxon-avg sampling rates by the earliest taxonavg sampling rates
		#i.e. the proportional increase over (a) the entire taxon's history or (b) the taxon duration (if specific)
		#either a single value for all taxa or taxon-specific values
	#all parameters can be given as single values or species-specific values
	#names<-paste("t",1:4,sep="");taxad<-cbind(c(250,230,210,200),c(240,215,205,150))
	#min.taxa=2;minInt=0.01;modern.samp.prob=1.0;plot=T;ranges.only=F;alt.method=F
	#r<-c(0.2,0.1,0.3,0.4);alpha<-4;beta<-4;rTimeRatio<-2
	#r<-c(0,0.1,0.3,0.4);alpha<-beta<-rTimeRatio<-2
	if(ncol(taxad)==6){				#also allow it to accept taxad objects
		living<-taxad[,5]
		cryptic<-sapply(1:nrow(taxad),function(x) if(taxad[x,1]!=taxad[x,6]){which(taxad[,1]==taxad[x,6])}else{NA})
		timeData<-taxad[,3:4,drop=FALSE]
		names<-if(is.null(rownames(taxad))){paste("t",taxad[,1],sep="")}else{rownames(taxad)}
		}
	if(ncol(taxad)==2){			#assumes it has two matrices
		living<-NULL
		cryptic<-rep(NA,nrow(taxad))
		timeData<-taxad
		names<-if(is.null(rownames(taxad))){paste("t",1:nrow(taxad),sep="")}else{rownames(taxad)}
		}	
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	#check input parameters for inconsistencies
	if(any(r<0)){stop("Error: some r values < 0")}
	if(any(alpha<=0) | any(beta<=0)){stop("Error: some shape parameters (alpha, beta) <= 0")}
	if(length(alpha)==1){alpha<-rep(alpha,length(names))
		}else{if(length(alpha)!=length(names)){
		stop("Error: Multiple alpha values input but not same length as number of taxa?")}}
	if(length(beta)==1){beta<-rep(beta,length(names))
		}else{if(length(beta)!=length(names)){
		stop("Error: Multiple beta values input but not same length as number of taxa?")}}
	if(length(rTimeRatio)!=1){if(length(rTimeRatio)!=length(names)){
		stop("Error: Multiple rTimeRatio values input but not same length as number of taxa?")}}
	midpoint<-(max(timeData)+min(timeData))/2
	taxa.midpoints<-((timeData[,1]+timeData[,2])/2)
	#need to calculate rTimeChange based on rTimeRatio
	r_start<-r/((rTimeRatio/2)+0.5)
	r_end<-2*r*(1-1/(rTimeRatio+1))
	if(length(r)==1 & length(rTimeRatio)==1){
		#thanks to Emily King and Brian Koch for helping me with this bit!
		rTimeChange<-(r_end-r_start)/(max(timeData)-min(timeData))
		rTimeChange<-rep(rTimeChange,length(names))
	}else{	#if there's multiple r values, multiple ratios or both...
		dur1<-timeData[,1]-timeData[,2]
		rTimeChange<-(r_end-r_start)/dur1
		}
	if(length(r)==1){		#get per-taxon r.avg for each taxon, if not input already
		taxa.timeChange<-midpoint-taxa.midpoints
		r<-r+(rTimeChange*taxa.timeChange)
		}
	if(length(r)!=length(names)){stop("Error: Multiple r values given but not same length as number of taxa?")}
	if(any(r<0)){stop("Error: rTimeRatio so low as to cause some species-specific avg r < 0")}
	if(any(alpha!=1) | any(beta!=1) | any(rTimeChange!=0) | alt.method){		#get per species time vectors of r
		rangesTimes<-apply(timeData,1,function(x) seq(x[1],x[2],by=-minInt))
		#get the time change component per minInt of the taxon ranges
		#first get the different of each min int from the midpoints, then calculate the change r accordint to rTimeChange
		taxa.shiftTimes<-lapply(1:length(r),function(x) taxa.midpoints[x]-rangesTimes[[x]])
		taxa.rShiftTime<-lapply(1:length(r),function(x) rTimeChange[x]*taxa.shiftTimes[[x]])
		#now get the hat component, combine with the time component
		scaledTimes<-lapply(rangesTimes,function(x) seq(0,1,length.out=length(x)))
		rHat<-lapply(1:length(r),function(x) r[x]*dbeta(scaledTimes[[x]],alpha[x],beta[x]))
		rHatTime<-lapply(1:length(r),function(x) rHat[[x]]+taxa.rShiftTime[[x]])
		#transform so that the hats always are above zero but have correct means and correct timechange slope (so complicated!)
		#rHatTime<-lapply(rHatTime,function(x) (x-min(x))/(1-(min(x)/mean(x))))
		rHatTime1<-lapply(rHatTime,function(x) ifelse(x<0,0,x))
		rHatTime<-lapply(1:length(rHatTime),function(x) if(mean(rHatTime1[[x]])>0){
				mean(rHatTime[[x]])*rHatTime1[[x]]/mean(rHatTime1[[x]])
			}else{mean(rHatTime[[x]])*rHatTime1[[x]]}
			)
	}else{	
		rHatTime<-NULL
		}
	if(plot){if(!is.null(rHatTime)){
		plot(0:1,c(min(unlist(rHatTime)),max(unlist(rHatTime))),type="n",xlim=c(max(timeData),min(timeData)),	
			xlab="Time (Time-Units before Present)",ylab=c("Instant. Sampling Rate","(per lineage time-units)"))
		for(i in 1:length(r)){lines(rangesTimes[[i]],rHatTime[[i]],lwd=2)}
	}else{
		plot(c(max(timeData),min(timeData)),c(min(r),max(r)),type="n",xlim=c(max(timeData),min(timeData)),
			xlab="Time (Time-Units before Present)",ylab=c("Instant. Sampling Rate","(per lineage time-units)"))
		for(i in 1:length(r)){lines(timeData[i,],c(r[i],r[i]),lwd=2)}
		}}
	#plot(rangesTimes[[2]],rHatTime[[2]],type="l",lwd=3,xlim=c(max(rangesTimes[[2]]),
	#	min(rangesTimes[[2]])),xlab="time",ylab="Instantaneous Sampling Rate (per Lmy)")
	#(x<-cbind(r,sapply(rHatTime,mean)));plot(x);abline(0,1)	#means should be close to actual r
	#okay, now to actually do the sampling
	redo<-TRUE
	while(redo){
		samp_occ<-list()	#sampled occurances
		for(i in 1:nrow(timeData)){
			if(is.null(rHatTime)){
				if(r[i]>0){
					samps<-timeData[i,1]	#the time of speciation is the lower bound
					while(min(samps)>timeData[i,2]){	#keep sampling until you go past the time of extinction
						samps<-c(samps,min(samps)-rexp(1,rate=r[i]))}
					samps<-samps[-c(1,length(samps))]
				}else{samps<-numeric(0)}
			}else{
				samps<-rangesTimes[[i]][sapply(rHatTime[[i]],function(x) runif(1)<=(x*minInt))]
				}
			if(!is.null(living) & modern.samp.prob>0){if(living[i]==1){		#rework modern.sampling as a probability
				if(runif(1)<=modern.samp.prob){samps<-c(samps,0)}
				}}
			if(length(samps)>0){samp_occ[[i]]<-samps}else{samp_occ[[i]]<-NA}
			}
		redo<-min.taxa>sum(sapply(samp_occ,function(x) !any(is.na(x))))
		}
	names(samp_occ)<-names
	if(any(!is.na(cryptic)) & merge.cryptic){for(i in which(!is.na(cryptic))){
		samp_occ[[cryptic[i]]]<-c(samp_occ[[cryptic[i]]],samp_occ[[i]])
		samp_occ[[i]]<-NA
		}}
	if(ranges.only){
		ranges<-cbind(sapply(samp_occ,max),sapply(samp_occ,min))
		rownames(ranges)<-names;colnames(ranges)<-c("FAD","LAD")
		res<-ranges
	}else{
		res<-samp_occ
		}
	return(res)
	}

taxa2cladogram<-function(taxad,drop.cryptic=FALSE,plot=FALSE){
	#take a taxad and turn it into an unscaled cladogram
		#do this by forming the tree as newick format first
		#essentially, this algorithm works by defining clades as only those sets of taxa which start with the FAD of the first taxon
		#This could be a 'clade' of a single OTU
		#or a clade of a bunch of taxa where one is a budding ancestor and the other are budding descendants or whatever
		#the important thing is that there's only so many nodes as there are instances where morphotaxa end/start allowing for homologies
	require(ape)
	if(any(taxad[,6]!=taxad[,1])){
		for(i in which(taxad[,1]!=taxad[,6])){
			#reset descendants of cryptic taxa so all bud off of first cryptic species		
			taxad[taxad[,2]==taxad[i,1],2]<-taxad[i,6]
			}
		}
	tlabs<-rownames(taxad)
	desc<-lapply(taxad[,1],function(x) (taxad[taxad[,2]==x,1])[!is.na(taxad[taxad[,2]==x,1])])
	ndesc<-sapply(desc,length)
	rank<-numeric(length(ndesc))
	rank[ndesc==0]<-1
	rank[rank==0]<-NA
	while(any(is.na(rank))){
		rank<-sapply(1:length(rank),function(x) ifelse(!is.na(rank[x]),rank[x],
				1+max(rank[sapply(desc[[x]],function(y) which(y==taxad[,1]))])))}
	#okay, now all taxa are ranked by their depth from the tips
	comp<-numeric(length(ndesc))
	lab<-list()
	lab[rank==1]<-tlabs[rank==1]
	comp[rank==1]<-1
	while(any(comp==0)){
		tpot<-comp==0
		tpot2<-rank==min(rank[tpot])
		tpick<-which(tpot & tpot2)[1]
		dlab<-paste(unlist(lab[desc[[tpick]]]),",",sep="",collapse="")
		lab[[tpick]]<-paste("(",dlab,tlabs[tpick],")",sep="")
		comp[tpick]<-1
		}
	tree1<-paste(lab[[1]],";",sep="")
	tree2<-read.tree(text=tree1)
	if(drop.cryptic & any(taxad[,6]!=taxad[,1])){
		tree2<-drop.tip(tree2,tlabs[taxad[,6]!=taxad[,1]])
		tree2<-collapse.singles(tree2)
		}
	if(plot){plot(ladderize(tree2),show.tip.label=FALSE)}
	return(tree2)
	}

taxa2phylo<-function(taxad,obs_time=NULL,plot=FALSE){
	#INPUT a taxad object and a vector of observation times for each species
		#if obs=NULL, LADs in taxad1 are used as observation times
		#all times must be in backwards format (zero is present)
	#root.time
		#ALL TREES ARE OUTPUT WITH ELEMENTs "$root.time"
		#this is the time of the root on the tree, which is important for comparing across trees
		#this must be calculated prior to adding anything to terminal branches
	#OUTPUT an ape phylo object with the tips at the times of observation
	require(ape)
	taxad1<-taxad[,1:4]
	if(any((taxad1[,4]-taxad1[,3])<0)){taxad1[,3:4]<-max(taxad1[,3:4])-taxad1[,3:4]}
	if(any((taxad1[,4]-taxad1[,3])<0)){stop("Error: Time Error! Check data in taxad")}
	if(is.null(obs_time)){obs<-taxad1[,4]}else{obs<-max(taxad[,3:4])-obs_time}
	if(nrow(taxad1)!=length(obs)){stop("#obs != # lineages !")}
	#make observations as fake taxa, assuming that observations are WITHIN actual taxon ranges
	fake_taxa<-matrix(sapply((1:nrow(taxad1))[!is.na(obs)],function(x) c(nrow(taxad1)+x,taxad1[x,1],obs[x],obs[x])),,4,byrow=TRUE)
	fake_taxa[,1]<-(1:nrow(fake_taxa))+nrow(taxad1)
	taxad2<-rbind(taxad1,fake_taxa)
	ntaxa<-nrow(taxad2)
	#MAKE IT INTO AN NODE/EDGE-BASED PHYLOGENY
	#find descendents of every taxon
	desc<-lapply(taxad2[,1],function(x) taxad2[c(FALSE,taxad2[-1,2]==x),1])						#desc of each taxon
	births2<-lapply(desc,function(x) if(length(x)>0){sapply(x,function(y) taxad2[y==taxad2[,1],3])})	#time of desc births
	desc<-lapply(1:length(desc),function(x)	#MORE STUPIDLY COMPLICATED CODE
		if(length(desc[[x]])>0){desc[[x]][match(1:length(births2[[x]]),rank(births2[[x]],ties.method="random"))]}else{numeric()})
	births<-lapply(desc,function(x) if(length(x)>0){sapply(x,function(y) taxad2[y==taxad2[,1],3])})	#time of desc births
	#make events list: first event is taxon birth, with that taxon as desc, next is desc births, extinction not recorded but implied
	events<-lapply(1:ntaxa,function(x) if(length(desc[[x]])>0){c(taxad2[x,1],desc[[x]])}else{c(taxad2[x,1])})
	#times of events: taxon birth, desc births, extinction
	times<-lapply(1:ntaxa,function(x) if(length(births[[x]])>0){c(taxad2[x,3],births[[x]],taxad2[x,4])}else{c(taxad2[x,3],taxad2[x,4])})
	#labels for each lineage segment between event times (these may as well represent the daughter node ID too)
		#use decimals to keep track of segments, set first segment as X.0
		#as long as a single taxon doesn't have millions of descendants, in which case the IDs may stop being unique...
	nseg<-sapply(times,length)-1
	seg_labs<-lapply(1:ntaxa,function(x) taxad2[x,1]+(1:nseg[x]/(nseg[x]+1)))
	seg_labs<-lapply(seg_labs,function(x) c(floor(x[1]),x[-1]))
	#now, find the mother segment for each taxon
	taxa_anc<-c(0,sapply(2:ntaxa,function(x) unlist(seg_labs[taxad2[x,2]==taxad2[,1]])[which(unlist(events[taxad2[x,2]==taxad2[,1]])==taxad2[x,1])-1]))
	moms2<-lapply(seg_labs,function(x) x[-length(x)])	#now make list of all seg ids for all mom segs
	moms<-lapply(1:ntaxa,function(x) c(taxa_anc[x],unlist(moms2[[x]])))
	lengths<-lapply(times,diff)	#get lengths of segments 
	term<-sapply(unlist(seg_labs),function(x) !any(unlist(moms)==x))	#which are terminal?
	#make edge data.frame, with id, anc-id, length and logical indicating terminal branches
	edgeD<-data.frame(id=unlist(seg_labs),anc=unlist(moms),brlen=unlist(lengths),term=term)
	MRCA<-min(edgeD[sapply(edgeD$id,function(x) 1<sum(x==edgeD$anc)),1])
	edgeD<-edgeD[-which(edgeD[,1]<=MRCA),]	#want to drop any extraneous root as well
	#add/leave only the terminals needed for the observations, remove others
	droppers<-which(edgeD[,4] & edgeD[,1]<(nrow(taxad1)+1))
	while(length(droppers)>0){
		edgeD<-edgeD[-droppers,]
		edgeD[,4]<-sapply(edgeD[,1],function(x) !any(edgeD[,2]==x))
		droppers<-which(as.logical(edgeD[,4]) & edgeD[,1]<(nrow(taxad1)+1))
		}
	edgeD[,4]<-sapply(edgeD[,1],function(x) !any(edgeD[,2]==x))
	#collapse single internal nodes
	ndesc<-sapply(edgeD[,1],function(x) sum(x==edgeD[,2]))	#ndesc from each node
	while(any(ndesc==1)){		#picks only internal branches with 1 desc
		epick<-edgeD[ndesc==1,];if(is.data.frame(epick)){epick<-epick[1,]}	#pick a single, if matrix, use first one
		edesc<-edgeD[edgeD$anc==epick$id,]
		#remove desc
		edgeD<-edgeD[-which(edgeD$id==edesc$id),]
		#replace picked edge w/new edge
		newe<-data.frame(id=edesc$id,anc=epick$anc,brlen=(epick$brlen+edesc$brlen),term=edesc$term)
		edgeD[edgeD$id==epick$id,]<-newe
		ndesc<-sapply(edgeD$id,function(x) sum(x==edgeD$anc))
		}
	ndesc<-sapply(edgeD$id,function(x) sum(x==edgeD$anc))
	#replace stupid decimal edge placeholders with clean numbers
	e_fix<-numeric()
	e_fix[edgeD$term]<-sapply(edgeD[edgeD$term,1],function(x) which(fake_taxa[,1]==x)) 	#the species in the original data
	e_fix[!edgeD$term]<-sum(edgeD$term)+1+(1:sum(!edgeD$term))
	ea_fix<-sapply(edgeD$anc,function(x) ifelse(x!=MRCA,e_fix[edgeD$id==x],sum(edgeD$term)+1))
	#NOW MAKE A TREE
	tlabs<-rownames(taxad1)[!is.na(obs)]
	edgf<-cbind(ea_fix,e_fix);colnames(edgf)<-NULL
	tree1<-list(edge=edgf,tip.label=tlabs,edge.length=edgeD[,3],Nnode=length(unique(edgf[,1])))
	class(tree1)<-"phylo"						#ITS A TREE!
	tree<-reorder(collapse.singles(tree1),"cladewise") 	#REORDER IT
	if(plot){plot(ladderize(tree),show.tip.label=FALSE);axisPhylo()}
	#now, root.time should be the time of the first obs PLUS the distance from the earliest tip to the root
	first_obs_time<-max(taxad[,3:4])-min(obs,na.rm=TRUE)
	tree$root.time<-first_obs_time+min(dist.nodes(tree)[1:Ntip(tree),Ntip(tree)+1])
	return(tree)
	}

dropExtant<-function(tree,tol=0.01){
	#drop all terminal taxa that are more than 0.001 from the modern
	require(ape)
	if(is.null(tree$root.time)){
		message("Warning: no tree$root.time! Assuming latest tip is at present (time=0)")
		}
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	dnode<-dist.nodes(tree)[1:Ntip(tree),Ntip(tree)+1]
	if(!is.null(tree$root.time)){if(tree$root.time>max(dnode)){stop("Error: all tips are extinct based on tree$root.time!")}}
	droppers<-which((dnode+tol)>max(dnode))
	if((Ntip(tree)-length(droppers))<2){stop("Error: Less than 2 tips extinct on the tree!")}
	stree<-drop.tip(tree,droppers)
	if(!is.null(tree$root.time)){
		#now need to add $root.time given the droppers
			#should be root.time MINUS distance from furthest tip in tree PLUS distance from latest tip to root of stree
		stree$root.time<-tree$root.time-max(dnode)+max(dist.nodes(stree)[1:Ntip(stree),Ntip(stree)+1])
		}
	return(stree)
	}

dropExtinct<-function(tree,tol=0.01,ignore.root.time=FALSE){
	#drop all terminal taxa that are less than 0.001 from the modern
	require(ape)
	if(is.null(tree$root.time)){
		message("No tree$root.time: Assuming latest tip is at present (time=0)")
		}
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	dnode<-dist.nodes(tree)[1:Ntip(tree),Ntip(tree)+1]
	if(!is.null(tree$root.time) & !ignore.root.time){if(tree$root.time>max(dnode)){stop("Error: all tips are extinct based on tree$root.time!")}}
	droppers<-which((dnode+tol)<max(dnode))
	if((Ntip(tree)-length(droppers))<2){stop("Error: Less than 2 tips are extant on the tree!")}
	stree<-drop.tip(tree,droppers)
	if(!is.null(tree$root.time)){
		#now need to add $root.time given the droppers
			#should be root.time MINUS distance from furthest tip in tree PLUS distance from latest tip to root of stree
		stree$root.time<-tree$root.time-max(dnode)+max(dist.nodes(stree)[1:Ntip(stree),Ntip(stree)+1])
		}
	return(stree)
	}

timeSliceTree<-function(ttree,sliceTime,drop.extinct=FALSE,plot=TRUE){
	#take a phylogeny and produce a phylogenetic 'slice' at time X (w/respect to root.time)
		#lineages extant at sliceTime sliced to end at that point
		#if no root.time, then it is assumed the tip furthest from the root is at 0 (present-day)
			#a warning will be issued if this is assumed
		#extinct tips will be dropped if drop.extinct=TRUE
	require(ape)
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
	stree$root.time<-ttree$root.time
	if(drop.extinct){
		stree1<-dropExtinct(stree,ignore.root.time=TRUE)
	}else{stree1<-stree}
	if(plot){layout(1:2);plot(ladderize(ttree),show.tip.label=FALSE);axisPhylo();
		plot(ladderize(stree1),show.tip.label=FALSE);layout(1)}
	return(stree1)
	}

plotTraitgram<-function(trait,tree,trait.name="'trait'",conf.int=TRUE,lwd=1.5){
	#traitgram plotted using ML ASR from geiger (or ace() from ape if ci=TRUE)
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
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
		#require(geiger)
		tr1<-c(trait,getAncStates(trait,tree));edges<-tree$edge
		plot(1,1,type="n",xlim=c(min(tr1)-0.1,max(tr1)+0.1),ylim=c(max(times),min(times)),
			xlab="Trait Values",ylab="Time  (Before Present)",
			main=paste("Traitgram of",trait.name))
		for(i in 1:nrow(edges)){anc<-edges[i,1];desc<-edges[i,2]
			lines(c(tr1[anc],tr1[desc]),c(times[anc],times[desc]),lwd=lwd)}
		}
	}

simFossilTaxa<-function(p,q,anag.rate=0,prop.bifurc=0,prop.cryptic=0,nruns=1,mintaxa=1,maxtaxa=1000,
	mintime=1,maxtime=1000,minExtant=0,maxExtant=NULL,min.cond=TRUE,count.cryptic=FALSE,print.runs=FALSE,plot=FALSE){
	#simulates taxon evolution as in a fossil record, birth, death and anagenesis as parameters
		#plot argument will produce a diversity curve everytime a new clade is made
		#Time-scale is backwards, as expected for paleo data (root is always expected to be at maxtime1
	#p is the rate of speciation as either budding cladogenesis or bifurcationg speciation
		#HOWEVER by default, w (the rate of anagensis) is zero and u (the proportion of speciation by bifurcation) is zero
		#these need to changed to allow other modes of morphological speciation in the fossil record
	#Simulations will not go longer than maxtime, period
		#when maxtaxa is hit, simulation will go up until FO of maxtaxa+1 taxon to avoid Hartmann et al. effect
		#if minExtant is set, simulation will end once minExtant is hit
			#unless maxExtant is zero, in which case 
		#if maxExtant is set, simulation will end once maxExtant is hit, if before maxtaxa or mintime is hit
			#if maxtaxa or mintime is hit, run is discarded if maxExtant is not satisfied
			#when maxExtant is hit, simulation will go up until FO of maxExtant+1 to avoid Hartmann et al. effect
	#
	#p=0.1;q=0.1;anag.rate=0.1;prop.bifurc=0.1;prop.cryptic=0;nruns=10;mintaxa=1;maxtaxa=200;mintime=1;maxtime=100;minExtant=0;maxExtant=0;plot=TRUE;print.runs=TRUE;min.cond=TRUE
	#
	#p=0.1;q=0.9;anag.rate=0;prop.bifurc=0;prop.cryptic=0;nruns=1;mintaxa=1;maxtaxa=100;mintime=1;maxtime=100;minExtant=0;maxExtant=0;plot=TRUE;print.runs=TRUE;min.cond=TRUE
	#min.cond example
	#set.seed(444);p=0.1;q=0.1;anag.rate=0;prop.bifurc=0;prop.cryptic=0;nruns=10;mintaxa=1;maxtaxa=1000;mintime=1;maxtime=100;minExtant=0;maxExtant=NULL;plot=TRUE;print.runs=TRUE;min.cond=FALSE
	#pure birth example
	#set.seed(444);p=0.1;q=0;anag.rate=0;prop.bifurc=0;prop.cryptic=0;nruns=1;mintaxa=10;maxtaxa=20;mintime=1;maxtime=10;minExtant=0;maxExtant=NULL;plot=TRUE;print.runs=TRUE;min.cond=TRUE
	#cryptic speciation
	#set.seed(444);p=0.1;q=0.1;anag.rate=0.1;prop.bifurc=0.5;prop.cryptic=0.5;nruns=1;mintaxa=10;maxtaxa=20;mintime=1;maxtime=10;minExtant=0;maxExtant=NULL;plot=TRUE;print.runs=TRUE;min.cond=TRUE;count.cryptic=FALSE
	#set.seed(444);p=0.1;q=0.1;anag.rate=0.1;prop.bifurc=0;prop.cryptic=1;nruns=1;mintaxa=10;maxtaxa=20;mintime=1;maxtime=10;minExtant=0;maxExtant=NULL;plot=TRUE;print.runs=TRUE;min.cond=TRUE;count.cryptic=TRUE
	#set.seed(444);p=0.1;q=0.1;anag.rate=0.1;prop.bifurc=0;prop.cryptic=1;nruns=1;mintaxa=10;maxtaxa=20;mintime=1;maxtime=10;minExtant=0;maxExtant=NULL;plot=TRUE;print.runs=TRUE;min.cond=TRUE;count.cryptic=FALSE
	#idiot proofing
	if(any(c(p,q,anag.rate,prop.bifurc,prop.cryptic)<0)){stop(
		"Error: bad parameters input, p, q, anag.rate, prop.bifurc or prop.cryptic are less than 0")}
	if(prop.bifurc>0 & prop.cryptic==1){stop()}
	if(nruns<1){stop("Error: nruns<1")}
	if(maxtaxa<0){stop("Error: maxtaxa<0")}
	if(mintaxa<1){stop("Error: mintaxa<1")}
	if(mintime<1){stop("Error: mintime<1")}
	if(maxtime<mintime){stop("Error: maxtime<mintime")}
	if(mintaxa>maxtaxa){stop("Error: mintaxa > maxtaxa")}
	if(maxtaxa>10000 & maxtime>10000){warning("Warning: Unrealistic limits for maxtaxa or maxtime")}
	if(minExtant<0){stop("Error: minExtant<0")}
	if(minExtant>mintaxa){mintaxa<-minExtant}
	if(!is.null(maxExtant)){
		if(maxExtant<0){stop("Error: maxExtant<0")}
		if(maxExtant>maxtaxa){maxtaxa<-maxExtant}
		if(minExtant>maxExtant){stop("Error: maxExtant is set higher than minExtant")}
		}
	if(!min.cond){message("No conditioning during simulation; run until max limits or total extinction")}
	#end idiot proofing
	#make new min/max extant
	minExtant1<-ifelse(minExtant==0,0,minExtant+1)
	results<-list()
	ntries<-0
	for(i in 1:nruns){
		ntries<-ntries+1
		taxad<-matrix(c(1,NA,0,NA,1),1,);pqw<-p+q+anag.rate
		maxtime1<-maxtime;continue<-TRUE;eval<-FALSE
		while(any(is.na(taxad[,4])) & continue){
			tpot<-is.na(taxad[,4])
			tpot2<-min(taxad[tpot,3])==taxad[,3]
			tpick<-which(tpot & tpot2)[1]
			tpick_FO<-taxad[tpick,3]
			wait<-0
			while(is.na(taxad[tpick,4]) & continue){
				wait<-rexp(1,rate=pqw)+wait
				type<-sample(1:3,1,prob=c(p/pqw,q/pqw,anag.rate/pqw))	#choose the event type: cladogenesis, extinction, anagenesis
				if(type==1){	#IF SPECIATION
					#now need to choose if cryptic, budding or bifurcation!
					type1<-sample(1:3,1,prob=c((1-prop.cryptic)-((1-prop.cryptic)*prop.bifurc),
						(1-prop.cryptic)*prop.bifurc,prop.cryptic))	
					if(type1==1){	#IF BUDDING
						taxad<-rbind(taxad,c(max(taxad[,1])+1,taxad[tpick,1],wait+tpick_FO,NA,max(taxad[,1])+1))
						}
					if(type1==2){	#IF BIFURCATION
						taxad[tpick,4]<-wait+tpick_FO
						taxad<-rbind(taxad,c(max(taxad[,1])+1,taxad[tpick,1],wait+tpick_FO,NA,max(taxad[,1])+1))
						taxad<-rbind(taxad,c(max(taxad[,1])+1,taxad[tpick,1],wait+tpick_FO,NA,max(taxad[,1])+1))
						}
					if(type1==3){	#IF CRYPTIC
						taxad<-rbind(taxad,c(max(taxad[,1])+1,taxad[tpick,1],wait+tpick_FO,NA,taxad[tpick,5]))
						}					
					}
				if(type==2){	#IF EXTINCTION
					taxad[tpick,4]<-wait+tpick_FO}
				if(type==3){	#IF ANAGENESIS
					taxad[tpick,4]<-wait+tpick_FO
					taxad<-rbind(taxad,c(max(taxad[,1])+1,tpick,wait+tpick_FO,NA,max(taxad[,1])+1))
					}
				#these loops ONLY end if maxtime1 is hit, so to kill a run, you need to change maxtime1
					#then you'll need to evaluate it again to make sure it meets criteria
				#count numtax and numext based on count.cryptic
				if(count.cryptic){numtax<-nrow(taxad)}else{numtax<-length(unique(taxad[,5]))}
				if(numtax>maxtaxa){
					maxtime1<-min(c(maxtime1,taxad[maxtaxa+1,3]))
					if(!min.cond){eval<-TRUE}
					}
				#under pure-birth, extinction will never happen: kill off lineage manually
				if(wait>maxtime1){taxad[tpick,4]<-wait}
				#have to kill this if the number of taxa just explodes
				if(sum(taxad[,3]<maxtime1)>(maxtaxa*2)){
					taxad[is.na(taxad[,4]),4]<-maxtime+1
					}
				#simulation MUST always terminate if maxtaxa or maxtime are hit (these are safety limits)
					#maxExtant is a similar soft bound; the real bounds that will determine the output are the mins...
				if(count.cryptic){
					numtax<-nrow(taxad)
					numext<-sum(is.na(taxad[,4]))+sum(taxad[!is.na(taxad[,4]),4]>=maxtime1)	#extant taxa w/cryptic
				}else{
					numtax<-length(unique(taxad[,5]))
					numext<-length(unique(taxad[(is.na(taxad[,4]) | taxad[,4]>=maxtime1),5]))
					}
				#want to end the function if >mintaxa,>mintime,>minExtant1 and >maxExtant
				if(max(taxad[,3:4],na.rm=TRUE)>=mintime & ifelse(numext>0,numtax>mintaxa,numtax>=mintaxa)
					& numext>=minExtant1 & ifelse(is.null(maxExtant),TRUE,numext<=maxExtant) & min.cond){
						#if conditions have been hit,reset maxtime1 to the FAD of the newest living taxa that broke conditions
					if(any(is.na(taxad[,4]))){		#if its dead, don't change maxtime1...
						maxtime2<-min(c(maxtime1,max(taxad[is.na(taxad[,4]) | taxad[,4]>=maxtime1,3])))
						if(maxtime2>mintime){maxtime1<-maxtime2}
						}
					eval<-TRUE
					}
				#are any "live" taxa below maxtime1? if so, continue
				continue<-ifelse(any(is.na(taxad[,4])),any(taxad[is.na(taxad[,4]),3]<=maxtime1),FALSE)
				if(!continue & !min.cond){eval<-TRUE}
				taxad_save<-taxad
				#print(c(nrow(taxad),sum(is.na(taxad[,4]))))
				}
			if(!continue & eval){
				#if continue is false (maxtime1 is hit!), evaluate!
				#don't just use one maxtime1, use a bunch 02-07-12: let's you use more runs!
				taxad<-matrix(taxad[taxad[,3]<maxtime1,],sum(taxad[,3]<maxtime1),)
				if(any(is.na(taxad[,4]))){stop("Error: Live creatures escaping simulation! Get out now while you still have time!")}
				posstimes<-sort(unique(c(taxad[,3:4],maxtime1)))
				maxtimes<-posstimes[posstimes>=mintime & posstimes<=maxtime1]				#make vector of maxtimes
				mtds<-lapply(maxtimes,function(x) matrix(taxad[taxad[,3]<x,],sum(taxad[,3]<x),)) 	#maxtime taxa datasets
				if(count.cryptic){
					numtaxa<-sapply(mtds,function(x) nrow(x))
					numexta<-sapply(1:length(mtds),function(x) sum(mtds[[x]][,4]>=maxtimes[x]))		#number of extant taxa
				}else{	#do not count cryptics
					numtaxa<-sapply(mtds,function(x) length(unique(x[,5])))
					numexta<-sapply(1:length(mtds),function(x) length(unique(mtds[[x]][mtds[[x]][,4]>=maxtimes[x],5])))
					}
				minta<-numtaxa>=mintaxa									#is the clade big enough, per mintaxa?
				maxta<-numtaxa<=maxtaxa									#is the clade small enough, per maxtaxa?
				minti<-maxtimes>=mintime							#is the simulation long enough, per mintime?
				maxti<-maxtimes<=maxtime							#is the simulation short enough, per maxtime?
				maxext<-if(!is.null(maxExtant)){maxExtant>=numexta}else{TRUE}	#is the number of extant taxa <= max?
				minext<-minExtant<=numexta						#is there the right number of extant taxa?
				evalcond<-maxext & minext & minta & maxta & minti & maxti	#evaluate conditions				
				#numext<-sum(is.na(taxad[,4]))+sum(taxad[!is.na(taxad[,4]),4]>=maxtime1)	
				if(any(evalcond)){
					chosen<-rev(which(evalcond))[1]	#choose the last time eval is good to avoid Gerehardt effect
					taxad<-mtds[[chosen]]
					maxtime1<-maxtimes[chosen]
				}else{eval<-FALSE}
				}
			if(!continue & !eval){
				#reset if the clade is done (continue=FALSE) but eval is FALSE (didn't hit conditions)
				taxad<-matrix(c(1,NA,0,NA,1),1,)
				ntries<-ntries+1
				continue<-TRUE;eval<-FALSE;evalcond<-NULL
				maxtime1<-maxtime
				}
			}
		taxad1<-cbind(taxad[,1:4,drop=FALSE],(taxad[,4]>=maxtime1),taxad[,5])	#make extinct/extant column
		#reorder time so that time is backwards
		taxad1[,3:4]<-maxtime1-taxad1[,3:4]
		taxad1[,3]<-round(taxad1[,3],digits=4)
		taxad1[,4]<-round(taxad1[,4],digits=4)
		taxad1[taxad1[,4]<0,4]<-0
		#change order so that taxa ids are sequential	
		taxad2<-cbind(matrix(1:nrow(taxad1),,1),matrix(match(taxad1[,2],taxad1[,1]),,1),
			matrix(taxad1[,3:5],,3),matrix(match(taxad1[,6],taxad1[,1]),,1))	
		#give column names to taxad2
		colnames(taxad2)<-c("taxon.id","ancestor.id","orig.time","ext.time","still.alive","looks.like")
		#give rownames
		names<-paste("t",taxad2[,1],sep="")
		if(any(taxad2[,6]!=taxad2[,1])){
			for(cry in which(sapply(taxad2[,6],function(x) sum(x==taxad2[,6])>1))){
				#name all cryptic taxa special
				names[cry]<-paste("t",taxad2[cry,6],".",sum(taxad2[1:cry,6]==taxad2[cry,6]),sep="")
			}}
		rownames(taxad2)<-names
		results[[i]]<-taxad2
		if(plot){
			taxicDivCont(results[[i]],int.length=0.2)
			if(nruns>1){title(paste("Run #",i," of ",nruns,sep=""))}
			}
		}
	if(print.runs){message(paste(nruns," runs accepted from ",ntries," total runs (",signif(nruns/ntries,2)," Acceptance Probability)",sep=""))}
	if(nruns==1){results<-results[[1]]}
	return(results)
	}

taxicDivCont<-function(timeData,int.length=1,int.times=NULL,plot=TRUE,plotLogRich=FALSE,timelims=NULL,drop.cryptic=FALSE){
	#This function estimates diversity for bins from continuous-time range data
	#input is a per-species matrix of backwards-time FADs and LADs in 2 columns (FADs first)
		#assumes time is in millions of years
	#time interval starts and ends can be pre-input as a 2 column matrix
		#int.length is ignored in this case
	#output (if TRUE) is matrix of bin-start, bit-end, div
	tblen<-int.length
	if(ncol(timeData)==6){	#also allow it to accept taxad objects
		if(!drop.cryptic){
			timeData<-timeData[,3:4,drop=FALSE]
		}else{
			timeDataF<-sapply(unique(timeData[,6]),function(x) max(timeData[x==timeData[,6],3]))
			timeDataL<-sapply(unique(timeData[,6]),function(x) min(timeData[x==timeData[,6],4]))
			timeData<-cbind(timeDataF,timeDataL)
			}
		}	
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	FAD<-as.numeric(timeData[,1]);LAD<-as.numeric(timeData[,2])
	if(is.null(int.times)){
		midtimes<-seq(max(FAD)+2*tblen,min(LAD)-2*tblen,by=-tblen)
		midtimes<-midtimes[midtimes>0]
		int.start<-midtimes+(tblen/2)
		int.end<-midtimes-(tblen/2)
	}else{
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	div<-sapply(1:length(midtimes),function(x) sum(FAD>=int.end[x])-sum(LAD>int.start[x]))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
				xlab="Time (Before Present)",ylab="taxic Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
				xlab="Time (Before Present)",ylab="Taxic Richness")
			}
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}

taxicDivDisc<-function(timeList,int.times=NULL,plot=TRUE,plotLogRich=FALSE,timelims=NULL,split.int=TRUE){
	#this function estimates diversity for binned intervals from discrete interval range data
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
	#time interval starts and ends can be pre-input as a 2 column matrix
		#HOWEVER this could be pretty misleading!
		#standing richness may never be high as the apparent richness of some bins
	#output (if TRUE) is 3 col matrix of int-start, int-end, div
	if(class(timeList[[1]])!="matrix"){if(class(timeList[[1]])=="data.frame"){timeList[[1]]<-as.matrix(timeList[[1]])
		}else{stop("Error: timeList[[1]] not of matrix or data.frame format")}}
	if(class(timeList[[2]])!="matrix"){if(class(timeList[[2]])=="data.frame"){timeList[[2]]<-as.matrix(timeList[[2]])
		}else{stop("Error: timeList[[2]] not of matrix or data.frame format")}}
	intMat<-timeList[[1]]	#the intervals the DATA is given in
	timeData<-timeList[[2]]
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(intMat,1,diff)>0)){stop("Error: timeList[[1]] not in intervals in time relative to modern")}
	if(any(intMat[,2]<0)){stop("Error: Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeData,1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	Fint<-as.numeric(timeData[,1]);Lint<-as.numeric(timeData[,2])
	FAD<-intMat[Fint,1];LAD<-intMat[Lint,2]
	if(is.null(int.times)){
		avg_dur<-abs(mean(apply(timeList[[1]],1,diff)))
		int.bounds<-unique(c(intMat,max(intMat)+avg_dur,min(intMat)-avg_dur))		#add a little space at start and end
		int.bounds<-int.bounds[order(-int.bounds)]
		intMat<-cbind(int.bounds[-length(int.bounds)],int.bounds[-1])
		int.start<-intMat[,1];int.end<-intMat[,2]
		midtimes<-apply(intMat,1,mean)
	}else{
		if(split.int){	#if split.int, then any interval times given are split at discrete time intervals
			splinters<-sort(unique(c(intMat)))
			mustSplit<-apply(int.times,1,function(x) any(sapply(splinters,function(y) x[1]>y & x[2]<y)))
			if(any(mustSplit)){
				for(i in which(mustSplit)){
					splitter<-splinters[sapply(splinters,function(y) int.times[i,1]>y & int.times[i,2]<y)]
					for(j in splitter){		#in case there is more than one splitter
						int.times<-rbind(int.times,c(int.times[i,1],j),c(j,int.times[i,2]))
						}
					}
				int.times<-int.times[-which(mustSplit),]
				int.times<-int.times[order(-int.times[,1]),]
				}
			}
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	div<-sapply(1:length(midtimes),function(x) sum(FAD>int.end[x])-sum(LAD>=int.start[x]))
	#div<-sapply(min(timeData):max(timeData),function(x) 	sum(FAD<=x & LAD>=x))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
				xlab="Time (Before Present)",ylab="Taxic Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
				xlab="Time (Before Present)",ylab="Taxic Richness")
			}
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}

phyloDiv<-function(tree,int.length=1,int.times=NULL,plot=TRUE,plotLogRich=FALSE,drop.ZLB=TRUE,timelims=NULL){
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
	require(ape)
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
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
				xlab="Time (Before Present)",ylab="Lineage Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
				ylim=c(0,max(div1)+1),
				xlab="Time (Before Present)",ylab="Lineage Richness")
			}
		par(parOrig);layout(1)
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}

multiDiv<-function(data,int.length=1,plot=TRUE,split.int=TRUE,drop.ZLB=TRUE,drop.cryptic=FALSE,plotLogRich=FALSE,timelims=NULL){
	#lines up a bunch of taxic or phylo objects and calculates diversity curves simulataneously
		#across all their objects; intuits the type of object without being told
		#it also calculates a "average" median curve and 95% quantile intervals
	#input is a list of dicrete interval or continuous time taxic data or a timetree
		#as in the respective functions
	#output is a list with third objects
		#the first object is a 2-column matrix with interval starts and ends
		#the second object is a matrix 
			#with the measured diversity for all the objects as columns, intervals as rows
	#3rd object consists of a 3col matrix of information related to median curve
		#first column is a per-interval median of the combined diversity curves
		#second and third columns are 95% quantile intervals on that median
	require(ape)
	tblen<-int.length
	dclass<-sapply(data,class)	#data classes
	dclass1<-numeric(length(dclass));dclass1[dclass=="matrix"]<-1;
		dclass1[dclass=="list"]<-2;dclass1[dclass=="phylo"]<-3
	if(any(dclass1==0)){stop("Error: Data of Unknown Type")}
	#get max and min times for each type
	if(any(dclass1==1)){
		lims1<-sapply(data[dclass1==1],function(x) 
			if(ncol(x)==6){
				c(min(x[,3:4],na.rm=T),max(x[,3:4],na.rm=T))	
			}else{
				c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))}
			)
	}else{lims1<-NA}
	if(any(dclass1==2)){
		lims2<-sapply(data[dclass1==2],function(x) 
			c(min(x[[1]][max(x[[2]]),]),max(x[[1]][min(x[[2]]),])))
	}else{lims2<-NA}
	if(any(dclass1==3)){
		lims3<-numeric()
		for(i in which(dclass1==3)){
			ttree<-data[[i]]
			if(is.null(ttree$root.time)){
				ntime<-dist.nodes(ttree)[,Ntip(ttree)+1]
				ntime<-max(ntime)-ntime
			}else{
				ntime<-dist.nodes(ttree)[,Ntip(ttree)+1]
				ntime<-ttree$root.time-ntime
				}
			lims3<-c(lims3,c(min(ntime),max(ntime)))
			}
	}else{lims3<-NA}
	end<-min(c(lims1,lims2,lims3),na.rm=TRUE)
	start<-max(c(lims1,lims2,lims3),na.rm=TRUE)
	midtimes<-seq(start+2*tblen,end-2*tblen,by=-tblen)
	midtimes<-midtimes[midtimes>0]
	int.start<-midtimes+(tblen/2);int.end<-midtimes-(tblen/2)
	int.times<-cbind(int.start,int.end)
	if(split.int & any(dclass1==2)){
		#for every single discrete time dataset, bins must be split at their boundaries
		splinters<-sapply(data[dclasss=2],function(x) x[[1]][,1])
		mustSplit<-apply(int.times,1,function(x) any(sapply(splinters,function(y) x[1]>y & x[2]<y)))
		if(any(mustSplit)){
				for(i in which(mustSplit)){
					splitter<-splinters[sapply(splinters[,1],function(y) int.times[i,1]>y & int.times[i,2]<y),1]
					int.times<-rbind(int.times,c(int.times[i,1],splitter),c(splitter,int.times[i,2]))
				}
			int.times<-int.times[-which(mustSplit),]
			int.times<-int.times[order(-int.times[,1]),]
			}
		midtimes<-(int.start+int.end)/2		
		}
	div<-matrix(,nrow(int.times),1)
	if(any(dclass1==1)){
		divs1<-sapply(data[dclass1==1],function(x) 
			taxicDivCont(timeData=x,int.times=int.times,plot=FALSE,drop.cryptic=drop.cryptic)[,3])
		div<-cbind(div,divs1)
		}
	if(any(dclass1==2)){
		divs2<-sapply(data[dclass1==2],function(x) 
			taxicDivDisc(timeList=x,int.times=int.times,split.int=F,plot=FALSE)[,3])
		div<-cbind(div,divs2)
		}
	if(any(dclass1==3)){
		divs3<-sapply(data[dclass1==3],function(x) 
			phyloDiv(tree=x,int.times=int.times,plot=FALSE,drop.ZLB=drop.ZLB)[,3])
		div<-cbind(div,divs3)
		}
	div<-div[,-1]
	colnames(div)<-paste("dataset",1:ncol(div),sep="") 
	#get median curve
	median<-apply(div,1,median)
	q1<-apply(div,1,quantile,probs=0.025)	#the low quantile
	q2<-apply(div,1,quantile,probs=0.975)	#the high quantile
	median.curve<-cbind(median=median,low.95.quantile=q1,high.95.quantile=q2)
	res<-list(int.times=int.times,div.counts=div,median.curve=median.curve)
	if(plot){plotMultiDiv(res,plotLogRich=plotLogRich,timelims=timelims)}
	return(invisible(res))
	}

plotMultiDiv<-function(results,plotLogRich=FALSE,timelims=NULL){
	#plots the median diversity curve for a multiDiv() result
	int.start<-results[[1]][,1]
	int.end<-results[[1]][,2]
	mdiv<-results[[3]]
	times1<-c(int.start,(int.end+((int.start-int.end)/100)))
	mdiv1<-rbind(mdiv,mdiv)[order(times1),]
	times1<-sort(times1)
	if(plotLogRich){
		mdiv1[mdiv1[,2]<1,2]<-1;mdiv1[mdiv1[,3]<1,3]<-1
		y_lim<-c(min(mdiv1[mdiv1>=1]),max(mdiv1[mdiv1>=1]))
		plot(times1[mdiv1[,3]>0],mdiv1[mdiv1[,3]>0,3],type="n",ylim=y_lim,log="y",
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
			xlab="Time (Before Present)",ylab="Log Lineage/Taxic Richness",
			main=paste("Median Diversity Curve"))
	}else{
		y_lim<-c(min(mdiv1),max(mdiv1))
		plot(times1,mdiv1[,3],type="n",ylim=y_lim,
				xlim=if(is.null(timelims)){c(max(times1),min(times1))}else{timelims},
			xlab="Time (Before Present)",ylab="Lineage/Taxic Richness",
			main=paste("Median Diversity Curve"))
		}
	polygon(c(times1,rev(times1)),c(mdiv1[,3],rev(mdiv1[,2])),col="gray",border=NA)
	lines(times1,mdiv1[,1],lwd=3)
	}

simPaleoTrees<-function(p,q,r,ntrees=1,all.extinct=FALSE,modern.samp.prob=1.0,mintime=1,maxtime=100,
	mintaxa=2,maxtaxa=500,drop.zlb=TRUE,print.runs=FALSE,plot=FALSE){
	#this is a wrapper which will create many paleo trees with at least two observed tips
		#uses simFossilTaxa, sampRanges,taxa2phylo, etc
		#good if you want to simulate many many trees with extinct taxa
		#divergence times for nodes will be perfectly known
	#by default: (minimal conditioning)
		#no conditioning on the number of extant taxa
		#living taxa are sampled perfectly at the present
		#zero-length branches are dropped
	#simPaleoTrees(p=0.1,q=0.1,r=0.1,ntrees=10)
	require(ape)
	if(mintaxa<2){stop("Error: Need at least two taxa per tree; increase mintaxa")}
	if(ntrees<1){stop("Error: number of trees to simulate is <1")}
	res<-rmtree(ntrees,2)
	ntries<-0
	for(i in 1:ntrees){
		rerun<-TRUE
		while(rerun){
			ntries<-ntries+1
			taxa<-suppressMessages(simFossilTaxa(p=p,q=q,anag.rate=0,prop.bifurc=0,prop.cryptic=0,nruns=1,mintaxa=mintaxa,
				maxtaxa=maxtaxa,maxtime=maxtime,maxExtant=ifelse(all.extinct,0,maxtaxa),min.cond=FALSE,plot=plot))
			ranges<-sampleRanges(taxa,r,min.taxa=0,modern.samp.prob=modern.samp.prob)
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

simFossilTaxa_SRCond<-function(r,avgtaxa,p,q,anag.rate=0,prop.bifurc=0,prop.cryptic=0,nruns=1,maxtime=1000,
	maxExtant=NULL,count.cryptic=FALSE,print.runs=FALSE,plot=FALSE){
	#wrapper for simulating clades large enough for 
		#getting some avg number of taxa under a given sampling parameter
	#for sampling: proper conditioning on number of taxa
		#how many taxa are needed in true clade to sample on average X taxa??
		#avgtaxa= average number of taxa you want to recover
		#mintaxa and maxtaxa will be set within 20% +/- of values necc for given avgtaxa
			#These should work well if avgtaxa is large (~50 or so)
	#simFossilTaxa_SRCond(r=0.1,p=0.1,q=0.1,nruns=10,avgtaxa=50,maxExtant=0)
	N<-avgtaxa/(1-exp(-r/(q+anag.rate+(prop.bifurc*p))))
	results<-simFossilTaxa(p=p,q=q,anag.rate=anag.rate,prop.bifurc=prop.bifurc,prop.cryptic=prop.cryptic,nruns=nruns,mintaxa=N,
			maxtaxa=2*N,maxtime=maxtime,maxExtant=maxExtant,count.cryptic=FALSE,print.runs=print.runs,plot=plot)
	#if(nruns==1){results<-results[[1]]}
	return(results)
	}

cladogeneticTraitCont<-function(taxa,rate=1,meanChange=0,rootTrait=0){
	#simulate speciational trait evolution for datasets from simFossilTaxa
	#idiot proofing
	taxa<-taxa[order(taxa[,1]),]
	if(any(taxa[-1,2]>=taxa[-1,1])){stop("Ancestors have higher IDs than Descendants?")}
	anctaxa<-sapply(taxa[-1,2],function(x) which(x==taxa[,1]))
	traits<-rootTrait
	for(i in 2:nrow(taxa)){
		if(taxa[i,1]==taxa[i,6]){
			traits[i]<-traits[anctaxa[i-1]]+rnorm(1,mean=meanChange,sd=sqrt(rate))
		}else{
			traits[i]<-traits[anctaxa[i-1]]
			}
		}
	names(traits)<-rownames(taxa)
	return(traits)
	}

compareNodeAges<-function(tree1,tree2){
	#output vector of shifts in node dates
	require(ape)
	if(class(tree1)!="phylo"){stop("Error: tree1 is not of class phylo")}
	if(class(tree2)!="phylo"){stop("Error: tree2 is not of class phylo")}
	matches1<-which(!is.na(match(tree1$tip.label,tree2$tip.label)))[1]
	tipmatch<-tree1$tip.label[matches1]
	if(length(matches1)<1){stop("Error: No shared taxa!")}
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
		if(min(ntime1)<0){stop("Error: tree1$root.time is less than total depth of tree1!")}
	}else{
		ntime1<-max(ntime1)-ntime1
		}
	if(!is.null(tree2$root.time)){
		tree2$root.time<-tree2$root.time-(mtimeB-mtime2)
		ntime2<-tree2$root.time-ntime2
		if(min(ntime2)<0){stop("Error: tree2$root.time is less than total depth of tree2!")}
	}else{
		ntime2<-max(ntime2)-ntime2
		}
	matches<-match(lapply(prop.part(tree1),function(x) sort(tree1$tip.label[x])),
		lapply(prop.part(tree2),function(x) sort(tree2$tip.label[x])))
	ages1<-ntime1[Ntip(tree1)+which(!is.na(matches))]
	ages2<-ntime2[Ntip(tree2)+matches[!is.na(matches)]]
	age_diff<-ages1-ages2
	names(age_diff)<-NULL
	return(age_diff)
	}

compareTermBranches<-function(tree1,tree2){
	#output vector of shifts in terminal branch lengths
	require(ape)
	if(class(tree1)!="phylo"){stop("Error: tree1 is not of class phylo")}
	if(class(tree2)!="phylo"){stop("Error: tree2 is not of class phylo")}
	tree1<-drop.tip(tree1,tree1$tip.label[is.na(match(tree1$tip.label,tree2$tip.label))])
	tree2<-drop.tip(tree2,tree2$tip.label[is.na(match(tree2$tip.label,tree1$tip.label))])
	term1<-tree1$edge.length[tree1$edge[,2]<=Ntip(tree1)]
	term1<-term1[order(tree1$edge[tree1$edge[,2]<=Ntip(tree1),2])]
	term2<-tree2$edge.length[tree2$edge[,2]<=Ntip(tree2)]
	term2<-term2[order(tree2$edge[tree2$edge[,2]<=Ntip(tree2),2])]
	term2<-term2[match(tree1$tip.label,tree2$tip.label)]
	term_diff<-term2-term1
	names(term_diff)<-tree1$tip.label
	return(term_diff)
	}

addTermBranchLength<-function(tree,addtime=0.001){
	require(ape)
	if(class(tree)!="phylo"){stop("Error: tree is not of class phylo")}
	tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]<-tree$edge.length[tree$edge[,2]<(Ntip(tree)+1)]+addtime
	if(any(tree$edge.length<0)){stop("Error: tree has negative branch lengths!")}
	if(!is.null(tree$root.time)){tree$root.time<-tree$root.time+addtime}
	return(tree)
	}

