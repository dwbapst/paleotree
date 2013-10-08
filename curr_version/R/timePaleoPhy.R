timePaleoPhy<-function(tree,timeData,type="basic",vartime=NULL,ntrees=1,randres=FALSE,timeres=FALSE,add.term=FALSE,
	inc.term.adj=FALSE,rand.obs=FALSE,node.mins=NULL,plot=FALSE){
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
	#tree<-rtree(10);tree$edge.length<-NULL;type="basic";vartime=NULL;add.term=FALSE;node.mins=NULL
	#timeData<-runif(10,30,200);timeData<-cbind(timeData,timeData-runif(10,1,20));rownames(timeData)<-tree$tip.label
	#node.mins<-runif(9,50,300)
	#require(ape)	
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	if(class(tree$tip.label)!="character"){stop("Error: tree tip labels are not a character vector")}
	if(ntrees<1){stop("Error: ntrees<1")}
	if(!add.term & rand.obs){stop("Error: Inconsistent arguments: add.term must be true for rand.obs to have any effect on output!")}
	if(ntrees>1 & !randres & !rand.obs){stop("Error: Time-scale more trees without randomly resolving or random obs?!")}
	if(ntrees==1 & randres){message("Warning: Do not interpret a single randomly-resolved tree")}
	if(ntrees==1 & rand.obs){message("Warning: Do not interpret a single tree with randomly-placed obs times")}
	if(randres & timeres){stop(
		"Error: Inconsistent arguments: You cannot randomly resolve polytomies and resolve with respect to time simultaneously!")}
	if(!add.term & inc.term.adj){stop(
		"Error: Inconsistent arguments: Terminal ranges cannot be used in adjustment of branch lengths if not added to tree!")}
	if(type=="basic" & inc.term.adj){stop(
		"Error: Inconsistent arguments: Terminal range adjustment of branch lengths does not affect basic time-scaling method!")}
	#remove taxa that are NA or missing in timeData
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		tree<-drop.tip(tree,droppers)
		if(is.null(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		timeData[which(!sapply(rownames(timeData),function(x) any(x==tree$tip.label))),1]<-NA
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	ttrees<-rmtree(ntrees,2)
	savetree<-tree			#save tree now so can replace with each loop for multi2di()
	saveTD<-timeData
	for(ntr in 1:ntrees){
		#resolve nodes, if tree is not binary
		tree<-savetree
		if(!is.binary.tree(savetree)){
			if(randres){tree<-multi2di(savetree)}
			if(timeres){tree<-timeLadderTree(savetree,timeData)}
			}
		if(rand.obs){timeData[,2]<-apply(saveTD,1,function(x) runif(1,x[2],x[1]))}else{timeData<-saveTD}
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
		#okay, now time to do add terminal branch lengths if inc.term.adj
		if(add.term & inc.term.adj){
			obs_ranges<-timeData[,1]-timeData[,2]
			term_id<-ttree$tip.label[ttree$edge[ttree$edge[,2]<=Ntip(ttree),2]]
			term_add<-sapply(term_id,function(x) obs_ranges[x])
			ttree$edge.length[ttree$edge[,2]<=Ntip(ttree)]<-ttree$edge.length[ttree$edge[,2]<=Ntip(ttree)]+term_add
			}
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
		#if add.term!=FALSE, then taxon observed ranges are added to the tree, with the LADs becoming the location of the tips
		if(add.term & !inc.term.adj){
			obs_ranges<-timeData[,1]-timeData[,2]
			term_id<-ttree$tip.label[ttree$edge[ttree$edge[,2]<=Ntip(ttree),2]]
			term_add<-sapply(term_id,function(x) obs_ranges[x])
			ttree$edge.length[ttree$edge[,2]<=Ntip(ttree)]<-ttree$edge.length[ttree$edge[,2]<=Ntip(ttree)]+term_add
			}
		#now add root.time: 
		if(add.term){
			#should be time of earliest LAD + distance of root from earliest tip
			ttree$root.time<-max(timeData[ttree$tip.label,2])+min(dist.nodes(ttree)[1:Ntip(ttree),Ntip(ttree)+1])	
		}else{
			#should be time of earliest FAD + distance of root from earliest tip
			ttree$root.time<-max(timeData[ttree$tip.label,1])+min(dist.nodes(ttree)[1:Ntip(ttree),Ntip(ttree)+1])	
			}
		if(plot){
			parOrig<-par(mar=c(2.5,1,1,0.5));layout(1:2)
			plot(ladderize(tree),show.tip.label=TRUE,use.edge.length=FALSE)
			plot(ladderize(ttree),show.tip.label=TRUE);axisPhylo()
			layout(1);par(parOrig)
			}
		names(ttree$edge.length)<-NULL
		ttrees[[ntr]]<-ttree
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}