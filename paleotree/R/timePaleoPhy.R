#' Timescaling of Paleo-Phylogenies
#' 
#' Timescales an unscaled cladogram of fossil taxa using information on their
#' temporal ranges, using various methods. Also can resolve polytomies randomly
#' and output samples of randomly-resolved trees.
#' 
#' @details These functions are an attempt to unify and collect previously used and
#' discussed methods for time-scaling phylogenies of fossil taxa.
#' Unfortunately, it is occasional difficult to attribute some time-scaling methods to
#' specific references in the literature.
#' 
#' There are five method types that can be used by timePaleoPhy. Four of these
#' use some value of absolute time, chosen a priori, to time-scale the tree.
#' This is handled by the argument vartime, which is NULL by default and unused
#' for type "basic".
#' 
#' \describe{

#'  \item{"basic"}{This most simple of methods ignores vartime and
#' scales nodes so they are as old as the first appearance of their oldest
#' descendant (Smith, 1994). This method produces many zero-length branches
#' (Hunt and Carrano, 2010).}

#'  \item{"equal"}{The 'equal' method defined by G. Lloyd and used in Brusatte
#' et al. (2008) and Lloyd et al. (2012). Originally usable in code supplied by
#' G. Lloyd, it is recreated here. This method works by increasing the time of
#' the root divergence by some amount (here set by vartime) and then adjusting
#' zero-length branches so that time on early branches is re-apportioned out
#' along those later branches equally.}

#' \item{"aba"}{All branches additive. This method takes the "basic" tree and
#' adds vartime to all branches.}

#' \item{"zlba"}{Zero-length branches additive. This method adds vartime to all
#' zero-length branches in the "basic" tree. Discussed (possibly?) by Hunt and Carrano,
#' 2010.} 

#' \item{"mbl"}{Minimum branch length. Scales all branches so they are
#' greater than or equal to vartime, and subtract time added to later branches
#' from earlier branches in order to maintain the temporal structure of events.
#' A version of this was first introduced by Laurin (2004).} }
#' 
#' These functions cannot time-scale branches relative to reconstructed
#' character changes along branches, as used by Lloyd et al. (2012).
#' 
#' These functions will intuitively drop taxa from the tree with NA for range
#' or are missing from timeData or timeList. Taxa dropped from the tree will be
#' will be listed in a 
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' timePaleoPhy is designed for direct application to datasets where taxon first 
#' and last appearances are precisely known in continuous time, with no stratigraphic
#' uncertainty. This is an uncommon form of data to have from the fossil record, 
#' although not an impossible form (micropaleontologists often have very precise 
#' range charts, for example). This means that most users SHOULD NOT use timePaleoPhy directly, 
#' unless they have written their own code to deal with stratigraphic uncertainty. For
#' some groups, the more typical 'first' and 'last' dates represent the minimum
#' and maximum absolute ages for the fossil collections that a taxon is known
#' is known from. Presumably, the first and last appearances of that taxon in
#' the fossil record is at unknown dates within these bounds. These should not
#' be mistaken as the FADs and LADs desired by timePaleoPhy, as timePaleoPhy 
#' will use the earliest dates provided to calibrate node ages, which is either
#' an overly conservative approach to time-scaling or fairly nonsensical.
#'
#' Alternatively to using timePaleoPhy, bin_timePaleoPhy is a wrapper of 
#' timePaleoPhy which produces timescaled trees for datasets which only have 
#' interval data available. For each output tree, taxon first and last appearance 
#' dates are placed within their listed intervals under a uniform distribution. 
#' Thus, a large sample of time-scaled trees will approximate the uncertainty in 
#' the actual timing of the FADs and LADs. 
#'
#' The input timeList object can have overlapping (i.e. non-sequential) intervals,
#' and intervals of uneven size. Taxa alive in the modern should be listed as last 
#' occurring in a time interval that begins at time 0 and ends at time 0. If taxa 
#' occur only in single collections (i.e. their first and last appearance in the 
#' fossil record is synchronous, the argument point.occur will force all taxa
#' to have instantaneous durations in the fossil record. Otherwise, by default,
#' taxa are assumed to first and last appear in the fossil record at different points
#' in time, with some positive duration. The sites matrix can be used to force
#' only a portion of taxa to have simultaneous first and last appearances.
#' 
#' By setting the argument nonstoch.bin to TRUE for bin_timePaleoPhy, the dates are NOT
#' stochastically pulled from uniform bins but instead FADs are assigned to the
#' earliest time of whichever interval they were placed in and LADs are placed
#' at the most recent time in their placed interval. This option may be useful
#' for plotting. The sites argument becomes arbitrary if nonstoch.bin is TRUE.
#' 
#' If timeData or the elements of timeList are actually data.frames (as output
#' by read.csv or read.table), these will be coerced to a matrix.
#' 
#' A tutorial for applying the time-scaling functions in paleotree, along with
#' an example using real (graptolite) data, can be found here:
#' http://nemagraptus.blogspot.com/2013/06/a-tutorial-to-cal3-time-scaling-using.html
#' 
#' @aliases timePaleoPhy bin_timePaleoPhy

#' @param tree An unscaled cladogram of fossil taxa, of class 'phylo'. Tip labels
#' must match the taxon labels in the respective temporal data.

#' @param timeData Two-column matrix of first and last occurrences in absolute
#' continous time, with row names as the taxon IDs used on the tree. This means the
#' first column is very precise FADs (first appearance dates) and the second 
#' column is very precise LADs (last appearance dates), reflect the precise points
#' in time when taxa first and last appear. If there is stratigraphic uncertainty in
#' when taxa appear in the fossil record, please use the 'bin' time-scaling functions. 

#' @param type Type of time-scaling method used. Can be "basic", "equal",
#' "aba", "zbla" or "mbl". Type="basic" by default. See the note below for more
#' details

#' @param vartime Time variable; usage depends on the method 'type' argument.
#' Ignored if type = "basic"

#' @param ntrees Number of time-scaled trees to output. If ntrees is greater
#' than one and both randres and rand.obs are false, the function will fail and
#' a warning is issued, as these arguments would simply produce multiple
#' identical time-scaled trees.

#' @param randres Should polytomies be randomly resolved? By default,
#' timePaleoPhy does not resolve polytomies, instead outputting a time-scaled
#' tree that is only as resolved as the input tree. If randres=T, then
#' polytomies will be randomly resolved using \code{\link{multi2di}} from the
#' package ape. If randres=T and ntrees=1, a warning is printed that users
#' should analyze multiple randomly-resolved trees, rather than a single such
#' tree, although a tree is still output.

#' @param timeres Should polytomies be resolved relative to the order of
#' appearance of lineages? By default, timePaleoPhy does not resolve
#' polytomies, instead outputting a time-scaled tree that is only as resolved
#' as the input tree. If timeres=T, then polytomies will be resolved with
#' respect to time using the paleotree function \code{\link{timeLadderTree}}.
#' See that functions help page for more information; the result of time-order
#' resolving of polytomies generally does not differ across multiple uses,
#' unlike use of multi2di.

#' @param add.term If true, adds terminal ranges.By default, this function will
#' not add the ranges of taxa when time-scaling a tree, so that the tips
#' correspond temporally to the first appearance datums of the given taxa. If
#' add.term=T, then the 'terminal ranges' of the taxa are added to the tips
#' after tree is time-scaled, such that the tips now correspond to the last
#' appearance datums.

#' @param inc.term.adj If true, includes terminal ranges in branch length
#' estimates for the various adjustment of branch lengths under all methods
#' except 'basic' (i.e. a terminal length branch will not be treated as zero
#' length is this argument is TRUE if the taxon at this tip has a non-zero
#' duration). By default, this argument is FALSE and this function will not
#' include the ranges of taxa when adjusting branch lengths, so that
#' zero-length branches before first appearance times will be extended. An
#' error is returned if this argument is true but type="basic" or
#' add.term=FALSE, as this argument is inconsistent with those argument
#' options.

#' @param rand.obs Should the tips represent observation times uniform
#' distributed within taxon ranges? This only impacts the location of tip-dates, 
#' i.e. the 'times of observation' for taxa, and does not impact the dates used to 
#' determine node ages. Thus, this is an alternative to using only the LADs or only the FADs
#' as the per-taxon times of observation. For these functions,rand.obs is TRUE can only impact
#' the result when add.term is FALSE (otherwise the time of observation can only be the FADs) 
#' and so the function fails and a warning is issued. If rand.obs is TRUE, then it is assumed
#' that users wish the tips to represent observations made with some temporal
#' uncertainty, such that they might have come from any point within a taxon's
#' known range.  This might be the case, for example, if a user is interested in
#' applying phylogeny-based approaches to studying trait evolution, but have
#' per-taxon measurements of traits that come from museum specimens with
#' uncertain temporal placement. When rand.obs is TRUE, the tips are placed randomly
#' within taxon ranges, as if uniformly distributed, and thus multiple trees should be 
#' created and analyzed.

#' @param node.mins Minimum ages of nodes on the tree. The minimum dates of
#' nodes can be set using node.mins; this argument takes a vector of the same
#' length as the number of nodes, with dates given in the same order as nodes
#' are they are numbered in the tree$edge matrix (note that in tree$edge, the
#' tips are given first Ntip numbers and these are ignored here). Not all nodes
#' need be set; those without minimum dates can be given as NA in node.mins.

#' @param noisyDrop If TRUE (the default), any taxa dropped from tree due to not
#' having a matching entry in the time data will be listed in a system message.

#' @param plot If true, plots the input and output phylogenies.

#' @param timeList A list composed of two matrices giving interval times and
#' taxon appearance dates. The rownames of the second matrix should be the taxon IDs,
#' identical to the tip.labels for tree. See details.

#' @param nonstoch.bin If true, dates are not stochastically pulled from
#' uniform distributions. See below for more details.

#' @param sites Optional two column matrix, composed of site IDs for taxon FADs
#' and LADs. The sites argument allows users to constrain the placement of
#' dates by restricting multiple fossil taxa whose FADs or LADs are from the
#' same very temporally restricted sites (such as fossil-rich Lagerstatten) to
#' always have the same date, across many iterations of time-scaled trees. To
#' do this, simply give a matrix where the "site" of each FAD and LAD for every
#' taxon is listed, as corresponding to the second matrix in timeList. If no
#' sites matrix is given (the default), then it is assumed all fossil come from
#' different "sites" and there is no shared temporal structure among the
#' events.

#' @param point.occur If true, will automatically produce a 'sites' matrix
#' which forces all FADs and LADs to equal each other. This should be used when
#' all taxa are only known from single 'point occurrences', i.e. each is only
#' recovered from a single bed/horizon, such as a Lagerstatten.

#' @return The output of these functions is a time-scaled tree or set of
#' time-scaled trees, of either class phylo or multiphylo, depending on the
#' argument ntrees. All trees are output with an element $root.time. This is
#' the time of the root on the tree and is important for comparing patterns
#' across trees.
#' 
#' Trees created with bin_timePaleoPhy will output with some additional
#' elements, in particular $ranges.used, a matrix which records the
#' continuous-time ranges generated for time-scaling each tree. (Essentially a
#' pseudo-timeData matrix.)

#' @note Please account for stratigraphic uncertainty in your analysis.
#' Unless you have exceptionally resolved data, use a wrapper with timePaleoPhy,
#' either the provided bin_timePaleoPhy or code a wrapper function of your
#' own that accounts for stratigraphic uncertainty in your dataset. Remember that
#' the FADs (earliest dates) given to timePaleoPhy will *always* be used to 
#' calibrate node ages!

#' @author David W. Bapst, heavily inspired by code supplied by Graeme Lloyd
#' and Gene Hunt.

#' @seealso \code{\link{cal3TimePaleoPhy}}, \code{\link{binTimeData}},
#' \code{\link{multi2di}}

#' @references 
#' Bapst, D. W. 2013. A stochastic rate-calibrated method for time-scaling
#' phylogenies of fossil taxa. \emph{Methods in Ecology and Evolution}.
#' 4(8):724-733.
#' 
#' Brusatte, S. L., M. J. Benton, M. Ruta, and G. T. Lloyd. 2008 Superiority,
#' Competition, and Opportunism in the Evolutionary Radiation of Dinosaurs.
#' \emph{Science} \bold{321}(5895):1485-91488.
#' 
#' Hunt, G., and M. T. Carrano. 2010 Models and methods for analyzing
#' phenotypic evolution in lineages and clades. In J. Alroy, and G. Hunt, eds.
#' Short Course on Quantitative Methods in Paleobiology. Paleontological
#' Society.
#' 
#' Laurin, M. 2004. The Evolution of Body Size, Cope's Rule and the Origin of
#' Amniotes. \emph{Systematic Biology} 53(4):594-622.
#' 
#' Lloyd, G. T., S. C. Wang, and S. L. Brusatte. 2012 Identifying Heterogeneity
#' in Rates of Morphological Evolutio: Discrete Character Change in the
#' Evolution of Lungfish(Sarcopterygii, Dipnoi). \emph{Evolution}
#' \bold{66}(2):330--348.
#' 
#' Smith, A. B. 1994 Systematics and the fossil record: documenting
#' evolutionary patterns. Blackwell Scientific, Oxford.

#' @examples
#' 
#' #Simulate some fossil ranges with simFossilTaxa
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' cladogram <- taxa2cladogram(taxa,plot=TRUE)
#' #Now let's try timePaleoPhy using the continuous range data
#' ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#' #plot diversity curve 
#' phyloDiv(ttree)
#' 
#' #that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#' #let's add those terminal ranges back on with add.term
#' ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#' #plot diversity curve 
#' phyloDiv(ttree)
#' 
#' #that tree didn't look very resolved, does it? (See Wagner and Erwin 1995 to see why)
#' #can randomly resolve trees using the argument randres
#' #each resulting tree will have polytomies randomly resolved in different ways using multi2di
#' ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=1,randres=TRUE,
#'     add.term=TRUE,plot=TRUE)
#' #notice well the warning it prints!
#' #we would need to set ntrees to a large number to get a fair sample of trees
#' 
#' #if we set ntrees>1, timePaleoPhy will make multiple time-trees
#' ttrees <- timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=9,randres=TRUE,
#'     add.term=TRUE,plot=TRUE)
#' #let's compare nine of them at once in a plot
#' layout(matrix(1:9,3,3));parOrig <- par(mar=c(1,1,1,1))
#' for(i in 1:9){plot(ladderize(ttrees[[i]]),show.tip.label=FALSE,no.margin=TRUE)}
#' #they are all a bit different!
#' 
#' #we can also resolve the polytomies in the tree according to time of first appearance
#' 	#via the function timeLadderTree, by setting the argument 'timeres' to TRUE
#' ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=1,timeres=TRUE,
#'     add.term=TRUE,plot=TRUE)
#' 
#' #can plot the median diversity curve with multiDiv
#' layout(1);par(parOrig)
#' multiDiv(ttrees)
#' 
#' #compare different methods of timePaleoPhy
#' layout(matrix(1:6,3,2));parOrig <- par(mar=c(3,2,1,2))
#' plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="basic",vartime=NULL,add.term=TRUE)))
#'     axisPhylo();text(x=50,y=23,"type=basic",adj=c(0,0.5),cex=1.2)
#' plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="equal",vartime=10,add.term=TRUE)))
#'     axisPhylo();text(x=55,y=23,"type=equal",adj=c(0,0.5),cex=1.2)
#' plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="aba",vartime=1,add.term=TRUE)))
#'     axisPhylo();text(x=55,y=23,"type=aba",adj=c(0,0.5),cex=1.2)
#' plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="zlba",vartime=1,add.term=TRUE)))
#'     axisPhylo();text(x=55,y=23,"type=zlba",adj=c(0,0.5),cex=1.2)
#' plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="mbl",vartime=1,add.term=TRUE)))
#'     axisPhylo();text(x=55,y=23,"type=mbl",adj=c(0,0.5),cex=1.2)
#' layout(1);par(parOrig)
#' 
#' #using node.mins
#' #let's say we have (molecular??) evidence that node #5 is at least 1200 time-units ago
#' #to use node.mins, first need to drop any unshared taxa
#' droppers <- cladogram$tip.label[is.na(
#'       match(cladogram$tip.label,names(which(!is.na(rangesCont[,1])))))]
#' cladoDrop <- drop.tip(cladogram, droppers)
#' # now make vector same length as number of nodes
#' nodeDates <- rep(NA, Nnode(cladoDrop))
#' nodeDates[5] <- 1200
#' ttree1 <- timePaleoPhy(cladoDrop,rangesCont,type="basic",
#'   	randres=FALSE,node.mins=nodeDates,plot=TRUE)
#' ttree2 <- timePaleoPhy(cladoDrop,rangesCont,type="basic",
#'    	randres=TRUE,node.mins=nodeDates,plot=TRUE)
#' 
#' #Using bin_timePaleoPhy to timescale with discrete interval data
#' #first let's use binTimeData() to bin in intervals of 1 time unit
#' rangesDisc <- binTimeData(rangesCont,int.length=1)
#' ttreeB1 <- bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,randres=TRUE,
#'     add.term=TRUE,plot=FALSE)
#' #notice the warning it prints!
#' phyloDiv(ttreeB1)
#' #with time-order resolving via timeLadderTree
#' ttreeB2 <- bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,timeres=TRUE,
#'     add.term=TRUE,plot=FALSE)
#' phyloDiv(ttreeB2)
#' #can also force the appearance timings not to be chosen stochastically
#' ttreeB3 <- bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,
#'     nonstoch.bin=TRUE,randres=TRUE,add.term=TRUE,plot=FALSE)
#' phyloDiv(ttreeB3)
#' 
#' \dontrun{
#' #simple three taxon example for testing inc.term.adj
#' ranges1<-cbind(c(3,4,5),c(2,3,1));rownames(ranges1)<-paste("t",1:3,sep="")
#' clado1<-read.tree(file=NA,text="(t1,(t2,t3));")
#' ttree1<-timePaleoPhy(clado1,ranges1,type="mbl",vartime=1)
#' ttree2<-timePaleoPhy(clado1,ranges1,type="mbl",vartime=1,add.term=TRUE)
#' ttree3<-timePaleoPhy(clado1,ranges1,type="mbl",vartime=1,add.term=TRUE,inc.term.adj=TRUE)
#' layout(1:3)
#' ttree1$root.time;plot(ttree1);axisPhylo()
#' ttree2$root.time;plot(ttree2);axisPhylo()
#' ttree3$root.time;plot(ttree3);axisPhylo()
#' -apply(ranges1,1,diff)
#' }
#'
#' @export
timePaleoPhy<-function(tree,timeData,type="basic",vartime=NULL,ntrees=1,randres=FALSE,timeres=FALSE,add.term=FALSE,
	inc.term.adj=FALSE,rand.obs=FALSE,node.mins=NULL,noisyDrop=TRUE,plot=FALSE){
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
	if(!add.term & rand.obs){stop(
		"Error: Inconsistent arguments: randomized observation times are treated as LAST appearance times, so add.term must be true for rand.obs to have any effect on output!"
		)}
	if(ntrees>1 & !randres & !rand.obs){stop("Error: Time-scale more trees without randomly resolving or random obs?!")}
	if(ntrees==1 & randres){message("Warning: Do not interpret a single randomly-resolved tree")}
	if(ntrees==1 & rand.obs){message("Warning: Do not interpret a single tree with randomly-placed obs times")}
	if(randres & timeres){stop(
		"Error: Inconsistent arguments: You cannot randomly resolve polytomies and resolve with respect to time simultaneously!")}
	if(!add.term & inc.term.adj){stop(
		"Error: Inconsistent arguments: Terminal ranges cannot be used in adjustment of branch lengths if not added to tree!")}
	if(type=="basic" & inc.term.adj){stop(
		"Error: Inconsistent arguments: Terminal range adjustment of branch lengths does not affect basic time-scaling method!")}
	originalInputTree<-tree
	#remove taxa that are NA or missing in timeData
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){stop("Error: Absolutely NO valid taxa shared between the tree and temporal data!")}
		if(noisyDrop){message(paste("Warning: Following taxa dropped from tree:",paste0(droppers,collapse=", ")))}
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
		if(!is.binary.tree(savetree)  | !is.rooted(savetree)){
			if(randres){tree<-multi2di(savetree)}
			if(timeres){tree<-timeLadderTree(savetree,timeData)}
			}
		if(rand.obs){timeData[,2]<-apply(saveTD,1,function(x) runif(1,x[2],x[1]))}else{timeData<-saveTD}
		ntime<-sapply(1:Nnode(tree),function(x) 
			max(timeData[tree$tip.label[unlist(prop.part(tree)[x])],1]))	#first, get node times
		ntime<-c(timeData[tree$tip.label,1],ntime)
		if(!is.null(node.mins)){	#if there are node.mins, alter ntime as necessary
			#needs to be same length as nodes in originalInputTree
			if(Nnode(savetree)!=length(node.mins)){
				stop("node.mins must be same length as number of nodes in the input tree!")}
			#of course, node.mins is referring to nodes in unresolved originalInputTree
			#need to figure out which nodes are which now if randres; remake node.mins
			if(length(droppers)>0){
				stop("node.mins not compatible with datasets where some taxa are dropped; drop before analysis instead")}
			if((!is.binary.tree(originalInputTree) | !is.rooted(tree)) & randres){
				node_changes<-match(prop.part(originalInputTree),prop.part(tree))
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

#' @rdname timePaleoPhy
#' @export
bin_timePaleoPhy<-function(tree,timeList,type="basic",vartime=NULL,ntrees=1,nonstoch.bin=FALSE,randres=FALSE,timeres=FALSE,
	sites=NULL,point.occur=FALSE,add.term=FALSE,inc.term.adj=FALSE,rand.obs=FALSE,node.mins=NULL,noisyDrop=TRUE,plot=FALSE){
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
	if(!is(tree, "phylo")){stop("Error: tree is not of class phylo")}
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
		if(noisyDrop){message(paste("Warning: Following taxa dropped from tree:",paste0(droppers,collapse=", ")))}
		tree<-drop.tip(tree,droppers)
		if(Ntip(tree)<2){stop("Error: Less than two valid taxa shared between the tree and temporal data!")}
		timeList[[2]][which(!sapply(rownames(timeList[[2]]),function(x) any(x==tree$tip.label))),1]<-NA
		}
	if(!is.null(node.mins)){
		if(length(droppers)>0){	#then... the tree has changed unpredictably, node.mins unusable
			stop("node.mins not compatible with datasets where some taxa are drop; drop before analysis instead")}
		if(Nnode(tree)!=length(node.mins)){
			stop("node.mins must be same length as number of nodes in the input tree!")}
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
		if(!is.binary.tree(tree) | !is.rooted(tree)){
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