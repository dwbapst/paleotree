#' Simulating Taxa in the Fossil Record
#' 
#' Functions for simulating taxon ranges and relationships under various models
#' of evolution
#' 
#' simFossilTaxa simulates a birth-death process (Kendall, 1948; Nee, 2006),
#' but unlike most functions for this implemented in R, this function enmeshes
#' the simulation of speciation and extinction with explicit models of how
#' lineages are morphologically differentiated, as morphotaxa are the basic
#' units of paleontological estimates of diversity and phylogenetics.
#' 
#' Any particular use of simFossilTaxa will probably involve iteratively
#' running many simulations of diversification. Simulation runs are only
#' accepted for output if and when they meet the conditioning criteria defined
#' in the arguments, both minima and maxima. If min.cond is true (the default),
#' simulations will be stopped and accepted when clades satisfy mintime,
#' mintaxa, minExtant and maxExtant (if the later is set). To reduce the effect
#' of one conditioning criterion, simply set that limit to an arbitrarily low
#' or high number (depending if a minimum or maximum constraint is involved).
#' If min.cond is false, simulation runs are not stopped and evaluated for
#' output acceptance until they (a) go completely extinct or (b) hit either
#' maxtaxa or maxtime. Whether the simulation runs are accepted or not for
#' output is still dependent on mintaxa, mintime, minExtant and maxExtant. Note
#' that some combinations of conditions, such as attempting to condition on a
#' specific non-zero value of minExtant and maxExtant, may take a long time to
#' find any acceptable simulation runs.
#' 
#' Hartmann et al. (2011) recently discovered a potential statistical artifact
#' when branching simulations are conditioned on some maximum number of taxa.
#' Thus, this function continues the simulation once mintaxa or minExtant is
#' hit, until the next taxon (limit +1) originates. Once the simulation
#' terminates, it is judged whether it is acceptable for all conditions given
#' and if so, the run is accepted as a dataset for output.
#' 
#' Please note that mintaxa and maxtaxa refer to the number of static
#' morphotaxa birthed over the entire evolutionary history of the simulated
#' clade, not the extant richness at the end of the simulation. Use minExtant
#' and maxExtant if you want to condition on the number of taxa living at some
#' time.
#' 
#' The simFossilTaxa function can effectively simulate clades evolving under
#' any combination of the three "modes" of speciation generally referred to by
#' paleontologists: budding cladogenesis, branching cladogenesis and anagenesis
#' (Foote, 1996). The first two are "speciation" in the typical sense used by
#' biologists, with the major distinction between these two modes being whether
#' the ancestral taxon shifts morphologically at the time of speciation. The
#' third is where a morphotaxon changes into another morphotaxon with no
#' branching, hence the use of the terms "pseudoextinction" and
#' "pseudospeciation". As bifurcation and budding are both branching events,
#' both are controlled by the p, the instantaneous rate, while the probability
#' of a branching event being either is set by u. By default, only budding
#' cladogenesis occurs To have these three modes occur in equal proportions,
#' set p to be twice the value of w and set u to 0.5.
#' 
#' This function also includes the ability to simulate cryptic cladogenesis.
#' The available patterns of morphological speciation thus form a gradient:
#' cryptic cladogenesis has no morphological shifts in either daughter branches
#' after a branching event, budding cladogenesis has one morphological shift in
#' the two daughter lineages and, in bifurcating cladogenesis, shifts occur in
#' both daughter lineages. The argument prop.cryptic dictates what proportion
#' of branching/cladogenesis events (the overall occurance of which with rate
#' p) are cryptic versus those that have some morphological divergence (either
#' budding of bifurcating. prop.bifurc controls the proportion of
#' morphologically divergent cladogenesis which is bifurcating relative to
#' budding. Thus, for example, the probability of a given cladogenesis event
#' being budding is (1-prop.cryptic)*prop.bifurc.
#' 
#' When there is cryptic speciation, by default, the conditioning arguments
#' involving numbers of taxa (mintaxa, maxtaxa, minExtant and maxExtant) count
#' the number of unique morphologically distinguishable taxa is checked (i.e.
#' the number of unique values in column 6 of the simulated data). This
#' behavior can be changed with the argument count.cryptic.See below about the
#' output data structure to see how information about cryptic cladogenesis is
#' recorded. The functions taxa2phylo, taxa2cladogram and taxicDivCont each
#' handle cryptic species in different ways, as described in their respective
#' help files.
#' 
#' If maxExtant is 0, then the function will be limited to only accepting
#' simulations that end in total clade extinction before maxtime.
#' 
#' If conditions are such that a clade survives to maxtime, then maxtime will
#' become the time of first appearance for the first taxa. Unless maxtime is
#' very low, however, it is more likely the maxtaxa limit will be reached
#' first, in which case the point in time at which maxtaxa is reached will
#' become the present data and the entire length of the simulation will be the
#' time of the first appearance of the first taxon.
#' 
#' simFossilTaxa simulates single taxa until they go extinct or exceed maxtime.
#' This means the function may have fully simulated some lineages for thousands
#' of time-steps while others are not yet simulated, and thus sometimes
#' overshoot constraints on the number of taxa. This function will
#' automatically discard any runs where the number of taxa exceeds 2 x maxtaxa
#' to avoid blowing up computation time. This is likely to happen under a
#' pure-birth scenario; I suggest using low maxtime settings if doing a
#' pure-birth simulation.
#' 
#' simFossilTaxa_SRCond is a wrapper for simFossilTaxa for when clades of a
#' particular size are desired, post-sampling. For more details, see the help
#' file at \code{\link{simFossilTaxa_SRCond}}.
#' 
#' More details on this function's design can be read here:
#' http://nemagraptus.blogspot.com/2012/04/simulating-fossil-record.html

#' @param p Instantaneous rate of origination/branching per lineage-time units.

#' @param q Instantaneous rate of extinction per lineage-time units.

#' @param anag.rate Instantaneous rate of anagenesis (i.e. pseudospeciation/pseudoextinction).

#' @param prop.bifurc Proportion of morphological branching by bifurcating
#' cladogenesis relative to budding cladogenesis.

#' @param prop.cryptic Proportion of branching events with no morphological
#' differentiation (i.e. cryptic speciation) relative to branching events
#' associated with morphological differentiation (budding and/or bifurcating 
#' cladogenesis).

#' @param nruns Number of datasets to accept, save and output.

#' @param mintaxa Minimum number of total taxa over the entire history of a
#' clade necessary for a dataset to be accepted.
#' @param maxtaxa Maximum number of total taxa over the entire history of a
#' clade necessary for a dataset to be accepted.

#' @param mintime Minimum time units to run any given simulation before
#' stopping.
#' @param maxtime Maximum time units to run any given simulation before
#' stopping.

#' @param minExtant Minimum number of living taxa allowed at end of
#' simulations.
#' @param maxExtant Maximum number of living taxa allowed at end of
#' simulations.

#' @param min.cond If TRUE, the default, simulations are stopped when they meet
#' all minimum conditions. If FALSE, simulations will continue until they hit
#' maximum conditions, but are only accepted as long as they still meet all
#' minimum conditions in addition.

#' @param count.cryptic If TRUE, cryptic taxa are counted as separate taxa for
#' conditioning limits, such as maxtaxa or maxExtant. If FALSE, then each cryptic
#' complex (i.e. each distinguishable morphotaxon) is treated as a single taxon

#' @param print.runs If TRUE, prints the proportion of simulations accepted for
#' output to the terminal.

#' @param sortNames If TRUE, output taxonomic matrices are sorted by the taxon
#' names (rownames; thus sorting cryptic taxa together) rather than by taxon id
#' (which is the order they were simulated in).

#' @param plot If TRUE, plots the diversity curves of accepted simulations.

#' @return This function gives back a list containing nruns number of taxa
#' datasets, where each element is a matrix. If nruns=1, the output is not a
#' list but just a single matrix. Sampling has not been simulated in the output
#' for either function; the output represents the 'true' history of the
#' simulated clade.
#' 
#' For each dataset, the output is a six column per-taxon matrix where all
#' entries are numbers, with the first column being the taxon ID, the second
#' being the ancestral taxon ID (the first taxon is NA for ancestor), the third
#' column is the first appearance date of a species in absolute time, the
#' fourth column is the last appearance data and the fifth column records
#' whether a species is still extant at the time the simulation terminated (a
#' value of 1 indicates a taxon is still alive, a value of 0 indicates the
#' taxon is extinct). The sixth column (named "looks.like") gives information
#' about the morphological distinguishability of taxa; if they match the taxon
#' ID, they are not cryptic. If they do not match, then this column identifies
#' which taxon id they would be identified as.
#' 
#' Each matrix of simulated data also has rownames, generally of the form "t1"
#' and "t2", where the number is the taxon id. Cryptic taxa are instead named
#' in the form of "t1.2" and "t5.3", where the first number is the taxon which
#' they are a cryptic descendant of (i.e. column 6 of the matrix,
#' "looks.like"). The second number, after the period, is the rank order of
#' taxa in that cryptic group of taxa. Taxa which are the common ancestor of a
#' cryptic lineage are also given a unique naming convention, of the form
#' "t1.1" and "t5.1", where the first number is the taxon id and the second
#' number communicates that this is the first species in a cryptic lineage.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' @author David W. Bapst
#' @seealso \code{\link{simFossilTaxa_SRCond}}, \code{\link{sampleRanges}},
#' \code{\link{simPaleoTrees}}, \code{\link{taxa2phylo}},
#' \code{\link{taxa2cladogram}}
#' @references Foote, M. 1996 On the Probability of Ancestors in the Fossil
#' Record. \emph{Paleobiology} \bold{22}(2):141--151.
#' 
#' Hartmann, K., D. Wong, and T. Stadler. 2010 Sampling Trees from Evolutionary
#' Models. \emph{Systematic Biology} \bold{59}(4):465--476.
#' 
#' Kendall, D. G. 1948 On the Generalized "Birth-and-Death" Process. \emph{The
#' Annals of Mathematical Statistics} \bold{19}(1):1--15.
#' 
#' Nee, S. 2006 Birth-Death Models in Macroevolution. \emph{Annual Review of
#' Ecology, Evolution, and Systematics} \bold{37}(1):1--17.
#' 
#' Solow, A. R., and W. Smith. 1997 On Fossil Preservation and the
#' Stratigraphic Ranges of Taxa. \emph{Paleobiology} \bold{23}(3):271--277.
#' @examples
#' 
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#' #let's see what the 'true' diversity curve looks like in this case
#' #plot the FADs and LADs with taxicDivCont
#' taxicDivCont(taxa[,3:4])
#' #can also see this by setting plot=TRUE in simFossilTaxa
#' 
#' #make datasets with multiple speciation modes
#' #following has anagenesis, budding cladogenesis and bifurcating cladogenesis
#'     #all set to 1/2 extinction rate
#' set.seed(444)
#' res <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0.05,prop.bifurc=0.5,mintaxa=30,maxtaxa=60,
#'     maxExtant=0,nruns=1,plot=TRUE)
#' #what does this mix of speciation modes look like as a phylogeny?
#' tree <- taxa2phylo(res,plot=TRUE)
#' 
#' \dontrun{
#' #some other options with cryptic speciation
#' taxaCrypt1 <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0,prop.bifurc=0,prop.crypt=0.5,mintaxa=30,
#'     maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
#' tree1 <- taxa2phylo(taxaCrypt1,plot=TRUE)
#' taxaCrypt2 <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0.05,prop.bifurc=0.5,prop.crypt=0.5,
#'     mintaxa=30,maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
#' tree2 <- taxa2phylo(taxaCrypt2,plot=TRUE)
#' taxaCrypt3 <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0.05,prop.bifurc=0,prop.crypt=1,
#'     mintaxa=30,maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
#' tree3 <- taxa2phylo(taxaCrypt2,plot=TRUE)
#' }
#' 
#' set.seed(444)
#' #can choose to condition on total morphologically-distinguishable taxa 
#'     #or total taxa including cryptic taxa with count.cryptic=FALSE
#' taxa<-simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,mintaxa=20,
#'     count.cryptic=FALSE,plot=TRUE)
#' nrow(taxa)                 #number of lineages (inc. cryptic)
#' length(unique(taxa[,6]))            #number of morph-distinguishable taxa
#' #now with count.cryptic=TRUE
#' taxa <- simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,mintaxa=20,
#'     count.cryptic=TRUE,plot=TRUE)
#' nrow(taxa)                  #number of lineages (inc. cryptic)
#' length(unique(taxa[,6]))              #number of morph-distinguishable taxa
#' 
#' #now let's look at extant numbers of taxa
#' #can choose to condition on total morphologically-distinguishable living taxa 
#'     #or total living taxa including cryptic taxa with count.cryptic=FALSE
#' taxa <- simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,minExtant=20,
#'     count.cryptic=FALSE,plot=TRUE)
#'     sum(taxa[,5])             #number of still-living lineages (inc. cryptic)
#'     length(unique(taxa[taxa[,5]==1,6]))	  #number of still-living morph-dist. taxa
#' #now with count.cryptic=TRUE
#' taxa <- simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,minExtant=20,
#'     count.cryptic=TRUE,plot=TRUE)
#'     sum(taxa[,5])              #number of still-living lineages (inc. cryptic)
#'     length(unique(taxa[taxa[,5]==1,6]))	   #number of still-living morph-dist, taxa
#' 
#' #can generate datasets that meet multiple conditions: time, # total taxa, # extant taxa
#' set.seed(444)
#' res <- simFossilTaxa(p=0.1,q=0.1,mintime=10,mintaxa=30,maxtaxa=40,minExtant=10,maxExtant=20,
#'     nruns=20,plot=FALSE,print.runs=TRUE)
#' #use print.run to know how many simulations were accepted of the total generated
#' layout(1:2)
#' #histogram of # taxa over evolutionary history
#' hist(sapply(res,nrow),main="#taxa")
#' #histogram of # extant taxa at end of simulation
#' hist(sapply(res,function(x) sum(x[,5])),main="#extant")
#' 
#' \dontrun{
#' #pure-birth example
#' #note that conditioning is tricky
#' layout(1)
#' taxa <- simFossilTaxa(p=0.1,q=0,mintime=10,mintaxa=100,maxtime=100,maxtaxa=100,
#'     nruns=10,plot=TRUE)
#' 
#' #can generate datasets where simulations go until extinction or max limits
#'      #and THEN are evaluated whether they meet min limits
#'      #good for producing unconditioned birth-death trees
#' set.seed(444)
#' res <- simFossilTaxa(p=0.1,q=0.1,maxtaxa=100,maxtime=100,nruns=10,plot=TRUE,
#'     print.runs=TRUE,min.cond=FALSE)
#' #hey, look, we accepted everything! (That's what we want.)
#' layout(1:2)
#' #histogram of # taxa over evolutionary history
#' hist(sapply(res,nrow),main="#taxa")
#' #histogram of # extant taxa at end of simulation
#' hist(sapply(res,function(x) sum(x[,5])),main="#extant")
#' }
#' 
#' layout(1)
#' 
#' 
#' @export simFossilTaxa
simFossilTaxa<-function(p,q,anag.rate=0,prop.bifurc=0,prop.cryptic=0,nruns=1,mintaxa=1,maxtaxa=1000,
	mintime=1,maxtime=1000,minExtant=0,maxExtant=NULL,min.cond=TRUE,count.cryptic=FALSE,print.runs=FALSE,
	sortNames=FALSE,plot=FALSE){
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
	#set.seed(444);p=0.1;q=0.1;anag.rate=0.1;prop.bifurc=0;prop.cryptic=1;nruns=1;mintaxa=10;maxtaxa=20;mintime=1;maxtime=10;minExtant=0;maxExtant=NULL;plot=TRUE;print.runs=TRUE;min.cond=TRUE;count.cryptic=FALSE
	#idiot proofing
	if(any(c(p,q,anag.rate,prop.bifurc,prop.cryptic)<0)){stop(
		"Error: bad parameters input, p, q, anag.rate, prop.bifurc or prop.cryptic are less than 0")}
	if(prop.bifurc>0 & prop.cryptic==1){stop("Error: Prop.bifurc greater than 0 even though cryptic cladogenesis = 1??")}
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
		taxad<-matrix(c(1,NA,0,NA,1),1,)
		pqw<-p+q+anag.rate
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
				if(length(maxtimes)==0){
					eval<-FALSE
				}else{
					mtds<-lapply(maxtimes,function(x) matrix(taxad[taxad[,3]<x,],sum(taxad[,3]<x),)) 	#maxtime taxa datasets
					if(count.cryptic){
						numtaxa<-sapply(mtds,function(x) nrow(x))
						numexta<-sapply(1:length(mtds),function(x) sum(mtds[[x]][,4]>=maxtimes[x]))		#number of extant taxa
					}else{	#do not count cryptics
						numtaxa<-sapply(mtds,function(x) length(unique(x[,5])))
						numexta<-sapply(1:length(mtds),function(x) length(unique(mtds[[x]][mtds[[x]][,4]>=maxtimes[x],5])))
						}
					minta<-numtaxa>=mintaxa								#is the clade big enough, per mintaxa?
					maxta<-numtaxa<=maxtaxa								#is the clade small enough, per maxtaxa?
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
		if(sortNames){
			taxad2<-taxad2[order(as.numeric(substring(rownames(taxad2),2))),]
			}
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
