#' Full-Scale Simulations of the Fossil Record with Birth, Death and Sampling of Morphotaxa
#'
#' 

#' @details

#' @inheritParams

#' @param

#' @return

 #taxon.id ancestor.id orig.time ext.time still.alive looks.like

#' @aliases

#' @seealso

#' @author 
#' David W. Bapst, inspired by code written by Peter Smits.

#' @references

#' @examples
#' 
#' set.seed(2)
#' 
#' #a quick birth-death-sampling simulation with 1 run, 100 taxa
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1,
#' 	nTotalTaxa=50, plot=TRUE)
#' 
#' ###################################################################
#' 
#' # examining multiple runs of simulations
#' 
#' #example of repeated pure birth simulations over 50 time-units
#' records <- simFossilRecord(p=0.1, q=0, r=0, nruns=10,
#' 	totalTime=50, plot=TRUE)
#' #plot multiple diversity curves on a log scale
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE,plotLogRich=TRUE)
#' #histogram of total number of taxa
#' hist(sapply(records,nrow))
#' 
#' #example of repeated birth-death-sampling simulations over 50 time-units
#' records <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=10,
#' 	totalTime=50, plot=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' 
#' #like above, but conditioned instead on having 10 extant taxa
#' 	# between 1 and 100 time-units
#' set.seed(4)
#' records <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=10,
#' 	totalTime=c(1,300), nExtant=10, plot=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' 
#' ################################################
#' 
#' # How probable were the runs I accepted?
#' 	# The effect of conditions
#' 
#' set.seed(1)
#' 
#' # Let's look at an example of a birth-death process
#' 	# with high extinction relative to branching
#' # use default run conditions (barely any conditioning)
#' # use print.runs to look at acceptance probability
#' records <- simFossilRecord(p=0.1, q=0.8, r=0, nruns=10,
#' 	print.runs=TRUE, plot=TRUE)
#' # 10 runs accepted from a total of 10 !
#' 
#' # now let's give much more stringent run conditions
#' 	# require 3 extant taxa at minimum, 5 taxa total minimum
#' records <- simFossilRecord(p=0.1, q=0.8, r=0, nruns=10,
#' 	nExtant=c(3,100), nTotalTaxa=c(5,100),
#' 	print.runs=TRUE, plot=TRUE)
#' # thousands of simulations to just obtail 10 accepable runs!
#' 	# most ended in extinction before minimums were hit
#' 
#' # beware analysis of simulated where acceptance conditions 
#' 	# are too stringent: your data will be a 'special case'
#' 	# of the simulation parameters
#' # it will also take you a long time to generate reasonable
#' 	# numbers of replicates for whatever analysis you are doing
#' 
#' # TLDR: You should look at print.runs=TRUE
#' 
#' ######################
#' 
#' # Using 
#' 
#' # Diversity Dependent Models!
#' 
#' # First, let's try Diversity-Dependent Branching over 50 time-units
#' 
#' # first, let's write the rate equation
#' # We'll use the diversity dependent rate equation model
#' 	# from Ettienne et al. 2012 as an example here
#' 	# Under this equation, p=q at carrying capacity K
#' # Many others are possible!
#' # Note that we don't need to use max(0,rate) as negative rates
#' 	# are converted to zero by default, as controlled by
#' 	# the argument negRatesAsZero
#' 
#' # From Ettiene et al.
#' #	lambda = lambda0 - (lambda0 - mu)*(n/K)
#' # lambda and mu are branching rate and extinction rate
#' 	# lambda and mu == p and q in paleotree (i.e. Foote convention)
#' # lambda0 is the branching rate at richness=0
#' # K is the carrying capacity
#' # n is the richness
#' 
#' # 'N' is the algebra symbol for standing taxonomic richness 
#' 	# for simFossilRecord's simulation capabilities
#' # also branching rate cannot reference extinction rate
#' # we'll have to set lambda0, mu and K in the rate equation directly
#' 
#' lambda0 <- 0.3	# branching rate at 0 richness in Ltu
#' K <- 40		# carrying capacity 
#' mu <- 0.1		# extinction rate will 0.1 Ltu (= 1/3 of lambda0)
#' 
#' # technically, mu here represents the lambda at richness=K
#' 	# i.e. lambdaK
#' # Ettienne et al. are just implicitly saying that the carrying capacity
#' 	# is the richness at which lambda==mu
#' 
#' # construct the equation programmatically using paste0
#' branchingRateEq<-paste0(lambda0,"-(",lambda0,"-",mu,")*(N/",K,")")
#' # and take a look at it...
#' branchingRateEq
#' # its a thing of beauty, folks
#' 
#' # now let's try it
#' records <- simFossilRecord(p=branchingRateEq, q=mu, r=0, nruns=3,
#' 	totalTime=100, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' # those are some happy little diversity plateaus!
#' 
#' 
#' # now let's do diversity-dependent extinction
#' 
#' # let's slightly modify the model from Ettiene et al.
#' #	mu = mu0 + (mu0 - muK)*(n/K)
#' 
#' mu0<-0.001		# mu at n=0
#' muK<-0.1		# mu at n=K (should be equal to lambda at K)
#' K<-40
#' lambda<-muK		# equal to muK
#' 
#' # construct the equation programmatically using paste0
#' extRateEq<-paste0(mu0,"-(",mu0,"-",muK,")*(N/",K,")")
#' extRateEq
#' 
#' # now let's try it
#' records <- simFossilRecord(p=lambda, q=extRateEq, r=0, nruns=3,
#' 	totalTime=100, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' 
#' # these plateaus looks a little more spiky 
#' 	#( maybe there is more turnover at K? )
#' # also, it took a longer for the rapid rise to occur
#' 
#' # Now let's try an example with time-dependent origination
#' 	# and extinction constrained to equal origination
#' 
#' # First, let's define a time-dependent rate equation
#' 	# "T" is the symbol for time passed
#' timeEquation<-"0.4-(0.007*T)"
#' 
#' #in this equation, 0.4 is the rate at time=0
#' 	# and it will decrease by 0.007 with every time-unit
#' 	# at time=50, the final rate will be 0.05
#' # We can easily make it so extinction is always equal to branching rate
#' # "P" is the algebraic equivalent for "branching rate" in simFossilRecord
#' 
#' # now let's try it
#' records <- simFossilRecord(p=timeEquation, q="P", r=0, nruns=3,
#' 	totalTime=50, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' # high variability that seems to then smooth out as turnover decreases
#' 
#' ##########################################################
#' 
#' # Speciation Modes
#' # Some examples of varying the 'speciation modes' in simFossilRecord
#' 
#' # The default is pure budding cladogenesis
#' 	# anag.rate = prop.bifurc = prop.cryptic = 0
#' # let's just set those for the moment anyway
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0, plot=TRUE)
#' 
#' #convert and plot phylogeny
#' 	# note this will not reflect the 'budding' pattern
#' 	# branching events will just appear like bifurcation
#' 	# its a typical convention for phylogeny plotting
#' converted<-fossilRecord2fossilTaxa(record)
#' tree<-taxa2phylo(converted,plot=TRUE)
#' 
#' #now, an example of pure bifurcation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=1, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' 
#' # all the short branches are due to ancestors that terminate
#' 	# via pseudoextinction at bifurcation events
#' 
#' # an example with anagenesis = branching
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=0, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' # lots of pseudoextinction
#' 
#' # an example with anagenesis, pure bifurcation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=1, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' # lots and lots of pseudoextinction
#' 
#' # an example with half cryptic speciation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=c(20,30), nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' 
#' # notice that the tree has many more than the maximum of 30 tips:
#' 	# that's because the cryptic taxa are not counted as
#' 	# separate taxa by default, as controlled by count.cryptic
#' 
#' # an example with anagenesis, bifurcation, cryptic speciation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=0.5, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=c(20,30), nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' # note in this case, 50% of branching is cryptic
#' 	# 25% is bifurcation, 25% is budding
#' 
#' # an example with anagenesis, pure cryptic speciation
#' 	# morphotaxon identity will thus be entirely indep of branching!
#' 	# I wonder if this is what is really going on, sometimes...
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=0, prop.cryptic=1,
#' 	nruns=1, nTotalTaxa=c(20,30), nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' 
#' #############
#' 
#' # playing with count.cryptic with simulations of pure cryptic speciation
#' 
#' #can choose to condition on total morphologically-distinguishable taxa
#'     #or total taxa including cryptic taxa with count.cryptic=FALSE
#' 
#' # an example with pure cryptic speciation with count.cryptic=TRUE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=1,
#' 	nruns=1, totalTime=50, nTotalTaxa=c(10,100), count.cryptic=TRUE)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record))
#' plot(tree);axisPhylo()
#' # notice how the tip labels indicate all are the same morphotaxon
#' 
#' # we'll replace the # of taxa constraints with a time constraint
#' 	# or else the count.cryptic=FALSE simulation will never end!
#' 
#' # an example with pure cryptic speciation with count.cryptic=FALSE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=1,
#' 	nruns=1, totalTime=50, count.cryptic=FALSE)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record))
#' plot(tree);axisPhylo()
#' 
#' #let's look at numbers of taxa returned when varying count.cryptic
#' 	# with prop.cryptic=0.5
#' 
#' #simple simulation going for 50 total taxa	
#' 
#' #first, count.cryptic=FALSE (default)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=50, count.cryptic=FALSE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' nrow(taxa)                 		#number of lineages (inc. cryptic)
#' length(unique(taxa[,6]))            #number of morph-distinguishable taxa
#' 
#' # and count.cryptic=TRUE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=50, count.cryptic=TRUE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' nrow(taxa)                 		#number of lineages (inc. cryptic)
#' length(unique(taxa[,6]))            #number of morph-distinguishable taxa
#' 
#' # okay...
#' # now let's try with 50 extant taxa
#' 
#' #first, count.cryptic=FALSE (default)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nExtant=10, totalTime=c(1,100), count.cryptic=FALSE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' sum(taxa[,5])             		  	#number of still-living lineages (inc. cryptic)
#' length(unique(taxa[taxa[,5]==1,6]))	   	#number of still-living morph-dist. taxa
#' 
#' # and count.cryptic=TRUE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nExtant=10, totalTime=c(1,100), count.cryptic=TRUE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' sum(taxa[,5])             		  	#number of still-living lineages (inc. cryptic)
#' length(unique(taxa[taxa[,5]==1,6]))	   	#number of still-living morph-dist. taxa
#' 
#' 




	
#' \donttest{
#'
#' res<-simFossilRecord(0.1,0.1,0.1,nTotalTaxa=10,nExtant=0,nruns=1000,plot=TRUE)
#' anyLive<-any(sapply(res,function(z) any(sapply(z,function(x) x[[1]][5]==1))))
#' if(anyLive){
#'	stop("Runs have extant taxa under conditioning for none?")
#'	}
#'
#' }

#' @name
#' @rdname
#' @export



#' Hartmann et al. (2011) recently discovered a potential statistical artifact
#' when branching simulations are conditioned on some maximum number of taxa.








#' @param totalTime,nTotalTaxa,nExtant,NSamp These arguments represent stopping and
#' acceptance conditions for simulation runs. They are respectively \code{totalTime}, the
#' total length of the simulation in time-units, \code{nTotalTaxa}, the total number of taxa
#' over evolutionary history in the clade, \code{nExtant}, the total number of extant taxa at
#' the end of the simulation and \code{nSamp} the total number of sampled taxa (not counting extant
#' taxa sampled at the modern day). These are used to determine when to end simulation runs, and whether to accept
#' or reject them as output. They can be input as a vector of two numbers, representing minimum
#' and maximum values of a range for accepted simulation runs (i.e. the simulation length can be between 0 and
#' 1000 time-steps, by default), or as a single number, representing a point condition (i.e. if
#' \code{nSamp = 100} then the only simulation runs with exactly 100 taxa sampled will be output).
#' Note that it is easy to set combinations of parameters and run conditions that are impossible
#' to produce satisfactory input under, in which case \code{simFossilRecord} would run in a nonstop loop.

#' @param negRatesAsZero A logical. Should rates calculated as a negative number cause the simulation to fail
#' with an error message (\code{ = FALSE}) or should these be treated as zero (\code{"= TRUE"}, the default).




#' @param p,q,r,anag.rate These parameters control the instantaneous ('per-capita') rates of branching, extinction,
#' sampling and anagenesis, respectively.
#' Rates set to \code{= Inf} are treated as if 0. When a rate is set to 0, this event type will not occur in the simulation.

#' @export
simFossilRecord<-function(

	# model parameters
	#
	p, q, r, anag.rate=0, prop.bifurc=0, prop.cryptic=0,
	modern.samp.prob=1, startTaxa=1, nruns=1,

	# run conditions can be given as vectors of length 1 or length 2 (= min,max)
	#
	totalTime = c(0, 1000), nTotalTaxa = c(1, 1000),
	nExtant = c(0, 1000), nSamp = c(0, 1000),

	#control parameters
	#
	tolerance=10^-4, shiftRoot4TimeSlice="withExtantOnly",
	count.cryptic=FALSE, negRatesAsZero=TRUE, print.runs=FALSE, sortNames=FALSE, plot=FALSE){

	# NOT USED (but in simFossilTaxa):	min.cond=TRUE

	#example parameter sets
	#
	# DEFAULTS (without rates set)
		# anag.rate=0;prop.bifurc=0;prop.cryptic=0;startTaxa=1;nruns=1;
		# nTotalTaxa=c(1,1000);totalTime=c(1,1000);nSamp=c(0,1000);nExtant=c(0,1000);
		# plot=TRUE;count.cryptic=FALSE;print.runs=TRUE;sortNames=FALSE	
	#
	# BASIC RUN	with diversity-dep extinction
	# p=0.1;q='0.01*N';r=0.1
		# anag.rate=0;prop.bifurc=0;prop.cryptic=0;startTaxa=1;nruns=1;
		# nTotalTaxa=c(10,200);totalTime=c(1,1000);nSamp=c(0,1000);nExtant=c(0,0);plot=TRUE;
		# count.cryptic=FALSE;print.runs=TRUE;sortNames=FALSE;set.seed(444)
	#
	
	

	#p=0.1;q=0.1;r=0.1,anag.rate=0.1;prop.bifurc=0.1;prop.cryptic=0;nruns=10;mintaxa=1;maxtaxa=200;mintime=1;maxtime=100;minExtant=0;maxExtant=0;plot=TRUE;print.runs=TRUE;min.cond=TRUE
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


	
	#functions that will live only in simFossilRecord's namespace
		#so not to crowd paleotree namespace
	
	makeParFunct<-function(par,isBranchRate){
		#things to look for
			# N = number of extant taxa
			# T = time
			# P = branching rate
		acceptedArg<-if(isBranchRate){c('N','T')
			}else{c('N','T','P')}
		if(is.numeric(par)){par<-as.character(par)}
		if(is.numeric(type.convert(par,as.is=TRUE))){
			res<-type.convert(par,as.is=TRUE)
			#convert any 'Inf' rates to 0: this event type cannot occur !
			if(is.infinite(res)){res<-0}
			#return error if any rate is negative
			if(res<0){
				stop("input rates must be at least 0 if input can be coerced to type numeric")}
			#now convert formula expression to function
			parFunct<-if(isBranchRate){
				function(N,T){}
			}else{
				function(N,T,P){}
				}
			body(parFunct)<-res
			res<-parFunct
		}else{
			#first convert par to a formula
			#check arguments match only accepted list
			args<-all.vars(as.formula(paste0('XXdave~',par)))
			args<-args[!(args=='XXdave')]
			if(!all(sapply(args,function(x) any(x==acceptedArg)))){
				if(isBranchRate){
					stop(paste0('Incorrect parameterization of branching rate formula, \n',
						'Only N and T allowed as variables'))
				}else{
					stop(paste0('Incorrect parameterization of a non-branching rate formula, \n',
						'Only P, N and T allowed as variables'))
					}
				}
			#now convert formula expression to function
			parFunct<-if(isBranchRate){
				function(N,T){}
			}else{
				function(N,T,P){}
				}
			body(parFunct)<-parse(text=par)
			res<-parFunct
			}
		return(res)
		}

	#RANDOM EXAMPLES TO TEST makeParFunct with...
	#makeParFunct(0.1,isBranchRate=TRUE)
	#
	#z<-0.1
	#makeParFunct(z^2,isBranchRate=TRUE)
	#
	#makeParFunct('P-0.1*N',isBranchRate=TRUE)
	#
	#makeExtRate<-makeParFunct('P-0.1*N',isBranchRate=FALSE)
	#makeExtRate(P=0.1,N=10,T=100)
	#
	#makeParFunct('0.1+T*0.2-0.1^N',isBranchRate=FALSE)

	# get rate vector

	getRateVector<-function(taxa,timePassed,
			getBranchRate,getExtRate,getSampRate,getAnagRate,
			prop.cryptic,prop.bifurc,negRatesAsZero){
		#
		#get some basic summary statistics first
		#vector of which taxa are still alive
		whichExtant<-whichLive(taxa)
		#standing number of extant lineages 
		# (i.e. number of lineages that stuff can happen to)
		nLive<-length(whichExtant)	
		#
		###########################################################
		# calculate rates (which may be time of diversity dependent)
			#use rate-getting functions from above
		#get the new branching rate, extinction rate, sampling rate, anagenesis rate
		branchRate<-getBranchRate(N=nLive,T=timePassed)
		extRate<-getExtRate(N=nLive,T=timePassed,P=branchRate)
		sampRate<-getSampRate(N=nLive,T=timePassed,P=branchRate)
		anagRate<-getAnagRate(N=nLive,T=timePassed,P=branchRate)
		##
		# now deal with proportional types of branching
		#get cryptic, budding and bifurcation components
		crypticRate<-branchRate*(prop.cryptic)
		#rate of morph differentiation per branching event
		morphRate<-branchRate*(1-prop.cryptic)
		buddRate<-morphRate*(1-prop.bifurc)
		bifurcRate<-morphRate*(prop.bifurc)		
		#
		#get probabilities of event types
		rateVector<-c(buddRate,bifurcRate,anagRate,crypticRate,extRate,sampRate)
		names(rateVector)<-c('budd','bifurc','anag','crypt','ext','samp')
		#check rates, make sure none are less than zero
		if(any(rateVector<0)){
			if(negRatesAsZero){
				rateVector[rateVector<0]<-0
			}else{
				stop(paste0(names(which(rateVector<0)),'rate calculated less than zero'))
				}
			}
		# check that not *all* rates are 0\
		if(!(sum(rateVector)>0)){
			stop("Simulation found scenario in which all rates are zero (?!)")
			}	
		#
		return(rateVector)
		}

	#internal functions for branching/extinction/anagenesis processes	

	initiateTaxa<-function(startTaxa,time){
		newTaxa<-lapply(1:startTaxa,function(x) 
			newTaxon(newID=x,ancID=NA,time=time,looksLike=x)
			)
		return(newTaxa)
		}
		
	newTaxon<-function(newID,ancID,time,looksLike){
		#creates an entirely new just-originated taxon
		#store taxa as a list structure
			# $taxa.data, exactly like output from simFossilTaxa
		taxaData<-c(newID,ancID,time,NA,1,looksLike)
		names(taxaData)<- c('taxon.id','ancestor.id','orig.time','ext.time','still.alive','looks.like')
		# $sampling.times = times of sampling events for this taxon
			#thus can come up with quick/simple ways of evaluating run conditions
			# e.g. evaluate number of sampled taxa by sum(length($sampling.times)>0)
		taxon<-list(taxa.data=taxaData, sampling.times=numeric())
		return(taxon)
		}

	origination<-function(taxa,ancID,time,looksLike=NULL){
		#adds a new taxon via branching or anagenesis
		newID<-max(sapply(taxa,function(x) x[[1]][1]))+1
		if(is.null(looksLike)){
			looksLike<-newID
		}else{
			#looksLike needs to be looksLike of parent
			looksLike<-taxa[[looksLike]][[1]][6]
			}
		newTaxonData<-newTaxon(newID=newID,ancID=ancID,time=time,looksLike=looksLike)	
		newTaxa<-c(taxa,list(newTaxonData))
		return(newTaxa)
		}

	termination<-function(taxa,target,time){
		#ends an existing taxon (from extinction or pseudoextinction)
		whichTarget<-which(sapply(taxa,function(x) x[[1]][1]==target))
		if(length(whichTarget)!=1){stop('taxon IDs repeated??')}
		taxa[[whichTarget]][[1]][4]<-time
		taxa[[whichTarget]][[1]][5]<-0
		return(taxa)
		}

	buddingEvent<-function(taxa,parent,time){
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
		return(taxa)
		}

	crypticEvent<-function(taxa,parent,time){
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=parent)
		return(taxa)
		}

	anagEvent<-function(taxa,parent,time){
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
		taxa<-termination(taxa=taxa,target=parent,time=time)
		return(taxa)
		}	
			
	bifurcEvent<-function(taxa,parent,time){
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
		taxa<-termination(taxa=taxa,target=parent,time=time)
		return(taxa)
		}	
		
	extEvent<-function(taxa,target,time){
		taxa<-termination(taxa=taxa,target=target,time=time)
		return(taxa)		
		}
		
	sampEvent<-function(taxa,target,time){
		whichTarget<-which(sapply(taxa,function(x) x[[1]][1]==target))
		taxa[[whichTarget]][[2]]<-c(taxa[[whichTarget]][[2]],time)
		return(taxa)		
		}
			
	eventOccurs<-function(taxa,target,type,time){
		#possible types : 'budd','bifurc','anag','crypt','ext','samp'
		if(type=="budd"){
			taxa<-buddingEvent(taxa=taxa,parent=target,time=time)
			}
		if(type=="bifurc"){
			taxa<-bifurcEvent(taxa=taxa,parent=target,time=time)		
			}
		if(type=="anag"){
			taxa<-anagEvent(taxa=taxa,parent=target,time=time)
			}
		if(type=="crypt"){
			taxa<-crypticEvent(taxa=taxa,parent=target,time=time)
			}
		if(type=="ext"){
			taxa<-extEvent(taxa=taxa,target=target,time=time)		
			}	
		if(type=="samp"){
			taxa<-sampEvent(taxa=taxa,target=target,time=time)
			}
		#check
		if(any(sapply(taxa,length)!=2)){
			stop("taxa object contains taxa with not 2 elements?")}
		return(taxa)
		}

	# functions for identifying live/sampled taxa and checking simulation runs
	
	whichLive<-function(taxa){
		res<-which(sapply(taxa,function(x) x[[1]][5]==1))
		res2<-which(sapply(taxa,function(x) is.na(x[[1]][4]) | identical(unname(x[[1]][4]),0)))
		if(!identical(unname(res),unname(res2))){
			#browser()
			stop("Disagreement on which taxa are extant")}
		return(res)
		}

	whichSampled<-function(taxa){
		res<-which(sapply(taxa,function(x) length(x[[2]])>0))
		return(res)
		}
	
	getRunVitals<-function(taxa,count.cryptic){
		#NOTE need to change vital measurement dependent on count.cryptic or not
		if(!is.list(taxa)){stop("handed getRunVitals a taxa object that isn't a list??")}
		whichExtant<-whichLive(taxa)
		whichSamp<-whichSampled(taxa)
		if(count.cryptic){
			nTaxa<-length(taxa)	#total number of taxa
			nLive<-length(whichExtant)
			nSampled<-length(whichSamp)	#total number of sampled taxa
		}else{
			looksLike<-sapply(taxa,function(x) x[[1]][6])
			#count number of unique taxa based on looksLike
			nTaxa<-length(unique(looksLike))
			#count number of unique extant taxa
			nLive<-length(unique(looksLike[whichExtant]))
			#count number of unique sampled taxa
			nSampled<-length(unique(looksLike[whichSamp]))
			}
		vitals<-c(nTotalTaxa=nTaxa,nExtant=nLive,nSamp=nSampled)
		return(vitals)
		}
	
	testContinue<-function(vitals,timePassed,runConditions){
		#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
			# none of these can ever REVERSE
		#
		#time passed
		timePassed<-runConditions$totalTime[2]-currentTime			
		#
		# test run conditions
		totalExtinction<-vitals[2]==0
		tooMuchTime<-timePassed>runConditions$totalTime[2]
		tooManyTaxa<-vitals[1]>runConditions$nTotalTaxa[2]
		tooManySamp<-vitals[3]>runConditions$nSamp[2]
		runStop<-totalExtinction | tooMuchTime | tooManyTaxa | tooManySamp			
		continue<-unname(!runStop)
		return(continue)
		}

	contiguousIntegerSeq<-function(vector){
		if(!is.vector(vector)){
			stop("contiguousIntegerSeq handed a non-vector?")}
		if(length(vector)<2){
			stop("contiguousIntegerSeq handed length 1 vector?")}
		#		
		vector<-as.integer(vector)
		#because unbelievably base R has no simple function for
			#pulling contiguous sequences of integers from 
		starts<-sapply(2:length(vector),function(i){
			(vector[i]-vector[i-1])>1
			})
		#browser()
		starts<-c(TRUE,starts)
		starts<-vector[starts]
		ends<-sapply(1:(length(vector)-1),function(i){
			(vector[i+1]-vector[i])>1
			})
		ends<-c(ends,TRUE)
		ends<-vector[ends]
		seqMat<-cbind(starts,ends)
		#
		# checks
		if(!is.matrix(seqMat)){stop("seqMat isn't a matrix")}
		if(ncol(seqMat)!=2){stop("seqMat doesn't have 2 columns")}
		return(seqMat)
		}

	insertRow<-function(table,row,rownum){
		#because unbelievably base R has no simple function for inserting a row
		# insert new immediately at this number, shifts row currently at that location DOWN
		table<-rbind(table[1:(rownum-1),],row,table[-(1:(rownum-1)),])
		return(table)
		}
		
	worthCheckingVitalsRecord<-function(vitalsRecord,runConditions){
		#
		# check that labels for vitalsRecord and runConditions match
		labMatch<-colnames(vitalsRecord)[2:4]==names(runConditions)[2:4]
		if(!all(labMatch)){
			stop("runConditions and vitalsRecord objects are mislabeled/misordered")
			}
		#
		lastVitals<-vitalsRecord[nrow(vitalsRecord),]
		reachMin<-sapply(c(1,2,4),function(i){
			lastVitals[i]>=runConditions[[i]][1]
			})
		worthChecking<-all(reachMin)
		#
		return(worthChecking)
		}

	testVitalsRecord<-function(vitalsRecord,runConditions,tolerance){
		#
		# check that labels for vitalsRecord and runConditions match
		labMatch<-colnames(vitalsRecord)[2:4]==names(runConditions)[2:4]
		if(!all(labMatch)){
			stop("runConditions and vitalsRecord objects are mislabeled/misordered")
			}
		#
		#first INSERT FAKE EVENTS INTO VITALS MAT
			# FOR MIN TIME AND MAX TIME
		#
		# for min time
		if(vitalsRecord[1,1]<runConditions$totalTime[1] 
			& vitalsRecord[nrow(vitalsRecord),1]>runConditions$totalTime[1]
			& all(vitalsRecord[,1]!=runConditions$totalTime[1])){
			#
			#what row to insert at
			whereInsert<-which(vitalsRecord[,1]>runConditions$totalTime[1])[1]
			newRow<-c(runConditions$totalTime[1],vitalsRecord[whereInsert-1,-1])
			#vitalsRecord<-insertRow(table=vitalsRecord,row=newRow,rownum=whereInsert-1)
			vitalsRecord<-rbind(vitalsRecord,newRow)
			}
		vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
		#
		# for max time
		if(vitalsRecord[1,1]<runConditions$totalTime[2] 
			& vitalsRecord[nrow(vitalsRecord),1]>runConditions$totalTime[2]
			& all(vitalsRecord[,1]!=runConditions$totalTime[2])){
			#
			#what row to insert at
			whereInsert<-rev(which(vitalsRecord[,1]<runConditions$totalTime[2]))[1]
			newRow<-c(runConditions$totalTime[2],vitalsRecord[whereInsert,-1])
			#vitalsRecord<-insertRow(table=vitalsRecord,row=newRow,rownum=whereInsert+1)
			vitalsRecord<-rbind(vitalsRecord,newRow)
			}
		vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
		#
		################################################################
		#
		# NOW need to essentially duplicate EVERY ROW with time-stamp of row immediately after it
		newVitalsRecord<-cbind(vitalsRecord[2:nrow(vitalsRecord),1,drop=FALSE],
			vitalsRecord[1:(nrow(vitalsRecord)-1),2:4,drop=FALSE])
		pastIncrement<-diff(vitalsRecord[,1])
		if(any(pastIncrement>0)){
			pastIncrement<-min(pastIncrement[pastIncrement>0])/1000
			pastIncrement<-min(c(tolerance,pastIncrement))
		}else{
			pastIncrement<-tolerance
			}
		newTimes<-newVitalsRecord[,1]-pastIncrement
		newTimes<-ifelse(newTimes>0,newTimes,0)
		newVitalsRecord[,1]<-newTimes
		vitalsRecord<-rbind(vitalsRecord,newVitalsRecord)
		vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
		#
		###########################################################################
		# identify all rows where timePassed, nTaxa, nExtant and nSamp are good
		#
		okayVitalsMat<-sapply(1:4,function(i){
			var<-vitalsRecord[,i]
			varRange<-runConditions[[i]]
				var>=varRange[1] & var<=varRange[2]
				})
		okayVitals<-apply(okayVitalsMat,1,all)				
		#
		#########################################################
		# Now test if there are any, if so, sequence
		#
		if(any(okayVitals)){
			# need to build a matrix of the paired-date sequences
			whichOkay<-which(okayVitals)
			if(length(whichOkay)>1){
				seqVitals<-contiguousIntegerSeq(whichOkay)
			}else{
				seqVitals<-matrix(whichOkay,1,2)
				}
			#replaced with the passedTime dates
			seqVitals<-apply(seqVitals,2,function(x) vitalsRecord[x,1])
			if(is.vector(seqVitals)){
				seqVitals<-matrix(seqVitals,,2)
				}
			#
			# checks
			if(!is.matrix(seqVitals)){stop("seqVitals isn't a matrix")}
			if(ncol(seqVitals)!=2){stop("seqVitals doesn't have 2 columns")}
		}else{
			seqVitals<-NA
			}
		#
		#check to make sure it makes sense
		if(length(seqVitals)==0){stop("seqVitals constructed incorrectly")}
		#
		return(seqVitals)
		}

	sampleSeqVitals<-function(seqVitals){
		cumSumSeq<-cumsum(apply(seqVitals,1,diff))
		totalSum<-rev(cumSumSeq)[1]
		if(totalSum>0){
			placedDate<-runif(1)*totalSum
			findRow<-which(cumSumSeq>=placedDate)[1]
			earlierRowCumSum<-ifelse(findRow==1,0,cumSumSeq[findRow-1])
			date<-seqVitals[findRow,1]+placedDate-earlierRowCumSum
		}else{
			# no probability density to sample
			# randomly pick a row
			if(length(seqVitals[,1])>1){
				date<-sample(seqVitals[,1],1)
			}else{
				date<-seqVitals[,1]
				}
			}					
		return(date)
		}

	getTaxaNames<-function(taxa){	
		#name each normal taxon as t + ID 
			#cryptic taxa are cryptic id + . taxon number within that complex
		cryptIDs<-sapply(taxa,function(x) x[[1]][6])
		taxonIDs<-sapply(taxa,function(x) x[[1]][1])
		whichCrypt<-sapply(cryptIDs,function(x) sum(x==cryptIDs)>1)
		newIDs<-sapply(1:length(taxonIDs),function(x){
			if(whichCrypt[x]){
				#find what number this taxon is, within the cryptic complex
				matchCrypt<-cryptIDs==cryptIDs[x]
				nCrypt<-which(taxonIDs[matchCrypt]==taxonIDs[x])-1
				res<-paste0(cryptIDs[x],".",nCrypt)
			}else{
				res<-taxonIDs[x]
				}
			return(res)
			})
		newIDs<-paste0("t",newIDs)
		return(newIDs)
		}

	testFinal<-function(taxa,timePassed,runConditions,count.cryptic){
		# test that the produced taxa object actually passed the runConditions
		finalVitals<-getRunVitals(taxa=taxa,count.cryptic=count.cryptic)
		finalVitals<-c(timePassed=timePassed,finalVitals)
			#time
			#okayTime<-(timePassed>=runConditions[[1]][1] & timePassed<=runConditions[[1]][2])
		#other vitals
		okayVitals<-sapply(1:4,function(i){
			var<-finalVitals[i]
			varRange<-runConditions[[i]]
			var>=varRange[1] & var<=varRange[2]
			})
		#finalCheck<-all(finalVitals)
		if(any(!okayVitals)){
			print(finalVitals)
			#browser()
			stop(paste0("Accepted run as outside of bounds set for conditions:",
				names(runConditions)[!okayVitals],collapse=", "))
			}
		finalCheck<-all(okayVitals)
		return(finalVitals)
		}		


	##################################################################################
	#
	#simplified birth-death-sampling simulator for fossil record
	#
	#p is the rate of branching
		#may be either budding (prob = 1-prop.bifurc) or bifurcation (prob = prop.bifurc)
	#anagenesis is a separate rate (anag.rate)
	#q is rate of extinction
	#r is rate of sampling
	#
	#ARGUMENT CHECKING
	#
	# number of starting taxa and runs must be at least 1
	if(nruns<1){
		stop("nruns must be at least 1")}
	if(startTaxa<1){
		stop("startTaxa must be at least 1")}
	#nruns, starting taxa must be integer values
	if(!all(sapply(c(nruns,startTaxa),function(x) x==round(x)))){
			stop("nruns and startTaxa must coercible to whole number integers")}
	# check that prop.bifurc, prop.cryptic, modern.samp.prob are greater than 0 and less than 1
	if(any(c(prop.bifurc, prop.cryptic, modern.samp.prob)<0) |
			any(c(prop.bifurc, prop.cryptic, modern.samp.prob)>1)){
		stop("bad parameters input: prop.bifurc, prop.cryptic and modern.samp.prob must be between 0 and 1")}	
	# is prop.bifurc and prop.cryptic consistent?
	if(prop.bifurc>0 & prop.cryptic==1){
		stop("Prop.bifurc greater than 0 when probability of branching being cryptic is 1")}
 	#check that min nSamp isn't higher that 0, if r = 0 or Inf
	if((r==0 | is.infinite(r)) & nSamp[1]>0){
		stop("Minimum number of required sampled taxa is >0 but sampling rate is zero (or infinite)")}
	# check count.cryptic
	if(!count.cryptic){ #if false
		#check that min nTotalTaxa, nExtant, nSamp isn't higher that 1, if prop.cryptic=1
		if((prop.cryptic==1 & anag.rate==0) & ( nTotalTaxa[1]>1 | nExtant[1]>1 | nSamp[1]>1 )){
			stop(paste0("Minimum number of required  nTotalTaxa, nExtant and/or nSamp is >1 but these",
				" constraints cannot be reached as count.cryptic=FALSE, prop.cryptic=1 and anag.rate=0)"))
			}
		}
	#check that count.cryptic,negRatesAsZero,print.runs,sortNames,plot are all logicals
	if(!all(sapply(c(count.cryptic,negRatesAsZero,print.runs,sortNames,plot),is.logical))){
		stop("count.cryptic, negRatesAsZero, print.runs, sortNames, and plot arguments must be logicals")}
	#
	##################################
	# CHECK RUN CONDITIONS
	#
	# nTotalTaxa, nExtant, nSamp must all be integer values
	if(!all(sapply(c(nTotalTaxa,nExtant,nSamp),function(x) x==round(x)))){
			stop("nTotalTaxa, nExtant, nSamp must coercible to whole number integers")}		
	#
	runConditions<-list(totalTime=totalTime,nTotalTaxa=nTotalTaxa,nExtant=nExtant,nSamp=nSamp)
	#check that all are numeric
	if(any(!sapply(runConditions,is.numeric))){
		stop("Run condition arguments must be all of type numeric")
		}
	#are length of 1 or 2
	if(any(sapply(runConditions,length)>2) | any(sapply(runConditions,length)<1)){
		stop("Run condition arguments must be of length 1 or 2")
		}
	# run conditions can be given as vectors of length 1 or 2
		# i.e. a point condition or range
	# turn run conditions of length 1 into vectors of length 2
	runConditions<-lapply(runConditions,function(x)
		if(length(x)==1){c(x,x)}else{x}
		)
	#all values are over or equal to zero
	if(any(!sapply(runConditions,function(x) all(x>=0)))){
		stop("Run Condition values must be equal to or greater than 0")
		}	
	#with minimums less than maximums
	if(any(!sapply(runConditions,function(x) x[1]<=x[2]))){
		stop("Run condition misordered: values given as a range must have the minimum before the maximum")
		}	
	###########################
	#get the basic rate functions
	getBranchRate<-makeParFunct(p,isBranchRate=TRUE)
	getExtRate<-makeParFunct(q,isBranchRate=FALSE)
	getSampRate<-makeParFunct(r,isBranchRate=FALSE)
	getAnagRate<-makeParFunct(anag.rate,isBranchRate=FALSE)
	#
	##############################################
	#now iterate for nruns
	results<-list()
	ntries<-0
	for(i in 1:nruns){
		accept<-FALSE
		while(!accept){
			ntries<-ntries+1
			#
			#initiate the taxa dataset
			timePassed<-0
			#currentTime is the max time from runConditions
			currentTime<-runConditions$totalTime[2]
			taxa<-initiateTaxa(startTaxa=startTaxa,time=currentTime)
			#
			#get vitals
			startVitals<-getRunVitals(taxa=taxa,count.cryptic=count.cryptic)
			#start vitals table		
			vitalsRecord<-cbind(timePassed=timePassed,t(as.matrix(startVitals)))
			#test to make sure run conditions aren't impossible
			continue<-testContinue(vitals=startVitals,timePassed=timePassed,
				runConditions=runConditions)
			if(!continue){
				stop("Initial starting point already matches given run conditions")
				}
			while(continue){
				#only as long as continue=TRUE
				#
				#timePassed from the initiation of the simulation
				timePassed<-runConditions$totalTime[2]-currentTime
				#
				# get rates, sample new event, have it occur
				#
				#get event probability vector
				rateVector<-getRateVector(taxa=taxa, timePassed=timePassed,
					getBranchRate=getBranchRate, getExtRate=getExtRate,
					getSampRate=getSampRate, getAnagRate=getAnagRate,
					prop.cryptic=prop.cryptic, prop.bifurc=prop.bifurc,
					negRatesAsZero=negRatesAsZero)
				#
				#sum the rates
				sumRates<-sum(rateVector)
				eventProb<-rateVector/sumRates
				#
				#pull type of event (from Peter Smits)
				event <- sample( names(eventProb), 1, prob = eventProb)
				#
				#vector of which taxa are still alive
				whichExtant<-whichLive(taxa)
				#select which lineage does it occur to
				if(length(whichExtant)>1){
					target<-sample(whichExtant,1)
				}else{
					target<-whichExtant
					}
				#
				#draw waiting time to an event (from Peter Smits)
				changeTime <- rexp(1, rate =sumRates*length(whichExtant))
				newTime<- currentTime - changeTime
				newTimePassed<-timePassed+changeTime
				#
				# make the new event so!
				taxa<-eventOccurs(taxa=taxa,target=target,type=event,time=newTime)
				#
				####################################################
				#
				#evaluate run conditions NOW for stopping
				#
				#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
					# none of these can REVERSE
				#
				#get vitals
				currentVitals<-getRunVitals(taxa=taxa,count.cryptic=count.cryptic)
				# continue ??
				continue<-testContinue(vitals=currentVitals,timePassed=newTimePassed,
					runConditions=runConditions)
				#
				# Updated vitals table
					#for (2), keep a table that records changes in nTotalTaxa, nExtant, nSamp with timePassed
				#then can quickly evaluate (2)
				currentVitals<-c(timePassed=newTimePassed,t(as.matrix(currentVitals)))
				vitalsRecord<-rbind(vitalsRecord,currentVitals)
				# set new current time
				currentTime<-newTime	
				#
				###############################################################
				# some archived debugging lines for posterity
				#if(newTimePassed>74.5){browser()}
				#if(newTimePassed>120){if(taxa[[4]][[1]][4]<120){browser()}}
				#
			###########################################
			#
			#accepting or rejecting runs
			#
			#discussion with Smits 05/11/15
				#real run condition is max limits / total ext for typical birth-death simulators
				#minimums are just for acceptability of runs when they hit run conditions
			#
			# NEED TO AVOID HARTMANN ET AL. EFFECT ---- simFossilTaxa did it wrong!!
				# sample simulation from intervals where it 'matched' run conditions
			#
			#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
				# none of these can REVERSE
			#(2) then go back, find all interval for which run conditions were met
				# if no acceptable intervals, reject run
			#(3) randomly sample within intervals for a single date, apply timeSliceFossilRecord
			#
			###########################################
			#
			# use vitalRecords to identify intervals of acceptable parameter values
			#		
			#is it even worth checking? (were mins reached)
			worthyVitals<-worthCheckingVitalsRecord(vitalsRecord=vitalsRecord,runConditions=runConditions)
			if(worthyVitals){
			#test with testVitalsRecord to get seqVitals
				seqVitals<-testVitalsRecord(vitalsRecord=vitalsRecord,runConditions=runConditions
					,tolerance=tolerance)
				if(all(!is.na(seqVitals))){
					#hey, if its an acceptable simulation!!!!!!
					accept<-TRUE
					}
				}
			}
		#sample the sequences for a date
		passedDate<-sampleSeqVitals(seqVitals=seqVitals)
		#this date is in timePassed units: convert to backwards currentTime
		currentDate<-runConditions$totalTime[2]-passedDate
		# now time slice
			# if stop and there are extant, evaluate if sampled at modern
			# 0< modern.samp.prob <1 need to randomly sample
		taxa<-timeSliceFossilRecord(fossilRecord=taxa,sliceTime=currentDate,
			shiftRoot4TimeSlice=shiftRoot4TimeSlice, modern.samp.prob=modern.samp.prob)
		#
		##############################################################################
		# FINAL CHECKS
		# test that the produced taxa object actually passed the runConditions
		finalTest<-testFinal(taxa=taxa,timePassed=passedDate,
			runConditions=runConditions,count.cryptic=count.cryptic)
		#are there any non-identical taxa in a simulation with pure cryptic speciation?
			if(anag.rate==0 & prop.cryptic==1 & startTaxa==1){
				taxaIDsTest<-sapply(taxa,function(x) x[[1]][6])
				if(any(!sapply(taxaIDsTest,function(x) all(x==taxaIDsTest)))){
					stop("non-cryptic taxa created in a simulation with pure cryptic speciation?!")
					}
				}
			}
		################################################################################
		#
		#name each normal taxon as t + ID 
			#cryptic taxa are cryptic id + . taxon number within that complex
		names(taxa)<-getTaxaNames(taxa=taxa)
		#sort if sortNames
		if(sortNames){
			taxa<-taxa[order(names(taxa))]
			}
		#
		results[[i]]<-taxa
		if(plot){
			taxaConvert<-fossilRecord2fossilTaxa(fossilRecord=taxa)
			#taxicDivCont(taxaConvert,int.length=0.2)
			fossilRanges<-fossilRecord2fossilRanges(fossilRecord=taxa, merge.cryptic=TRUE, ranges.only = TRUE)
			curveList<-list(taxaConvert,fossilRanges)
			multiDiv(curveList,plotMultCurves=TRUE,
				divPalette=c("black","red"),divLineType=c(1,2),main="")
			legend("topleft",legend=c("True Richness", "Sampled Richness"),
				col=c("black","red"),lty=c(1,2))
			if(nruns>1){
				title(paste0("Run Number ",i," of ",nruns))
				}
			}
		}
	if(print.runs){
		message(paste(
			nruns," runs accepted from ",ntries," total runs ("
			,signif(nruns/ntries,2)," Acceptance Probability)",sep=""))
		}
	if(nruns==1){results<-results[[1]]}
	return(results)	
	}	
