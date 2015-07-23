"

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

#' @name
#' @rdname
#' @export











#' @param totalTime,nTotalTaxa,nExtant,NSamp These arguments represent stopping conditions,
#' respectively the total length of the simulation in time-units, the total number of taxa
#' over evolutionary history in the clade, the total number of extant taxa at the end of the
#' simulation and the total number of sampled taxa (not counting extant taxa sampled at the
#' modern day). These are used to determine when to end simulation runs, and whether to accept
#' or reject them as output. They can be input as a vector of two numbers, representing minimum
#' and maximum values of a range for accepted simulation runs (i.e. the simulation length can be between 0 and
#' 1000 time-steps, by default), or as a single number, representing a point condition (i.e. if
#' \code{nSamp = 100} then the only simulation runs with exactly 100 taxa sampled will be output).
#' Note that it is easy to set combinations of parameters and stopping conditions that are impossible
#' to produce satisfactory input under, in which case \code{simFossilRecord} would run in a nonstop loop.

#' @param negRatesAsZero A logical. Should rates calculated as a negative number cause the simulation to fail
#' with an error message (\code{ = FALSE}) or should these be treated as zero (\code{"= TRUE"}, the default).




#' @param p,q,r,anag.rate These parameters control the instantaneous ('per-capita') rates of branching, extinction,
#' sampling and anagenesis, respectively.

# infinite rates are treated as if 0 : this event type cannot occur!

#WRITE CHECKS

simFossilRecord<-function(

	# model parameters
	p, q, r, anag.rate=0, prop.bifurc=0, prop.cryptic=0,
	modern.samp.prob=1, startTaxa=1,  nruns=1,

	# stopping conditions can be given as vectors of length 1 or length 2 (= min,max)

	totalTime = c(0, 1000), nTotalTaxa = c(0, 1000), nExtant = c(0, 1000), nSamp = c(0, 1000))

	#control parameters
	
	count.cryptic=FALSE, negRatesAsZero=TRUE, print.runs=FALSE, sortNames=FALSE, plot=FALSE){

	min.cond=TRUE

	#example parameter sets
	#
	# DEFAULTS
	# anag.rate=0.1;prop.bifurc=0.1;prop.cryptic=0;startTaxa=1;nruns=1;
		# nTotalTaxa=c(1,1000);totalTime=c(1,1000);nSamp=c(0,1000);nExtant=c(0,1000);
		# plot=TRUE;count.cryptic=FALSE;print.runs=TRUE;sortNames=FALSE	
	#
	# BASIC RUN	with diversity-dep extinction
	# p=0.1;q='0.01*N';r=0.1;anag.rate=0.1;prop.bifurc=0.1;prop.cryptic=0;startTaxa=1;nruns=1;
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

	#FUNCTIONALIZE EACH EVENT TYPE
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
			#thus can come up with quick/simple ways of evaluating stopping conditions
			# e.g. evaluate number of sampled taxa by sum(length($sampling.times)>0)
		taxon<-list(taxa.data=taxaData, sampling.times=numeric())
		return(taxon)
		}
	
	origination<-function(taxa,ancID,time,looksLike=NULL){
		#adds a new taxon via branching or anagenesis
		newID<-max(sapply(taxa,function(x) x[[1]][1]))+1
		if(is.null(looksLike)){
			looksLike<-newID
			}
		newTaxonData<-newTaxon(newID=newID,ancID=ancID,time=time,looksLike=looksLike)	
		newTaxa<-c(taxa,newTaxonData)
		return(newTaxa)
		}

	termination<-function(taxa,targetID,time){
		#ends an existing taxon (from extinction or pseudoextinction)
		whichTarget<-which(sapply(taxa,function(x) x[[1]][1]==targetID))
		if(length(whichTarget)!=1){stop('taxon IDs repeated??')}
		taxa[[whichTarget]][[1]][4]<-time
		taxa[[whichTarget]][[1]][5]<-0
		return(taxa)
		}

	#newTaxaID<-function(taxa){
	#	newID<-max(sapply(taxa,function(x) x[[1]][1]))+1
	#	return(newID)
	#	}

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
		taxa<-termination(taxa=taxa,targetID=parent,time=time)
		return(taxa)
		}	
			
	bifurcEvent<-function(taxa,parent,time){
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
		taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
		taxa<-termination(taxa=taxa,targetID=parent,time=time)
		return(taxa)
		}	
		
	extEvent<-function(taxa,target,time){
		taxa<-termination(taxa=taxa,targetID=target,time=time)
		return(taxa)		
		}
		
	sampEvent<-function(taxa,target,time){
		whichTarget<-which(sapply(taxa,function(x) x[[1]][1]==targetID))
		taxa[[whichTarget]][[2]]<-c(taxa[[whichTarget]][[2]],time)
		return(taxa)		
		}
			
	eventOccurs<-(taxa,target,type)
		#possible types : 'budd','bifurc','anag','crypt','ext','samp'
		if(type=="budd"){
			taxa<-buddEvent(taxa=taxa,parent=target,time=time)
			}
		if(type=="bifurc"){
			taxa<-bifurcEvent(taxa=taxa,parent=target,time=time)		
			}
		if(type=="anag"){
			taxa<-anagEvent(taxa=taxa,parent=target,time=time)
			}
		if(type=="crypt"){
			taxa<-cryptEvent(taxa=taxa,parent=target,time=time)
			}
		if(type=="ext"){
			taxa<-extEvent(taxa=taxa,target=target,time=time)		
			}	
		if(type=="samp"){
			taxa<-sampEvent(taxa=taxa,target=target,time=time)
			}
		return(taxa)
		}
	


	




	#simplified birth-death-sampling simulator for fossil record
	#
	#p is the rate of branching
		#may be either budding (prob = 1-prop.bifurc) or bifurcation (prob = prop.bifurc)
	#anagenesis is a separate rate (anag.rate)
	#q is rate of extinction
	#r is rate of sampling
	#
	#ARGUMENT CHECKING

	


	

#OLD SIM FOSSIL TAXA CODE
'
	
	#CHECKING
	if(any(c(p,q,anag.rate,prop.bifurc,prop.cryptic)<0)){stop(
		"bad parameters input, p, q, anag.rate, prop.bifurc or prop.cryptic are less than 0")}
	if(prop.bifurc>0 & prop.cryptic==1){stop("Prop.bifurc greater than 0 even though cryptic cladogenesis = 1??")}
	if(nruns<1){stop("nruns<1")}
	if(maxtaxa<0){stop("maxtaxa<0")}
	if(mintaxa<1){stop("mintaxa<1")}
	if(mintime<1){stop("mintime<1")}
	if(maxtime<mintime){stop("maxtime<mintime")}
	if(mintaxa>maxtaxa){stop("mintaxa > maxtaxa")}
	if(maxtaxa>10000 & maxtime>10000){warning("Warning: Unrealistic limits for maxtaxa or maxtime")}
	if(minExtant<0){stop("minExtant<0")}
	if(minExtant>mintaxa){mintaxa<-minExtant}
	if(!is.null(maxExtant)){
		if(maxExtant<0){stop("maxExtant<0")}
		if(maxExtant>maxtaxa){maxtaxa<-maxExtant}
		if(minExtant>maxExtant){stop("maxExtant is set higher than minExtant")}
		}
	if(!min.cond){message("No conditioning during simulation; run until max limits or total extinction")}
	#end idiot proofing

'		

 

		

	
	#
	##################################
	# CHECK STOPPING CONDITIONS
	stoppingConditions<-list(totalTime=totalTime,nTotalTaxa=nTotalTaxa,nExtant=nExtant,nSamp=nSamp)
	#check that all are numeric
	if(any(!sapply(stoppingConditions,is.numeric))){
		stop("Stopping condition arguments must be all of type numeric")
		}
	#are length of 1 or 2
	if(any(sapply(stoppingConditions,length)>2) | any(sapply(stoppingConditions,length)<1)){
		stop("Stopping condition arguments must be of length 1 or 2")
		}
	# stopping conditions can be given as vectors of length 1 or 2
		# i.e. a point condition or range
	# turn stopping conditions of length 1 into vectors of length 2
	stoppingConditions<-lapply(stoppingConditions,function(x)
		if(length(x)==1){c(x,x)}else{x}
		)
	#all values are over zero
	if(any(!sapply(stoppingConditions,function(x) all(x>=0)))){
		stop("Stopping Condition values must be equal to or greater than 0")
		}	
	#with minimums less than maximums
	if(any(!sapply(stoppingConditions,function(x) x[1]<=x[2]))){
		stop("Stopping condition misordered: values given as a range must have the minimum before the maximum")
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
		ntries<-ntries+1
		#initiate the taxa dataset
			#currentTime is the max time from stoppingConditions
		currentTime<-stoppingConditions$totalTime[2]
		taxa<-initiateTaxa(startTaxa=startTaxa,time=currentTime)
		continue<-TRUE
		while(continue){
			#get some basic summary statistics first
			#vector of which taxa are still alive
			whichLive<-which(sapply(taxa,function(x) x[[1]][5]==1))
			#standing number of extant lineages 
			# (i.e. number of lineages that stuff can happen to)
			nLive<-length(whichLive)
			#timePassed from the initiation of the simulation
			timePassed<-totalTime[2]-currentTime	
			#
			#########################
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
			#
			#sum the rates
			sumRates<-sum(rateVector)
			eventProb<-rateVector/sumRates
			##################################
			#
			#pull type of event (from Peter Smits)
			event <- sample( names(rateVector), 1, prob = eventProb)
			#select which lineage does it occur to
			target<-sample(whichLive,1)
			#
			#draw waiting time to an event (from Peter Smits)
			changeTime <- rexp(1, rate =sumRates*nLive)
			newTime<- currentTime + changeTime
			#
			###########################################
			#
			#evalutate stopping conditions NOW
			#
			# NEED TO AVOID HARMANN ET AL. EFFECT
				# run
			#
			#discussion with Smits 05/11/15
				#real stopping condition is max limits for typical birth-death simulators
					#except total extinction...
				#minimums are just for acceptability of runs when they hit stopping conditions
			#
			#how should treat min/max bounds?
				#
				
				#how simFossiltaxa does it:
				#Simulations will not go longer than maxtime, period
				#when maxtaxa is hit, simulation will go up until FO of maxtaxa+1 taxon to avoid Hartmann et al. effect
				#if minExtant is set, simulation will end once minExtant is hit
					#unless maxExtant is zero, in which case 
				#if maxExtant is set, simulation will end once maxExtant is hit, if before maxtaxa or mintime is hit
					#if maxtaxa or mintime is hit, run is discarded if maxExtant is not satisfied
					#when maxExtant is hit, simulation will go up until FO of maxExtant+1 to avoid Hartmann et al. effect
		
		if(newTime>	





		
	# if stop and there are extant, evaluate if sampled at modern
		# 0< modern.samp.prob <1 need to randomly sample
		
			

			#
			############################################
			# additional tasks
			# only add new event if continue = TRUE
			if(continue = TRUE){
				#update taxa with the new event
				taxa<-eventOccurs(taxa=taxa,target=target,type=event)
				}
			# set new current time
			currentTime<-newTime	

		

		
		
		






################################################



	#make new min/max extant
	minExtant1<-ifelse(minExtant==0,0,minExtant+1)


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
					#then youll need to evaluate it again to make sure it meets criteria
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
					if(any(is.na(taxad[,4]))){		#if its dead, dont change maxtime1...
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
				#dont just use one maxtime1, use a bunch 02-07-12: lets you use more runs!
				taxad<-matrix(taxad[taxad[,3]<maxtime1,],sum(taxad[,3]<maxtime1),)
				if(any(is.na(taxad[,4]))){stop("Live creatures escaping simulation! Get out now while you still have time!")}
				posstimes<-sort(unique(c(taxad[,3:4],maxtime1)))
				maxtimes<-posstimes[posstimes>=mintime & posstimes<=maxtime1]				#make vector of maxtimes
				if(length(maxtimes)==0){
					eval<-FALSE
				}else{
					mtds<-lapply(maxtimes,function(x) matrix(taxad[taxad[,3]<x,],sum(taxad[,3]<x),)) 	#maxtime taxa datasets
					if(count.cryptic){
						numtaxa<-sapply(mtds,function(x) nrow(x))
						numexta<-sapply(1:length(mtds),function(x) sum(mtds[[x]][4]>=maxtimes[x]))		#number of extant taxa
					}else{	#do not count cryptics
						numtaxa<-sapply(mtds,function(x) length(unique(x[,5])))
						numexta<-sapply(1:length(mtds),function(x) length(unique(mtds[[x]][mtds[[x]][4]>=maxtimes[x],5])))
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
				#reset if the clade is done (continue=FALSE) but eval is FALSE (didnt hit conditions)
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
			if(nruns>1){title(paste("Run Num.",i," of ",nruns,sep=""))}
			}
		}
	if(print.runs){message(paste(nruns," runs accepted from ",ntries," total runs (",signif(nruns/ntries,2)," Acceptance Probability)",sep=""))}
	if(nruns==1){results<-results[[1]]}
	return(results)	
	
}

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


# a function that transforms a simfossilrecord to a taxa object

# a function that transforms a simfossilrecord to a set of ranges (like from sampleRanges)
	# merge.cryptic = TRUE or FALSE
	# ranges.only or sampling times?

# a function that wraps taxa2phylo for simFossilRecord, providing time-scaled tree of sampled taxa
	# merge.cryptic = TRUE or FALSE
	#ala simPaleoTrees:
		# tree<-taxa2phylo(taxa,obs_time=ranges1[,2],plot=plot)	

