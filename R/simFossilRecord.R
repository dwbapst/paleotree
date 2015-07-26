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



#' Hartmann et al. (2011) recently discovered a potential statistical artifact
#' when branching simulations are conditioned on some maximum number of taxa.




#maybe export?

#' @export
timeSliceFossilRecord<-function(fossilRecord,sliceTime,modern.samp.prob=1){
	#take a fossilRecord data object and cut it at some specific date
	#
	#drop all taxa that originate after the sliceTime
	droppers<-which(sapply(fossilRecord,function(x) x[[1]][3]<sliceTime))
	fossilRecord<-fossilRecord[-droppers]
	#
	#turn all taxa that went extinct after sliceTime so they are still alive
	stillAlive<-which(sapply(fossilRecord,function(x) x[[1]][4]<sliceTime))
	fossilRecord[stillAlive][[1]][4:5]<-c(NA,1)
	#
	#remove all sampling events after sliceTime
	#adjust all dates so cutdate becomes 0
	for(i in 1:length(fossilRecord)){
		#remove all sampling events after sliceTime
		fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]][fossilRecord[[i]][[2]]>=sliceTime]
		#adjust all dates so cutdate becomes 0
		fossilRecord[[i]][[1]][3:4]<-fossilRecord[[i]][[1]][3:4]-sliceTime
		fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]]-sliceTime
		}
	#
	# sample at modern based on modern.samp.prob
	whichExtant<-which(sapply(fossilRecord,function(x) x[[1]][5]==1))
	nLive<-length(whichExtant)
	liveSampled<-as.logical(rbinom(n=nLive, size=1, prob=modern.samp.prob))
	whichSampled<-whichExtant[liveSampled]
	#add sampling event at modern
	for(i in whichSampled){
		fossilRecord[[i]][[2]]<-c(fossilRecord[[i]][[2]],0)
		}
	#
	return(fossilRecord)
	}





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

# infinite rates are treated as if 0 : this event type cannot occur!

#WRITE CHECKS

	min.cond=TRUE

#' @export
simFossilRecord<-function(

	# model parameters
	#
	p, q, r, anag.rate=0, prop.bifurc=0, prop.cryptic=0,
	modern.samp.prob=1, startTaxa=1,  nruns=1,

	# run conditions can be given as vectors of length 1 or length 2 (= min,max)
	#
	totalTime = c(0, 1000), nTotalTaxa = c(1, 1000),
	nExtant = c(0, 1000), nSamp = c(0, 1000))

	#control parameters
	#
	count.cryptic=FALSE, negRatesAsZero=TRUE, print.runs=FALSE, sortNames=FALSE, plot=FALSE){



	#example parameter sets
	#
	# DEFAULTS (without rates set)
		# anag.rate=0.1;prop.bifurc=0.1;prop.cryptic=0;startTaxa=1;nruns=1;
		# nTotalTaxa=c(1,1000);totalTime=c(1,1000);nSamp=c(0,1000);nExtant=c(0,1000);
		# plot=TRUE;count.cryptic=FALSE;print.runs=TRUE;sortNames=FALSE	
	#
	# BASIC RUN	with diversity-dep extinction
	# p=0.1;q='0.01*N';r=0.1
		# anag.rate=0.1;prop.bifurc=0.1;prop.cryptic=0;startTaxa=1;nruns=1;
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

	# get rate vector

	getRateVector<-function(taxa,timePassed,
			getBranchRate,getExtRate,getSampRate,getAnagRate,
			prop.cryptic,prop.bifurc,negRatesAsZero){
		#
		#get some basic summary statistics first
		#vector of which taxa are still alive
		whichLive<-whichLive(taxa)
		#standing number of extant lineages 
		# (i.e. number of lineages that stuff can happen to)
		nLive<-length(whichLive)	
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
			
	eventOccurs<-(taxa,target,type,time)
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

	# functions for identifying live/sampled taxa and checking simulation runs
	
	whichLive<-function(taxa){
		res<-which(sapply(taxa,function(x) x[[1]][5]==1))
		return(res)
		}

	whichSampled<-function(taxa){
		res<-which(sapply(taxa,function(x) length(x[[2]])>0))
		return(res)
		}
	
	getRunVitals<-function(taxa,count.cryptic){
		#NOTE need to change vital measurement dependent on count.cryptic or not
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
		vector<-as.integer(vector)
		#because unbelievably base R has no simple function for
			#pulling contiguous sequences of integers from 
		starts<-sapply(2:length(vector),function(i){
			vector[i]-vector[i-1]>1
			})
		starts<-c(TRUE,starts)
		starts<-vector[starts]
		ends<-sapply(1:(length(vector)-1),function(i){
			vector[i+1]-vector[i]>1
			})
		ends<-c(ends,TRUE)
		ends<-vector[ends]
		seqMat<-cbind(starts,ends)
		return(seqMat)
		}

	insertRow<-function(table,row,rownum){
		#because unbelievably base R has no simple function for inserting a row
		# insert new immediately at this number, shifts row currently at that location DOWN
		table<-rbind(table[1:(rownum-1),],row,table[-(1:(rownum-1)),])
		return(table)
		}

	testVitalsRecord<-function(vitalsRecord,runConditions){
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
		# NOW need to essentially duplicated EVERY ROW with time-stamp of row immediately after it
		newVitalsRecord<-cbind(vitalsRecord[2:nrow(vitalsRecord),1],
			vitalsRecord[1:(nrow(vitalsRecord)-1),2:4])
		pastIncrement<-diff(vitalsRecord[,1])
		pastIncrement<-min(pastIncrement[pastIncrement>0])/1000
		newTimes<-newVitalsRecord[,1]-pastIncrement
		newTimes<-ifelse(newTimes>0,newTimes,0)
		newVitalsRecord[,1]<-newTimes
		vitalsRecord<-rbind(vitalsRecord,newVitalsRecord)
		vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
		#
		###########################################################################
		# identify all rows where nTaxa, nExtant and nSamp are good
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
			seqVitals<-contiguousIntegerSeq(which(okayVitals))
			#replaced with the passedTime dates
			seqVitals<-apply(seqVitals,2,function(x) vitalsRecord[x,1])
		}else{
			seqVitals<-NA
			}
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
			date<-sample(seqVitals[,1],1)
			}					
		return(date)
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

	
# number of starting taxa must be at least 1

	# check that prop.bifurc, prop.cryptic are greater than 0
	if(any(c(prop.bifurc,prop.cryptic)<0)){stop(
		"bad parameters input: prop.bifurc or prop.cryptic are less than 0")}	

#OLD SIM FOSSIL TAXA CODE
'
	
	#CHECKING

	if(prop.bifurc>0 & prop.cryptic==1){stop("Prop.bifurc greater than 0 even though cryptic cladogenesis = 1??")}
	if(nruns<1){stop("nruns<1")}
	if(maxtaxa<0){stop("maxtaxa<0")}
	if(mintaxa<1){stop("mintaxa<1")}
	if(mintime<1){stop("mintime<1")}

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

 	#check that min nSamp isn't higher that 0, if r = 0 or Inf

# nTotalTaxa, nExtant, nSamp must all be integer values
		

	
	#
	##################################
	# CHECK RUN CONDITIONS
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
	#all values are over zero
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
			startVitals<-getRunVitals(taxa,count.cryptic=count.cryptic)
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
				rateVector<-getRateVector(taxa=taxa,timePassed=timePassed,
					getBranchRate=getBranchRate,getExtRate=getExtRate,
					getSampRate=getSampRate,getAnagRate=getAnagRate,
					prop.cryptic=prop.cryptic,prop.bifurc=prop.bifurc,
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
				whichLive<-whichLive(taxa)
				#select which lineage does it occur to
				target<-sample(whichLive,1)
				#
				#draw waiting time to an event (from Peter Smits)
				changeTime <- rexp(1, rate =sumRates*nLive)
				newTime<- currentTime - changeTime
				newTimePassed<-timePassed+changeTime
				#
				# make the new event so!
				taxa<-eventOccurs(taxa=taxa,target=target,type=event,time=newTime)
				#
				####################################################
				#
				#evalutate run conditions NOW
				#
				#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
					# none of these can REVERSE
				#
				#get vitals
				currentVitals<-getRunVitals(taxa,count.cryptic=count.cryptic)
				# continue ??
				continue<-testContinue(vitals=currentVitals,timePassed=newTimePassed,
					runConditions=runConditions)
				#
				# Updated vitals table
					#for (2), keep a table that records changes in nTotalTaxa, nExtant, nSamp with timePassed
				#then can quickly evaluate (2)
				currentVitals<-c(timePassed=timePassed,t(as.matrix(currentVitals)))
				vitalsRecord<-rbind(vitalsRecord,currentVitals)
				# set new current time
				currentTime<-newTime	
				}
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
			#test with testVitalsRecord to get seqVitals
			seqVitals<-testVitalsRecord(vitalsRecord,runConditions)
			if(!is.na(seqVitals)){
				#hey, if its an acceptable simulation!!!!!!
				accept<-TRUE
				}
			}
		#sample the sequences for a date
		passedDate<-sampleSeqVitals(seqVitals=seqVitals)
		#this date is in timePassed units: convert to currentTime
		currentDate<-runConditions$totalTime[2]-passedDate
		timeSliceFossilRecord(fossilRecord=taxa,sliceTime=currentDate)
		# if stop and there are extant, evaluate if sampled at modern
		# 0< modern.samp.prob <1 need to randomly sample
	

		}










		

		
			



		

		
		
		






################################################



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

