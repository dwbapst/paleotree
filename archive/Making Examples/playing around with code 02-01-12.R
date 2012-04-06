source(file.choose())#source functions
#library(phangorn);library(geiger)

#timeSliceTree
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
taxicDivCont(taxa)
#that's the whole ("true") diversity curve
#with timeSliceTree we could look at the lineage accumulation curve we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=T,drop.extinct=F)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=T,drop.extinct=T)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=T)

#example using dropExtinct and dropExtant
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=40,maxtime=1000,maxExtant=20)
tree<-taxa2phylo(taxa)
phyloDiv(tree)
tree1<-dropExtinct(tree)
phyloDiv(tree1)
tree2<-dropExtant(tree)
phyloDiv(tree2)

taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)

res<-simFossilTaxa(p=0.1,q=0.1,mintime=10,mintaxa=20,nruns=100,plot=T,print.runs=T)
hist(sapply(res,nrow))

#simFossilTaxa

set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont
taxicDivCont(taxa[,3:4])

#make datasets with multiple speciation modes
#following has anagenesis, budding cladogenesis and bifurcating cladogenesis
	#all set to 1/2 extinction rate
set.seed(444)
res<-simFossilTaxa(p=0.1,q=0.1,w=0.05,u=0.5,mintaxa=30,maxtaxa=60,maxExtant=0,nruns=1,plot=T)
#what does this mix of speciation modes look like as a phylogeny? 
	#(note this ignores static morphotaxon identity)
tree<-taxa2phylo(res,plot=T)

#can generate datasets that meet multiple conditions: time, # total taxa, # extant taxa
set.seed(444)
res<-simFossilTaxa(p=0.1,q=0.1,mintime=10,mintaxa=30,maxtaxa=40,minExtant=10,maxExtant=20,nruns=20,plot=F,print.runs=T)
#use print.run to know how many simulations were accepted of the total generated
layout(matrix(1:2,2,))
#histogram of # taxa over evolutionary history
hist(sapply(res,nrow),main="#taxa")
#histogram of # extant taxa at end of simulation
hist(sapply(res,function(x) sum(x[,5])),main="#extant")

#can generate datasets where simulations go until extinction or max limits
	#and THEN are evaluated whether they meet min limits
	#good for producing unconditioned birth-death trees
set.seed(444)
res<-simFossilTaxa(p=0.1,q=0.1,maxtaxa=100,maxtime=100,nruns=10,plot=T,print.runs=T,no.cond=T)
#hey, look, we accepted everything!
layout(matrix(1:2,2,))
#histogram of # taxa over evolutionary history
hist(sapply(res,nrow),main="#taxa")
#histogram of # extant taxa at end of simulation
hist(sapply(res,function(x) sum(x[,5])),main="#extant")

#using the SRcond version
set.seed(444)
avgtaxa<-50
r<-0.5
taxa<-simFossilTaxa_SRCond(r=r,p=0.1,q=0.1,nruns=20,avgtaxa=avgtaxa)
#now let's use sampleRanges and count number of sampled taxa
ranges<-lapply(taxa,sampleRanges,r=r)
ntaxa<-sapply(ranges,function(x) sum(!is.na(x[,1])))
hist(ntaxa);mean(ntaxa)
#works okay... some parameter combinations are difficult to get right number of taxa

######################


#Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
#for discrete data we can estimate the sampling probability per interval (R)
	#i.e. this is not the same thing as the instantaneous sampling rate (r)
#can use sRate2sProb to see what we would expect
sRate2sProb(r)
#expect R = ~0.39
#now we can use maximum likelihood to taxon ranges to get sampling probability
SPres1<-getSampProbDisc(rangesDisc)
sProb<-SPres1[[2]][2]
print(sProb)
#est. R = ~0.42; not too off what we would expect!
#for the src based timescaling methods, we want an estimate of the instantaneous samp rate
#we can use sProb2sRate() to get the rate. We will also need to also tell it the int.length
sRate<-sProb2sRate(sProb,int.length=1)
print(sRate)
#estimates that r=0.54... Not bad!
#Note: for real data, you may need to use an average int.length (no constant length)

\dontrun{
#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebinsAIC should be higher (less informative models)
randomgroup<-sample(1:2,nrow(rangesDisc[[2]]),replace=TRUE)
SPres2<-getSampProbDisc(rangesDisc,grp1=randomgroup)
SPres3<-getSampProbDisc(rangesDisc,n_tbins=2)
print(c(SPres1$AICc,SPres2$AICc,SPres3$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)
}

##############







