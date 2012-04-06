pkgname <- "paleotree"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('paleotree')

assign(".oldSearch", search(), pos = 'CheckExEnv')
assign(".ExTimings", "paleotree-Ex.timings", pos = 'CheckExEnv')
cat("name\tuser\tsystem\telapsed\n", file=get(".ExTimings", pos = 'CheckExEnv'))
assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  format(x[1L:3L])
},
pos = 'CheckExEnv')

cleanEx()
nameEx("DiversityCurves")
### * DiversityCurves

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: DiversityCurves
### Title: Diversity Curves
### Aliases: taxicDivCont taxicDivDisc phyloDiv
### Keywords: datagen

### ** Examples

set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
layout(matrix(1:2,,2))
taxicDivCont(rangesCont)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)
#compare to the continuous time diversity curve
#Now let's make a tree using taxa2phylo()
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
phyloDiv(tree)

#a simple example with phyloDiv
  #using a tree from rtree in ape
set.seed(444)
tree<-rtree(100)
phyloDiv(tree)

#a neat example of using phyDiv with timeSliceTree 
 #to simulate doing molecular-phylogeny studies 
 #of diverification...in the past
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
taxicDivCont(taxa)
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve 
 #we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=TRUE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("DiversityCurves", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SamplingConv")
### * SamplingConv

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SamplingConv
### Title: Converting Sampling Estimates
### Aliases: sProb2sRate sRate2sProb pqsRate2sProb qsProb2Comp qsRate2Comp
###   probAnc
### Keywords: datagen

### ** Examples

sRate2sProb(r=0.5)
sProb2sRate(R=0.1)
pqsRate2sProb(r=0.5,p=0.1,q=0.1)
qsProb2Comp(R=0.1,q=0.1)
qsRate2Comp(r=0.1,q=0.1)
probAnc(p=0.1,q=0.1,R=0.5)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("SamplingConv", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("binTimeData")
### * binTimeData

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: binTimeData
### Title: Bin Temporal Ranges in Discrete Intervals
### Aliases: binTimeData
### Keywords: datagen

### ** Examples

#Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("binTimeData", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("degradeTree")
### * degradeTree

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: degradeTree
### Title: Randomly collapse nodes on a phylogeny
### Aliases: degradeTree
### Keywords: datagen

### ** Examples

set.seed(444)
tree<-rtree(100)
tree1<-degradeTree(tree,0.5)
#let's compare the input and output
layout(matrix(1:2,,2))
plot(tree,show.tip.label=FALSE,use.edge.length=FALSE)
plot(tree1,show.tip.label=FALSE,use.edge.length=FALSE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("degradeTree", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("depthRainbow")
### * depthRainbow

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: depthRainbow
### Title: Paint Tree Branch Depth by Color
### Aliases: depthRainbow
### Keywords: hplot

### ** Examples

set.seed(444)
tree<-rtree(500)
depthRainbow(tree)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("depthRainbow", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("droppingBranches")
### * droppingBranches

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: droppingBranches
### Title: Drop Terminal Zero-Length Branches
### Aliases: dropZLB dropExtinct dropExtant
### Keywords: datagen

### ** Examples

##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's make a tree using taxa2phylo()
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
#compare the two trees
layout(matrix(1:2,,2))
plot(ladderize(tree))
plot(ladderize(dropZLB(tree)))


#example using dropExtinct and dropExtant
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=40,maxtime=1000,maxExtant=20)
tree<-taxa2phylo(taxa)
phyloDiv(tree)
tree1<-dropExtinct(tree)
phyloDiv(tree1)
tree2<-dropExtant(tree)
phyloDiv(tree2)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("droppingBranches", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("expandTaxonTree")
### * expandTaxonTree

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: expandTaxonTree
### Title: Extrapolating Lower Taxon Phylogenies from Higher Taxon Trees
### Aliases: expandTaxonTree
### Keywords: datagen

### ** Examples

set.seed(444)
#lets make our hypothetical simulated tree of higher taxa
taxtr<-rtree(10)
taxd<-sample(taxtr$tip.label,30,replace=TRUE)	#taxa to place within higher taxa
names(taxd)<-paste(taxd,"_x",1:30,sep="")
coll<-sample(taxtr$tip.label,3)	#what to collapse?
expandTaxonTree(taxonTree=taxtr,taxaData=taxd,collapse=coll,plot=TRUE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("expandTaxonTree", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getSampProbDisc")
### * getSampProbDisc

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getSampProbDisc
### Title: Fit Models of Sampling Probability to Discrete-Interval Taxon
###   Ranges
### Aliases: getSampProbDisc
### Keywords: datagen

### ** Examples

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
sRate2sProb(r=0.5)
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

## Not run: 
##D #this data was simulated under homogenous sampling rates, extinction rates
##D #if we fit a model with random groups and allow for multiple timebins
##D 	#AIC should be higher (less informative models)
##D randomgroup<-sample(1:2,nrow(rangesDisc[[2]]),replace=TRUE)
##D SPres2<-getSampProbDisc(rangesDisc,grp1=randomgroup)
##D SPres3<-getSampProbDisc(rangesDisc,n_tbins=2)
##D print(c(SPres1$AICc,SPres2$AICc,SPres3$AICc))
##D #and we can see the most simple model has the lowest AICc (most informative model)
##D 
##D #testing temporal change in sampling rate
##D set.seed(444)
##D taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=100,maxtaxa=125,maxtime=1000,maxExtant=0,plot=T)
##D #let's see what the 'true' diversity curve looks like in this case
##D #simulate two sets of ranges at r=0.7 and r=0.1
##D rangesCont<-sampleRanges(taxa,r=1.1)
##D rangesCont2<-sampleRanges(taxa,r=0.2)
##D #now make it so that taxa which originated after 850 have r=0.1
##D rangesCont[taxa[,3]<850,]<-rangesCont2[taxa[,3]<850,]
##D rangesDisc<-binTimeData(rangesCont)
##D #lets plot the diversity curve
##D taxicDivDisc(rangesDisc)
##D SPres1<-getSampProbDisc(rangesDisc)
##D SPres2<-getSampProbDisc(rangesDisc,n_tbins=2)
##D #compare the AICc of the models
##D print(c(SPres1$AICc,SPres2$AICc)) #model 2 looks pretty good
##D #when does it find the break in time intervals?
##D print(rangesDisc[[1]][SPres2$t_ends[2],1])
##D #not so great: estimates 940, not 850 
##D 	#but look at the diversity curve: most richness in bin 1 is before 940
##D 	#might have found the right break time otherwise...
##D #the parameter values it found are less great. Finds variation in q	
## End(Not run)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("getSampProbDisc", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getSampRateCont")
### * getSampRateCont

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getSampRateCont
### Title: Fit Models of Sampling Rates to Continuous-Time Taxon Ranges
### Aliases: getSampRateCont
### Keywords: datagen

### ** Examples

##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
(SRres1<-getSampRateCont(rangesCont))
#that's all the results...
sRate<-SRres1$pars[2]
print(sRate)	#estimates that sRate=~0.4 (not too bad...)

#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins
	#AIC should be higher (less informative)
randomgroup<-sample(1:2,nrow(rangesCont),replace=TRUE)
SRres2<-getSampRateCont(rangesCont,grp1=randomgroup)
SRres3<-getSampRateCont(rangesCont,n_tbins=2)
SRres4<-getSampRateCont(rangesCont,n_tbins=3,grp1=randomgroup)
print(c(SRres1$AICc,SRres2$AICc,SRres3$AICc,SRres4$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("getSampRateCont", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multiDiv")
### * multiDiv

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multiDiv
### Title: Calculating Diversity Curves Across Multiple Datasets
### Aliases: multiDiv plotMultiDiv
### Keywords: datagen

### ** Examples
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
rangesCont<-sampleRanges(taxa,r=0.5)
rangesDisc<-binTimeData(rangesCont,int.length=1)
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#using multiDiv with very different data types
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=FALSE)
input<-list(rangesCont,rangesDisc,ttree)
multiDiv(input,plot=TRUE)

#using multiDiv with samples of trees
ttrees<-timePaleoPhy(cladogram,rangesCont,type="basic",randres=TRUE,ntrees=10,add.term=TRUE)
multiDiv(ttrees)
#uncertainty in diversity history is solely due to 
 #the random resolution of polytomies

#multiDiv can also take output from simFossilTaxa
#what do many simulations run under some conditions 'look' like on average?
set.seed(444)
taxa<-simFossilTaxa(p=0.3,q=0.1,nruns=20,maxtime=20,maxtaxa=100,plot=TRUE,min.cond=FALSE)
multiDiv(taxa)
#increasing cone of diversity! Even better on a log scale:
multiDiv(taxa,plotLogRich=TRUE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("multiDiv", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("paleotree-package")
### * paleotree-package

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: paleotree-package
### Title: paleotree: Paleontological and Phylogenetic Analyses of
###   Evolution
### Aliases: paleotree-package paleotree
### Keywords: datagen

### ** Examples

##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa)

#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
layout(matrix(1:2,,2))
taxicDivCont(rangesCont)

#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)
#compare to the continuous time diversity curve

#taxa2phylo assumes we know speciation events perfectly... what if we don't?
#first, let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy() using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#plot diversity curve 
phyloDiv(ttree,drop.ZLB=TRUE)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("paleotree-package", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotTraitgram")
### * plotTraitgram

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotTraitgram
### Title: Plot a Traitgram for Continuous Traits
### Aliases: plotTraitgram
### Keywords: hplot

### ** Examples

require(geiger)
set.seed(444)
tree<-rtree(10)
trait<-rTraitCont(tree)
#first, traitgram without conf intervals
plotTraitgram(trait,tree,conf.int=FALSE)

#now, with
plotTraitgram(trait,tree)
#not much confidence, eh?




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("plotTraitgram", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sampleRanges")
### * sampleRanges

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sampleRanges
### Title: Sampling Taxon Ranges
### Aliases: sampleRanges
### Keywords: datagen

### ** Examples

set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
layout(matrix(1:2,2,))
#let's see what the 'true' diversity curve looks like in this case
taxicDivCont(taxa)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
taxicDivCont(rangesCont)
#compare the true history to what we would observe!




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("sampleRanges", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simFossilTaxa")
### * simFossilTaxa

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simFossilTaxa
### Title: Simulating Taxa in the Fossil Record
### Aliases: simFossilTaxa simFossilTaxa_SRCond
### Keywords: datagen

### ** Examples

set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont
taxicDivCont(taxa[,3:4])
#can also see this by setting plot=TRUE in simFossilTaxa

#make datasets with multiple speciation modes
#following has anagenesis, budding cladogenesis and bifurcating cladogenesis
	#all set to 1/2 extinction rate
set.seed(444)
res<-simFossilTaxa(p=0.1,q=0.1,w=0.05,u=0.5,mintaxa=30,maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
#what does this mix of speciation modes look like as a phylogeny?
tree<-taxa2phylo(res,plot=TRUE)

#can generate datasets that meet multiple conditions: time, # total taxa, # extant taxa
set.seed(444)
res<-simFossilTaxa(p=0.1,q=0.1,mintime=10,mintaxa=30,maxtaxa=40,minExtant=10,maxExtant=20,nruns=20,plot=FALSE,print.runs=TRUE)
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
res<-simFossilTaxa(p=0.1,q=0.1,maxtaxa=100,maxtime=100,nruns=10,plot=TRUE,print.runs=TRUE,min.cond=FALSE)
#hey, look, we accepted everything! (That's what we want.)
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




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("simFossilTaxa", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simPaleoTrees")
### * simPaleoTrees

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simPaleoTrees
### Title: Simulating Un-Conditioned Trees of Fossil Taxa
### Aliases: simPaleoTrees
### Keywords: datagen

### ** Examples

set.seed(444)
#simulate trees conditioned to have no living descendants
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=10,all.extinct=TRUE,maxtime=100,plot=TRUE)
#number of tips
sapply(trees,Ntip)

#simulate trees conditioned to (possibly) have living descendants and perfect sampling at modern
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=10,all.extinct=FALSE,maxtime=100,modern.samp=TRUE,plot=TRUE)
#number of tips
sapply(trees,Ntip)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("simPaleoTrees", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("srcTimescaling")
### * srcTimescaling

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: srcTimescaling
### Title: SampRate-Calibrated Timescaling of Paleo-Phylogenies
### Aliases: srcTimePaleoPhy bin_srcTimePaleoPhy
### Keywords: datagen

### ** Examples

##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#this library allows one to use SRC type time-scaling methods (Bapst, in prep.)
#to use these, we need an estimate of the sampling rate (we set it to 0.5 above)
SRres<-getSampRateCont(rangesCont)
sRate<-SRres$pars[2]
#now let's try srcTimePaleoPhy(), which timescales using a sampling rate to calibrate
#This can also resolve polytomies based on sampling rates, with some stochastic decisions
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,plot=TRUE)
#notice the warning it gives!
phyloDiv(ttree)

#by default, srcTimePaleoPhy() is allowed to predict indirect ancestor-descendant relationships
#can turn this off by setting anc.wt=0
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,anc.wt=0,plot=TRUE)

#to get a fair sample of trees, let's increse ntrees
#can do an example of such an analysis via multDiv
ttrees<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=10,plot=FALSE)
multiDiv(ttrees)

#example with time in discrete intervals
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#we can do something very similar for the discrete time data (can be a bit slow)
SPres<-getSampProbDisc(rangesDisc)
sProb<-SPres$pars[2]
#but that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
sRate1<-sProb2sRate(sProb,int.length=1)
#estimates that r=0.3... kind of low (simulated sampling rate is 0.5)
#Note: for real data, you may need to use an average int.length (no constant length)
ttree<-bin_srcTimePaleoPhy(cladogram,rangesDisc,sampRate=sRate1,ntrees=1,plot=TRUE)
phyloDiv(ttree)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("srcTimescaling", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("taxa2cladogram")
### * taxa2cladogram

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: taxa2cladogram
### Title: Converting taxon data into cladogram
### Aliases: taxa2cladogram
### Keywords: datagen

### ** Examples

set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
layout(matrix(1:2,,2))
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#compare the real relationships (taxa2phylo) to the ideal cladogram
tree<-taxa2phylo(taxa,plot=TRUE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("taxa2cladogram", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("taxa2phylo")
### * taxa2phylo

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: taxa2phylo
### Title: Convert taxon data into Phylogeny
### Aliases: taxa2phylo
### Keywords: datagen

### ** Examples

set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
tree<-taxa2phylo(taxa)
phyloDiv(tree)
#now a phylogeny with tips placed at the apparent time of extinction
rangesCont<-sampleRanges(taxa,r=0.5)
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
phyloDiv(tree,drop.ZLB=FALSE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("taxa2phylo", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("timePaleoPhy")
### * timePaleoPhy

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: timePaleoPhy
### Title: Timescaling of Paleo-Phylogenies
### Aliases: timePaleoPhy bin_timePaleoPhy
### Keywords: datagen

### ** Examples

##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)

#that tree didn't look very resolved, does it? (The curse of sampled ancestry!)
#if we set ntrees>1, timePaleoPhy() will make multiple time-trees
#each resulting tree will have polytomies randomly resolved in different ways using multi2di()
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=1,randres=TRUE,add.term=TRUE,plot=TRUE)
#notice well the warning it prints!
#now let's plot the first tree (both trees will be identical because we used set.seed)
phyloDiv(ttree)
#we would need to set ntrees to a large number to get a fair sample of trees

#compare different methods of timePaleoPhy
layout(matrix(1:6,3,2));par(mar=c(3,2,1,2))
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="basic",vartime=NULL,add.term=TRUE)))
	axisPhylo();text(x=50,y=23,"type=basic",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="equal",vartime=10,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=equal",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="aba",vartime=1,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=aba",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="zlba",vartime=1,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=zlba",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="mbl",vartime=1,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=mbl",adj=c(0,0.5),cex=1.2)

#Using bin_timePaleoPhy to timescale with discrete interval data
#first let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
ttree<-bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,randres=TRUE,add.term=TRUE,plot=FALSE)
#notice the warning it prints!
phyloDiv(ttree)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("timePaleoPhy", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("timeSliceTree")
### * timeSliceTree

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: timeSliceTree
### Title: Taking a Timeslice of a Tree
### Aliases: timeSliceTree
### Keywords: datagen

### ** Examples

#a neat example of using phyloDiv with timeSliceTree 
 #to simulate doing molecular-phylogeny studies 
 #of diverification...in the past
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
taxicDivCont(taxa)
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve 
 #we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=TRUE)




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("timeSliceTree", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("unitLengthTree")
### * unitLengthTree

flush(stderr()); flush(stdout())

assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: unitLengthTree
### Title: Scale Tree to Unit-Length
### Aliases: unitLengthTree
### Keywords: datagen

### ** Examples

set.seed(444)
tree<-rtree(10)
layout(matrix(1:2,,2))
plot(tree)
plot(unitLengthTree(tree))




assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
cat("unitLengthTree", get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
