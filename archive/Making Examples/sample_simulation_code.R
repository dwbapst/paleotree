source(file.choose())#source functions

##Example using:
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa[,3:4])

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

#taxa2phylo assumes we know speciation events perfectly... what if we don't?
#first, let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=T)
#Now let's try timePaleoPhy() using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=T)
#plot diversity curve 
phyloDiv(ttree)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=T,plot=T)
#plot diversity curve 
phyloDiv(ttree)
#that tree didn't look very resolved, does it? (The curse of sampled ancestry!)


#if we set ntrees>1, timePaleoPhy() will make multiple time-trees
#each resulting tree will have polytomies randomly resolved in different ways using multi2di()
ttrees<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=2,add.term=T,plot=T)
#now let's plot the first tree (both trees will be identical because we used set.seed)
phyloDiv(ttrees[[1]]);is.binary.tree(ttrees[[1]])
#we would need to set ntrees to a large number to get a fair sample of trees

#this library allows one to use SRC type time-scaling methods (Bapst, in prep.)
#to use these, we need an estimate of the sampling rate (we set it to 0.5 above)
(SRres<-getSampRateCont(rangesCont))
sRate<-SRres[[2]][2]
print(sRate)	#estimates that sRate=~0.4 (not too bad...)
#now let's try srcTimePaleoPhy(), which timescales using a sampling rate to calibrate
#This can also resolve polytomies based on sampling rates, with some stochastic decisions
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,plot=T)
phyloDiv(ttree);is.binary.tree(ttree)
#Again, we would need to set ntrees to a large number to get a fair sample of trees

#by default, srcTimePaleoPhy() is allowed to predict indirect ancestor-descendant relationships
#can turn this off by setting anc.wt=0
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,anc.wt=0,plot=T)
phyloDiv(ttree);is.binary.tree(ttree)

#we can do something very similar for the discrete time data (can be a bit slow)
(SPres<-getSampProbDisc(rangesDisc))
sProb<-SPres[[2]][2]
#but that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
(sRate1<-sProb2sRate(sProb,int.length=1))
#estimates that r=0.3... kind of low (simulated sampling rate is 0.5)
#Note: for real data, you may need to use an average int.length (no constant length)
ttree<-bin_srcTimePaleoPhy(cladogram,rangesDisc,sampRate=sRate1,ntrees=1,plot=T)
phyloDiv(ttree);is.binary.tree(ttree)

#as with srcTimePaleoPhy(), set anc.wt=0 to not allow for ancestral relationships
ttree<-bin_srcTimePaleoPhy(cladogram,rangesDisc,sampRate=sRate1,ntrees=1,anc.wt=0,plot=T)
phyloDiv(ttree);is.binary.tree(ttree)

#can also calculate median diversity curves across multiple different data types
multiDiv(list(rangesCont,rangesDisc,ttree))

#examples of using conversion metrics
sRate2sProb(r=0.5)
sProb2sRate(R=0.1)
pqsRate2sProb(r=0.5,p=0.1,q=0.1)
qsProb2Comp(R=0.1,q=0.1)
qsRate2Comp(r=0.1,q=0.1)
probAnc(p=0.1,q=0.1,R=0.5)

#traitgram
require(geiger)
tree<-rtree(10)
trait<-rTraitCont(tree)
plotTraitgram(trait,tree)

#depthRainbow
set.seed(444)
tree<-rtree(100)
depthRainbow(tree)

#unitLengthTree
set.seed(444)
tree<-rtree(10)
layout(matrix(1:2,,2))
plot(tree);plot(unitLengthTree(tree))

#timeSliceTree
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
taxicDivCont(taxa[,3:4])
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=T,drop.extinct=F)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=T,drop.extinct=T)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=T)

#expandTaxonTree
taxtr<-rtree(10)
coll<-sample(taxtr$tip.label,3)
taxd<-sample(taxtr$tip.label,30,replace=T)
names(taxd)<-paste(taxd,"_x",1:30,sep="")
expandTaxonTree(taxonTree=taxtr,taxaData=taxd,collapse=coll,plot=T)

#degradeTree
set.seed(444)
tree<-rtree(100)
tree1<-degradeTree(tree,0.5,node.depth=NA)
layout(matrix(1:2,,2))
plot(tree,show.tip.label=F,use.edge.length=F)
plot(tree1,show.tip.label=F,use.edge.length=F)

#multiDiv
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
rangesCont<-sampleRanges(taxa,r=0.5)
rangesDisc<-binTimeData(rangesCont,int.length=1)
cladogram<-taxa2cladogram(taxa,plot=T)
#using multiDiv with very different data types
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=T,plot=F)
multiDiv(list(rangesCont,rangesDisc,ttree),plot=T,output=F)
#using multiDiv with samples of trees
ttrees<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=10,add.term=T,plot=F)
multiDiv(ttrees)

#simFossilTaxa_SRcond
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa[,3:4])

set.seed(444)
avgtaxa=50
r<-0.5
taxa<-simFossilTaxa_SRCond(r=r,p=0.1,q=0.1,nruns=20,avgtaxa=avgtaxa)
#sample and count number of taxa
ranges<-lapply(taxa,sampleRanges,r=r)
ntaxa<-sapply(ranges,function(x) sum(!is.na(x[,1])))
#works okay... some parameter combinations are difficult to get right number of taxa
hist(ntaxa);mean(ntaxa)

#simPaleoTrees
#simulate trees with no living descendants
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=2,nExtant=0,modern.samp=T,drop.zlb=T,plot=T)
sapply(trees,Ntip)
#simulate trees with living descendants
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=2,nExtant=10,modern.samp=T,drop.zlb=T,plot=T)
sapply(trees,Ntip)

#getSampRateCont
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
(SRres1<-getSampRateCont(rangesCont))
#that's all the results...
sRate<-SRres$pars[2]
print(sRate)	#estimates that sRate=~0.4 (not too bad...)
#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins, AIC should be higher (less informative)
randomgroup<-sample(1:2,nrow(rangesCont),replace=T)
SRres2<-getSampRateCont(rangesCont,grp1=randomgroup)
SRres3<-getSampRateCont(rangesCont,n_tbins=2)
SRres4<-getSampRateCont(rangesCont,n_tbins=3,grp1=randomgroup)
print(c(SRres1$AICc,SRres2$AICc,SRres3$AICc,SRres4$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)

#getSampProbDisc
#Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
#we can do something very similar for the discrete time data (can be a bit slow)
SPres1<-getSampProbDisc(rangesDisc)
sProb<-SPres1[[2]][2]
print(sProb)
#note that that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we want the rate for srcTimePaleoPhy()
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
sRate<-sProb2sRate(sProb,int.length=1)
print(sRate)
#estimates that r=0.3... kind of low... (simulated sampling rate is 0.5)
#Note: for real data, you may need to use an average int.length (no constant length)
#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins, AIC should be higher (less informative)
randomgroup<-sample(1:2,nrow(rangesDisc[[2]]),replace=T)
SPres2<-getSampProbDisc(rangesDisc,grp1=randomgroup)
SPres3<-getSampProbDisc(rangesDisc,n_tbins=2)
print(c(SPres1$AICc,SPres2$AICc,SPres3$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)

########################

#timePaleoPhy
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=T)
#Now let's try timePaleoPhy() using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=T)
#plot diversity curve 
phyloDiv(ttree)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=T,plot=T)
#plot diversity curve 
phyloDiv(ttree)

#that tree didn't look very resolved, does it? (The curse of sampled ancestry!)
#if we set ntrees>1, timePaleoPhy() will make multiple time-trees
#each resulting tree will have polytomies randomly resolved in different ways using multi2di()
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=1,randres=T,add.term=T,plot=T)
#notice that the warning it prints!
#now let's plot the first tree (both trees will be identical because we used set.seed)
phyloDiv(ttree)
#we would need to set ntrees to a large number to get a fair sample of trees

#compare different methods of timePaleoPhy
layout(matrix(1:6,3,2));par(mar=c(3,2,1,2))
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="basic",vartime=NULL,add.term=T)))
	axisPhylo();text(x=50,y=23,"type=basic",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="equal",vartime=10,add.term=T)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=equal",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="aba",vartime=1,add.term=T)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=aba",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="zlba",vartime=1,add.term=T)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=zlba",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="mbl",vartime=1,add.term=T)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=mbl",adj=c(0,0.5),cex=1.2)

#Using bin_timePaleoPhy to timescale with discrete interval data
#first let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
ttree<-bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,randres=T,add.term=T,plot=T)
#notice that the warning it prints!
phyloDiv(ttree)

################

#SRC methods

##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=T)
#this library allows one to use SRC type time-scaling methods (Bapst, in prep.)
#to use these, we need an estimate of the sampling rate (we set it to 0.5 above)
SRres<-getSampRateCont(rangesCont)
sRate<-SRres$pars[2]
#now let's try srcTimePaleoPhy(), which timescales using a sampling rate to calibrate
#This can also resolve polytomies based on sampling rates, with some stochastic decisions
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,plot=T)
#notice the warning it gives!
phyloDiv(ttree)

#Again, we would need to set ntrees to a large number to get a fair sample of trees
#can do an example of such an analysis via multDiv
ttrees<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=10,plot=F)
multiDiv(ttrees)

#by default, srcTimePaleoPhy() is allowed to predict indirect ancestor-descendant relationships
#can turn this off by setting anc.wt=0
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,anc.wt=0,plot=T)

#we can do something very similar for the discrete time data (can be a bit slow)
SPres<-getSampProbDisc(rangesDisc)
sProb<-SPres$pars[2]
#but that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
sRate1<-sProb2sRate(sProb,int.length=1)
#estimates that r=0.3... kind of low (simulated sampling rate is 0.5)
#Note: for real data, you may need to use an average int.length (no constant length)
ttree<-bin_srcTimePaleoPhy(cladogram,rangesDisc,sampRate=sRate1,ntrees=1,plot=T)
phyloDiv(ttree)


