\name{simFossilTaxa}
\alias{simFossilTaxa}
\alias{simFossilTaxa_SRCond}
\title{Simulating Taxa in the Fossil Record}
\description{Functions for simulating taxon ranges and relationships under various models of evolution}
\usage{
simFossilTaxa(p, q, w = 0, u = 0, nruns = 1, mintaxa = 1, 
   maxtaxa = 1000, mintime = 1,maxtime = 1000, minExtant = 0, 
   maxExtant = NULL, min.cond = TRUE, print.runs = FALSE, plot = FALSE)

simFossilTaxa_SRCond(r, avgtaxa, p, q, w = 0, u = 0, nruns = 1, 
   maxtime = 1000, maxExtant = NULL, print.runs = FALSE, plot = FALSE)
}
\arguments{
  \item{p}{Instantaneous rate of speciation/branching}
  \item{q}{Instantaneous rate of extinction}
  \item{w}{Instantaneous rate of pseudoextinction/anagensis}
  \item{u}{Proportion of branching by bifurcating cladogenesis relative to budding cladogenesis}
  \item{nruns}{Number of datasets to be output}
  \item{mintaxa}{Minimum number of total taxa over the entire history of a clade necessary for a dataset to be accepted}
  \item{maxtaxa}{Maximum number of total taxa over the entire history of a clade necessary for a dataset to be accepted}
  \item{mintime}{Minimum time units to run any given simulation before stopping it}
  \item{maxtime}{Maximum time units to run any given simulation before stopping it}
  \item{minExtant}{Minimum number of living taxa allowed at end of simulations}
  \item{maxExtant}{Maximum number of living taxa allowed at end of simulations}
  \item{min.cond}{Stop simulations when they hit minimum conditions or go until they hit maximum conditions?}
  \item{print.runs}{Print the proportion of simulations accepted for output?}
  \item{plot}{Plot the diversity curves of the accepted datasets as they are simulated?}
  \item{r}{Instantaneous sampling rate per time unit}
  \item{avgtaxa}{Desired average number of taxa}
}
\details{
simFossilTaxa simulates a birth-death process (Kendall, 1948; Nee, 2006), but unlike most functions for this implemented in R, this function simulates the diversification of clades where taxa are relatively morphologically static over long time intervals. The output is a description of the temporal and phylogenetic relationships of those morphotaxa. This is meant to emulate the sort of data that paleobiologists often work with, especially in well-sampled groups.

If min.cond=T (the default), simulations will stop when clades satisfy mintime, mintaxa, minExtant and maxExtant (if the later is set). To reduce the effect of one condition, simply set the limit to an arbitrarily low number. If min.cond=F, then the simulations are not stopped until they (a) go extinct or (b) hit either maxtaxa or maxtime. Whether they are accepted or not for output is still dependent on mintaxa, mintime, minExtant and maxExtant. Note that some combinations of conditions, such as setting minExtant=maxExtant>0 

Hartmann et al. (2011) recently discovered a potential statistical artifact when branching simulations are conditioned on some maximum number of taxa. Thus, this function continues the simulation once mintaxa or minExtant is hit, until the next taxon (limit +1) originates. Once the simulation terminates, it is judged whether it is acceptable for all conditions given and if so, it is accepted as a dataset to be output.

Please note that mintaxa and maxtaxa refer to the number of static morphotaxa that were birthed over the entire evolutionary history of the simulated clade, not the extant richness at the end of the simulation. Use minExtant and maxExtant if you want to condition on the number of taxa living at some time.

The simFossilTaxa function can effectively simulate clades evolving any combination of the three "modes" of speciation generally referred to by paleontologists: budding cladogenesis, branching cladogenesis and anagenesis (Foote, 1996). The first two are "speciation" in the typical sense used by biologists, with the major distinction between these two modes being whether the ancestral taxon shifts morphologically at the time of speciation. The third is where a morphotaxon changes into another morphotaxon with no branching, hence the use of the terms "pseudoextinction" and "pseudospeciation". As bifurcation and budding are both branching events, both are controlled by the p, the instantaneous rate, while the probability of a branching event being either is set by u. By default, only budding cladogenesis occurs To have these three modes occur in equal proportions, set p to be twice the value of w and set u to 0.5. There is no option for cryptic speciation in this function.

If maxExtant is 0, then the function will be limited to only accepting simulations that end in total clade extinction before maxtime.

If conditions are such that a clade survives to maxtime, then maxtime will become the time of first appearance for the first taxa. Unless maxtime is very low, however, it is more likely the maxtaxa limit will be reached first, in which case the point in time at which maxtaxa is reached will become the present data and the entire length of the simulation will be the time of the first appearance of the first taxon.

simFossilTaxa simulates single taxa until they go extinct or exceed maxtime. This means the function may have fully simulated some lineages for thousands of time-steps while others are not yet simulated, meaning it can sometimes overshoot constraints on the number of taxa. This function will automatically discard any runs where the number of taxa exceeds 2 x maxtaxa to avoid blowing up the computation time. This is likely to happen under a pure-birth scenario; I suggest using low maxtime settings if doing a pure-birth simulation.

simFossilTaxa_SRCond is a wrapper for simFossilTaxa for when you want clades of a particular size, post-sampling. This function accomplishes this task by first calculating the expected proportion of taxa sampled, given the sampling rate and the rates which control lineage termination: extinction, anagenesis and bifurcation. The average original clade size needed to produce the number of sampled taxa given by avgtaxa is calculated with the following equation:

\eqn{N = (average number of taxa desired) / ( 1 - exp( -r / (q + w + (u * p)) ))}

We will call that quantity N. Note that u * p describes the rate of bifurcation, as u is the ratio of budding to bifurcating cladogenesis.

This value is used with simFossilTaxa, where mintaxa is set to N and maxtaxa set to 2*N. simFossilTaxa_SRcond will generally produce simulated datasets that are generally of that size or larger post-sampling (although there can be some variance). Some combinations of p, q, r and avgtaxa may take an extremely long time to find large enough datasets. Some combinations may produce very strange datasets that may have weird structure that is only a result of the conditioning (for example, the only clades that have many taxa when net diversification is low or negative will have lots of very early divergences, which could impact analyses). Needless to say, conditioning can be very difficult.
}
\value{
Both of these functions give back a list containing nruns number of taxa datasets, where each element is a matrix. If nruns=1, the output is  not a list but just a single matrix. Sampling has not been simulated in the output for either function; the output represents the 'true' history of the simualted clade.

For each dataset, the output is a five column per-taxon matrix where all entries are numbers, with the first column being the taxon ID, the second being the ancestral taxon ID (the first taxon is NA for ancestor), the third column is the first appearance date of a species in absolute time, the fourth column is the last appearance data and the fifth column records whether a species is still extant at the time the simulation terminated (a value of 1 indicates a taxon is still alive, a value of 0 indicates the taxon is extinct).

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.
}
\references{
Foote, M. 1996. On the Probability of Ancestors in the Fossil Record. Paleobiology 22(2):141-151.

Hartmann, K., D. Wong, and T. Stadler. 2010. Sampling Trees from Evolutionary Models. Systematic Biology 59(4):465-476.

Kendall, D. G. 1948. On the Generalized "Birth-and-Death" Process. The Annals of Mathematical Statistics 19(1):1-15.

Nee, S. 2006. Birth-Death Models in Macroevolution. Annual Review of Ecology, Evolution, and Systematics 37(1):1-17.
}
\author{David W. Bapst}
\seealso{\code{\link{sampleRanges}},\code{\link{simPaleoTrees}},\code{\link{taxa2phylo}},\code{\link{taxa2cladogram}},
}
\examples{
dev.new()
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

#pure-birth example
#note that conditioning is tricky
simFossilTaxa(p=0.1,q=0,mintime=10,mintaxa=100,maxtime=100,maxtaxa=100,nruns=10,plot=T)

#can generate datasets where simulations go until extinction or max limits
	#and THEN are evaluated whether they meet min limits
	#good for producing unconditioned birth-death trees
dev.new()
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
}
