\name{simFossilTaxa}
\alias{simFossilTaxa}
\alias{simFossilTaxa_SRCond}
\title{Simulating Taxa in the Fossil Record}
\description{Functions for simulating taxon ranges and relationships under various models of evolution}
\usage{
simFossilTaxa(p, q, w = 0, u = 0, nruns = 1, mintaxa = 1, maxtaxa = 1000, maxtime = 100, nExtant = NULL, plot = F)
simFossilTaxa_SRCond(r, avgtaxa, p, q, w = 0, u = 0, nruns = 1, maxtime = 100, nExtant = NULL, plot = F)
}
\arguments{
  \item{p}{Instantaneous rate of speciation/branching}
  \item{q}{Instantaneous rate of extinction}
  \item{w}{Instantaneous rate of pseudoextinction/anagensis}
  \item{u}{Proportion of branching by bifurcating cladogenesis relative to budding cladogenesis}
  \item{nruns}{Number of datasets to be output}
  \item{mintaxa}{Minimum number of total taxa over the entire history of a clade necessary for a dataset to be accepted}
  \item{maxtaxa}{Maximum number of total taxa over the entire history of a clade necessary for a dataset to be accepted}
  \item{maxtime}{Maximum time units to run any given simulation before stopping it}
  \item{nExtant}{Maximum number of living taxa allowed in simulations}
  \item{plot}{Plot the diversity curves of the accepted datasets as they are simulated?}
  \item{r}{Instantaneous sampling rate per time unit}
  \item{avgtaxa}{Desired average number of taxa}
}
\details{
This function will simulate the diversification of clades where taxa are relatively morphologically static over long time intervals. 

Simulations will stop when various conditions are met (generally maxtaxa or maxtime). To reduce the effect of one condition, simply set the limit to an arbitrarily high number. Warning, this may cause the function to enter an unending loop. Hartmann et al. (2011) recently discovered a potential statistical artifact when branching simulations are conditioned on some maximum number of taxa. Thus, this function continues the simulation once maxtaxa is hit, until the next taxon (maxtaxa+1) originates. Once the simulation terminates, it is judged whether it is acceptable for all conditions given and if so, it is accepted as a dataset to be output.

Please note that mintaxa and maxtaxa refer to the number of static morphotaxa that were birthed over the entire history of the simulated clade. Use nExtant if you want to condition on the maximum number of taxa living at some time.

The simFossilTaxa function can effectively simulate clades evolving any combination of the three "modes" of speciation generally referred to by paleontologists: budding cladogenesis, branching cladogenesis and anagenesis (Foote, 1996). The first two are "speciation" in the typical sense used by biologists, with the major distinction between these two modes being whether the ancestral taxon shifts morphologically at the time of speciation. The third is where a morphotaxon changes into another morphotaxon with no branching, hence the use of the terms "pseudoextinction" and "pseudospeciation". As bifurcation and budding are both branching events, both are controlled by the p, the instantaneous rate, while the probability of a branching event being either is set by u. By default, only budding cladogenesis occurs To have these three modes occur in equal proportions, set p to be twice the value of w and set u to 0.5. There is no option for cryptic speciation in this function.

If nExtant is 0, then the function will be limited to only accepting simulations that end in total clade extinction before maxtime.

If conditions are such that a clade survives to maxtime, then maxtime will become the time of first appearance for the first taxa. Unless maxtime is very low, however, it is more likely the maxtaxa limit will be reached first, in which case the point in time at which maxtaxa is reached will become the present data and the entire length of the simulation will be the time of the first appearance of the first taxon.

simFossilTaxa_SRCond is a wrapper for simFossilTaxa for when you want clades of a particular size, post-sampling. This function accomplishes this task by calculating the probability of sampling per-taxon and calculating the average clade size needed to produce the number of sampled taxa given by avgtaxa. We will call that quantity N. Then, it uses simFossilTaxa, with mintaxa set to N and maxtaxa set to 2*N. It will generally produce simulated datasets that are generally of that size or larger post-sampling (although there can be some variance). Some combinations of p, q, r and avgtaxa may take an extremely long time to find large enough datasets. Some combinations may produce very strange datasets that may have weird structure that is only a result of the conditioning (for example, the only clades that have many taxa when net diversification is low or negative will have lots of very early divergences, which could impact analyses). Needless to say, conditioning can be very difficult.
}
\value{
Both of these functions give back a list containing nruns number of taxa datasets. Sampling has not been simulated in the output for either function; the output represents the 'true' history of the simualted clade.

For each dataset, the output is a five column per-taxon matrix where all entries are numbers, with the first column being the taxon ID, the second being the ancestral taxon ID (the first taxon is NA for ancestor), the third column is the first appearance date of a species in absolute time, the fourth column is the last appearance data and the fifth column records whether a species is still extant at the time the simulation terminated (1 means a taxon is still alive).

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.
}
\references{
Foote, M. 1996. On the Probability of Ancestors in the Fossil Record. Paleobiology 22(2):141-151.
}
\author{David W. Bapst}
\seealso{\code{\link{sampleRanges}},\code{\link{simPaleoTrees}},\code{\link{taxa2phylo}},\code{\link{taxa2cladogram}},
}
\examples{
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont
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
}
\keyword{datagen}
