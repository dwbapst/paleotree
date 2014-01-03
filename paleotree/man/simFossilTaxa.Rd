\name{simFossilTaxa}
\alias{simFossilTaxa}
\title{Simulating Taxa in the Fossil Record}
\usage{
simFossilTaxa(p, q, anag.rate = 0, prop.bifurc = 0, prop.cryptic = 0,
  nruns = 1, mintaxa = 1, maxtaxa = 1000, mintime = 1, maxtime = 1000,
  minExtant = 0, maxExtant = NULL, min.cond = TRUE,
  count.cryptic = FALSE, print.runs = FALSE, sortNames = FALSE,
  plot = FALSE)
}
\arguments{
  \item{p}{Instantaneous rate of origination/branching per
  lineage-time units.}

  \item{q}{Instantaneous rate of extinction per
  lineage-time units.}

  \item{anag.rate}{Instantaneous rate of anagenesis (i.e.
  pseudospeciation/pseudoextinction).}

  \item{prop.bifurc}{Proportion of morphological branching
  by bifurcating cladogenesis relative to budding
  cladogenesis.}

  \item{prop.cryptic}{Proportion of branching events with
  no morphological differentiation (i.e. cryptic
  speciation) relative to branching events associated with
  morphological differentiation (budding and/or bifurcating
  cladogenesis).}

  \item{nruns}{Number of datasets to accept, save and
  output.}

  \item{mintaxa}{Minimum number of total taxa over the
  entire history of a clade necessary for a dataset to be
  accepted.}

  \item{maxtaxa}{Maximum number of total taxa over the
  entire history of a clade necessary for a dataset to be
  accepted.}

  \item{mintime}{Minimum time units to run any given
  simulation before stopping.}

  \item{maxtime}{Maximum time units to run any given
  simulation before stopping.}

  \item{minExtant}{Minimum number of living taxa allowed at
  end of simulations.}

  \item{maxExtant}{Maximum number of living taxa allowed at
  end of simulations.}

  \item{min.cond}{If TRUE, the default, simulations are
  stopped when they meet all minimum conditions. If FALSE,
  simulations will continue until they hit maximum
  conditions, but are only accepted as long as they still
  meet all minimum conditions in addition.}

  \item{count.cryptic}{If TRUE, cryptic taxa are counted as
  separate taxa for conditioning limits, such as maxtaxa or
  maxExtant. If FALSE, then each cryptic complex (i.e. each
  distinguishable morphotaxon) is treated as a single
  taxon}

  \item{print.runs}{If TRUE, prints the proportion of
  simulations accepted for output to the terminal.}

  \item{sortNames}{If TRUE, output taxonomic matrices are
  sorted by the taxon names (rownames; thus sorting cryptic
  taxa together) rather than by taxon id (which is the
  order they were simulated in).}

  \item{plot}{If TRUE, plots the diversity curves of
  accepted simulations.}
}
\value{
This function gives back a list containing nruns number of
taxa datasets, where each element is a matrix. If nruns=1,
the output is not a list but just a single matrix. Sampling
has not been simulated in the output for either function;
the output represents the 'true' history of the simulated
clade.

For each dataset, the output is a six column per-taxon
matrix where all entries are numbers, with the first column
being the taxon ID, the second being the ancestral taxon ID
(the first taxon is NA for ancestor), the third column is
the first appearance date of a species in absolute time,
the fourth column is the last appearance data and the fifth
column records whether a species is still extant at the
time the simulation terminated (a value of 1 indicates a
taxon is still alive, a value of 0 indicates the taxon is
extinct). The sixth column (named "looks.like") gives
information about the morphological distinguishability of
taxa; if they match the taxon ID, they are not cryptic. If
they do not match, then this column identifies which taxon
id they would be identified as.

Each matrix of simulated data also has rownames, generally
of the form "t1" and "t2", where the number is the taxon
id. Cryptic taxa are instead named in the form of "t1.2"
and "t5.3", where the first number is the taxon which they
are a cryptic descendant of (i.e. column 6 of the matrix,
"looks.like"). The second number, after the period, is the
rank order of taxa in that cryptic group of taxa. Taxa
which are the common ancestor of a cryptic lineage are also
given a unique naming convention, of the form "t1.1" and
"t5.1", where the first number is the taxon id and the
second number communicates that this is the first species
in a cryptic lineage.

As with many functions in the paleotree library, absolute
time is always decreasing, i.e. the present day is zero.
}
\description{
Functions for simulating taxon ranges and relationships
under various models of evolution
}
\details{
simFossilTaxa simulates a birth-death process (Kendall,
1948; Nee, 2006), but unlike most functions for this
implemented in R, this function enmeshes the simulation of
speciation and extinction with explicit models of how
lineages are morphologically differentiated, as morphotaxa
are the basic units of paleontological estimates of
diversity and phylogenetics.

Any particular use of simFossilTaxa will probably involve
iteratively running many simulations of diversification.
Simulation runs are only accepted for output if and when
they meet the conditioning criteria defined in the
arguments, both minima and maxima. If min.cond is true (the
default), simulations will be stopped and accepted when
clades satisfy mintime, mintaxa, minExtant and maxExtant
(if the later is set). To reduce the effect of one
conditioning criterion, simply set that limit to an
arbitrarily low or high number (depending if a minimum or
maximum constraint is involved). If min.cond is false,
simulation runs are not stopped and evaluated for output
acceptance until they (a) go completely extinct or (b) hit
either maxtaxa or maxtime. Whether the simulation runs are
accepted or not for output is still dependent on mintaxa,
mintime, minExtant and maxExtant. Note that some
combinations of conditions, such as attempting to condition
on a specific non-zero value of minExtant and maxExtant,
may take a long time to find any acceptable simulation
runs.

Hartmann et al. (2011) recently discovered a potential
statistical artifact when branching simulations are
conditioned on some maximum number of taxa. Thus, this
function continues the simulation once mintaxa or minExtant
is hit, until the next taxon (limit +1) originates. Once
the simulation terminates, it is judged whether it is
acceptable for all conditions given and if so, the run is
accepted as a dataset for output.

Please note that mintaxa and maxtaxa refer to the number of
static morphotaxa birthed over the entire evolutionary
history of the simulated clade, not the extant richness at
the end of the simulation. Use minExtant and maxExtant if
you want to condition on the number of taxa living at some
time.

The simFossilTaxa function can effectively simulate clades
evolving under any combination of the three "modes" of
speciation generally referred to by paleontologists:
budding cladogenesis, branching cladogenesis and anagenesis
(Foote, 1996). The first two are "speciation" in the
typical sense used by biologists, with the major
distinction between these two modes being whether the
ancestral taxon shifts morphologically at the time of
speciation. The third is where a morphotaxon changes into
another morphotaxon with no branching, hence the use of the
terms "pseudoextinction" and "pseudospeciation". As
bifurcation and budding are both branching events, both are
controlled by the p, the instantaneous rate, while the
probability of a branching event being either is set by u.
By default, only budding cladogenesis occurs To have these
three modes occur in equal proportions, set p to be twice
the value of w and set u to 0.5.

This function also includes the ability to simulate cryptic
cladogenesis. The available patterns of morphological
speciation thus form a gradient: cryptic cladogenesis has
no morphological shifts in either daughter branches after a
branching event, budding cladogenesis has one morphological
shift in the two daughter lineages and, in bifurcating
cladogenesis, shifts occur in both daughter lineages. The
argument prop.cryptic dictates what proportion of
branching/cladogenesis events (the overall occurance of
which with rate p) are cryptic versus those that have some
morphological divergence (either budding of bifurcating.
prop.bifurc controls the proportion of morphologically
divergent cladogenesis which is bifurcating relative to
budding. Thus, for example, the probability of a given
cladogenesis event being budding is
(1-prop.cryptic)*prop.bifurc.

When there is cryptic speciation, by default, the
conditioning arguments involving numbers of taxa (mintaxa,
maxtaxa, minExtant and maxExtant) count the number of
unique morphologically distinguishable taxa is checked
(i.e. the number of unique values in column 6 of the
simulated data). This behavior can be changed with the
argument count.cryptic.See below about the output data
structure to see how information about cryptic cladogenesis
is recorded. The functions taxa2phylo, taxa2cladogram and
taxicDivCont each handle cryptic species in different ways,
as described in their respective help files.

If maxExtant is 0, then the function will be limited to
only accepting simulations that end in total clade
extinction before maxtime.

If conditions are such that a clade survives to maxtime,
then maxtime will become the time of first appearance for
the first taxa. Unless maxtime is very low, however, it is
more likely the maxtaxa limit will be reached first, in
which case the point in time at which maxtaxa is reached
will become the present data and the entire length of the
simulation will be the time of the first appearance of the
first taxon.

simFossilTaxa simulates single taxa until they go extinct
or exceed maxtime. This means the function may have fully
simulated some lineages for thousands of time-steps while
others are not yet simulated, and thus sometimes overshoot
constraints on the number of taxa. This function will
automatically discard any runs where the number of taxa
exceeds 2 x maxtaxa to avoid blowing up computation time.
This is likely to happen under a pure-birth scenario; I
suggest using low maxtime settings if doing a pure-birth
simulation.

simFossilTaxa_SRCond is a wrapper for simFossilTaxa for
when clades of a particular size are desired,
post-sampling. For more details, see the help file at
\code{\link{simFossilTaxa_SRCond}}.

More details on this function's design can be read here:

http://nemagraptus.blogspot.com/2012/04/simulating-fossil-record.html
}
\examples{
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont
taxicDivCont(taxa[,3:4])
#can also see this by setting plot=TRUE in simFossilTaxa

#make datasets with multiple speciation modes
#following has anagenesis, budding cladogenesis and bifurcating cladogenesis
    #all set to 1/2 extinction rate
set.seed(444)
res <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0.05,prop.bifurc=0.5,mintaxa=30,maxtaxa=60,
    maxExtant=0,nruns=1,plot=TRUE)
#what does this mix of speciation modes look like as a phylogeny?
tree <- taxa2phylo(res,plot=TRUE)

\dontrun{
#some other options with cryptic speciation
taxaCrypt1 <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0,prop.bifurc=0,prop.crypt=0.5,mintaxa=30,
    maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
tree1 <- taxa2phylo(taxaCrypt1,plot=TRUE)
taxaCrypt2 <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0.05,prop.bifurc=0.5,prop.crypt=0.5,
    mintaxa=30,maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
tree2 <- taxa2phylo(taxaCrypt2,plot=TRUE)
taxaCrypt3 <- simFossilTaxa(p=0.1,q=0.1,anag.rate=0.05,prop.bifurc=0,prop.crypt=1,
    mintaxa=30,maxtaxa=60,maxExtant=0,nruns=1,plot=TRUE)
tree3 <- taxa2phylo(taxaCrypt2,plot=TRUE)
}

set.seed(444)
#can choose to condition on total morphologically-distinguishable taxa
    #or total taxa including cryptic taxa with count.cryptic=FALSE
taxa<-simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,mintaxa=20,
    count.cryptic=FALSE,plot=TRUE)
nrow(taxa)                 #number of lineages (inc. cryptic)
length(unique(taxa[,6]))            #number of morph-distinguishable taxa
#now with count.cryptic=TRUE
taxa <- simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,mintaxa=20,
    count.cryptic=TRUE,plot=TRUE)
nrow(taxa)                  #number of lineages (inc. cryptic)
length(unique(taxa[,6]))              #number of morph-distinguishable taxa

#now let's look at extant numbers of taxa
#can choose to condition on total morphologically-distinguishable living taxa
    #or total living taxa including cryptic taxa with count.cryptic=FALSE
taxa <- simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,minExtant=20,
    count.cryptic=FALSE,plot=TRUE)
    sum(taxa[,5])             #number of still-living lineages (inc. cryptic)
    length(unique(taxa[taxa[,5]==1,6]))	  #number of still-living morph-dist. taxa
#now with count.cryptic=TRUE
taxa <- simFossilTaxa(0.1,0.1,prop.cryptic=1,anag.rate=0.05,minExtant=20,
    count.cryptic=TRUE,plot=TRUE)
    sum(taxa[,5])              #number of still-living lineages (inc. cryptic)
    length(unique(taxa[taxa[,5]==1,6]))	   #number of still-living morph-dist, taxa

#can generate datasets that meet multiple conditions: time, # total taxa, # extant taxa
set.seed(444)
res <- simFossilTaxa(p=0.1,q=0.1,mintime=10,mintaxa=30,maxtaxa=40,minExtant=10,maxExtant=20,
    nruns=20,plot=FALSE,print.runs=TRUE)
#use print.run to know how many simulations were accepted of the total generated
layout(1:2)
#histogram of # taxa over evolutionary history
hist(sapply(res,nrow),main="#taxa")
#histogram of # extant taxa at end of simulation
hist(sapply(res,function(x) sum(x[,5])),main="#extant")

\dontrun{
#pure-birth example
#note that conditioning is tricky
layout(1)
taxa <- simFossilTaxa(p=0.1,q=0,mintime=10,mintaxa=100,maxtime=100,maxtaxa=100,
    nruns=10,plot=TRUE)

#can generate datasets where simulations go until extinction or max limits
     #and THEN are evaluated whether they meet min limits
     #good for producing unconditioned birth-death trees
set.seed(444)
res <- simFossilTaxa(p=0.1,q=0.1,maxtaxa=100,maxtime=100,nruns=10,plot=TRUE,
    print.runs=TRUE,min.cond=FALSE)
#hey, look, we accepted everything! (That's what we want.)
layout(1:2)
#histogram of # taxa over evolutionary history
hist(sapply(res,nrow),main="#taxa")
#histogram of # extant taxa at end of simulation
hist(sapply(res,function(x) sum(x[,5])),main="#extant")
}

layout(1)
}
\author{
David W. Bapst
}
\references{
Foote, M. 1996 On the Probability of Ancestors in the
Fossil Record. \emph{Paleobiology} \bold{22}(2):141--151.

Hartmann, K., D. Wong, and T. Stadler. 2010 Sampling Trees
from Evolutionary Models. \emph{Systematic Biology}
\bold{59}(4):465--476.

Kendall, D. G. 1948 On the Generalized "Birth-and-Death"
Process. \emph{The Annals of Mathematical Statistics}
\bold{19}(1):1--15.

Nee, S. 2006 Birth-Death Models in Macroevolution.
\emph{Annual Review of Ecology, Evolution, and Systematics}
\bold{37}(1):1--17.

Solow, A. R., and W. Smith. 1997 On Fossil Preservation and
the Stratigraphic Ranges of Taxa. \emph{Paleobiology}
\bold{23}(3):271--277.
}
\seealso{
\code{\link{simFossilTaxa_SRCond}},
\code{\link{sampleRanges}}, \code{\link{simPaleoTrees}},
\code{\link{taxa2phylo}}, \code{\link{taxa2cladogram}}
}

