#' Simulating Taxa in the Fossil Record, Conditioned on Sampling Rate
#' 
#' Function for simulating taxon ranges and relationships under various models
#' of evolution, conditioned on creating clades large enough to sample some
#' average number of taxa at some sampling rate
#' 
#' simFossilTaxa_SRCond is a wrapper for simFossilTaxa for when clades of a
#' particular size are desired, post-sampling. simFossilTaxa simulates a
#' birth-death process (Kendall, 1948; Nee, 2006), but unlike most functions
#' for this implemented in R, this function enmeshes the simulation of
#' speciation and extinction with explicit models of how lineages are
#' morphologically differentiated, as morphotaxa are the basic units of
#' paleontological estimates of diversity and phylogenetics. For more details
#' on the workings of simFossilTaxa and many of the arguments involved, please
#' see the help file for simFossilTaxa (\code{\link{simFossilTaxa}}).
#' simFossilTaxa_SRCond first calculates the expected proportion of taxa
#' sampled, given the sampling rate and the rates which control lineage
#' termination: extinction, anagenesis and bifurcation. The average original
#' clade size needed to produce, on average, a given number of sampled taxa
#' (the argument 'avgtaxa') is calculated with the following equation:
#' 
#' \eqn{N = avgtaxa / (r / (r + (q + anag.rate + (prop.bifurc * p)) ))}
#' 
#' We will call that quantity N. Note that the quantity (prop.bifurc * p)
#' describes the rate of bifurcation when there is no cryptic cladogenesis, as
#' prop.bifurc is the ratio of budding to bifurcating cladogenesis. This
#' equation will diverge in ways that are not easily predicted as the rate of
#' cryptic speciation increases. Note, as of version 1.5, this equation was
#' altered to the form above. The previous form was similar and at values of
#' avgtaxa greater than 10 or so, produces almost identical values. The above
#' is preferred for its relationship to taxonomic completeness (Solow and
#' Smith, 1997).
#' 
#' Next, this value is used with simFossilTaxa, where mintaxa is set to N and
#' maxtaxa set to 2*N. simFossilTaxa_SRcond will generally produce simulated
#' datasets that are generally of that size or larger post-sampling (although
#' there can be some variance). Some combinations of parameters may take an
#' extremely long time to find large enough datasets. Some combinations may
#' produce very strange datasets that may have weird structure that is only a
#' result of the conditioning (for example, the only clades that have many taxa
#' when net diversification is low or negative will have lots of very early
#' divergences, which could impact analyses). Needless to say, conditioning can
#' be very difficult.
#' 
#' @param r Instantaneous sampling rate per time unit.
#' @param avgtaxa Desired average number of taxa.
#' @param p Instantaneous rate of speciation/branching.
#' @param q Instantaneous rate of extinction.
#' @param anag.rate Instantaneous rate of pseudoextinction/anagensis.
#' @param prop.bifurc Proportion of morphological branching by bifurcating
#' cladogenesis relative to budding cladogenesis.
#' @param prop.cryptic Proportion of cryptic speciation by relative to
#' morphological branching, such as bifurcating and budding.
#' @param nruns Number of datasets to be output.
#' @param mintime Minimum time units to run any given simulation before
#' stopping.
#' @param maxtime Maximum time units to run any given simulation before
#' stopping.
#' @param minExtant Minimum number of living taxa allowed at end of
#' simulations.
#' @param maxExtant Maximum number of living taxa allowed at end of
#' simulations.
#' @param count.cryptic If TRUE, cryptic taxa are counted as seperate taxa for
#' conditioning.
#' @param print.runs If TRUE, prints the proportion of simulations accepted for
#' output to the terminal.
#' @param plot If TRUE, plots the diversity curves of accepted simulations.
#' @return This function gives back a list containing nruns number of taxa
#' datasets, where each element is a matrix. If nruns=1, the output is not a
#' list but just a single matrix. Sampling has not been simulated in the output
#' for either function; the output represents the 'true' history of the
#' simualted clade.
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
#' @seealso \code{\link{simFossilTaxa}}, \code{\link{sampleRanges}},
#' \code{\link{simPaleoTrees}}, \code{\link{taxa2phylo}},
#' \code{\link{taxa2cladogram}}
#' @references Kendall, D. G. 1948 On the Generalized "Birth-and-Death"
#' Process. \emph{The Annals of Mathematical Statistics} \bold{19}(1):1--15.
#' 
#' Nee, S. 2006 Birth-Death Models in Macroevolution. \emph{Annual Review of
#' Ecology, Evolution, and Systematics} \bold{37}(1):1--17.
#' 
#' Solow, A. R., and W. Smith. 1997 On Fossil Preservation and the
#' Stratigraphic Ranges of Taxa. \emph{Paleobiology} \bold{23}(3):271--277.
#' @examples
#' 
#' set.seed(444)
#' 
#' avgtaxa <- 50
#' r <- 0.5
#' 
#' #using the SRcond version
#' taxa <- simFossilTaxa_SRCond(r=r,p=0.1,q=0.1,nruns=20,avgtaxa=avgtaxa)
#' #now let's use sampleRanges and count number of sampled taxa
#' ranges <- lapply(taxa,sampleRanges,r=r)
#' ntaxa <- sapply(ranges,function(x) sum(!is.na(x[,1])))
#' hist(ntaxa)
#' mean(ntaxa)
#' #works okay... some parameter combinations are difficult to get right number of taxa
#' 
#' layout(1)
#' 
#' @export simFossilTaxa_SRCond
simFossilTaxa_SRCond<-function(r,avgtaxa,p,q,anag.rate=0,prop.bifurc=0,prop.cryptic=0,nruns=1,mintime=1,maxtime=1000,
	minExtant=0,maxExtant=NULL,count.cryptic=FALSE,print.runs=FALSE,plot=FALSE){
	#wrapper for simulating clades large enough for 
		#getting some avg number of taxa under a given sampling parameter
	#for sampling: proper conditioning on number of taxa
		#how many taxa are needed in true clade to sample on average X taxa??
		#avgtaxa= average number of taxa you want to recover
		#mintaxa and maxtaxa will be set within 20% +/- of values necc for given avgtaxa
			#These should work well if avgtaxa is large (~50 or so)
	#simFossilTaxa_SRCond(r=0.1,p=0.1,q=0.1,nruns=10,avgtaxa=50,maxExtant=0)
	N<-avgtaxa/(r/(r+(q+anag.rate+(prop.bifurc*p))))
	results<-simFossilTaxa(p=p,q=q,anag.rate=anag.rate,prop.bifurc=prop.bifurc,prop.cryptic=prop.cryptic,nruns=nruns,mintaxa=N,
			maxtaxa=2*N,mintime=mintime,maxtime=maxtime,minExtant=minExtant,maxExtant=maxExtant,
			count.cryptic=FALSE,print.runs=print.runs,plot=plot)
	#if(nruns==1){results<-results[[1]]}
	return(results)
	}
