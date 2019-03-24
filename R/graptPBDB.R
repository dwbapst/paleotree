#' Example Occurrence and Taxonomic Datasets of the Graptolithina from the Paleobiology Database
#' 
#' Example datasets consisting of (a) occurrence data and (b) taxonomic data
#' downloaded from the Paleobiology Database API for the Graptolithina.
#' In order to make sure to catch anything that might be considered a graptolite,
#' the actual taxon searched for was the Pterobranchia, the larger clade
#' that includes graptolites within it (Mitchell et al., 2013).

#' @name graptPBDB
#' @rdname graptPBDB
#' @aliases graptPBDB graptOccPBDB graptTaxaPBDB

#' @details
#' This example PBDB data is included here for testing
#' functions involving occurrence data and taxonomy
#' in \code{paleotree}.

#' @format 
#' The example occurrence dataset (\code{graptOccPBDB}) is a
#' \code{data.frame} consisting of 5900 occurrences (rows) and 35 variables (columns).
#' The example taxonomy dataset (\code{graptTaxaPBDB}) is a
#' \code{data.frame} consisting of 364 formal taxa (rows) and 53 variables (columns).
#' Variables are coded in the 'pbdb' vocabulary of the PBDB API v1.2.

#' @seealso
#' \code{\link{taxonSortPBDBocc}}, \code{\link{occData2timeList}},
#' \code{\link{makePBDBtaxonTree}}, \code{\link{plotOccData}}

#' @references
#' Mitchell, C. E., M. J. Melchin, C. B. Cameron, and J. Maletz.
#' 2013. Phylogenetic analysis reveals that Rhabdopleura
#' is an extant graptolite. \emph{Lethaia} 46(1):34-56.
#' 
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.

#' @source 
#' See examples for the full R code used to obtain the data from the API.
#' You can find the Paleobiology Database at \url{http://paleobiodb.org}
#' 
#' The occurrence data was entered by (in order of relative portion) P. Novack-Gottshall, M. Krause, M. Foote,
#' A. Hendy, T. Hanson, M. Sommers and others. This same data was authorized mainly by A. Miller,
#' W. Kiessling, M. Foote, A. Hendy, S. Holland, J. Sepkoski and others.

#' @keywords datasets

#' @docType data

#' @examples
#' 
#' 
#' # let's look for pterobranch genera
#'    # pterobranchs are the larger group containing graptolites
#' 
#' taxon <- "Pterobranchia"
#' selectRank <- "genus"
#' 
#' \dontrun{
#' 
#' library(paleotree)
#' 
#' # get taxon data
#' 	# default variables
#' 
#' graptTaxaPBDB<-getCladeTaxaPBDB(taxon)
#' 
#' graptTree <- makePBDBtaxontree(graptTaxaPBDB,
#' 				rank = selectRank,
#' 				cleanDuplicate = TRUE)
#' 
#' # date the tree using the ranges
#' 	# provided directly by the PBDB
#' graptTimeTree <- dateTaxonTreePBDB(graptTree)
#' 
#' library(strap)
#' geoscalePhylo(graptTimeTree, 
#' 	ages=graptTimeTree$ranges.used)
#' nodelabels(graptTimeTree$node.label,
#' 	cex=0.7,
#' 	adj=c(0.3,0)
#' 	)
#' 
#' # we could also date the tree using the occurrence data
#' 	# default variables
#' 
#' graptOccPBDB <- getPBDBocc(taxon)
#' 
#' graptOccSort <- taxonSortPBDBocc(graptOccPBDB, 
#' 			rank = selectRank,
#' 			onlyFormal = FALSE, 
#' 			cleanUncertain = FALSE)
#' 
#' graptTimeList <- occData2timeList(occList = graptOccSort)
#' 
#' graptTimeTreeFromOcc <- bin_timePaleoPhy(
#' 	graptTree,
#' 	timeList = graptTimeList,
#' 	nonstoch.bin = TRUE,
#' 	type = "mbl",
#' 	vartime = 3)
#'   
#' library(strap)
#' geoscalePhylo(
#' 	graptTimeTreeFromOcc, 
#' 	ages = graptTimeTreeFromOcc$ranges.used)
#' nodelabels(
#' 	graptTimeTreeFromOcc$node.label,
#' 	cex = 0.7,
#' 	adj = c(0.3,0)
#' 	)
#' 
#' save(graptOccPBDB,
#' 		graptTaxaPBDB,
#' 		graptTree,
#'		graptTimeTree,
#'		file = "graptPBDB.rdata")
#' 
#' }
#' 
#' # load archived example data
#' data(graptPBDB)
#' 
#' # let's visualize who entered the majority of the occurrence data
#' pie(sort(table(graptOccPBDB$enterer)))
#' # and now who authorized it
#' pie(sort(table(graptOccPBDB$authorizer)))
#' 
#' # I *sort of* apologize for using pie charts.
#' 
#' # Let's look at age resolution of these occurrences
#' hist(graptOccPBDB$early_age - graptOccPBDB$late_age,
#' 		main = "Age Resolution of Occurrences",
#'		xlab = "Ma")
#' 
#' #distribution of taxa among taxonomic ranks
#' table(graptTaxaPBDB$taxon_rank)
#' barplot(table(graptTaxaPBDB$taxon_rank))
#' 
NULL