#' Obtaining Data for Sets of Taxa From Paleobiology Database API
#' 
#' The Paleobiology Database API () is very easy to use, and generally any data one wishes to collect

#' @details

#' @param

#' @return

#' @name getTaxaDataPBDB

#' @aliases getCladeTaxaPBDB getSpecificTaxonTreePBDB

#' @seealso
#' See \code{\link{getTaxaDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhylopicTreePBDB}}.


#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.


#' @examples
#' 
#' #an example here




# status -> all, accepted, valid
	# accepted -> only senior synonyms
	# valid -> snior synonyms + valid subjective synonyms
	# all -> valid taxa + repressed invalid taxa


#' @param taxon A single name of a of a higher taxon which you wish to catch
#' all taxonomic 'children' (included members - i.e. subtaxa) of, from 
#' within the Paleobiology Database.
	

#' @rdname getTaxaDataPBDB
#' @export
getCladeTaxaPBDB <- function(taxon,
		show = c("class","img","app","parent"),
		status = "accepted"){
	##########################################
	# check that only a single taxon is given
	if(length(taxon) != 1){
		stop("taxon should only be a single name of a higher taxon which you wish to catch all children of")
		}
	###################################
	# 12-30-18: modified for API version 1.2
	#let's get some taxonomic data
	taxaData <- read.csv(
		paste0("http://paleobiodb.org/",
			"data1.2/taxa/list.txt?base_name=", taxon,
			"&show=",paste0(show,collapse = ","),
			# status -> all, accepted, valid
				# accepted -> only senior synonyms
				# valid -> snior synonyms + valid subjective synonyms
				# all -> valid taxa + repressed invalid taxa
			"&taxon_status=",status
			),
	    stringsAsFactors = FALSE)
	#######################################
	return(taxaData)
	}

#' @rdname getTaxaDataPBDB
#' @export
getSpecificTaxonTreePBDB <- function(taxa,
		PDFfile = "myPBDBphylogeny.pdf",
		show = c("class","img","app","parent"),
		status = "accepted",
		addressOnly = FALSE){
	#####################################
	if(length(taxa)>1){
		# collapse taxa to a vector
		taxa <- paste0(taxa,
			collapse=",")		
		}
	apiAddressTaxa <- paste0(
		"http://paleobiodb.org/data1.2/taxa/list.txt?name=",taxa,
		"&show=",paste0(show,collapse = ","),
		# status -> all, accepted, valid
			# accepted -> only senior synonyms
			# valid -> snior synonyms + valid subjective synonyms
			# all -> valid taxa + repressed invalid taxa
		"&taxon_status=",status
		)
	# browseURL(apiAddressTaxa)
	#####################################
	taxaData <- read.csv(apiAddressTaxa,
		stringsAsFactors = FALSE)
	##################################
	# return
	return(taxaData)
	}	




