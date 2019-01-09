#' @details

#' @inheritParams

#' @param

#' @return

#' @name getTaxaDataPBDB

#' @aliases getCladeTaxaPBDB getSpecificTaxonTreePBDB

#' @seealso

#' @author David W. Bapst

#' @references

#' @examples




# status -> all, accepted, valid
	# accepted -> only senior synonyms
	# valid -> snior synonyms + valid subjective synonyms
	# all -> valid taxa + repressed invalid taxa




#' @rdname getTaxaDataPBDB
#' @export
getCladeTaxaPBDB <- function(taxon,
		show = c("class","img","app","parent"),
		status = "accepted"){
	##########################################
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
		plot = FALSE){
	#####################################
	taxa <- paste0(taxa,
		collapse=",")
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




