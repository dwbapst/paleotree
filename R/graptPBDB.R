#' Example Occurrence Dataset of the Graptolithina from the Paleobiology Database

#' An example dataset consisting of occurrence data downloaded from the Paleobiology Database API
#' for the Graptolithina.

#' @name graptPBDB
#' @rdname graptPBDB
#' @aliases graptPBDB graptOccPBDB

#' @details
#' This example PBDB data is included here for testing functions involving occurrence data and taxonomy
#' in \code{paleotree}.

#' @format 
#' The example dataset is a data.frame consisting of 5900 occurrences (rows) and 35 variables (columns). Variables are
#' coded in the 'pbdb' vocabulary of the PBDB API v1.1.




# @seealso


#' @source 
#' This dataset was downloaded on March 19th, 2015, from the Paleobiology Database
#' using the API version 1.1, using the following API command:
#' 
#' \code{http://paleobiodb.org/data1.1/occs/list.txt?base_name=Graptolithina&show=ident,phylo,entname&limit=all}
#' And was read into R using the command \code{read.csv}.
#'
#' You can find the Paleobiology Database at http://paleobiodb.org
#' 
#' This data was entered by (in order of relative portion) P. Novack-Gottshall, M. Krause, M. Foote,
#' A. Hendy, T. Hanson, M. Sommers and others. This data was authorized mainly by A. Miller,
#' W. Kiessling, M. Foote, A. Hendy, S. Holland, J. Sepkoski and others.
#'
#' See examples for the full R code used.

#' @keywords datasets

#' @docType data

#' @examples
#'
#' \dontrun{
#' 
#' #original code used to obtain this dataset on March 19th, 2015
#' graptOccPBDB<-read.csv(paste0(
#' 	"http://paleobiodb.org/data1.1/occs/list.txt?",
#' 		"base_name=Graptolithina&show=ident,phylo,entname&limit=all"))	
#' save(graptOccPBDB,file="graptPBDB.Rdata")
#'
#' }
#'
#' # load archived example data
#' data(graptPBDB)
#'
#' # let's visualize who entered the majority of this data
#' pie(sort(table(graptOccPBDB$authorizer)))
#' # and now who authorized it
#' pie(sort(table(graptOccPBDB$enterer)))
#' # I apologize for using pie charts.
#' 
#' # Let's look at age resolution of these occurrences
#' hist(graptOccPBDB$early_age-graptOccPBDB$late_age,
#'		main="Age Resolution of Occurrences", xlab="Ma")
#' 
#' 
#'
NULL