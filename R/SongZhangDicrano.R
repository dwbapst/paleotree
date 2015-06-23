#' Cladistic Data for Dicranograptid Graptolites from Song and Zhang (2014)
#'
#' Character matrix and majority-rule cladogram for 12 dicranograptid
#' (and outgroup) graptoloids, taken from Song and Zhang (2014). Included
#' here for use with functions related to character change.

#' @name SongZhangDicrano
#' @rdname SongZhangDicrano
#' @aliases charMatDicrano cladogramDicrano

#' @details
#' This example dataset is composed of a small cladistic character data for 12 taxa and 24 characters,
#' taken from Song and Zhang (2014). Note that character 22 is a biostratigraphic character, which was
#' not included in all analyses by Song and Zhang.
#'
#' The included cladogram is the majority-rule consensus of a maximum-parsimony analysis on 12 taxa with
#' 24 characters, including a biostratigraphic character. This tree is included (among the four depicted)
#' as it appeared to be the basis for the majority of Song and Zhang's discussion of dicranograptid
#' systematics.
#'
#' Both the matrix and the tree were entered by hand from their flat graphic depiction in Song and Zhang's
#' manuscript.

#' @format 
#' Loading this dataset adds two objects to the R environment.
#' \code{charMatDicrano} is a \code{data.frame} object composed of multiple factors, with \code{NA} values
#' representing missing values (states coded as '?'), read in with \code{readNexus} from package
#' \code{phylobase}. \code{cladogramDicrano} is a cladogram, formatted as a \code{phylo} class object
#' for use with package \code{ape}, without branch-lengths (as this was a consensus tree from a
#' maximum-parsimony analysis).

#' @source 
#' Song, Y., and Y. Zhang. 2014. A preliminary study on the relationship of the
#' early dicranograptids based on cladistic analysis. \emph{GFF} 136(1):243-248.

#' @keywords datasets

#' @docType data

#' @examples
#'
#' \dontrun{
#'
#'  
#'
#' require(ape)
#' require(phylobase)
#' 
#' charMatDicrano<-readNexus(file.choose(),type="data",SYMBOLS = " 0 1 2")
#' 
#' cladogramDicrano<-read.nexus(file.choose())
#' 
#' save(charMatDicrano,cladogramDicrano,file="SongZhangDicrano.rdata")
#' 
#' }
#'
#' data(SongZhangDicrano)
#'
#' # plot majority rule tree from Song and Zhang
#' plot(cladogramDicrano,
#'	main="MajRule_24charX12Taxa_wBiostratChar")
#'
#'
NULL