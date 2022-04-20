#' Obtaining Data for Taxa or Occurrences From Paleobiology Database API
#' 
#' The Paleobiology Database API (\href{http://paleobiodb.org/data1.2}{link})
#' is very easy to use, and generally any data one wishes to collect can be obtained
#' in R through a variety of ways - the simplest being to wrap a data retrieval request
#' to the API, specified for CSV output, with R function \code{read.csv}. The functions
#' listed here, however, are some simple helper functions for doing tasks common to
#' users of this package - downloading occurrence data, or taxonomic information,
#' for particular clades, or for a list of specific taxa.

#' @details
#' In many cases, it might be easier to write your own query - these
#' functions are only made to make getting data for some very specific
#' applications in \code{paleotree} easier.

#' @param taxon A single name of a of a higher taxon which you wish to catch
#' all taxonomic 'children' (included members - i.e. subtaxa) of, from 
#' within the Paleobiology Database.

#' @param taxa A character vector listing taxa of interest that the user
#' wishes to download information on from the Paleobiology Database.
#' Multiple taxa can be listed as a single character string, with desired taxa
#' separated by a comma with no whitespace (ex.
#' \code{"Homo,Pongo,Gorilla"}) or as a vector of character strings
#' (ex. \code{c("Homo", "Pongo", "Gorilla")}),
#' which will then formatted for use in the API call.
    
#' @param showTaxa Which variables for taxonomic data should be requested
#' from the Paleobiology Database? The default is to include classification (\code{"class"}),
#' parent-child taxon information (\code{"parent"}), information on each
#' taxon's first and last appearance (\code{"app"}), information on the
#' PhyloPic silhouette images assigned to that taxon (\code{"img"}), and
#' the names of those who entered and authorized the taxonomic data you
#' have downloaded (\code{"entname"}). Multiple variable blocks can be
#' given as a single character string, with desired variable selections
#' separated by a comma with no whitespace (ex. \code{"class,img,app"}) or 
#' as a vector of character strings (ex. \code{c("class", "img", "app")}),
#' which will then formatted for use in the API call. Other options that
#' you might want to include, such as information on ecospace or taphonomy,
#' can be included: please refer to the full list at
#' the \href{http://paleobiodb.org/data1.2/taxa/list_doc.htm}{documentation for the API}.

#' @param showOccs Which variables for occurrence data should be requested
#' from the Paleobiology Database? The default is to include classification (\code{"class"}),
#' classification identifiers (\code{"classext"}), genus and subgenus
#' identifiers (\code{"subgenus"}), and species-level identifiers (\code{"ident"}).
#' Multiple variable blocks can be given as a single character string, with desired
#' variable selections separated by a comma with no whitespace (ex.
#' \code{"class,subgenus,ident"}) or  as a vector of character strings
#' (ex. \code{c("class", "subgenus", "ident")}), which will then formatted for
#' use in the API call. For full list of other options that you might want to include, please refer
#' to \href{https://paleobiodb.org/data1.2/occs/list_doc.html}{documentation for the API}.

#' @param status What taxonomic status should the pull taxa have? 
#' The default is \code{status = "accepted"}, which means 
#' only those taxa that are both valid taxa and 
#' \emph{the accepted senior homonym}. Other typical statuses
#' to consider are \code{"valid"}, which is all valid taxa:
#' senior homonyms and valid subjective synonyms, and \code{"all"},
#' which will return all valid taxa and all otherwise repressed invalid taxa.
#' For additional statuses that you can request, please see the documentation at  
#' the \href{http://paleobiodb.org/data1.2/taxa/list_doc.htm}{documentation for the API}.
    
    # status -> all, accepted, valid
    # accepted -> only senior synonyms
    # valid -> snior synonyms + valid subjective synonyms
    # all -> valid taxa + repressed invalid taxa

#' @param urlOnly If \code{FALSE} (the default), then the
#' function behaves as expected, the API is called and a
#' data table pulled from the Paleobiology Database is returned.
#' If \code{urlOnly = TRUE}, the URL of the API call is returned
#' instead as a character string. 

#' @param stopIfMissing If some taxa within the requested set appear
#' to be missing from the Paleobiology Database's taxonomy table, should
#' the function halt with an error?

#' @param failIfNoInternet If the Paleobiology Database or another 
#' needed internet resource cannot be accessed, perhaps because of
#' no internet connection, should the function fail (with an error)
#' or should the function return \code{NULL} and return an
#' informative message instead, thus meeting the CRAN policy
#' that such functionalities must 'fail gracefully'?
#' The default is \code{TRUE} but all examples that might be auto-run
#' use \code{FALSE} so they do not fail during R CHECK.

#' @return 
#' These functions return a \code{data.frame} containing
#' variables pulled for the requested taxon selection.
#' This behavior can be modified by argument \code{urlOnly}.

#' @name getDataPBDB

#' @aliases getCladeTaxaPBDB getSpecificTaxaPBDB

#' @seealso 
#' See \code{\link{makePBDBtaxonTree}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhyloPicTree}} for functions that use taxonomic data.
#' Occurrence data is sorted by taxon via \code{\link{taxonSortPBDBocc}},
#' and further utilized \code{\link{occData2timeList}} and  \code{\link{plotOccData}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.


#' @examples
#' \donttest{
#' # Note that all examples here use argument 
#'     # failIfNoInternet = FALSE so that functions do
#'     # not error out but simply return NULL if internet
#'     # connection is not available, and thus
#'     # fail gracefully rather than error out (required by CRAN).
#' # Remove this argument or set to TRUE so functions fail
#'     # when internet resources (paleobiodb) is not available.
#' 
#' #graptolites
#' graptData <- getCladeTaxaPBDB("Graptolithina", 
#'     failIfNoInternet = FALSE)
#' dim(graptData)
#' sum(graptData$taxon_rank == "genus")
#' 
#' # so we can see that our call for graptolithina returned 
#'     # a large number of taxa, a large portion of which are
#'     # individual genera
#' # (554 and 318 respectively, as of 03-18-19)
#' 
#' tetrapodList<-c("Archaeopteryx", "Columba", "Ectopistes",
#'    "Corvus", "Velociraptor", "Baryonyx", "Bufo",
#'    "Rhamphorhynchus", "Quetzalcoatlus", "Natator",
#'    "Tyrannosaurus", "Triceratops", "Gavialis",
#'    "Brachiosaurus", "Pteranodon", "Crocodylus",
#'    "Alligator", "Giraffa", "Felis", "Ambystoma",
#'     "Homo", "Dimetrodon", "Coleonyx", "Equus",
#'    "Sphenodon", "Amblyrhynchus")
#' 
#' tetrapodData <-getSpecificTaxaPBDB(tetrapodList, 
#'     failIfNoInternet = FALSE)
#' dim(tetrapodData)
#' sum(tetrapodData$taxon_rank == "genus")
#' # should be 26, with all 26 as genera
#' 
#' #############################################
#' # Now let's try getting occurrence data
#' 
#' # getting occurrence data for a genus, sorting it
#' # Dicellograptus
#' dicelloData <- getPBDBocc("Dicellograptus", 
#'     failIfNoInternet = FALSE)
#' dicelloOcc2 <- taxonSortPBDBocc(dicelloData, 
#'     rank = "species", onlyFormal = FALSE, 
#'     failIfNoInternet = FALSE)
#' names(dicelloOcc2)
#' 
#' }
#' 


# getPBDBocc 
    # occs documentation
    # https://paleobiodb.org/data1.2/occs/list_doc.html
    #
    # simple function for getting occurrence data from API v1.2 
             #cleans PBDB occurrence downloads of warnings


#' @rdname getDataPBDB 
#' @export 
getCladeTaxaPBDB <- function(
        taxon,
        showTaxa = c("class", "parent", "app", "img", "entname"),
        status = "accepted", 
        urlOnly = FALSE, 
        stopIfMissing = FALSE,
        failIfNoInternet = TRUE
        ){
    ##########################################
    # check that only a single taxon is given
    if(length(taxon) != 1){
        stop(paste0(
            "Input 'taxon' should only be a character string of length one:\n",
            "   A single name of a higher taxon which you wish to catch all children of"
            ))
        }
    ###################################
    # 12-30-18: modified for API version 1.2
    #let's get some taxonomic data
    requestURLPBDB <- paste0("http://paleobiodb.org/",
            "data1.2/taxa/list.txt?base_name=", taxon,
            "&show=",paste0(showTaxa,collapse = ","),
            # status -> all, accepted, valid
                # accepted -> only senior synonyms
                # valid -> snior synonyms + valid subjective synonyms
                # all -> valid taxa + repressed invalid taxa
            "&taxon_status=",status
            )
    if(urlOnly){
        res <- requestURLPBDB
    }else{
        res <- getPBDBtaxaCSV(
            requestURL = requestURLPBDB, 
            stopIfMissing = stopIfMissing,
            failIfNoInternet = failIfNoInternet
            )
        if(is.null(res)){return(NULL)}
        }
    #######################################
    return(res)
    }

#' @rdname getDataPBDB 
#' @export 
getSpecificTaxaPBDB <- function(
        taxa,
        showTaxa = c("class", "parent", "app", "img", "entname"),
        status = "accepted",
        urlOnly = FALSE, 
        stopIfMissing = FALSE,
        failIfNoInternet = TRUE
        ){
    #####################################
    # test for duplicates, remove and report them
    if(length(taxa) > length(unique(taxa))){
        duplicates <- taxa[duplicated(taxa)]
        duplicates <- paste0(duplicates, collapse = ", ")
        warning("Duplicated taxa found in input:\n   ", duplicates)
        taxa <- unique(taxa)
        }
    #
    if(length(taxa)>1){
        # collapse taxa to a vector
        taxaMerged <- paste0(taxa,
            collapse=",")        
    }else{
        taxaMerged <- taxa
        }
    #
    requestURLPBDB <- paste0(
        "http://paleobiodb.org/data1.2/taxa/list.txt?name=",
            taxaMerged,
        "&show=",
            paste0(showTaxa,collapse = ","),
        # status -> all, accepted, valid
            # accepted -> only senior synonyms
            # valid -> snior synonyms + valid subjective synonyms
            # all -> valid taxa + repressed invalid taxa
        "&taxon_status=", status
        )
    if(urlOnly){
        res <- requestURLPBDB
        # browseURL(apiAddressTaxa)
    }else{
        res <- getPBDBtaxaCSV(
            requestURL = requestURLPBDB, 
            stopIfMissing = stopIfMissing,
            failIfNoInternet = failIfNoInternet
            )
        if(is.null(res)){return(NULL)}
        }
    ################
    # report taxa not found in records 
    if(nrow(res) < length(taxa)){
        warning(paste0(
            "Fewer taxon records returned (",
                nrow(res),
                ") than the number of unique taxon names from input (",
                length(taxa),
                ").\n",
            "  This may reflect missing taxa, or empty records."
            ))
        if(stopIfMissing){
            stop(
                "Terminated as stopIfMissing = TRUE, and number of unique input taxa differs from the number of taxon records output."
                )
            }
        }
    ##################################
    # return
    return(res)
    }    
    

getPBDBtaxaCSV <- function(
        requestURL, 
        stopIfMissing = FALSE,
        failIfNoInternet = TRUE
        ){
    ##########################################
    # first test internet
    testConnect <- canConnectPBDB(fail = failIfNoInternet)    
    if(is.null(testConnect)){
        return(NULL)
        }
    #
    linesOut <- readLines(requestURL)
    ###############
    # stop if contains no records returned
    if(any(linesOut == "\"THIS REQUEST RETURNED NO RECORDS\"")){
        stop(paste0(
            "For this query, the PBDB API returned:\n    ",
            "\"THIS REQUEST RETURNED NO RECORDS\"",
            "\n Presumably this means a taxon is missing/orphaned?"
            ))
        }
    #
    # catch any warnings from the PBDB
    if(any(linesOut == "\"Records:\"")){
        # then there must be errors, report them
        isWarning <- startsWith(linesOut, 
            prefix = "\"Warning:\",")
        for(i in which(isWarning)){
            errorMsg <- linesOut[i]
            errorMsg <- gsub("\"Warning:\",\"",
                "", errorMsg)
            warning(errorMsg)
            }
        if(stopIfMissing){
            stop(paste0(
                "Some taxa not found in PBDB taxonomy table,",
                " see warnings and check API request URL: \n",
                requestURL
                ))
            }
        lineRemove <- c(which(isWarning),which(linesOut == "\"Records:\""))
        linesOut<-linesOut[-lineRemove]
        }
    # first test internet
    testConnect <- canConnectPBDB(fail = failIfNoInternet)    
    if(is.null(testConnect)){
        return(NULL)
        }
    #
    res <- read.csv(text = linesOut,
        stringsAsFactors = FALSE)
    ###############
    return(res)
    }    


#' @rdname getDataPBDB 
#' @export 
getPBDBocc <- function(
        taxa, 
        showOccs = c("class", "classext", "subgenus", "ident", "entname"),
        failIfNoInternet = TRUE
        ){
    #############
    # occs documentation
    # https://paleobiodb.org/data1.2/occs/list_doc.html
    #
    # simple function for getting occurrence data from API v1.2 
             #cleans PBDB occurrence downloads of warnings
    taxa <- paste(taxa,collapse = ",")
    taxa <- paste(unlist(strsplit(taxa,"_")),collapse = "%20")
    showOccs <- paste(showOccs,collapse = ",")
    command <- paste0(
        "http://paleobiodb.org/data1.2/occs/list.txt?base_name=",
        taxa,"&show=",showOccs,        
        # "&limit=all", # not a command in 1.2
        collapse = "")
    command <- paste(unlist(strsplit(command,split = " ")),
        collapse = "%20")
    # first test internet
    testConnect <- canConnectPBDB(fail = failIfNoInternet)
    if(is.null(testConnect)){
        return(NULL)
        }
    #    
    downData <- readLines(command)
    if(length(grep("Warning",downData)) != 0){
        start <- grep("Records",downData)
        warn <- downData[1:(start-1)]
        warn <- sapply(warn, function(x) 
            paste0(unlist(strsplit(unlist(strsplit(x,'"')),",")),
                 collapse = ""))
        warn <- paste0(warn,collapse = "\n")
        names(warn) <- NULL
        mat <- downData[-(1:start)]
        # first test internet
        #testConnect <- canConnectPBDB(fail = failIfNoInternet)
        #if(is.null(testConnect)){
        #    return(NULL)
        #    }
        #
        mat <- read.csv(textConnection(mat))
        message(warn)
    }else{
        mat <- downData
        # first test internet
        #testConnect <- canConnectPBDB(fail = failIfNoInternet)
        #if(is.null(testConnect)){
        #    return(NULL)
        #    }
        #
        mat <- read.csv(textConnection(mat))
        }
    return(mat)
    }
