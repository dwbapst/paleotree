# hidden utility functions for getting phylopics

getPhyloPicPNG<-function(
		picID_PBDB, 
		cacheDir = NULL,
		cacheImage = FALSE,
        failIfNoInternet = TRUE
		){
	####################################################
	# first try to find and load a cached version
	# if that doesn't work
		# try to load from phylopic using PBDB UID
	# if that doesn't work
		# try to load the image from PBDB		
	notCached <- FALSE
	picPNG <- NULL
	###########################################
	# First try to get a cached version
	if(!is.null(cacheDir)) {
		cachePath <- file.path(cacheDir,
					paste0(picID_PBDB, ".png")
					)
		if(file.exists(cachePath)){
			picPNG <- png::readPNG(cachePath)
		}else{
			notCached <- TRUE
			}
		}
	##################################################
	# if that doesn't work
		# try to load from phylopic using PBDB UID
	if(is.null(picPNG)){
		picUID <- paste0("https://paleobiodb.org/data1.2/taxa/thumb.txt?id="
							,picID_PBDB)
	    # first test internet
        testConnect <- canConnectPBDB(fail = failIfNoInternet)
        if(is.null(testConnect)){
            return(NULL)
            }
        #
		picUID <- read.csv(picUID, stringsAsFactors = FALSE)
		picUID <- picUID$uid
		picPNG <- getPhyloPicFromPhyloPic(picUID)
		}
	################################################
	# if that doesn't work
		# try to load the image from PBDB		
	if(is.null(picPNG)){
		picPNG <- getPhyloPicPNG_PBDB(
		        picID_PBDB = picID_PBDB,
		        failIfNoInternet = failIfNoInternet
            	)
        if(is.null(picPNG) & failIfNoInternet){
            return(NULL)
            }
		}
	#########################
	if(cacheImage & notCached){
		if(!dir.exists(cacheDir)){
			dir.create(
				cacheDir,
				showWarnings = FALSE
				)
			}
		#
		png::writePNG(picPNG,
			target=file.path(
				cacheDir,
				paste0(picID_PBDB, ".png")
				)
			)
		}
	#########
	return(picPNG)	
	}		

getPhyloPicUIDsTableFromPBDB <- function(
        picIDs,
        failIfNoInternet = TRUE
        ){
    # first test internet
    testConnect <- canConnectPBDB(fail = failIfNoInternet)
    if(is.null(testConnect)){
        return(NULL)
        }
    #
	imgNumAPIurl <- "https://paleobiodb.org/data1.2/taxa/thumb.txt?id="
	URLwithNums <- paste0(imgNumAPIurl, picIDs)	
	names(URLwithNums) <- names(picIDs)
	res <- do.call(rbind,
          lapply(URLwithNums, read.csv, stringsAsFactors=FALSE)
		  )
	return(res)
	}
	
getPhyloPicFromPhyloPic <- function(picUID){
	# get image info
	picInfoURL <- paste0("https://phylopic.org/api/a/image/",
		picUID,"?options=credit+licenseURL+pngFiles")
	if(RCurl::url.exists(picInfoURL)){
		picInfo <- jsonlite::fromJSON(picInfoURL,  
			simplifyVector = FALSE)
		#################
		#
		if(length(picInfo$result$pngFiles)>0){
			#need to check if there is a URL at all
			picPNGurl <- picInfo$result$pngFiles[[
					length(picInfo$result$pngFiles)
				]]$url
			picPNGurl <- paste0("https://phylopic.org",
				picPNGurl)
			##############
			if(RCurl::url.exists(picPNGurl)){
				#need to check if there is a PNG at that URL
				pngContent <- RCurl::getURLContent(picPNGurl)
				pngContent  <-  png::readPNG(pngContent)
			}else{
				pngContent <- NULL
				}
		}else{
			pngContent <- NULL
			}
	}else{
		pngContent <- NULL
		}
	return(pngContent)
	}

	
getPhyloPicPNG_PBDB<-function(
		picID_PBDB,
        failIfNoInternet = TRUE
        ){
	############################################
	#	
	# require(png);require(RCurl)
	# png::readPNG RCurl::getURLContent
	#
	# pause 1 second
	Sys.sleep(runif(1,1,2))
	#
	# GET IMAGE
	# get the URL address for the pic via API
	apiPicURL <- paste0(
		"https://paleobiodb.org/data1.2/taxa/thumb.png?id=",
		picID_PBDB)
	#
	# first test internet
    testConnect <- canConnectPBDB(fail = failIfNoInternet)
    if(is.null(testConnect)){
        return(NULL)
        }
	# get picPNG		
	picPNG <-  png::readPNG(RCurl::getURLContent(apiPicURL))
	############
	# phylopic PNG is on the fourth slice
	#  image(picPNG [,,4])
	########################
	return(picPNG)
	}
	
getPhyloPicIDNumFromPBDB <- function(
        taxaData, 
        tree, 
        failIfNoInternet = TRUE
        ){
	# check or obtain the PBDB phylopic IDs from PBDB
		# get the phylopic-specific IDs  as well
	#
	if(is.null(taxaData)){
		# get image ID numbers using PBDB API calls for each
			# tip taxon in the tree using the tip labels
		tiptaxa <- paste0(tree$tip.label, 
			collapse = ",")
		apiAddressTaxa <- paste0(
			"https://paleobiodb.org/data1.2/taxa/list.txt?name=",
			tiptaxa, "&rel=exact&show=img"
			)	
		# call PBDB API
        # first test internet
        testConnect <- canConnectPBDB(fail = failIfNoInternet)
        if(is.null(testConnect)){
            return(NULL)
            }
		taxaData <- read.csv(apiAddressTaxa,
			stringsAsFactors = FALSE)
		# get the PBDB image IDs and label with tip labels
		phylopicIDsPBDB<- taxaData$image_no[
			match(tree$tip.label, taxaData$taxon_name)
			]
		names(phylopicIDsPBDB) <- tree$tip.label
	}else{
		#
		# get the PBDB image IDs and label with tip labels
		phylopicIDsPBDB <- taxaData$image_no[
			match(tree$tip.label, taxaData$taxon_name)
			]
		phylopicIDsPBDB <- as.character(phylopicIDsPBDB)
		names(phylopicIDsPBDB) <- tree$tip.label
		#
		# CHECKS
		# does the provided number of IDs equal the number
			# of tips?
		if(length(phylopicIDsPBDB) != Ntip(tree)){
			#print(tree)
			#print(phylopicIDsPBDB)
			stop(paste0(
				"If provided, phylopicIDsPBDB must be have the same",
				" length as the number of tips on the input tree"))
			}
		# make sure its a character vector
		if(!is.character(phylopicIDsPBDB)){
			phylopicIDsPBDB <- as.character(phylopicIDsPBDB)
			if(!is.character(phylopicIDsPBDB)){		
				stop("Cannot coerce phylopicIDsPBDB to be type character")
				}
			}
		# Are there NAs or blanks?
		if(any(is.na(phylopicIDsPBDB))){
			stop("If provided, phylopicIDsPBDB cannot contain any NA values")
			}
		if(any(phylopicIDsPBDB == "")){
			stop("If provided, phylopicIDsPBDB cannot contain any blank values")
			}
		# make sure user provided IDs are labeled
		if(length(names(phylopicIDsPBDB)) != Ntip(tree)){
			#print(length(names(phylopicIDsPBDB)))
			#print(phylopicIDsPBDB)
			stop(paste0("If provided, phylopicIDsPBDB must have",
				" labels matching the tree's tip labels"))	
			}
		}
	# make sure phylopicIDsPBDB is sorted relative to tip order
	phylopicIDsPBDB <- phylopicIDsPBDB[
		match(tree$tip.label, names(phylopicIDsPBDB))
		]
	return(phylopicIDsPBDB)
	}

# test if there is internet
canConnectPBDB <- function(fail = TRUE){
    # ip = "8.8.8.8"
    # is PBDB up
    res <- RCurl::url.exists("https://paleobiodb.org/")
    if(!res){ 
        connectMessage <- "Cannot connect to Paleobiology Database at https://paleobiodb.org/"
        if(fail){
            stop(connectMessage)
            }
        message(connectMessage)
        res <- NULL
        return(NULL)
        }
    #
    # is PBDB data service up
    res <- RCurl::url.exists(
        "https://paleobiodb.org/data1.2/taxa/single.txt?name=Dicellograptus/"
        )
    if(!res){
        connectMessage <- "Cannot connect to Paleobiology Database API at https://paleobiodb.org/data1.2/"
        if(fail){
            stop(connectMessage)
            }
        message(connectMessage)
        return(NULL)
        }
    
    # res <- try(read.csv(
    #    "https://paleobiodb.org/data1.2/taxa/single.txt?name=Dicellograptus",
    #    stringsAsFactors = TRUE
    #    ))
    return(TRUE)
    }
	
	