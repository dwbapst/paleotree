#' Sorting Unique Taxa of a Given Rank from Paleobiology Database Occurrence Data

#' Functions for sorting out unique taxa from Paleobiology Database occurrence downloads,
#' which should accept several different formats resulting from different versions of the
#' PBDB API and different vocabularies available from the API.



#the following nonsense is stuff from a blog post that I'm going to turn into documentation

# Alright, now that we picked a consistent vocabulary, we need to write a function that sorts occurrence data into the unique taxa of a user-selected taxonomic rank. I think the preferable way to do this is to make a list where each element is a 'unique' taxon and all the occurrences attributed to it are under that element as a table. Additionally, we'll have to allow the function to pull either just the 'formal' identified taxa or 'all' the identified taxa. This is necessary if we want species that are listed under their genus and never got assigned a species-level taxonomic ID in the PBDB. 

# Its important we consider the order of data manipulations we need to perform in this process. Formal identified taxa need to be pulled first, and then informal taxa pulled from whatever occurrences are left over. Furthermore, its possible (but hopefully rare) that the same taxon might be listed both formally and informally in a set of PBDB occurrences, and it would be preferable that the 'informal' occurrences for a taxon be concatenated to the 'formal' occurrences for that same taxon, rather than be treated as separate taxon.

# Now, the functions I used as references when writing this function were `paleobioDB`'s `pbdb_temp_range` function ([code here](https://github.com/ropensci/paleobioDB/blob/master/R/pbdb_temporal_functions.R#L64-178)), and [Matthew Clapham](http://people.ucsc.edu/~mclapham)'s `taxonClean` function ([code here](https://github.com/mclapham/PBDB-R-scripts/blob/master/taxonClean.R)). I took a somewhat different approach though. Both of these write seperate code for the different taxonomic ranks, which is vaguely unsatisfying from a programming perspective, as it could lead to inconsistencies where you fail to update all parts. Thankfully, except for species, one can get the formal taxonomic name for a given taxonomic rank by looking for the variable with the same name as the taxonomic rank. Additionally, we can look for informal occurrences of a taxon by then checking for occurrences with `taxon_rank` identical to the required rank (although it seems this is rarely important).

#And, as noted above, we need a check to convert empty `""` values to `NA`s.

#(I apologize for the terribly uneven indenting above. A strange copy/paste error involving tabs in Rstudio, Notepad++ and R GUI, and the more I tab, the worse it gets. I promise a clean-looking version in the `paleotree` github repo soon.)

# Now, in the following, I have dropped or revised some of the checks from `taxonClean` and `pbdb_temp_ranges`. 

# pbdb_temp_ranges originally (before my issue tickets) checked to see if the same taxon was listed under multiple `taxon_no` ID numbers, and returned an error if this was true. This resulted in some interesting cases where the error arose due to subtle differences in taxonomic name, like *Pseudoclimacograptus (Metaclimacograptus) hughesi* and *Pseudoclimacograptus hughesi* being listed (for whatever reason) under different `taxon_no` ID numbers. I'm blanket presuming that if a taxon's listed taxon name (e.g. formal names `matched_name`, `genus`, etc. as well as informal like `genus_name + species_name` and `taxon_name`.) are the same, its the same taxon. This assumption is tempered by the order of operations I discussed above: formal taxa are identified first, and informal taxa are identified from the occurrences that are 'leftover', and informal occurrences assigned to a taxon with a formal ID are concatenated to the 'formal' occurrences. Otherwise, we could mistakingly link the same occurrences to multiple taxa.

# For what its worth, with respect to my above example, it looks like all taxa with `taxon_name` of *Pseudoclimacograptus (Metaclimacograptus) hughesi* or *Pseudoclimacograptus hughesi* have identical `genus_name` and `species_name` variables... but not identical `matched_name` variables, as only some of these occurrences have been assigned to formal taxon *Diplograptus_hughesi*, while others are assigned to genus-level formal taxa. But that's a whole other can of worms: the inconsistent taxonomy of graptolites in the PBDB.

# taxonClean` removes taxa with question marks, *cf.* or other taxonomic-uncertainty flim-flam from the names it uses to identify unique taxa. Now there are two seperate concerns here:

# Occurrences with taxonomic uncertainy in our data. This really only pertains when `rank` is `"species"` or `"genus"`, and can be easily taken care of by removing occurrences where `species_reso` or `genus_reso` does not equal a 'safe' value, a functionality controlled by the argument `cleanUncertain`, which is `TRUE` by default. The 'safe' values for `species_reso` or `genus_reso`  are listed in the `cleanResoValues` argument, and by default includes common notifiers like *'n. sp.'* which indicate something other than taxonomic uncertainty.
	
# The 'cleaning' of names so they don't have unwanted taxonomic cruff on them, like that one sock in your laundry that attracts all the distasteful lint. This is particularly problematic as I use names, and not taxon ID numbers, to match occurrences as having the same tazon. However, my use-case is somewhat different from Clapham's, as `cleanTaxon` uses `taxon_name`, while I use `matched_name` for formal species and the appropriate taxonomic variable for higher-taxa (i.e. `genus` for genera) and, for informal taxa, combining `genus_name` and `species_name` (for species), with `taxon_name` only being called for referencing informal supraspecific taxa. My working assumption is that taxonomic name rubbish has been kept out of `matched_name` and the various 'formal' supraspecific taxon variables, as well as `genus_name` and `species_name` for the taxonomic names given in the reference for that occurrence. That could be a bad assumption, but so far I haven't found any. I know this rubbish exists in `taxon_name`, but (a) its hard to identify every use-case of such name-trash, and plus the case where an occurrence isn't given a formal higher-taxon is very rare (in the data I've looked at). 
	
#' @details
#' Data input for \code{taxonSortPBDBocc} are expected to be from version 1.2 API
#' with the 'pbdb' vocabulary. However, datasets are passed to internal function \code{translatePBDBocc},
#' which attempts to correct any necessary field names and field contents used by
#' \code{taxonSortPBDBocc}.

#' @param data A table of occurrence data collected from the Paleobiology Database. 

#' @param rank The selected taxon rank; must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'.

#' @param onlyFormal If TRUE (the default) only taxa formally accepted by the Paleobiology Database are returned. If FALSE, then the identified name fields are searched for any additional 'informal' taxa with the proper taxon. If their taxon name happens to match any formal taxa, their occurrences are merged onto the formal taxa. This argument generally has any appreciable effect when rank=species.

#' @param cleanUncertain If TRUE (the default) any occurrences with an entry in the respective 'resolution' field that is *not* found in the argument cleanResoValue will be removed from the dataset. These are assumed to be values indicating taxonomic uncertainty, i.e. 'cf.' or '?'.

#' @param cleanResoValues The set of values that can be found in a 'resolution' field that do not cause a taxon to be removed, as they do not seem to indicate taxonomic uncertainty.

#' @return
#' Returns a list where each element is different unique taxon obtained by the sorting function, and named with that taxon name. Each element is composed of a table containing all the same occurrence data fields as the input (potentially with some fields renamed and some field contents change, due to vocabulary translation).
	
#' @seealso
#' Example dataset of a PBDB download can be found at \code{\link{graptPBDB}}

#' @author 
#' David W. Bapst, but partly inspired by Matthew Clapham's \code{cleanTaxon} (found at https://github.com/mclapham/PBDB-R-scripts/blob/master/taxonClean.R on github) and R package paleobioDB's \code{pbdb_temp_range} function (found at https://github.com/ropensci/paleobioDB/blob/master/R/pbdb_temporal_functions.R#L64-178 on github.

#' @examples
#' #load example graptolite PBDB occ dataset
#' data(graptPBDB)
#' 
#' #get formal genera
#' occGenus<-taxonSortPBDBocc(graptOccPBDB, rank="genus")
#' length(occGenus)
#' 
#' #get formal species
#' occSpeciesFormal<-taxonSortPBDBocc(graptOccPBDB, rank="species")
#' length(occSpeciesFormal)
#' 
#' #yes, there are fewer 'formal' graptolite species in the PBDB then genera
#' 
#' #get formal and informal species
#' occSpeciesInformal<-taxonSortPBDBocc(graptOccPBDB, rank="species",
#' 	 onlyFormal=FALSE)
#' length(occSpeciesInformal)
#' 
#' #way more graptolite species are 'informal' in the PBDB
#' 
#' #get formal and informal species 
#' 	#including from occurrences with uncertain taxonomy
#' 	#basically everything and the kitchen sink
#' occSpeciesEverything<-taxonSortPBDBocc(graptOccPBDB, rank="species",
#' 		onlyFormal=FALSE, cleanUncertain=FALSE)
#' length(occSpeciesEverything)
#' 
#' \dontrun{
#' 
#' #try a PBDB API download with lots of synonymization
#' 	#this should have only 1 species
#' acoData<-read.csv(paste0("http://paleobiodb.org/data1.1/occs/list.txt?",
#' 	"base_name=Acosarina%20minuta&show=ident,phylo&limit=all"))
#' x<-taxonSortPBDBocc(acoData, rank="species", onlyFormal=FALSE)
#' names(x)
#'
#' #make sure works with API v1.2
#' 		#won't work until v1.2 goes live at the regular server
#' dicelloDataCom<-read.csv(paste0("http://paleobiodb.org",
#' 	"/data1.2/occs/list.txt?base_name=Dicellograptus",
#' 	"&show=ident,phylo&limit=all"))
#' dicelloOcc2<-taxonSortPBDBocc(dicelloDataCom, rank="species", onlyFormal=FALSE)
#' names(dicelloOcc2)
#' 
#' #make sure works with compact vocab v1.1
#' dicelloData<-read.csv(paste0("http://paleobiodb.org",
#' 	"/data1.1/occs/list.txt?base_name=Dicellograptus",
#' 	"&show=ident,phylo&limit=all&vocab=com"))
#' dicelloOccCom1<-taxonSortPBDBocc(dicelloData, rank="species", onlyFormal=FALSE)
#' names(dicelloOccCom1)
#' head(dicelloOccCom1[[1]])[,1:7]
#' 
#' 
#'
#' }
#' 
#' 


#' @name taxonSortPBDBocc
#' @rdname taxonSortPBDBocc
#' @export
taxonSortPBDBocc<-function(data,rank, onlyFormal=TRUE, cleanUncertain=TRUE, 
								cleanResoValues=c(NA, '"', "", "n. sp.", "n. gen.")){
	#this function inspired by Matt Clapham's taxonClean and paleobioDB's pbdb_temp_range
		#onlyFormal=FALSE;rank="species"
		#onlyFormal=FALSE;rank="genus"
	#pull occurrences out of a data table and sort by unique taxa into a list structure
	#translated vocabs!
	data<-translatePBDBocc(data)
	#Second, some warning checks to see if occurrences downloaded correctly: 
	# need phylo and ident data
		#the following is taken with minor modification from code in paleobioDB package
	if (!any("primary_name"==colnames(data)) | !any("genus"==colnames(data))){	
		stop("need to add 'show=c('phylo', 'ident')' to pbdb_occurrences query\n *or*\n  'show=ident,phylo' to PBDB API query")
		}
	#additional checks for rank
	if(length(rank)!=1){stop("length of rank must be 1")}
	if(!any(sapply(c("species","genus","family","order","class","phylum"),function(x) x==rank))){
		stop("rank must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")}
	#for inconsistencies between rank and onlyFormal - NOT TRUE, IGNORE THIS CHECK
	#if(!onlyFormal & (rank!="species" | rank!="genus")){
	#  stop("Informal taxon does not exist for above genus level, please use 'onlyFormal=TRUE'")}
	#need to replace any empty string values with NAs (due perhaps to use of read.csv with the API)
	data[data==""]<-NA
	#now,  pull taxa using the relevant variable from 'phylo' for formal ID
		#this matches what paleobioDB now does
	#sort occurrences by unique taxa in each level and then append valid names
	if(rank=="species"){
		if(cleanUncertain){ #remove uncertain species taxonomy
		data<-data[data[,"species_reso"] %in% cleanResoValues,,drop=FALSE]
		}
	#first formal taxa
    #get species names for taxa formally recognized as species level
    whichFormal<-(data[,"accepted_rank"]==rank)
    taxonVar<-as.character(data[,"accepted_name"])
    taxaNames<-as.character(unique(taxonVar[whichFormal]))  #get unique taxa
    #drop empty entries (occurrences listed at a higher taxonomic level, formally)
    taxaNames<-taxaNames[!is.na(taxaNames)]
    #take rows these taxa are in and drop them into elements of a list
		sortedOcc<-lapply(taxaNames,function(x) data[which(x==taxonVar),,drop=FALSE])
    names(sortedOcc)<-taxaNames
    if(!onlyFormal){
		#now use taxon_rank to identify useful informal occurrence
		taxonVar2<-data[,c("primary_name","species_name")]
		taxonVar2<-apply(taxonVar2,1,paste,collapse=" ")
		#which of these occs are useful as informal occs
		stillUseful<-which(!whichFormal & (data[,"identified_rank"]==rank))
		taxaNames2<-as.character(unique(taxonVar2[stillUseful]))
		#drop any weird empties
		taxaNames2<-taxaNames2[taxaNames2!=" " & !is.na(taxaNames2)]  
		sortedOcc2<-lapply(taxaNames2,function(x) data[which(x==taxonVar2),,drop=FALSE])
		names(sortedOcc2)<-taxaNames2
		# merge sortedocc2 with sortedOcc
		share1<-sapply(taxaNames,function(x) any(x==taxaNames2))
		share2<-sapply(taxaNames2,function(x) any(x==taxaNames))
		if(sum(!share1)>0){sortedOccU<-sortedOcc[!share1]}else{sortedOccU<-list()}
		if(sum(!share2)>0){sortedOcc2U<-sortedOcc2[!share2]}else{sortedOcc2U<-list()}
		#and for shared names
		if(sum(share1)>0){
			shared<-lapply(taxaNames[share1],function(x) 
				cbind(sortedOcc[[taxaNames==x]],sortedOcc2[[taxaNames2==x]]))
			names(shared)<-names(taxaNames[share1])
		}else{shared<-list()}
		sortedOcc<-c(sortedOccU,sortedOcc2U,shared)
		}
	}else{   #if not at the species rank
		if(cleanUncertain){  #removing uncertain taxonomy if appended to primary_name
			data<-data[data[,"primary_reso"] %in% cleanResoValues,,drop=FALSE]
			}
		taxonVar<-data[,rank] #then our taxonomic variable of interest is just so!
		taxaNames<-as.character(unique(taxonVar))  #get unique taxa
		taxaNames<-taxaNames[!is.na(taxaNames)]
		#take rows these taxa are in and drop them into elements of a list
		sortedOcc<-lapply(taxaNames,function(x) data[which(x==taxonVar),,drop=FALSE])
		names(sortedOcc)<-taxaNames
		if(!onlyFormal){
			#now use taxon_rank to identify useful informal occurrence
			taxonVar2<-data[,"primary_name"]
			#which of these occs are useful as informal occs
			stillUseful<-which(is.na(taxonVar) & (data[,"identified_rank"]==rank))
			taxaNames2<-as.character(unique(taxonVar2[stillUseful]))
			taxaNames2<-taxaNames2[!is.na(taxaNames2)]
			sortedOcc2<-lapply(taxaNames2,function(x) data[which(x==taxonVar2),,drop=FALSE])
			names(sortedOcc2)<-taxaNames2
			# merge sortedocc2 with sortedOcc
			share1<-sapply(taxaNames,function(x) any(x==taxaNames2))
			share2<-sapply(taxaNames2,function(x) any(x==taxaNames))
			if(sum(!share1)>0){sortedOccU<-sortedOcc[!share1]}else{sortedOccU<-list()}
			if(sum(!share2)>0){sortedOcc2U<-sortedOcc2[!share2]}else{sortedOcc2U<-list()}
			#and for shared names
			if(sum(share1)>0){
				shared<-lapply(taxaNames[share1],function(x)
					cbind(sortedOcc[[taxaNames==x]],sortedOcc2[[taxaNames2==x]]))
				names(shared)<-names(taxaNames[share1])
			}else{shared<-list()}
			sortedOcc<-c(sortedOccU,sortedOcc2U,shared)
			}
		}
	# sort occurrences by taxon name
	sortedOcc<-sortedOcc[order(names(sortedOcc))]
	return(sortedOcc)
	}

translatePBDBocc<-function(data){
	#translate PBDB occ data	
	if(any(colnames(data)=="taxon_name")){
		#from PBDB API version 1.1 to 1.2
		colnames(data)[colnames(data)=="taxon_name"]<-"identified_name"	
		colnames(data)[colnames(data)=="taxon_no"]<-"identified_no"
		colnames(data)[colnames(data)=="taxon_rank"]<-"identified_rank"	
		colnames(data)[colnames(data)=="matched_name"]<-"accepted_name"
		colnames(data)[colnames(data)=="matched_no"]<-"accepted_no"
		colnames(data)[colnames(data)=="matched_rank"]<-"accepted_rank"
		colnames(data)[colnames(data)=="genus_name"]<-"primary_name"
		colnames(data)[colnames(data)=="genus_reso"]<-"primary_reso"
		}
	if(any("tna"==colnames(data))){
		#need to translate data from 1 vocab to another
			#do on by-function basis - this is for taxonSort
		#for this just translate colname of relevant taxon variables
			#also will need to translate 'taxon_rank' and 'matched_rank' contents
		if(any("mna"==colnames(data))){
			#then compact vocab v1.1
			if(!all(c("tna","rnk","tid","mna","mra","mid","idt","ids",
				"gnl","fml","odl","cll","phl","rst","rss") %in% colnames(data))){
					stop("Not all fields founds for compact vocab v1.1")}	
			colnames(data)[colnames(data)=="tna"]<-"identified_name"	
			colnames(data)[colnames(data)=="tid"]<-"identified_no"
			colnames(data)[colnames(data)=="rnk"]<-"identified_rank"	
			colnames(data)[colnames(data)=="mna"]<-"accepted_name"
			colnames(data)[colnames(data)=="mid"]<-"accepted_no"
			colnames(data)[colnames(data)=="mra"]<-"accepted_rank"
			colnames(data)[colnames(data)=="idt"]<-"primary_name"
			colnames(data)[colnames(data)=="ids"]<-"species_name"
			colnames(data)[colnames(data)=="rss"]<-"primary_reso"
			colnames(data)[colnames(data)=="rst"]<-"species_reso"
			colnames(data)[colnames(data)=="gnl"]<-"genus"
			colnames(data)[colnames(data)=="fml"]<-"family"
			colnames(data)[colnames(data)=="odl"]<-"order"
			colnames(data)[colnames(data)=="cll"]<-"class"	
			colnames(data)[colnames(data)=="phl"]<-"phylum"			
			}
		if(any("idn"==colnames(data))){		
			#then compact vocab v1.2
			if(!all(c(	"idn","idr","iid","tna","rnk","tid","idg","ids",
				"gnl","fml","odl","cll","phl","rsg","rss") %in% colnames(data))){
					stop("Not all fields founds for compact vocab v1.1")}	
			colnames(data)[colnames(data)=="idn"]<-"identified_name"	
			colnames(data)[colnames(data)=="iid"]<-"identified_no"
			colnames(data)[colnames(data)=="idr"]<-"identified_rank"	
			colnames(data)[colnames(data)=="tna"]<-"accepted_name"
			colnames(data)[colnames(data)=="tid"]<-"accepted_no"
			colnames(data)[colnames(data)=="rnk"]<-"accepted_rank"
			colnames(data)[colnames(data)=="idg"]<-"primary_name"
			colnames(data)[colnames(data)=="ids"]<-"species_name"
			colnames(data)[colnames(data)=="rss"]<-"primary_reso"
			colnames(data)[colnames(data)=="rsg"]<-"species_reso"
			colnames(data)[colnames(data)=="gnl"]<-"genus"
			colnames(data)[colnames(data)=="fml"]<-"family"
			colnames(data)[colnames(data)=="odl"]<-"order"
			colnames(data)[colnames(data)=="cll"]<-"class"	
			colnames(data)[colnames(data)=="phl"]<-"phylum"			
			#stop("need to add 'vocab='pbdbd'' to pbdb_occurrences query\n  *or*\n 'vocab=pbdb' to your PBDB API query")
			}
		# taxon rank translation vectors for compact vocab
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")
		taxRankCOM<-2:26
		#change contents of "identified_rank" and "accepted_rank"
		data$identified_rank<-sapply(data$identified_rank,function(x) taxRankPBDB[x==taxRankCOM])
		data$accepted_rank<-sapply(data$accepted_rank,function(x) taxRankPBDB[x==taxRankCOM])
		}
	return(data)
	}
	



	
	