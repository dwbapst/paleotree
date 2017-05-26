#' Construct a Fully Formatted NEXUS Script for Performing Tip-Dating Analyses With MrBayes
#' 
#' This function is meant to expedite the creation of NEXUS files formatted
#' for performing tip-dating analyses in the popular phylogenetics software \emph{MrBayes},
#' particularly clock-less tip-dating analyses executed with 'empty' morphological matrices
#' (i.e. where all taxa ae coded for a single missing character), although a pre-existing
#' morphological matrix can also be input by the user (see argument \code{origNexusFile}).
#' 
#' The minimum required input for this function is a data set of tip ages (in various formats),
#' which are used to construct age calibrations commands on the tip taxa 
#' (via paleotree function \code{\link{createMrBayesTipCalibrations}}), and a set
#' of taxa designated as the outgroup, which is then converted into a command constraining
#' th monophyly on the ingroup taxa, which is presumed to be all taxa \emph{not} listed in the outgroup.
#' Many of the options available with \code{\link{createMrBayesTipCalibrations}} are available with this function,
#' allowing users to choose between fixed calibrations or uniform priors that approximate stratigraphic uncertainty.
#' In addition, a user may supply a tree which is then converted into a series of hard topological
#' constraints (via function \code{\link{createMrBayesConstraints}}, or a path to a text file
#' presumed to be a NEXUS file containing character data formatted for use with \emph{MrBayes}.
#' 
#' The resulting full NEXUS script is output as a set of character strings either
#' printed to the R console, or output to file which is then overwritten.

#' @details
#' The taxa listed in \code{tipTimes} must match the taxa in 
#' \code{treeConstraints}, if such is supplied. The taxa in \code{outgroupTaxa}
#' must be contained within this same set of taxa, and should not contradict
#' relationships on \code{treeConstraint}. The taxa in any
#' character matrix given in \code{origNexusFile} is \code{not} checked against
#' these two sources: it is up to the user to ensure the same
#' taxa are found in all three.
#' 
#' Note that because the same set of taxa must be contained in both inputs, 
#' relationships are constrained as 'hard' constraints, rather than 'partial' constraints,
#' which allows some taxa to float across a partially fixed topology. 
#' See the documentation for \code{\link{createMrBayesConstraints}},
#' for more details.

#' @inheritParams createMrBayesTipCalibrations

#' @param outgroupTaxa A vector of type 'character', containing taxon names designating the outgroup.
#'  All taxa not listed in the outgroup will be constrained to be a monophyletic ingroup, for sake of rooting
#' the resulting dated tree.  This outgroup selection should not disagree with any of the
#' relationships on \code{treeConstraint}, if such is supplied.

#' @param origNexusFile Filename (possibly with path) as a character
#' string leading to a NEXUS text file, presumably containing a matrix
#' of character date formateed for \emph{MrBayes}. If supplied
#' (it does not need to be supplied), the listed file is read as a text file, and
#' concatenated with the \emph{MrBayes} script produced by this function, so as to reproduce. 
#' Note that the taxa in this NEXUS file are \emph{NOT} checked against the user
#' input \code{tipTimes} and \code{treeConstraints}, so it is up to the user to
#' ensure the taxa are the same across the three data sources.

#' @param newFile Filename (possibly with path) as a character string
#' leading to a file which will be overwritten with the output tip age calibrations.
#' If \code{NULL}, tip calibration commands are output to the console.

#' @param createEmptyMorphMat If \code{origNexusFile} is not specified (implying there is no
#' pre-existing morphological character matrix for this dataset), then an 'empty' NEXUS-formatted matrix will be
#' appended to the set of \emph{MrBayes} commands if this command is \code{TRUE} (the default). This
#' 'empty' matrix will have each taxon in \code{tipTimes} coded for a single missing character
#' (i.e., '?'). This allows tip-dating analyses with hard topological constraints, and ages
#' determined entirely by the fossilized birth-death prior, with no impact from a
#' presupposed morphological clock (thus a 'clock-less analysis').

#' @param runName The name of the run, used for naming the log files. 
#' If not set, the name will be taken from the name given for outputting
#' the NEXUS script (\code{newFile}). If \code{newFile} is not given, and
#' \code{runName} is not set by the user, the default run name will be  "new_run_paleotree".

#' @param doNotRun If \code{TRUE}, the commands that cause a script to automatically begin running in 
#' \emph{MrBayes} will be left out. Useful for troubleshooting initial runs of scripts for non-fatal errors and
#' warnings (such as ignored constraints). Default for this argument is \code{FALSE}.

#' @param treeConstraints An object of class \code{phylo}, 
#' from which (if \code{treeConstraints} is supplied) the set topological constraints are derived, as
#' as described for argument \code{tree} for function \code{createMrBayesConstraints}.

#' @param morphModel This argument can be used to switch between two end-member models of 
#' morphological evolution in MrBayes, here named 'strong' and 'relaxed', for the 'strong assumptions'
#' and 'relaxed assumptions' models described by Bapst, Schreiber and Carlson (Systematic Biology).
#' The default is a model which makes very 'strong' assumptions about the process of morphological evolution,
#' while the 'relaxed' alternative allows for considerably more heterogeneity in the rate
#' of morphological evolution across characters, and in the forward and reverse transition
#' rates between states. Note that in both cases, the character data is assumed to be filtered
#' to only parsimony-informative characters, without autapomorphies.

#' @param cleanNames If \code{TRUE} (the default), then special characters
#' (currently, thsi only contains the forward-slashes: '/') are removed from
#' taxon names before construction of the NEXUS file.

#' @param printExecute If \code{TRUE} (the default) and if output is directed to a \code{newFile}
#' (i.e. a \code{newFile) is specified), a line for pasting into MrBayes for executing the newly created file
#' will be messaged to the terminal.

#' @note 
#' This function allows a user to take an undated phylogenetic tree in R, and a set of age estimates
#' for the taxa on that tree, and produce a posterior sample of dated trees using the MCMCMC in \emph{MrBayes},
#' while treating an 'empty' morphological matrix as an uninformative set of missing characters.
#' This 'clock-less tip-dating' approach is essentially an alternative to the \emph{cal3} method in paleotree, 
#' sharing the same fundamental theoretical model (a version of the fossilized birth-death model), but
#' with a better algorithm that considers the whole tree simultaneously, rather than evaluating each node
#' individually, from the root up to the tips (as \emph{cal3} does it, and which may cause artifacts). 
#' That said, \emph{cal3} still has a few advantages: tip-dating as of April 2017 still only treats OTUs as
#' point observations, contained in a single time-point, while cal3 can consider taxa as having durations with
#' first and last occurrences. This means it may be more straightforward to assess the extent of budding cladogenesis
#' patterns of ancestor-descendant relationships in \emph{cal3}, than in tip-dating.

#' @return
#' If argument \code{newFile} is \code{NULL}, then the text of the 
#' generated NEXUS script is ouput to the console as a series of character strings.

#' @seealso
#' \code{\link{createMrBayesConstraints}}, \code{\link{createMrBayesTipCalibrations}}, , \code{\link{cal3}}

#' @author
#' David W. Bapst. This code was produced as part of a project 
#' funded by National Science Foundation grant EAR-1147537 to S. J. Carlson.
#' 
#' The basic \emph{MrBayes} commands utilized in the output script are a collection
#' of best practices taken from studying NEXUS files supplied by April Wright,
#' William Gearty, Graham Slater, Davey Wright, and
#' guided by the recommendations of Matzke and Wright, 2016 in Biol. Lett.

#' @references
#' The basic fundamentals of tip-dating, and tip-dating with the fossilized
#' birth-death model are introduced in these two papers: 
#' 
#' Ronquist, F., S. Klopfstein, L. Vilhelmsen, S. Schulmeister, D. L. Murray,
#' and A. P. Rasnitsyn. 2012. A Total-Evidence Approach to Dating with Fossils,
#' Applied to the Early Radiation of the Hymenoptera. \emph{Systematic Biology} 61(6):973-999.
#' 
#' Zhang, C., T. Stadler, S. Klopfstein, T. A. Heath, and F. Ronquist. 2016. 
#' Total-Evidence Dating under the Fossilized Birth-Death Process.
#' \emph{Systematic Biology} 65(2):228-249. 
#' 
#' For recommended best practices in tip-dating analyses, please see:
#' 
#' Matzke, N. J., and A. Wright. 2016. Inferring node dates from tip dates
#' in fossil Canidae: the importance of tree priors. \emph{Biology Letters} 12(8).
#'
#' The rationale behind the two alternative morphological models are described in more detail here:
#' 
#' Bapst, D. W., H. A. Schreiber, and S. J. Carlson. In press. Combined analysis of extant Rhynchonellida
#' (Brachiopoda) using morphological and molecular data. \emph{Systematic Biology} doi: 10.1093/sysbio/syx049
#' 



#' @examples
#'
#' # let's do some examples
#' 
#' # load retiolitid dataset
#' data(retiolitinae)
#' 
#' # let's try making a NEXUS file!
#' 
#' # Use a uniform prior, with a 10 million year offset for
#' 	# the expected tree age from the earliest first appearance
#' # set average tree age to be 10 Ma earlier than first FAD
#' 
#' outgroupRetio<-"Rotaretiolites" # sister to all other included taxa
#' 
#' # the following will create a NEXUS file with an 'empty' morph matrix
#' 	# with the only topological constraint on ingroup monophyly
#' 	# Probably shouldn't do this: leaves too much to the FBD prior
#'  
#' # with doNotRun set to TRUE for troubleshooting
#' 
#' createMrBayesTipDatingNexus(tipTimes=retioRanges,
#' 		outgroupTaxa=outgroupRetio,treeConstraints=NULL,
#' 		ageCalibrationType="uniformRange",whichAppearance="first",
#' 		treeAgeOffset=10,	newFile=NULL,	
#' 		origNexusFile=NULL,createEmptyMorphMat=TRUE,
#' 		runName="retio_dating",doNotRun=TRUE)
#' 
#' # let's try it with a tree for topological constraints
#' 
#' # let's set doNotRun to FALSE
#' 
#' createMrBayesTipDatingNexus(tipTimes=retioRanges,
#' 		outgroupTaxa=outgroupRetio,treeConstraints=retioTree,
#' 		ageCalibrationType="uniformRange",whichAppearance="first",
#' 		treeAgeOffset=10,	newFile=NULL,	
#' 		origNexusFile=NULL,createEmptyMorphMat=TRUE,
#' 		runName="retio_dating",doNotRun=FALSE)
#' 
#' # the above is essentially cal3 with a better algorithm,
#' 		# and no need for a priori rate estimates
#' # just need a tree and age estimates for the tips!
#' 
#' #############################################################################
#' # some more variations for testing purposes
#' 
#' # no morph matrix supplied or generated
#' 	# you'll need to manually append to an existing NEXUS file
#' createMrBayesTipDatingNexus(tipTimes=retioRanges,
#' 		outgroupTaxa=outgroupRetio,treeConstraints=retioTree,
#' 		ageCalibrationType="uniformRange",whichAppearance="first",
#' 		treeAgeOffset=10,	newFile=NULL,	
#' 		origNexusFile=NULL,createEmptyMorphMat=FALSE,
#' 		runName="retio_dating",doNotRun=TRUE)
#' 
#' \dontrun{
#' 
#' # let's actually try writing an example with topological constraints
#' 	# to file and see what happens
#' 
#' # here's my super secret MrBayes directory
#' file<-"D:\\dave\\workspace\\mrbayes\\exampleRetio.nex"
#' 
#' createMrBayesTipDatingNexus(tipTimes=retioRanges,
#' 		outgroupTaxa=outgroupRetio,treeConstraints=retioTree,
#' 		ageCalibrationType="uniformRange",whichAppearance="first",
#' 		treeAgeOffset=10,	newFile=file,	
#' 		origNexusFile=NULL,createEmptyMorphMat=TRUE,
#' 		runName="retio_dating",doNotRun=FALSE)
#' 
#' }
#' 

## TO DO !!
#
# -Need to check that ingroup constraint isn't on treeConstraints, and if so, delete it
	# actually maybe just make it so no ingroup constraint is defined if treeConstraint is defined - presumably its already part of provided tree!!
# need to read the NEXUS file so can duplicate taxa for FAD/LAD/etc




#' @aliases tipdating
#' @name createMrBayesTipDatingNexus
#' @rdname createMrBayesTipDatingNexus
#' @export
createMrBayesTipDatingNexus<-function(tipTimes,outgroupTaxa,treeConstraints=NULL,morphModel="strong",
							ageCalibrationType,whichAppearance="first",treeAgeOffset,minTreeAge=NULL,
							collapseUniform=TRUE,origNexusFile=NULL,createEmptyMorphMat=TRUE,newFile=NULL,
							runName="new_run_paleotree",doNotRun=FALSE,cleanNames=TRUE,printExecute=TRUE){
	################################################################################################
	#         # a wooper of a function ... here's some ASCII from 'Psyduck'
	#
	#                       =@@@@@@                         
	#                     =@@.......@@@                      
	#                   =@..   ........@@    =@              
	#             =    @...    ..........@   =@  @           
	#          @  =@  @....   ............@  =@@ @           
	#          @  =@ =.....................@  =@ @==@        
	#         @@  =@ @...@.........@.......@  =@==@@         
	#         @====@@.......................@===@@=          
	#          @@@==@.=..............=....@===@@  @@         
	#         @@ @@@=.@..............@....@=@@=@   @         
	#         @   @@=..@=..........@=......@@ @=@            
	#             @ @....=@@@@@@@@=.........@ @=@            
	#            @@ @.......................@  @=            
	#            @   @.....................@                 
	#                =.....................@                 
	#                 @...................@                  
	#                  @.................@                   
	#                   @@.............@@                    
	#                     @@=.......=@@                      
	#                       @@.......@          @@@          
	#                       @.=====...@       @@   @         
	#                      @...===.....@     @      @        
	#                      @...........@    @.      @        
	#                     @..==...==....@ @@. .     @        
	#                     @...=====.....@@.... . .  @        
	#                     @.............@...... . .@         
	#                      @.==...==....@.........@          
	#                       @.=====....@=@@.....@@           
	#                     =@@=@......@@    @@@@@             
	#                    =...==@@@@@@==@@                    
	#                    @....@     =....@                   
	#                     @@@@      @....@                   
	#                                @@@@
	#
	#################################################################################################
	# CHECK TAXON NAMES
	# taxa in tipTimes is king
		# all taxa in input treeConstraints must be in tipTimes (and vice versa)
		# origNexusFile is *not* checked
	if(is.data.frame(tipTimes)){
		tipTimes<-as.matrix(tipTimes)}
	if(is.list(tipTimes)){
		taxaTipTimes<-rownames(tipTimes[[2]])
		}else{
		taxaTipTimes<-rownames(tipTimes)
		}
	if(!is.null(treeConstraints)){
		if(!is(treeConstraints,"phylo")){
			stop("treeConstraints must be a phylogeny object of type 'phylo'")
			}
		taxaTree<-treeConstraints$tip.label
		# check both
		missingTip<-taxaTipTimes[sapply(taxaTipTimes,function(x) all(x!=taxaTree))]
		missingTree<-taxaTree[sapply(taxaTree,function(x) all(x!=taxaTipTimes))]
		# stop if length>0
		if(length(c(missingTip,missingTree))>0){
			stop(paste0("Following taxa in tipTimes not found on treeConstraints: ",paste0(missingTip,collapse=" "),
				"\nFollowing taxa on treeConstraint not found in tipTimes: ",paste0(missingTree,collapse=" ")))
			}
		if(length(taxaTree)!=length(taxaTipTimes)){
			stop("Somehow have taxa missing from either tipTimes or treeConstraint")
			}
		if(!identical(sort(taxaTipTimes),sort(taxaTree))){
			stop("Nope, taxa in tipTimes and on treeConstraint STILL not identical!!")
			}
		}
	# create 'clean' version of taxon names - remove all '/' 
	cleanTaxonNames<-taxaTipTimes
	if(length(unique(cleanTaxonNames))!=length(cleanTaxonNames)){
		stop("Some taxon names were identical duplicates")
		}
	if(cleanNames){
		cleanTaxonNames<-gsub("/","",cleanTaxonNames)
		# check that unique
		if(length(unique(cleanTaxonNames))!=length(cleanTaxonNames)){
			stop("Some taxon names were identical duplicates when special character (/) was removed")
			}
		# gsub("[^A-Za-z0-9]", "", a) # removes every character except letters and numbers - only use if dire (i.e. backslashes)
		#
		if(is.list(tipTimes)){
			rownames(tipTimes[[2]])<-gsub("/","",rownames(tipTimes[[2]]))
		}else{
			rownames(tipTimes)<-gsub("/","",rownames(tipTimes))
			}	
		#
		if(!is.null(treeConstraints)){
			treeConstraints$tip.label<-gsub("/","",treeConstraints$tip.label)
			}
		}
	#
	taxonnames<-cleanTaxonNames
	#
	# check other arguments
	if(length(morphModel)!=1){
		stop("morphModel must be of length 1")
		}
	if(all(morphModel!=c("relaxed","strong"))){
		stop("morphModel must be one 'relaxed' or 'strong'")
		}
	if(morphModel=="relaxed" & is.null(origNexusFile)){
		warning("Why are you relaxing the morphological model without supplying an original matrix with origNexusFile? I hope you know what you are doing.")
		}
	# eh, origNexusFile might be a connection...
	#if(length(origNexusFile)!=1){
	#	stop("origNexusFile must be of length 1") 
	#	}
	if(length(createEmptyMorphMat)!=1 | !is.logical(createEmptyMorphMat)){
		stop("createEmpthyMorphMat must be  of length 1, and either TRUE or FALSE")
		}
	#####################################################################
	# get original Nexus file (or create an empty morph matrix
	if(!is.null(origNexusFile)){
		morphNexus<-readLines(con=origNexusFile,warn=FALSE)
	}else{
		if(createEmptyMorphMat){
			morphNexus<-makeEmptyMorphNexusMrB(taxonNames=taxonnames)
		}else{
			morphNexus<-" "
			}
		}
	########################################################
	# get runName if not supplied
	if(is.null(runName)){
		if(is.null(newFile)){
			stop("runName must be supplied if name of new file is not designated")
			}
		# get run name - everything after the last / or \\
		runName<-rev(strsplit(newFile,split="\\\\")[[1]])[1]
		runName<-rev(strsplit(runName,split="/")[[1]])[1]
		}
	# use run name as log file name
	logfileline<-paste0('log start filename="',runName,'.out" replace;')
	########################################################
	#	
	# get the ingroup constraint
	ingroupConstraint<-makeIngroupConstraintMrB(outgroupTaxa=outgroupTaxa,allTaxa=taxonnames)	
	#
	# get topological constraints, if given in input
	if(!is.null(treeConstraints)){
		topologicalConstraints<-createMrBayesConstraints(tree=treeConstraints,partial=FALSE,
			file=NULL,includeIngroupConstraint=TRUE)
	}else{
		topologicalConstraints<-" "
		}
	#
	# get age calibration block
	ageCalibrations<-createMrBayesTipCalibrations(tipTimes=tipTimes,
			ageCalibrationType=ageCalibrationType,whichAppearance=whichAppearance,
			treeAgeOffset=treeAgeOffset,minTreeAge=minTreeAge,
			collapseUniform=collapseUniform,file=NULL)
	#
	# make the final MrBayes Block
	MrBayesBlock<-makeMrBayesBlock(logBlock=logfileline,ingroupBlock=ingroupConstraint,
		ageBlock=ageCalibrations,constraintBlock=topologicalConstraints,
		morphModel=morphModel,doNotRun=doNotRun)
	###############################################
	# combine morph matrix block with MrBayes command block
	finalText<-c(morphNexus,MrBayesBlock)
	if(!is.null(newFile)){
		write(x=finalText,file=newFile)
		if(printExecute){
			# have function print command for pasting into MrBayes to execute: 
				# e.g. ' Execute "C://fossil data/myNexus.nex" '
			message("Now go to MrBayes and paste in this line:")
			message(paste0('Execute "',normalizePath(newFile),'";'))
			}
	}else{
		return(finalText)
		}
	}




#################################################################################

# Supplemental functions

###################################################################################

makeIngroupConstraintMrB<-function(outgroupTaxa,allTaxa){
	# use the list of outgroup taxa to define the in-group taxa
	# write a MrBayes line that defines all ingroup taxa
		# as being monophyletic
	################################################
	# checks
	if(length(outgroupTaxa)<1){
		stop("outgroupTaxa must be of length 1 or greater")
		}
	if(!is.character(outgroupTaxa) | !is.vector(outgroupTaxa)){
		stop("outgroupTaxa must be a vector of type character")
		}	
	if(!is.character(allTaxa) | !is.vector(allTaxa)){
		stop("taxon labels cannot be turned into a vector of strings")
		}	
	# are there any outgroup not in all?
	notPresent<-sapply(outgroupTaxa,function(x) all(x!=allTaxa))
	if(any(notPresent)){
		stop(paste("Following outgroup taxa not present in input tipTimes",
			paste(outgroupTaxa[notPresent],collapse=" ")))
		}
	#############################################
	# get the ingroup taxa
	inGroup<-allTaxa[sapply(allTaxa,function(x) all(x!=outgroupTaxa))]
	# check that there are SOME ingroup
	if(length(inGroup)<2){
		stop("Less than two taxa are in the ingroup??")
		}
	# make ingroup block
	ingroupBlock<-paste0("constraint ingroup = ",
		paste(inGroup,collapse=" "),";")
	ingroupBlock<-c(ingroupBlock,
		"prset topologypr = constraints(ingroup);")
	return(ingroupBlock)
	}

makeEmptyMorphNexusMrB<-function(taxonNames){
	#given a vector of taxonNames, make NEXUS character block
		# formatted for MrBayes
		# each character is a single ?
	# check
	if(!is.character(taxonNames) | !is.vector(taxonNames)){
		stop("taxon names not formatted correctly for making empty NEXUS matrix")
		}
	################
	# get standard header and 	
	nTaxa<-length(taxonNames)
	if(nTaxa<1){stop("vector of taxon names is zero-length")}
	# build header
	headMat<-c("#NEXUS","BEGIN DATA;",
		paste0("DIMENSIONS  NTAX=",nTaxa," NCHAR=1;"),
		"FORMAT DATATYPE=STANDARD GAP=- MISSING=?;",
		"MATRIX")
	endMat<-c("",";","END;")
	###############################
	# build body
	bodyMat<-sapply(taxonNames,function(x)
		paste0(x,"   ?")
		)
	#########################
	finalMatBlock<-c(headMat,bodyMat,endMat)
	names(finalMatBlock)<-NULL
	return(finalMatBlock)
	}


##########################################################################################################################################



makeMrBayesBlock<-function(logBlock,ingroupBlock,ageBlock,
							constraintBlock,morphModel="strong",doNotRun=FALSE){
#########################################################################################						
	###
	### # what follows will look very messy due to its untabbed nature
	### # I'm so sorry - dwb
	###
##################################################
block1<-"
begin mrbayes;

[This block is a combination of best practices taken from NEXUS files from April Wright,
     William Gearty, Graham Slater, Davey Wright, and guided by the 
	 recommendations of Matzke and Wright, 2016, Biol. Lett.]

[no autoclose, I like it to ask, but you might want it]
	[set autoclose=yes;]
"
##########################################################################
block2<-"
[DATA]

[edits to taxon or character selection go here]
			
[CONSTRAIN MONOPHYLY OF INGROUP]

	[constrain all members of the ingroup to be monophyletic]
"
#######################################################################
block3<-"	
[TOPOLOGICALCONSTRAINTS FOR ADDITIONAL NODES]

	[[EXAMPLE]
    constraint node1 = t1 t2 t3;
    constraint node2 = t4 t5;
	prset topologypr = constraints(ingroup,node1,node2); [need to include ingroup]]
"
##############################################################################
morphModelBlock_Strong<-"

[CHARACTER MODELS]
	
[morphology model settings]
	[make the number of beta categories for stationary frequencies 4 (the default)]
	[default: use pars-informative coding]
	
	[set coding and rates - default below maximizes information content]
		lset  nbetacat=5 rates=equal Coding=informative; [equal rate variation]
		[lset  nbetacat=5 rates=gamma Coding=informative;   [gamma distributed rate variation]]
			[gamma distributed rates may cause divide by zero problems with non-fixed symdiri]
	[symdirhyperpr prior, fixed vs. variable]	
		prset symdirihyperpr=fixed(infinity);		
		[prset symdirihyperpr=uniform(1,10);      [this range seems to avoid divide by zero error]]

"

morphModelBlock_Relaxed<-"

[CHARACTER MODELS]
	
[morphology model settings]
	[make the number of beta categories for stationary frequencies 4 (the default)]
	[default: use pars-informative coding]
	
	[set coding and rates - default below maximizes information content]
		[lset  nbetacat=5 rates=equal Coding=informative; [equal rate variation]]
		lset  nbetacat=5 rates=gamma Coding=informative;   [gamma distributed rate variation]
			[gamma distributed rates may cause divide by zero problems with non-fixed symdiri]
	[symdirhyperpr prior, fixed vs. variable]	
		[prset symdirihyperpr=fixed(infinity);]		
		prset symdirihyperpr=uniform(1,10);      [this range seems to avoid divide by zero error]

"


#############################################################################
block4<-"
[TIP AND TREE AGE CALIBRATIONS] 
	
		[[EXAMPLE WITH BOUNDS: min age, max age]	
		calibrate t1 = uniform(105, 113);
		calibrate t2 = uniform(97, 115);
		[EXAMPLE WITH FIXED AGE]
		calibrate t3 = fixed (100);]
	
	[prior on tree age]
	[with offset exponential in MrB, first par is min age, second is expected mean age]
		[mean date must be greater than min date - duh]
		
		[EXAMPLE]
		[prset treeagepr = offsetexp(115,120); ]
			
"
##########################################################################
block5<-"	
[FOSSILIZED BIRTH DEATH MODEL & ASSOCIATED PARAMETERS]

    prset brlenspr=clock:fossilization;
	prset fossilizationpr = beta(1,1); [flat, sampling is psi/(mu+psi), 0-1 ]

    prset speciationpr = uniform(0,10);
	prset extinctionpr = beta(1,1); [flat, extinction is relative to speciation between 0-1]

	prset samplestrat = random; [default: bdss prior ]
	[prset samplestrat = fossiltip; [this would mean no sampled ancestors]]
	
	prset sampleprob = 1; [rho, sampling in the present, default fixed to 1]
	[Unclear what to set this to if all taxa are extinct...]
	
[CLOCK PRIORS]
	
	[clock rate: truncated normal, mean=0.08, but very flat; see Supp Table 2 of Matzke & Wright 2016]
    prset clockratepr = normal(0.0025,0.1);     
	prset clockvarpr=igr;
    prset igrvarpr = uniform(0.0001, 200); [vague prior that is actually vague; Matzke & Wright 2016]

	prset nodeagepr = calibrated; [set node age priors]
	
[SETTINGS FOR PARTITIONED ANALYSES]

	[Necessary if a combined analysis, I guess]

[RUN SETTINGS]

[FYI burnin settings for one command sets burnin settings for all]
	mcmcp ngen=100000000 relburnin=yes burninfrac=0.5 printfreq=1000 samplefreq=1000 nchains=4 nruns=2 savebrlens=yes;    
	set seed=4;
"
###################################################################
runBlock<-"	

[RUN]
	mcmc;
	sump;
	
	[Given this is a dating analysis, should really use MCCT,not summary trees]
		[I guess best summary tree though is allcompat?]
	[sumt contype=allcompat; ]    [allcompat has issues with sampled ancestors]
	
	log stop;
"
######################################################
	if(doNotRun){
		runBlock<-" "
		}
	if(morphModel=="strong"){
		morphModelBlock<-morphModelBlock_Strong
		}
	if(morphModel=="relaxed"){
		morphModelBlock<-morphModelBlock_Relaxed
		}
	####################
	finalBlock<-c(block1,logBlock,block2,ingroupBlock,block3,
		constraintBlock,morphModelBlock,block4,
		ageBlock,block5,runBlock," ","end;")
	return(finalBlock)
	}



