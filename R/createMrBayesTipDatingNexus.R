#' Construct a Fully Formatted NEXUS Script for Performing Tip-Dating Analyses With MrBayes
#'
#' 






#' @details
#' The taxa listed in \code{tipTimes} must match the taxa in \code{treeConstraints}.

# taxa in tipTimes is king
	# all taxa in input treeConstraints must be in tipTimes (and vice versa)
	# origNexusFile is *not* checked


#' @inheritParams createMrBayesConstraints createMrBayesTipCalibrations

#' @param outgroupTaxa

#' @param origNexusFile Filename (possibly with path) as a character
#' string leading to a NEXUS text file, presumably containing a matrix
#' of character date formateed for MrBayes. If supplied
#' (it does not need to be supplied), the listed file is read as a text file, and
#' concatenated with the MrBayes script produced by this function, so as to reproduce. 

#' @param newFile Filename (possibly with path) as a character string
#' leading to a file which will be overwritten with the output tip age calibrations.
#' If \code{NULL}, tip calibration commands are output to the console.

#' @param createEmptyMorphMat

#' @param runName The name of the run, used for naming the log files. 
#' If not set, the name will be taken from the name given for outputting
#' the NEXUS script (\code{newFile}). If \code{newFile} is not given, and
#' \code{runName} is not set by the user, the default run name will be  "new_run_paleotree".

#' @param doNotRun=FALSE





#' @param treeConstraints An object of class \code{phylo}, 
#' from which (if \code{treeConstraints} is supplied) the set topological constraints are derived, as
#' as described for argument \code{tree} for function \code{createMrBayesConstraints}.



#' @param treeConstraints 

#' @param createEmptyMorphMat If origNexusFile is not supplied, should an empty NEXUS block formatted
#' for use with MrBayes, be created using the list of taxa in tipTimes? Default is TRUE.


#' @return
#' If argument \code{newFile} is \code{NULL}, then the text of the 
#' generated NEXUS script is ouput to the console as a series of character strings.

#' @seealso
#' \code{\link{createMrBayesConstraints}}, \code{\link{createMrBayesConstraints}}

#' @author David W. Bapst

#' @references
#' The basic fundamentals of tip-dating, and tip-dating with the fossilized
#' birth-death model are introduced in these two papers:
#'
#' Ronquist, F., S. Klopfstein, L. Vilhelmsen, S. Schulmeister, D. L. Murray,
#' and A. P. Rasnitsyn. 2012. A Total-Evidence Approach to Dating with Fossils,
#' Applied to the Early Radiation of the Hymenoptera. \emph{Systematic Biology} 61(6):973-999.
#'
#' Zhang, C., T. Stadler, S. Klopfstein, T. A. Heath, and F. Ronquist. 2016. 
#' Total-Evidence Dating under the Fossilized Birth–Death Process.
#' \emph{Systematic Biology} 65(2):228-249. 

#' @examples




#' @name createMrBayesTipDatingNexus
#' @rdname createMrBayesTipDatingNexus
#' @export
createMrBayesTipDatingNexus<-function(tipTimes,outgroupTaxa,treeConstraints=NULL,
							ageCalibrationType,whichAppearance="first",treeAgeOffset,minTreeAge=NULL,
							origNexusFile=NULL,newFile=NULL,createEmptyMorphMat=TRUE,
							runName="new_run_paleotree",doNotRun=FALSE){
	################################################################################################
	#         # a whopper of a function
	#
	#################################################################################################
	# CHECK TAXON NAMES
	# taxa in tipTimes is king
		# all taxa in input treeConstraints must be in tipTimes (and vice versa)
		# origNexusFile is *not* checked
	if(is.list(tipTimes)){
		taxaTipTimes<-rownames(tipTimes[[1]])
		}else{
		taxaTipTimes<-rownames(tipTimes)
		}
	if(!is.null(treeConstraints)){
		taxaTree<-treeConstraints$tip.labels
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
	taxonnames<-taxaTipTimes
	#
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
	ageCalibrations<-createMrBayesTipCalibrations<-function(tipTimes=tipTimes,
			ageCalibrationType=ageCalibrationType,whichAppearance=whichAppearance,
			treeAgeOffset=treeAgeOffset,minTreeAge=minTreeAge,file=NULL)
	#
	# make the final MrBayes Block
	MrBayesBlock<-makeMrBayesBlock(logBlock=logfileline,ingroupBlock=ingroupConstraint,
		ageBlock=ageCalibrations,constraintBlock=topologicalConstraints,
		doNotRun=doNotRun)
	###############################################
	# combine morph matrix block with MrBayes command block
	finalText<-c(morphNexus,MrBayesBlock)
	if(!is.null(newFile)){
		write(finalText,file)
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

makeMrBayesBlock<-function(logBlock,ingroupBlock,ageBlock,constraintBlock,doNotRun=FALSE){
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
	sumt contype=allcompat; 
	log stop;
"
######################################################
	if(doNotRun){
		runBlock<-" "
		}
	finalBlock<-c(block1,logBlock,block2,ingroupBlock,block3,
		constraintBlock,block4,ageBlock,block5,runBlock," ","end;")
	return(finalBlock)
	}



