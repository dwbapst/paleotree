# write simple sumt

# The standard MrBayes output format for summary trees is difficult to read into and depict in R with support values at nodes. 
# This script creates a NEXUS file with a MrBayes block that would load the output of existing MrBayes analyses, 
# then run the 'sumt' command with the correct formatting for use in R, thus quickly generating all desired summary trees in the correct format.


writeSimpleSumtNexus<-function(directory,origNexus=NULL,outputNexus=NULL){
	#
	# The standard MrBayes output format for summary trees is difficult to read into and depict in R with support values at nodes. 
	# This script creates a NEXUS file with a MrBayes block that would load the output of existing MrBayes analyses, 
	# then run the 'sumt' command with the correct formatting for use in R, thus quickly generating all desired summary trees in the correct format.
	#
	#change to dir for MB files
	#setwd(directory)
	#get files
	files<-list.files(directory)
	paths<-paste0(directory,"//",files)
	#get tree files
	treeFiles<-files[grep(x=files,pattern="run..t")]  
	#get the unique identifiers
	inputs<-sort(unique(sub(".run..t","",treeFiles)))
	#get the nexus files
	inputs<-sort(unique(sub(".tree.","",inputs)))
	#
	# this should catch partitioned unlinked analyses 
	# (with more than one tree) as separate instances
	# for calculating consensus trees as the rest of the nexus file is read...
	#
	#now start constructing sumt lines
	bayesBlock<-"#NEXUS"
	for(i in (1:length(inputs))){
		if(is.null(origNexus)){
			origNex<-scan(inputs[i],what="character",sep="\n",blank.lines.skip=FALSE)
		}else{
			origNex<-scan(origNexus,what="character",sep="\n",blank.lines.skip=FALSE)
			}
		cropNex<-origNex[2:grep(origNex,pattern="mcmcp")[1]]
		#remove log line
		cropNex<-cropNex[-grep(cropNex,pattern="log start")]
		sumtLine<-paste0("sumt filename=",inputs[i]," outputname=",
			paste0(inputs[i],".simple")," conformat=simple;")
		bayesBlock<-c(bayesBlock,cropNex,sumtLine,"","end;","")
		}
	#remove log line
	cropNex<-cropNex[-grep(cropNex,pattern="log start")]
	#now generate lines for all inputs
	for(i in 1:length(inputs)){
		sumtLine<-paste0("sumt filename=",inputs[i]," outputname=",
			paste0(inputs[i],".simple")," conformat=simple;")
		bayesBlock<-c(bayesBlock,cropNex,sumtLine,"","end;","")
		}
	if(!is.null(outputNexus)){
		write(bayesBlock,file=outputNexus)
	}else{
		return(bayesBlock)
		}
	}


