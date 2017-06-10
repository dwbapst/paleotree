



# function for scaling posterior trees to absolute time (get root.time)
#


kiko<-'
absTimescaleFossilMrB<-function(trees,fixedAges=attr(trees,"fixedTable"){
	



	
	#For all trees to be comparable, we will use the $root.time convention from paleotree (Bapst, 2012)

#add $root.time to BEAST and MrBayes trees
	#do all have Zanz. jr at 0?

whatIsYoungest<-function(tree){
  nodeDates<-suppressMessages(dateNodes(tree=tree))
  youngest<-which(nodeDates==min(nodeDates))[1]
  youngest<-c(tree$tip.label[[youngest]],min(nodeDates))
  youngest
}


#check to see if the youngest taxon in every tree is "Zanabazar_junior" with an age of 0
youngTruth<-c("Zanabazar_junior",0)
youngTest<-c(sapply(treesB2BDSS,function(tree) identical(whatIsYoungest(tree=tree),youngTruth)),
             sapply(treesB2SABD,function(tree) identical(whatIsYoungest(tree=tree),youngTruth)),
             sapply(treesMrB,function(tree) identical(whatIsYoungest(tree=tree),youngTruth)))
if(!all(youngTest)){
  stop(paste0("Zanabazar junior is not youngest taxon on trees ",
              paste0(which(!youngTest),collapse = " ")))
  }

#set root age
setRootAgeZanz<-function(tree){
  tree$root.time<-max(node.depth.edgelength(tree))+zanaDate
  return(tree)
  }

treesB2BDSStime<-lapply(treesB2BDSS,setRootAgeZanz)
treesB2SABDtime<-lapply(treesB2SABD,setRootAgeZanz)
treesMrBtime<-lapply(treesMrB,setRootAgeZanz)

	



	
	}

'