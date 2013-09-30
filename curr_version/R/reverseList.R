reverseList<-function(list,simplify=FALSE){
	#reverses primary and secondary list structure
		#i.e. if a list is 10 elements each 50 long, get 50 elements 10 long
	if(length(unique(sapply(list,length)))!=1){
		stop("Error: Not all lists equally long")}
	list1<-list()
	for(i in 1:length(list[[1]])){
		list1[[i]]<-sapply(list,function(x) x[[i]],simplify=simplify)
		}
	names(list1)<-names(list[[1]])
	return(list1)
	}