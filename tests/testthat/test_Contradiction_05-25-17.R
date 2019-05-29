test_that("contradiction distance tests work as expected", {

# test for hidden function testContradiction

#library(paleotree)

test<-c(
  paleotree:::testContradiction(c("A","B"),c("A","B")), #should be FALSE (no contradiction)
  paleotree:::testContradiction(c("A","B"),c("C","D")), #should be FALSE
  paleotree:::testContradiction(c("A","B"),c("A","C")), #should be TRUE
  paleotree:::testContradiction(c("A","B"),c("B","C")), #should be TRUE
  paleotree:::testContradiction(c("A","B"),c("A","D")), #should be TRUE
  paleotree:::testContradiction(c("A","B"),c("B","D"))  #should be TRUE
)
expected<-c(FALSE,FALSE,TRUE,TRUE,TRUE,TRUE)
if(!identical(test,expected)){stop("testContradiction failed")}

# test with three labels
test<-c(
  paleotree:::testContradiction(c("A","B","C"),c("A","B")), #should be FALSE
  paleotree:::testContradiction(c("A","B","C"),c("C","D")), #should be TRUE
  paleotree:::testContradiction(c("A","B","C"),c("A","C")), #should be FALSE
  paleotree:::testContradiction(c("A","B","C"),c("B","C")), #should be FALSE
  paleotree:::testContradiction(c("A","B","C"),c("A","D")), #should be TRUE
  paleotree:::testContradiction(c("A","B","C"),c("B","D"))  #should be TRUE
)
expected<-c(FALSE,TRUE,FALSE,FALSE,TRUE,TRUE)
if(!identical(test,expected)){stop("testContradiction failed")}

# test treeContradiction

# test with three labels
test<-c(
  paleotree:::testContradiction(c("A","B","C"),c("A","B")), #should be FALSE
  paleotree:::testContradiction(c("A","B","C"),c("C","D")), #should be TRUE
  paleotree:::testContradiction(c("A","B","C"),c("A","C")), #should be FALSE
  paleotree:::testContradiction(c("A","B","C"),c("B","C")), #should be FALSE
  paleotree:::testContradiction(c("A","B","C"),c("A","D")), #should be TRUE
  paleotree:::testContradiction(c("A","B","C"),c("B","D"))  #should be TRUE
)
expected<-c(FALSE,TRUE,FALSE,FALSE,TRUE,TRUE)
if(!identical(test,expected)){stop("testContradiction failed")}

set.seed(1)
treeA<-rtree(30,br=NULL)
treeB<-rtree(30,br=NULL)
# let's simulate the worst resolved tree possible: a star tree
treeC<-stree(30)

# Tree C (the star tree) has zero CD between it and trees A and B
test<-c(
	identical(treeContradiction(treeA,treeC),0),  # should be zero distance
	identical(treeContradiction(treeB,treeC),0),  # should be zero distance
	# two identical trees also have zero CD between them (as you'd hope) 
	identical(treeContradiction(treeA,treeA),0)  # should be zero distanceH
)
if(!all(test)){stop("test of treeContradiction failed!")}

})