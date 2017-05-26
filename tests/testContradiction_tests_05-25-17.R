# test for hidden function testContradiction
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