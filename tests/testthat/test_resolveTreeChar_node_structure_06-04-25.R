# test to make sure resolveTreeChar isn't numbering nodes wrong
    # based on examples from Zena Lapp 06-01-25

test_that("resolveTreeChar function is fine", {

set.seed(4)    
    
# example tree that doesn't work for orderedChar = TRUE
example_tree <- structure(list(
    edge = structure(
        c(6L, 6L, 6L, 6L, 6L, 1L, 2L, 3L, 4L, 5L), 
        dim = c(5L, 2L)),
    
    edge.length = c(0.0019837977, 0.012936646, 0.0265849463,
                    0.0033293995, 0.0053513114),
    Nnode = 1L, node.label = "", 
    tip.label = c("74", "458", "807", "74", "262")),
    class = "phylo", order = "cladewise")

#checkValidPhylo(example_tree)
#testEdgeMat(example_tree)

example_tree <- cleanNewPhylo(example_tree)

example_trait <- factor(example_tree$tip.label, 
                        levels = sort(unique(
                            as.numeric(example_tree$tip.label))))
names(example_trait) <- example_tree$tip.label

# works with orderedChar = FALSE
example_tree_resolved_unordered <- resolveTreeChar(
    example_tree, example_trait, orderedChar = FALSE)


# doesn't work with orderedChar = TRUE regardless of stateBias. 
example_tree_resolved_ordered <- resolveTreeChar(
    example_tree, example_trait, orderedChar = TRUE)


expect_true(testEdgeMat(example_tree_resolved_unordered))

expect_true(testEdgeMat(example_tree_resolved_ordered))


# changing the "74" tip labels to "31" fixes the problem:
example_tree2 <- example_tree
example_tree2$tip.label[example_tree2$tip.label == "74"] <- "31"

example_trait2 <- factor(example_tree2$tip.label, 
                         levels = sort(unique(as.numeric(example_tree2$tip.label))))
names(example_trait2) <- example_tree2$tip.label

example_tree_resolved_unordered2 <- resolveTreeChar(example_tree2, 
                                                    example_trait2, orderedChar = FALSE)
example_tree_resolved_ordered2 <- resolveTreeChar(example_tree2,
                                                  example_trait2, orderedChar = TRUE)

expect_true(testEdgeMat(example_tree_resolved_unordered))

expect_true(testEdgeMat(example_tree_resolved_ordered))

#####
# partially resolved example

example_tree <- structure(list(
    edge = structure(c(12L, 12L, 13L, 14L, 14L, 15L,
                       15L, 13L, 13L, 16L, 17L, 18L, 
                       18L, 18L, 18L, 17L, 16L, 1L, 13L,
                       14L, 2L, 15L, 3L, 4L, 5L, 16L, 
                       17L, 18L, 6L, 7L, 8L, 9L, 10L,
                       11L), 
                     dim = c(17L, 2L)), 
    edge.length = c(0.0104310095, 0.0068884473, 0.0027000338, 0.0103750076, 0.0091539735,
                    0.0064466731, 0.0152757103, 0.0319091095, 0.0706365098, 0.0083015911, 0.0016927818,
                    0.0171124659, 0.0207471466, 0.0085753685, 0.0208127876, 0.003382258, 0.0207352662 ), 
    Nnode = 7L, 
    node.label = c("0", "75.3", "75.6", "77.3", "75.1", "81.7", "77.5"), 
    tip.label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")), 
    class = "phylo", 
    order = "cladewise")

#str(example_tree)

#checkValidPhylo(example_tree)
#testEdgeMat(example_tree)
#example_tree <- cleanNewPhylo(example_tree)

#plot(example_tree)

example_trait <- structure(
    c('1' = 1L, '2' = 2L, '3' = 3L, '4' = 4L, '5' = 6L, '6' = 5L, 
      '7' = 7L, '8' = 8L, '9' = 10L, '10' = 7L, '11' = 9L), 
    levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
    class = "factor")

# works with orderedChar = FALSE
example_tree_resolved_unordered <- resolveTreeChar(
    example_tree, example_trait, orderedChar = FALSE)

# doesn't work with orderedChar = TRUE
example_tree_resolved_ordered <- resolveTreeChar(
    example_tree, example_trait, orderedChar = TRUE)

expect_true(testEdgeMat(example_tree_resolved_unordered))

expect_true(testEdgeMat(example_tree_resolved_ordered))

})
