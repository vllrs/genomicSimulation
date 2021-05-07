test_that("After deleting a group, it is no longer in memory", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #setup worked as expected
  expect_identical(length(see.group.data(g, "Indexes")), 6L)
  
  #testing function
  expect_output(delete.group(g), "6 genotypes were deleted")
  
  expect_identical(length(see.group.data(g, "Indexes")), 0L)
          })

test_that("After deleting a group, other data is shuffled correctly in memory", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  capture_output(g2 <- load.more.genotypes("helper_genotypes.txt"), print=F)
  capture_output(g3 <- load.more.genotypes("helper_genotypes.txt"), print=F)
  
  expect_output(delete.group(g2), "6 genotypes were deleted")
  expect_identical(see.group.data(g, "Indexes"), c(0L, 1L, 2L, 3L, 4L, 5L))
  expect_identical(see.group.data(g3, "Indexes"), c(6L, 7L, 8L, 9L, 10L, 11L))
  
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  capture_output(g2 <- load.more.genotypes("helper_genotypes_long.txt"), print=F)
  capture_output(g3 <- load.more.genotypes("helper_genotypes_long.txt"), print=F)
  
  expect_output(delete.group(g2), "1003 genotypes were deleted")
  expect_identical(see.group.data(g, "Indexes"), c(0L, 1L, 2L, 3L, 4L, 5L))
  expect_identical(see.group.data(g3, "Indexes"), c(6L:1008L))
  
  expect_output(delete.group(g), "6 genotypes were deleted")
  expect_identical(see.group.data(g3, "Indexes"), c(0L:1002L))
  
})