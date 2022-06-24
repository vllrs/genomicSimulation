#At the moment, all these test is that the functions run without crashing and make correctly-sized groups.
#More in-depth testing and testing of how GenOptions options work is not yet implemented.

test_that("cross.randomly works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- cross.randomly(g, n.crosses=3, offspring=2, give.names=TRUE, name.prefix="cr")
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  
  #This may fail just because it picked the same random cross multiple times ...
  #break.group.into.families(g2)
  #expect_identical(see.existing.groups()$GroupSize, "GroupSize"=c(6L,2L,2L,2L))
  
  # and check cross.randomly.between works (simple check, C version has more)
  g3 <- cross.randomly.between(g, g2, n.crosses=3, offspring=2)
  g4 <- cross.randomly.between(g2, g3, cap1=5, n.crosses=30)

  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3,g4),"GroupSize"=c(6L,6L,6L,30L)))
  
  clear.simdata()
})



test_that("cross.combinations works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- cross.combinations(first.parents=c(0L,1L), second.parents=c(3L,0L))
  g3 <- cross.combinations(first.parents=c("G01", "G02", "G03"), second.parents=c("G06","G05","G04"))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3),"GroupSize"=c(6L,2L,3L)))
  
  clear.simdata()
})

test_that("cross.combinations.file works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- cross.combinations.file("helper_crosses.txt")
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,2L)))
  
  clear.simdata()
})

test_that("cross.all.pairs works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- cross.all.pairs(g)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,15L)))
  
  clear.simdata()
})

test_that("cross.dc.combinations.file works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  f <- cross.all.pairs(g)
  g2 <- cross.dc.combinations.file("helper_dcrosses.txt")
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,f,g2),"GroupSize"=c(6L,15L,2L)))
  
  clear.simdata()
})

test_that("self.n.times works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- self.n.times(g, 2)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  
  clear.simdata()
})

test_that("make.doubled.haploids works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.doubled.haploids(g)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  
  clear.simdata()
})
