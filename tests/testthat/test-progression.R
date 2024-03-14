#At the moment, all these test is that the functions run without crashing and make correctly-sized groups.
#More in-depth testing and testing of how GenOptions options work is not yet implemented.

test_that("make.random.crosses works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.random.crosses(g, n.crosses=3, offspring=2, give.names=TRUE, name.prefix="cr")
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  
  # and check cross.randomly.between works (simple check, C version has more)
  g3 <- make.random.crosses.between(g, g2, n.crosses=3, offspring=2)
  g4 <- make.random.crosses.between(g2, g3, cap1=5, n.crosses=30)

  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3,g4),"GroupSize"=c(6L,6L,6L,30L)))
  
  clear.simdata()
})

test_that("make.targeted.crosses works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.targeted.crosses(first.parents=c(0L,1L), second.parents=c(3L,0L))
  g3 <- make.targeted.crosses(first.parents=c("G01", "G02", "G03"), second.parents=c("G06","G05","G04"))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3),"GroupSize"=c(6L,2L,3L)))
  
  clear.simdata()
})

test_that("make.crosses.from.file works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.crosses.from.file("helper_crosses.txt")
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,2L)))
  
  clear.simdata()
})

test_that("make.all.unidirectional.crosses works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.all.unidirectional.crosses(g)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,15L)))
  
  clear.simdata()
})

test_that("make.double.crosses.from.file works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  f <- make.all.unidirectional.crosses(g)
  g2 <- make.double.crosses.from.file("helper_dcrosses.txt")
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,f,g2),"GroupSize"=c(6L,15L,2L)))
  
  clear.simdata()
})

test_that("self.n.times works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- self.n.times(g, 2)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  
  clear.simdata()
})

test_that("make.doubled.haploids works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.doubled.haploids(g)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  
  # Check homozygous
  a <- see.group.data(g2,"Genotypes")
  expect_identical(paste(strsplit(a,"")[[1]][c(T,F)]),paste(strsplit(a,"")[[1]][c(F,T)]))
  
  clear.simdata()
})

test_that("make.clones works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #Make sure it doesn't crash and creates right-size group
  g2 <- make.clones(g)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,6L)))
  expect_identical(see.group.data(g2,"Genotypes"),see.group.data(g,"genotypes"))
  
  clear.simdata()
})
