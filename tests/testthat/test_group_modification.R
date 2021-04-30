test_that("make.group moves the correct genotypes to a new group", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g2 <- cross.randomly(g, n.crosses=5, offspring=1, give.names=TRUE, name.prefix="cr")
  
  #setup worked as expected
  expect_identical(see.group.data(g, "Indexes"), c(0L, 1L, 2L, 3L, 4L, 5L))
  expect_identical(see.group.data(g, "Names"), c("G01", "G02", "G03", "G04", "G05", "G06"))
  expect_identical(see.group.data(g2, "Indexes"), c(6L, 7L, 8L, 9L, 10L))
  expect_identical(see.group.data(g2, "Names"), c("cr7", "cr8", "cr9", "cr10", "cr11"))
  
  #testing make.group worked
  mg <- make.group(c(1L,4L, 9L, 10L))
  expect_identical(see.group.data(mg, "Indexes"), c(1L, 4L, 9L, 10L))
  expect_identical(see.group.data(mg, "Names"), c("G02", "G05", "cr10", "cr11"))
  expect_identical(see.group.data(g, "Indexes"), c(0L, 2L, 3L, 5L))
  expect_identical(see.group.data(g, "Names"), c("G01", "G03", "G04", "G06"))
  expect_identical(see.group.data(g2, "Indexes"), c(6L, 7L, 8L))
  expect_identical(see.group.data(g2, "Names"), c("cr7", "cr8", "cr9"))
  
})

test_that("combine.groups successfully merges 2+ groups", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g2 <- cross.randomly(g, n.crosses=5, offspring=1)
  g3 <- cross.randomly(g, n.crosses=5, offspring=1)
  g4 <- cross.randomly(g, n.crosses=7, offspring=1)
  
  #setup worked as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3,g4),"GroupSize"=c(6L,5L,5L,7L)))
  
  #Test combine groups worked
  g5 <- combine.groups(c(g2,g3,g4))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g5),"GroupSize"=c(6L,17L)))
  
  g6 <- combine.groups(c(g,g5))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g6),"GroupSize"=c(23L)))
  
})

test_that("break.group.into.individuals runs successfully", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #setup worked as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g),"GroupSize"=c(6L)))
  
  #function works
  g2 <- break.group.into.individuals(g)
  expect_identical(see.existing.groups(), data.frame("Group"=g2,"GroupSize"=rep(1L,times=6)))
  
})

test_that("break.group.into.families runs successfully", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  #setup
  g2 <- cross(0L,1L, give.names=T, name.prefix="two") #two7
  g3 <- cross(0L,4L, offspring=4, give.names=T, name.prefix="three") #three8, three9, three10, three11
  g4 <- cross(3L,6L, offspring=3, give.names=T, name.prefix="four") #four12, four13, four14
  gcom <- combine.groups(c(g2,g3,g4))
  
  #setup works as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,gcom),"GroupSize"=c(6L,8L)))
  
  #function works
  fs <- break.group.into.families(gcom)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,fs),"GroupSize"=c(6L,1L,4L,3L)))
  
  expect_identical(see.group.data(fs[1], "Name"), c("two7"))
  expect_identical(see.group.data(fs[2], "Name"), c("three8", "three9", "three10", "three11"))
  expect_identical(see.group.data(fs[3], "Name"), c("four12", "four13", "four14"))
})