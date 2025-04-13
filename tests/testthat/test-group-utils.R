test_that("make.group moves the correct genotypes to a new group", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  g2 <- make.random.crosses(g, n.crosses=5, offspring=1, give.names=TRUE, name.prefix="cr")
  
  #setup worked as expected
  expect_identical(see.group.data(g, "XIndexes"), c(0L, 1L, 2L, 3L, 4L, 5L))
  expect_identical(see.group.data(g, "Names"), c("G01", "G02", "G03", "G04", "G05", "G06"))
  expect_identical(see.group.data(g2, "XIndexes"), c(6L, 7L, 8L, 9L, 10L))
  expect_identical(see.group.data(g2, "Names"), c("cr7", "cr8", "cr9", "cr10", "cr11"))
  
  #testing make.group worked
  mg <- make.group(c(1L,4L, 9L, 10L))
  expect_identical(see.group.data(mg, "XIndexes"), c(1L, 4L, 9L, 10L))
  expect_identical(see.group.data(mg, "Names"), c("G02", "G05", "cr10", "cr11"))
  expect_identical(see.group.data(g, "XIndexes"), c(0L, 2L, 3L, 5L))
  expect_identical(see.group.data(g, "Names"), c("G01", "G03", "G04", "G06"))
  expect_identical(see.group.data(g2, "XIndexes"), c(6L, 7L, 8L))
  expect_identical(see.group.data(g2, "Names"), c("cr7", "cr8", "cr9"))
  
})

test_that("labels can be used to split groups", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  indivs <- break.group.into.individuals(g)
  
  expect_identical(make.label(1L), 1L)
  
  expect_identical(change.label.by.amount(1L, -2, indivs[1]), 0L)
  expect_identical(change.label.to.this(1L, 5L, indivs[3]), 0L)
  expect_identical(change.label.to.values(1L, c(2L,3L), skip=4), 0L)
  
  gnew <- combine.groups(indivs)
  expect_identical(see.group.data(gnew,"L"),c(-1L,1L,5L,1L,2L,3L))

  g1 <- break.group.by.label.value(1L, -1)
  g2 <- break.group.by.label.range(1L, 3, 10)
  
  expect_identical(see.group.data(g1, "Names"), c("G01"))
  expect_identical(see.group.data(g2, "Names"), c("G03", "G06"))
})

test_that("combine.groups successfully merges 2+ groups", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  g2 <- make.random.crosses(g, n.crosses=5, offspring=1)
  g3 <- make.random.crosses(g, n.crosses=5, offspring=1)
  g4 <- make.random.crosses(g, n.crosses=7, offspring=1)
  
  #setup worked as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3,g4),"GroupSize"=c(6L,5L,5L,7L)))
  
  #Test combine groups worked
  g5 <- combine.groups(c(g2,g3,g4))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g5),"GroupSize"=c(6L,17L)))
  
  g6 <- combine.groups(c(g,g5))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g6),"GroupSize"=c(23L)))
  
})

test_that("break.group.into.individuals runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup worked as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g),"GroupSize"=c(6L)))
  
  #function works
  g2 <- break.group.into.individuals(g)
  expect_identical(see.existing.groups(), data.frame("Group"=g2,"GroupSize"=rep(1L,times=6)))
  
})

test_that("break.group.into.families runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup
  g2 <- make.targeted.crosses(0L,1L, give.names=T, name.prefix="two") #two7
  g3 <- make.targeted.crosses(0L,4L, offspring=4, give.names=T, name.prefix="three") #three8, three9, three10, three11
  g4 <- make.targeted.crosses(3L,6L, offspring=3, give.names=T, name.prefix="four") #four12, four13, four14
  gcom <- combine.groups(c(g2,g3,g4))
  
  #setup works as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,gcom),"GroupSize"=c(6L,8L)))
  
  #function works
  fs <- break.group.into.families(gcom)
  df <- see.existing.groups()
  expect_identical(df$Group, c(g,fs))
  expect_identical(df$GroupSize, c(6L,1L,4L,3L))
  
  expect_identical(see.group.data(fs[1], "Name"), c("two7"))
  expect_identical(see.group.data(fs[2], "Name"), c("three8", "three9", "three10", "three11"))
  expect_identical(see.group.data(fs[3], "Name"), c("four12", "four13", "four14"))
})


test_that("break.group.into.halfsib.families runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup
  g2 <- make.targeted.crosses(0L,1L, give.names=T, name.prefix="two") #two7
  g3 <- make.targeted.crosses(0L,4L, offspring=4, give.names=T, name.prefix="three") #three8, three9, three10, three11
  g4 <- make.targeted.crosses(3L,6L, offspring=3, give.names=T, name.prefix="four") #four12, four13, four14
  gcom <- combine.groups(c(g2,g3,g4))
  
  #setup works as expected
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,gcom),"GroupSize"=c(6L,8L)))
  
  #function works
  fs <- break.group.into.halfsib.families(gcom,1L)
  df <- see.existing.groups()
  expect_identical(df$Group, c(g,fs))
  expect_identical(df$GroupSize, c(6L,5L,3L))
  
  expect_identical(see.group.data(fs[1], "Name"), c("two7", "three8", "three9", "three10", "three11"))
  expect_identical(see.group.data(fs[2], "Name"), c("four12", "four13", "four14"))
})

test_that("break.group.randomly runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup
  g2 <- make.random.crosses(g, 2000)
  
  fs <- break.group.randomly(g2)
  
  df <- see.existing.groups()
  
  expect_identical(df$Group, c(g,fs))
  expect_gt(df$GroupSize[[2]], 800)
  expect_lt(df$GroupSize[[2]], 1200)
  expect_identical(df$GroupSize[[2]] + df$GroupSize[[3]], 2000L)
  #does not check that the groups are shuffled properly
  
})


test_that("break.group.with.probabilities runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup
  g2 <- make.random.crosses(g, 2000)
  
  fs <- break.group.with.probabilities(g2, c(0.2,0.5))
  
  df <- see.existing.groups()
  
  expect_identical(df$Group, c(g,fs))
  expect_gt(df$GroupSize[[2]], 200)
  expect_lt(df$GroupSize[[2]], 600)
  expect_gt(df$GroupSize[[3]], 800)
  expect_lt(df$GroupSize[[3]], 1200)
  expect_identical(df$GroupSize[[2]] + df$GroupSize[[3]] + df$GroupSize[[4]], 2000L)
  #does not check that the groups are shuffled properly
  
  # does-not-crash test:
  expect_output(break.group.with.probabilities(df$Group[[3]], c(0.1,0.6,0.7)),
                 "NOTE! Provided probabilities add up to 1 or more: some buckets will not be filled")
  
})


test_that("break.group.evenly runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup
  g2 <- make.random.crosses(g, 20)
  
  fs <- break.group.evenly(g)
  f2s <- break.group.evenly(g2, 3)
  
  truth = data.frame("Group"=c(fs,f2s),"GroupSize"=c(3L,3L,7L,7L,6L))
  expect_identical(see.existing.groups(), truth[order(truth$Group),])
  #does not check that the groups are shuffled properly
  
})

test_that("break.group.into.buckets runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  #setup
  g2 <- make.random.crosses(g, 20)
  
  fs <- break.group.into.buckets(g2, c(3L,13L))
  
  truth <- data.frame("Group"=c(g,fs),"GroupSize"=c(6L,3L,13L,4L))
  expect_identical(see.existing.groups(), truth[order(truth$Group),])
  #does not check that the groups are shuffled properly
  
  #does-not-crash test
  expect_output(break.group.into.buckets(fs[[2]], c(2L,2L,2L,20L)), 
                 "NOTE! Provided capacities are larger than actual group: some buckets will not be filled")
})


test_that("break.group.by.GEBV runs successfully", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  expect_identical(see.group.data(break.group.by.GEBV(init$groupNum, number=5), "X"), c(0L,1L,2L,3L,4L))
  
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  expect_identical(see.group.data(break.group.by.GEBV(init$groupNum, low.score.best=T, number=2), "X"), c(3L,5L))
  
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  expect_identical(see.group.data(break.group.by.GEBV(init$groupNum, percentage=50), "X"), c(0L,1L,2L))
  
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  expect_identical(see.group.data(break.group.by.GEBV(init$groupNum, low.score.best=T, percentage=20), "X"), c(5L))
  
  expect_error(break.group.by.GEBV(g, low.score.best=T),"Exactly one of parameters `percentage` and `number` must be set.")
  expect_error(break.group.by.GEBV(g, percentage=20, number=2),"Exactly one of parameters `percentage` and `number` must be set.")
  
  clear.simdata()
})
