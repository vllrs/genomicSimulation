#At the moment, all these test is that the functions run without crashing and make correctly-sized groups.
#We're reliant on the testing of the underlying C library for everything else about correct functionality

test_that("GenOptions works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  #current_id <- 6
  
  # "track.pedigree" setting
  g2 <- make.targeted.crosses(c("G01","G02","G03"), c(2,3,4), track.pedigree=TRUE)
  expect_identical(see.group.data(g2,"P1"), c("G01","G02","G03"))
  expect_identical(see.group.data(g2,"P2D"), c(2L,3L,4L)+1L)
  capture_output(delete.group(g2))
  g2 <- make.targeted.crosses(c("G01","G02","G03"), c(2,3,4), track.pedigree=FALSE)
  expect_identical(see.group.data(g2,"P1"), c("0","0","0"))
  expect_identical(see.group.data(g2,"P2D"), c(0L,0L,0L))
  capture_output(delete.group(g2))
  #current_id <- 6 + 6
  
  # "retain" setting
  g2 <- make.targeted.crosses(c("G01","G02","G03"), c(2,3,4), retain=FALSE)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g),"GroupSize"=c(6L)))
  # current_id <- 6 + 6 + 3
  
  # "offspring" setting
  g2 <- make.targeted.crosses(c("G01","G02","G03"), c(2,3,4))
  g3 <- make.targeted.crosses(c("G01","G02","G03"), c(2,3,4), offspring=10)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3),"GroupSize"=c(6L,3L,30L)))
  expect_identical(unique(see.group.data(g3,"P1")), see.group.data(g2,"P1"))
  expect_identical(unique(see.group.data(g3,"P2")), see.group.data(g2,"P2"))
  capture_output(delete.group(c(g2,g3)))
  # current_id <- 6 + 6 + 3 + 33
  
  # "give.ids" setting
  g2 <- make.random.crosses(g, 10, give.ids = TRUE)
  g3 <- make.random.crosses(g, 10, give.ids = FALSE)
  expect_identical(see.group.data(g2,"D"), 49L:58L)
  expect_identical(see.group.data(g3,"D"), rep(0L,10))
  capture_output(delete.group(c(g2,g3)))
  
  # "give.names" and "name.prefix" setting
  g2 <- make.random.crosses(g, 10, give.names = FALSE)
  g3 <- make.random.crosses(g, 10, give.names = TRUE)
  g4 <- make.random.crosses(g, 10, give.names = TRUE, name.prefix = "p")
  g5 <- make.random.crosses(g, 10, give.names = FALSE, name.prefix = "unused")
  # should modify see.group.data("N") so we can distinguish between nameless g2 and number-named g3
  expect_identical(see.group.data(g2,"N"), as.character(59L:68L)) 
  expect_identical(see.group.data(g3,"N"), as.character(69L:78L))
  expect_identical(see.group.data(g4,"N"), paste0("p",79L:88L))
  expect_identical(see.group.data(g5,"N"), as.character(89L:98L))
  capture_output(delete.group(c(g2,g3,g4,g5)))
  
  # save-as-you-go file settings
  g2 <- make.clones(g,file.prefix="sayg1",save.genotype=TRUE, retain=FALSE)
  f_out <- readLines("sayg1-genotype.txt")
  expect_identical(length(f_out), 7L)
  expect_identical(f_out[1], "\tm1\tm2\tm3")
  expect_identical(f_out[2], "G01\tTT\tAA\tTT")
  expect_identical(f_out[3], "G02\tTT\tAA\tTT")
  expect_identical(f_out[4], "G03\tTT\tAA\tTA")
  expect_identical(f_out[5], "G04\tTA\tAA\tTA")
  expect_identical(f_out[6], "G05\tTT\tTT\tTT")
  expect_identical(f_out[7], "G06\tAT\tAA\tTT")
  file.remove("sayg1-genotype.txt")
  
  # current_id <- 98 + 6
  g2 <- make.clones(g,file.prefix="sayg2",save.gebv=TRUE, retain=FALSE)
  f_out <- readLines("sayg2-bv.txt")
  expect_identical(length(f_out), 6L)
  expect_identical(f_out[1], "105\tG01\t1.400000")
  expect_identical(f_out[2], "106\tG02\t1.400000")
  expect_identical(f_out[3], "107\tG03\t1.600000")
  expect_identical(f_out[4], "108\tG04\t-0.100000")
  expect_identical(f_out[5], "109\tG05\t0.600000")
  expect_identical(f_out[6], "110\tG06\t-0.300000")
  file.remove("sayg2-bv.txt")
  
  g2 <- make.doubled.haploids(g,file.prefix="sayg3",save.pedigree=TRUE, retain=FALSE)
  f_out <- readLines("sayg3-pedigree.txt")
  expect_identical(length(f_out), 6L)
  expect_identical(f_out[1], "111\t=(G01)")
  expect_identical(f_out[2], "112\t=(G02)")
  expect_identical(f_out[3], "113\t=(G03)")
  expect_identical(f_out[4], "114\t=(G04)")
  expect_identical(f_out[5], "115\t=(G05)")
  expect_identical(f_out[6], "116\t=(G06)")
  file.remove("sayg3-pedigree.txt")
  
  # Checking we can do multiple save-as-you-gos at once and they deal with nulls correctly (in names, in file prefixes, in ids)
  g2 <- make.clones(g, save.genotype=TRUE, save.gebv=TRUE, save.pedigree=TRUE, give.ids=FALSE, inherit.names=FALSE, give.names=FALSE)
  f_out <- readLines("out-genotype.txt")
  expect_identical(length(f_out), 7L)
  expect_identical(f_out[1], "\tm1\tm2\tm3")
  expect_identical(f_out[2], "0\tTT\tAA\tTT")
  expect_identical(f_out[3], "0\tTT\tAA\tTT")
  expect_identical(f_out[4], "0\tTT\tAA\tTA")
  expect_identical(f_out[5], "0\tTA\tAA\tTA")
  expect_identical(f_out[6], "0\tTT\tTT\tTT")
  expect_identical(f_out[7], "0\tAT\tAA\tTT")
  file.remove("out-genotype.txt")
  f_out <- readLines("out-bv.txt")
  expect_identical(length(f_out), 6L)
  expect_identical(f_out[1], "0\t\t1.400000")
  expect_identical(f_out[2], "0\t\t1.400000")
  expect_identical(f_out[3], "0\t\t1.600000")
  expect_identical(f_out[4], "0\t\t-0.100000")
  expect_identical(f_out[5], "0\t\t0.600000")
  expect_identical(f_out[6], "0\t\t-0.300000")
  file.remove("out-bv.txt")
  f_out <- readLines("out-pedigree.txt")
  expect_identical(length(f_out), 6L)
  expect_identical(f_out[1], "0\t=(G01)")
  expect_identical(f_out[2], "0\t=(G02)")
  expect_identical(f_out[3], "0\t=(G03)")
  expect_identical(f_out[4], "0\t=(G04)")
  expect_identical(f_out[5], "0\t=(G05)")
  expect_identical(f_out[6], "0\t=(G06)")
  file.remove("out-pedigree.txt")
  
  clear.simdata()
})

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
  g2 <- make.targeted.crosses(first.parents=c(0,1), second.parents=c(3,0))
  g3 <- make.targeted.crosses(first.parents=c("G01", "G02", "G03"), second.parents=c("G06","G05","G04"))
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3),"GroupSize"=c(6L,2L,3L)))
  expect_identical(see.group.data(c(g2,g3),"P1D"),c(1L,2L,1L,2L,3L))
  expect_identical(see.group.data(c(g2,g3),"P2"),c("G04","G01","G06","G05","G04"))
  
  # And test what happens if we give it bad input
  expect_snapshot(g4 <- make.targeted.crosses(first.parents=c(0,"abc"), second.parents=c(3,0)))
  expect_snapshot(g5 <- make.targeted.crosses(first.parents=c(0,10,3), second.parents=c(30,100,4)))
  expect_identical(g4, 0L)
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2,g3,g5),"GroupSize"=c(6L,2L,3L,1L)))
  expect_identical(see.group.data(g5,"P1D"),c(4L))
  expect_identical(see.group.data(g5,"P2D"),c(5L))
  
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
