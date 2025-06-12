test_that("see.existing.groups works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  g2 <- make.random.crosses(g, n.crosses=5, offspring=1, give.names=TRUE, name.prefix="cr")
  
  expect_identical(see.existing.groups(), data.frame("Group"=c(g,g2),"GroupSize"=c(6L,5L)))
  
  clear.simdata()
})

test_that("see.group.data works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  expect_identical(see.group.data(g, "Names"), c("G01","G02","G03","G04","G05","G06"))
  
  expect_identical(see.group.data(g, "D (that's ID)"), c(1L,2L,3L,4L,5L,6L))
  
  expect_identical(see.group.data(g, "X for Index"), c(0L,1L,2L,3L,4L,5L))
  
  gs <- see.group.data(g, "Genotypes")
  expect_identical(substr(gs[1],0,6), "TTAATT")
  expect_identical(substr(gs[3],0,6), "TTAATA")
  expect_identical(substr(gs[5],0,6), "TTTTTT")
  expect_identical(substr(gs[6],0,6), "ATAATT")
  
  expect_equal(see.group.data(g, "BVs"), c(1.4,1.4,1.6,-0.1,0.6,-0.3))
  
  g2 <- make.targeted.crosses(c(0,1), c(2,3))
  expect_identical(see.group.data(g2, "P1"), c("G01", "G02"))
  expect_identical(see.group.data(g2, "P2"), c("G03", "G04"))
  expect_identical(see.group.data(g2, "P2D"), as.integer(c(3,4)))
  g3 <- self.n.times(g2, 1)
  expect_identical(see.group.data(g3,"P2D"), c(7L,8L))
  expect_identical(see.group.data(g3,"P1D"), c(7L,8L))
  expect_identical(see.group.data(g3,"pedigree"), c("9\t=(7=(G01,G03))", "10\t=(8=(G02,G04))"))
  
  # little bit of excessive label testing snuck in here
  expect_error(see.group.data(g,"Label"),"Need to create at least one custom label before requesting custom label values")
  l1 <- create.new.label(10)
  expect_error(see.group.data(g,"L",label=(l1+5)),"`label` parameter does not match a current existing custom label id")
  expect_identical(see.group.data(g,"L",label=l1),rep(10L,6))
  l2 <- create.new.label(-3)
  expect_identical(see.group.data(g,"L",label=l1),rep(10L,6))
  expect_identical(see.group.data(g,"L",label=l2),rep(-3L,6))
  change.label.by.amount(l2,2)
  expect_identical(see.group.data(g,"L",label=l2),rep(-1L,6))
  expect_identical(see.group.data(g,"L",label=l1),rep(10L,6))
  change.label.to.values(l1,1:5,group=g,startIndex=2)
  expect_identical(see.group.data(g,"L",label=l1),c(10L,1L,2L,3L,4L,5L))
  expect_identical(see.group.data(g,"L",label=l2),rep(-1L,6))
  
  clear.simdata()
})


test_that("see.group.data works with multiple groups", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  expect_snapshot(ntmp <- see.group.data(c(g,50), "Names"))
  expect_identical(ntmp, c("G01","G02","G03","G04","G05","G06"))
  
  g2 <- make.random.crosses(g,n.crosses=10)
  g3 <- make.random.crosses(g,n.crosses=2)
  
  expect_identical(see.group.data(c(g3,g2,g),"X"), c(16:17, 6:15, 0:5))
  
  capture_output(delete.group(g))
  
  expect_snapshot(ntmp <- see.group.data(c(g3,g,g2),"D"))
  expect_identical(ntmp, c(17:18, 7:16))
  
  clear.simdata()
})

test_that("See.group.gene.data works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  mt <- read.table("helper_genotypes.txt", header=TRUE)
  rownames(mt) <- mt[,1]
  mt <- mt[,-1]
  mt <- as.matrix(mt)
  
  mt2 <- see.group.gene.data(1L)
  
  expect_true(all.equal(mt["m1",],mt2[1,]))
  expect_true(all.equal(mt["m2",],mt2[2,]))
  expect_true(all.equal(mt["m3",],mt2[3,]))
  
  # Check works when there's missing alleles
  capture_output(g2 <- load.genotypes("helper_genotypes_noheader_nullable_counts.txt",
                                      define.matrix.format.details(has.header=FALSE)))
  
  mt2b <- see.group.gene.data(g2,unknown.allele = '?')
  expect_equal(rownames(mt2b),rownames(mt2))
  expect_equal(which(mt2b=="TT"),which(mt2=="TT"))
  expect_equal(which(mt2b=="AA"),which(mt2=="AA"))
  expect_equal(mt2b[1,4],"??")
  expect_equal(sum(mt2b=="??"),1)
  expect_equal(which(mt2b=="AT" | mt2b=="TA" | mt2b=="??"),which(mt2=="AT"|mt2=="TA"))
  
  expect_snapshot(tmp <- see.group.gene.data(g2,unknown.allele='\n'))
  
  # Check counts of one allele
  mt3 <- see.group.gene.data(1,"A")
  
  mtc <- matrix(0,nrow=dim(mt2)[1],ncol=dim(mt2)[2],dimnames=list(rownames(mt2),colnames(mt2)))
  mtc[mt2 == "AA"] <- 2
  mtc[mt2 == "AT" | mt2 == "TA"] <- 1
  
  expect_true(all.equal(mtc,mt3))
  
  # Check counts of both alleles always add up to 2
  expect_equal(see.group.gene.data(g,"A") + see.group.gene.data(g,"T"), 
               matrix(2,nrow=3,ncol=6,dimnames=list(rownames(mt2),colnames(mt2))))
  
  g2 <- make.random.crosses(g,n.crosses=13)
  
  # check counts of neither allele
  expect_equal(see.group.gene.data(g,"C"), matrix(0,nrow=3,ncol=6),ignore_attr=TRUE)
  
  # Check other and multi-group calls seem to act as expected
  expect_equal(dim(see.group.gene.data(g2)), c(3,13))
  expect_equal(dim(see.group.gene.data(c(g,g2))), c(3,13+6))
  
})

test_that("see.genetic.map works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  truth <- data.frame("marker"=c("m1","m2","m3"), "chr"=c('1','1','3'), "pos"=c(0,8.3-5.2,NaN))
  expect_equal(see.genetic.map(), truth)
  
  clear.simdata()
  capture_output(init2 <- load.data(map.file="helper_map2.txt"))
  truth2 <- data.frame("marker"=c("a1","and&a2"), "chr"=c('1B','1a'), "pos"=c(NaN,NaN))
  expect_equal(see.genetic.map(init2$map), truth2)
  
})

test_that("see.marker.effects works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  truth <- data.frame("marker"=c("m1","m1","m2","m2","m3","m3"), 
                      "allele"=c("T","A","T","A","T","A"), 
                      "eff"=c(0.9,-0.8,-0.5,-0.1,-0.1,0.1))
  expect_equal(see.marker.effects(), truth)
  expect_equal(see.marker.effects(init$effectID), truth)
  
})

