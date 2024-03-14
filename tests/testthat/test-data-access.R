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
  
  g2 <- make.targeted.crosses(c(0L,1L), c(2L,3L))
  expect_identical(see.group.data(g2, "P1"), c("G01", "G02"))
  expect_identical(see.group.data(g2, "P2"), c("G03", "G04"))
  g3 <- self.n.times(g2, 1)
  expect_identical(see.group.data(g3,"P2"), c("7", "8"))
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

