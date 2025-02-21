# How many of these can be tested without going into the C input?

test_that("package can load files", {
  expect_snapshot(load.data("helper_genotypes_long.txt", "helper_map.txt", "helper_eff_2.txt"))

  expect_snapshot(load.data("helper_genotypes.txt", "helper_map.txt"))
  expect_snapshot(load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"))
  
  expect_snapshot(load.genotypes("helper_genotypes.txt"))
  expect_snapshot(load.genotypes("helper_genotypes_long.txt"))
  
  expect_snapshot(load.effects("helper_eff_2.txt"))
  clear.simdata()
})

#test_that("package is loading genotypes correctly", {})

#test_that("package is loading genetic map correctly", {})

#test_that("package is loading allele effects correctly", {})

#test_that("package can be set up without an effect file", {})

#test_that("package effect file can be swapped out for another", {})

#test_that("data can be deleted on completion", {
#  load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt")
#  clear_simdata()
#})

test_that("package can load files with manually-specified formats", {
  # Goal is not to test that automatic file format detection works as expected. That's covered in C tests
  # Just want to show that, if something did go wrong in auto format detection,
  # the manual format specification would work. 
  
  expect_identical(define.matrix.format.details(), list())
  
  expect_snapshot(load.data("helper_genotypes_long.txt", "helper_map.txt", 
                            format=define.matrix.format.details(has.header=TRUE,
                                                                markers.as.rows=TRUE,
                                                                cell.style="pairs")))
  
  expect_snapshot(load.data("helper_genotypes_inverted_counts.txt", "helper_map.txt", 
                            format=define.matrix.format.details(has.header=TRUE,
                                                                markers.as.rows=FALSE)))
  
  expect_snapshot(load.data("helper_genotypes_noheader_nullable_counts.txt", "helper_map.txt", 
                            format=define.matrix.format.details(has.header=FALSE,
                                                                markers.as.rows=TRUE,
                                                                cell.style="C")))
  
  expect_snapshot(load.genotypes("helper_genotypes_slash_nocorner.txt",
                                 define.matrix.format.details(cell.style="/",markers.as.rows=TRUE)))
  
  clear.simdata()
  
})

test_that("package can create and delete labels", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  expect_identical(create.new.label(7L), 1L)
  expect_identical(create.new.label(-2), 2L)
  
  expect_output(change.label.default(c(1L,2L),c(2L,-4L)),
                "Set the defaults of 2 labels.")
  
  expect_identical(delete.label(1L), 0L)
  expect_identical(create.new.label(0L), 1L)
  
  clear.simdata()
})

test_that("package can load multiple effect sets", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  expect_identical(init$effectID, 1L)
  capture_output(eff2 <- load.effects("helper_eff_2.txt"), print=F)
  expect_identical(eff2, 2L)
  
  delete.effect.set(2L)
  
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt"), print=F)
  capture_output(eff3 <- load.effects("helper_eff.txt"), print=F)
  expect_identical(eff3, 1L)
  
  clear.simdata()
})

test_that("Other functions don't run without data being loaded", {
  clear.simdata()
  expect_error(break.group.into.families(1L),"Please load.data first")
  expect_error(break.group.into.individuals(1L),"Please load.data first")
  expect_error(combine.groups(c(1L,2L)),"Please load.data first")
  expect_error(make.all.unidirectional.crosses(1L),"Please load.data first")
  expect_error(make.targeted.crosses(1L,2L),"Please load.data first")
  expect_error(make.crosses.from.file("imaginary"),"Please load.data first")
  expect_error(make.double.crosses.from.file("imaginary"),"Please load.data first")
  expect_error(make.random.crosses(1L),"Please load.data first")
  expect_error(make.random.crosses.between(1L,2L),"Please load.data first")
  expect_error(delete.group(1L),"Please load.data first")
  expect_error(delete.effect.set(1L),"Please load.data first")
  expect_error(delete.label(1L),"Please load.data first")
  expect_error(delete.recombination.map(1L),"Please load.data first")
  expect_error(find.crossovers("imaginary", "imaginary2"),"Please load.data first")
  expect_error(find.plot.crossovers("imaginary", "imaginary2"),"Please load.data first")
  expect_error(load.effects("imaginary"),"Please load.data first")
  expect_error(load.genotypes("imaginary"),"Please load.data first")
  expect_error(load.map("imaginary"),"Please load.data first")
  expect_error(make.doubled.haploids(1L),"Please load.data first")
  expect_error(make.group(c(3L,4L,5L)),"Please load.data first")
  expect_error(save.allele.counts("imaginary", allele="A"),"Please load.data first")
  expect_error(save.GEBVs("imaginary"),"Please load.data first")
  expect_error(save.genotypes("imaginary"),"Please load.data first")
  expect_error(save.local.GEBVs.blocks.from.file("imaginary", "imaginary2"),"Please load.data first")
  expect_error(save.local.GEBVs.blocks.from.chrsplit("imaginary", 2),"Please load.data first")
  expect_error(save.pedigrees("imaginary"),"Please load.data first")
  expect_error(see.existing.groups(),"Please load.data first")
  expect_error(see.group.data(1L, "X"),"Please load.data first")
  expect_error(see.group.data(1L,"Bv"),"Please load.data first")
  expect_error(see.optimal.haplotype(),"Please load.data first")
  expect_error(see.optimal.GEBV(),"Please load.data first")
  expect_error(see.minimal.GEBV(),"Please load.data first")
  expect_error(break.group.by.GEBV(1L, number=2),"Please load.data first")
  expect_error(self.n.times(2L,1L),"Please load.data first")
  expect_error(self.n.times(2L,1L),"Please load.data first")
  expect_error(create.new.label(1L),"Please load.data first")
  expect_error(delete.label(1L),"Please load.data first")
  expect_error(change.label.to.values(1L,c(1,2,3),0),"Please load.data first")
  expect_error(change.label.to.this(1L,3L),"Please load.data first")
  expect_error(change.label.by.amount(2L, -2),"Please load.data first")
  expect_error(change.label.default(1L,0L),"Please load.data first")
  expect_error(change.names.to.values(c("A","B")),"Please load.data first")
  expect_error(change.allele.symbol("A","B"),"Please load.data first")
  expect_error(break.group.by.label.value(1L,3L),"Please load.data first")
  expect_error(break.group.by.label.range(1L,2L,4L),"Please load.data first")
  
  capture_output(g0 <- load.data("helper_genotypes.txt", "helper_map.txt"), print=F)
  expect_error(see.group.data(g0$groupNum,"BV"),"Need to load at least one set of marker effects before requesting breeding values")
  expect_error(save.GEBVs("imaginary"),"Need to load at least one set of marker effects before requesting breeding values")
  #expect_error(save.local.GEBVs.blocks.from.chrsplit("imaginary",3),"Need to load at least one set of marker effects before requesting breeding values")
  #expect_error(save.local.GEBVs.blocks.from.file("imaginary","im2"),"Need to load at least one set of marker effects before requesting breeding values")
  expect_error(break.group.by.GEBV(g0$groupNum, number=2),"Need to load effect values before running this function")
  expect_error(see.optimal.GEBV(),"Need to load effect values before running this function")
  expect_error(see.optimal.haplotype(),"Need to load effect values before running this function")
  expect_error(see.optimal.possible.GEBV(g0$groupNum),"Need to load effect values before running this function")
  expect_error(see.optimal.possible.haplotype(g0$groupNum),"Need to load effect values before running this function")
  expect_error(see.minimal.GEBV(),"Need to load effect values before running this function")
  
  clear.simdata()
})