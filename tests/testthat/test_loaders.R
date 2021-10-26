# How many of these can be tested without going into the C input?

test_that("package can load files", {
  expect_output(load.data("helper_genotypes_long.txt", "helper_map.txt", "helper_eff_2.txt"),
                "1003 genotypes of 3 markers were loaded. 0 pairs of alleles could not be loaded\n3 markers with map positions. 0 markers remain unmapped.\n6 effect values spanning 2 alleles loaded.")
  expect_output(load.data("helper_genotypes.txt", "helper_map.txt"),
                "6 genotypes of 3 markers were loaded. 0 pairs of alleles could not be loaded\n3 markers with map positions. 0 markers remain unmapped." )
  expect_output(load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"),
                "6 genotypes of 3 markers were loaded. 0 pairs of alleles could not be loaded\n3 markers with map positions. 0 markers remain unmapped.\n6 effect values spanning 2 alleles loaded." )
  
  expect_output(load.more.genotypes("helper_genotypes.txt"),
                "6 genotypes were loaded.")
  expect_output(load.more.genotypes("helper_genotypes_long.txt"),
                "1003 genotypes were loaded")
  
  expect_output(load.different.effects("helper_eff_2.txt"),
                "6 effect values spanning 2 alleles loaded.")
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

test_that("Other functions don't run without data being loaded", {
  clear.simdata()
  expect_error(break.group.into.families(1L),"Please load.data first")
  expect_error(break.group.into.individuals(1L),"Please load.data first")
  expect_error(combine.groups(c(1L,2L)),"Please load.data first")
  expect_error(cross.all.pairs(1L),"Please load.data first")
  expect_error(cross.combinations(1L,2L),"Please load.data first")
  expect_error(cross.combinations.file("imaginary"),"Please load.data first")
  expect_error(cross.dc.combinations.file("imaginary"),"Please load.data first")
  expect_error(cross.randomly(1L),"Please load.data first")
  expect_error(delete.group(1L),"Please load.data first")
  expect_error(find.crossovers("imaginary", "imaginary2"),"Please load.data first")
  expect_error(find.plot.crossovers("imaginary", "imaginary2"),"Please load.data first")
  expect_error(load.different.effects("imaginary"),"Please load.data first")
  expect_error(load.more.genotypes("imaginary"),"Please load.data first")
  expect_error(make.doubled.haploids(1L),"Please load.data first")
  expect_error(make.group(c(3L,4L,5L)),"Please load.data first")
  expect_error(save.allele.counts("imaginary", allele="A"),"Please load.data first")
  expect_error(save.GEBVs("imaginary"),"Please load.data first")
  expect_error(save.genome.model("imaginary"),"Please load.data first")
  expect_error(save.genotypes("imaginary"),"Please load.data first")
  expect_error(save.local.GEBVs.by.file("imaginary", "imaginary2"),"Please load.data first")
  expect_error(save.local.GEBVs.by.chr("imaginary", 2),"Please load.data first")
  expect_error(save.pedigrees("imaginary"),"Please load.data first")
  expect_error(see.existing.groups(),"Please load.data first")
  expect_error(see.group.data(1L, "X"),"Please load.data first")
  expect_error(see.group.data(1L,"Bv"),"Please load.data first")
  expect_error(see.optimal.haplotype(),"Please load.data first")
  expect_error(see.optimal.GEBV(),"Please load.data first")
  expect_error(select.by.gebv(1L, number=2),"Please load.data first")
  expect_error(self.n.times(2L,1L),"Please load.data first")
  
  capture_output(g0 <- load.data("helper_genotypes.txt", "helper_map.txt"), print=F)
  expect_error(see.group.data(g0,"BV"),"Need to load effect values before running this function")
  expect_error(save.GEBVs("imaginary"),"Need to load effect values before running this function")
  expect_error(select.by.gebv(g0, number=2),"Need to load effect values before running this function")
  
  clear.simdata()
})