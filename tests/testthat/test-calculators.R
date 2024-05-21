test_that("GEBVs are correctly calculated and shared with the see function", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  expect_equal(see.group.data(g, "bv", 1L), c(1.4,1.4,1.6,-0.1,0.6,-0.3))
  
  capture_output(eff2 <- load.effects("helper_eff_2.txt"), print=F)
  
  expect_equal(see.group.data(g,"bv", eff2), c(0.804,0.804,1.404,2.502,-0.696,1.902))
})

test_that("GEBVs are correctly calculated and shared with the save function", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  expect_equal(save.GEBVs("imaginar", group=g),0)
  f_out <- readLines("imaginar")
  expect_identical(length(f_out), 6L)
  
  f_out_split <- scan(text=f_out[1], what=" ", quiet=TRUE)
  expect_identical(as.integer(f_out_split[1]), 1L)
  expect_identical(f_out_split[2], "G01")
  expect_equal(as.numeric(f_out_split[3]), 1.4)
  
  f_out_split <- scan(text=f_out[4], what=" ", quiet=TRUE)
  expect_identical(as.integer(f_out_split[1]), 4L)
  expect_identical(f_out_split[2], "G04")
  expect_equal(as.numeric(f_out_split[3]), -0.1)
  
  f_out_split <- scan(text=f_out[6], what=" ", quiet=TRUE)
  expect_identical(as.integer(f_out_split[1]), 6L)
  expect_identical(f_out_split[2], "G06")
  expect_equal(as.numeric(f_out_split[3]), -0.3)
  
  file.remove("imaginar")
  clear.simdata()
})


test_that("Local GEBVs using blocks from file are correctly calculated and saved", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  expect_equal(save.local.GEBVs.blocks.from.file("imaginar", "helper_blocks.txt", group=g),0)
  f_out <- readLines("imaginar")
  expect_identical(length(f_out), 12L)
  
  f_out_split <- scan(text=f_out[1], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G01_1")
  expect_equal(as.numeric(f_out_split[2]), 0.8)
  expect_equal(as.numeric(f_out_split[3]), -0.1)
  expect_identical(length(f_out_split), 3L)
  
  f_out_split <- scan(text=f_out[2], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G01_2")
  expect_equal(as.numeric(f_out_split[2]), 0.8)
  expect_equal(as.numeric(f_out_split[3]), -0.1)
  expect_identical(length(f_out_split), 3L)
  
  f_out_split <- scan(text=f_out[8], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G04_2")
  expect_equal(as.numeric(f_out_split[2]), -0.9)
  expect_equal(as.numeric(f_out_split[3]), 0.1)
  expect_identical(length(f_out_split), 3L)
  
  f_out_split <- scan(text=f_out[9], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G05_1")
  expect_equal(as.numeric(f_out_split[2]), 0.4)
  expect_equal(as.numeric(f_out_split[3]), -0.1)
  expect_identical(length(f_out_split), 3L)
  
  file.remove("imaginar")
  clear.simdata()
})

test_that("Local GEBVs using blocks from slicing are correctly calculated and saved", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff_2.txt"), print=F)
  g <- init$groupNum
  
  expect_equal(save.local.GEBVs.blocks.from.chrsplit("imagina", 2, group=g),0)
  f_out <- readLines("imagina")
  expect_identical(length(f_out), 12L)
  
  f_out_split <- scan(text=f_out[1], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G01_1")
  expect_equal(as.numeric(f_out_split[2]), 0.002)
  expect_equal(as.numeric(f_out_split[3]), 0.7)
  expect_equal(as.numeric(f_out_split[4]), -0.3)
  expect_equal(as.numeric(f_out_split[5]), 0)
  expect_identical(length(f_out_split), 5L)
  
  f_out_split <- scan(text=f_out[8], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G04_2")
  expect_equal(as.numeric(f_out_split[2]), 1.1)
  expect_equal(as.numeric(f_out_split[3]), 0.7)
  expect_equal(as.numeric(f_out_split[4]), 0.3)
  expect_equal(as.numeric(f_out_split[5]), 0)
  expect_identical(length(f_out_split), 5L)
  
  f_out_split <- scan(text=f_out[9], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G05_1")
  expect_equal(as.numeric(f_out_split[2]), 0.002)
  expect_equal(as.numeric(f_out_split[3]), -0.05)
  expect_equal(as.numeric(f_out_split[4]), -0.3)
  expect_equal(as.numeric(f_out_split[5]), 0)
  expect_identical(length(f_out_split), 5L)
  
  f_out_split <- scan(text=f_out[12], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "G06_2")
  expect_equal(as.numeric(f_out_split[2]), 0.002)
  expect_equal(as.numeric(f_out_split[3]), 0.7)
  expect_equal(as.numeric(f_out_split[4]), -0.3)
  expect_equal(as.numeric(f_out_split[5]), 0)
  expect_identical(length(f_out_split), 5L)
  
  file.remove("imagina")
  clear.simdata()
})

test_that("Functions to see optimal genotype and GEBV work", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  capture_output(eff2 <- load.effects("helper_eff_2.txt"), print=F)
  
  expect_identical(see.optimal.haplotype(), "TAA")
  expect_equal(see.optimal.GEBV(), 1.8)
  expect_equal(see.minimal.GEBV(),-2.8)
  
  expect_identical(see.optimal.haplotype(eff2), "AAA")
  expect_equal(see.optimal.GEBV(eff2), 4.2)
  expect_equal(see.minimal.GEBV(eff2),-0.696)
  
  delete.effect.set(c(1L,2L))
  expect_error(see.optimal.GEBV()) # there are no effect sets left
  
  clear.simdata()
})

test_that("Functions to see optimal genotype and GEBV of a group work", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  expect_identical(see.optimal.haplotype(), "TAA")
  expect_equal(see.optimal.GEBV(), 1.8)
  expect_equal(see.minimal.GEBV(),-2.8)
  
  expect_identical(see.optimal.possible.haplotype(g), "TAA")
  expect_equal(see.optimal.possible.GEBV(g), 1.8)
  
  g2 <- make.group(c(4L,5L))
  expect_identical(see.optimal.possible.haplotype(g2), "TAT")
  expect_equal(see.optimal.possible.GEBV(g2), 1.4)
  
  clear.simdata()
})