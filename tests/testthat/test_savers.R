test_that("save.genome.model works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genome.model("imaginary")
  f_out <- readLines("imaginary")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "name\tchr\tpos\tA\tT")
  
  f_out_split <- scan(text=f_out[2], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m1")
  expect_identical(as.integer(f_out_split[2]), 1L)
  expect_identical(as.numeric(f_out_split[3]), 5.2)
  expect_identical(as.numeric(f_out_split[4]), -0.8)
  expect_identical(as.numeric(f_out_split[5]), 0.9)
  
  f_out_split <- scan(text=f_out[4], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m3")
  expect_identical(as.integer(f_out_split[2]), 3L)
  expect_identical(as.numeric(f_out_split[3]), 15)
  expect_identical(as.numeric(f_out_split[4]), 0.1)
  expect_identical(as.numeric(f_out_split[5]), -0.1)
  clear.simdata()
})

test_that("save.genotypes in regular format with group works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genotypes("imaginary2", group=g, type="R")
  f_out <- readLines("imaginary2")
  expect_identical(length(f_out), 7L)
  
  expect_identical(f_out[1], "1\tm1\tm2\tm3")
  expect_identical(f_out[2], "G01\tTT\tAA\tTT")
  expect_identical(f_out[5], "G04\tTA\tAA\tTA")
  expect_identical(f_out[7], "G06\tAT\tAA\tTT")
  clear.simdata()
})

test_that("save.genotypes in regular format without group works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genotypes("imaginary2", type="R")
  f_out <- readLines("imaginary2")
  expect_identical(length(f_out), 7L)
  
  expect_identical(f_out[1], "\tm1\tm2\tm3")
  expect_identical(f_out[2], "G01\tTT\tAA\tTT")
  expect_identical(f_out[5], "G04\tTA\tAA\tTA")
  expect_identical(f_out[7], "G06\tAT\tAA\tTT")
  clear.simdata()
})

test_that("save.genotypes in transposed format with group works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genotypes("imaginary3", group=g, type="T")
  f_out <- readLines("imaginary3")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "1\tG01\tG02\tG03\tG04\tG05\tG06")
  expect_identical(f_out[2], "m1\tTT\tTT\tTT\tTA\tTT\tAT")
  expect_identical(f_out[3], "m2\tAA\tAA\tAA\tAA\tTT\tAA")
  expect_identical(f_out[4], "m3\tTT\tTT\tTA\tTA\tTT\tTT")
  clear.simdata()
})

test_that("save.genotypes in transposed format without group works", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genotypes("imaginary3", type="T")
  f_out <- readLines("imaginary3")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "\tG01\tG02\tG03\tG04\tG05\tG06")
  expect_identical(f_out[2], "m1\tTT\tTT\tTT\tTA\tTT\tAT")
  expect_identical(f_out[3], "m2\tAA\tAA\tAA\tAA\tTT\tAA")
  expect_identical(f_out[4], "m3\tTT\tTT\tTA\tTA\tTT\tTT")
  clear.simdata()
})

test_that("save.allele.counts works with group", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.allele.counts("imaginary4", group=g, allele="T")
  f_out <- readLines("imaginary4")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "1\tG01\tG02\tG03\tG04\tG05\tG06")
  f_out_split <- scan(text=f_out[2], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m1")
  expect_equal(as.numeric(f_out_split[-1]), c(2,2,2,1,2,1))
  f_out_split <- scan(text=f_out[3], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m2")
  expect_equal(as.numeric(f_out_split[-1]), c(0,0,0,0,2,0))
  f_out_split <- scan(text=f_out[4], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m3")
  expect_equal(as.numeric(f_out_split[-1]), c(2,2,1,1,2,2))
  
  clear.simdata()
})

test_that("save.allele.counts works without group", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.allele.counts("imaginary4", allele="A")
  f_out <- readLines("imaginary4")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "\tG01\tG02\tG03\tG04\tG05\tG06")
  f_out_split <- scan(text=f_out[2], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m1")
  expect_equal(as.numeric(f_out_split[-1]), c(0,0,0,1,0,1))
  f_out_split <- scan(text=f_out[3], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m2")
  expect_equal(as.numeric(f_out_split[-1]), c(2,2,2,2,0,2))
  f_out_split <- scan(text=f_out[4], what=" ", quiet=TRUE)
  expect_identical(f_out_split[1], "m3")
  expect_equal(as.numeric(f_out_split[-1]), c(0,0,1,1,0,0))
  clear.simdata()
})

test_that("save.pedigrees in one-step format works", {})

test_that("save.pedigrees in recursive format works", {})