test_that("save.genotypes in regular format with group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  save.genotypes("imaginary2", group=g)
  f_out <- readLines("imaginary2")
  expect_identical(length(f_out), 7L)
  
  expect_identical(f_out[1], "1\tm1\tm2\tm3")
  expect_identical(f_out[2], "G01\tTT\tAA\tTT")
  expect_identical(f_out[5], "G04\tTA\tAA\tTA")
  expect_identical(f_out[7], "G06\tAT\tAA\tTT")
  
  file.remove("imaginary2")
  clear.simdata()
})

test_that("save.genotypes in regular format without group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genotypes("imaginary2", markers.as.rows=FALSE)
  f_out <- readLines("imaginary2")
  expect_identical(length(f_out), 7L)
  
  expect_identical(f_out[1], "\tm1\tm2\tm3")
  expect_identical(f_out[2], "G01\tTT\tAA\tTT")
  expect_identical(f_out[5], "G04\tTA\tAA\tTA")
  expect_identical(f_out[7], "G06\tAT\tAA\tTT")
  
  file.remove("imaginary2")
  clear.simdata()
})

test_that("save.genotypes in transposed format with group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  save.genotypes("imaginary3", group=g, TRUE)
  f_out <- readLines("imaginary3")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "1\tG01\tG02\tG03\tG04\tG05\tG06")
  expect_identical(f_out[2], "m1\tTT\tTT\tTT\tTA\tTT\tAT")
  expect_identical(f_out[3], "m2\tAA\tAA\tAA\tAA\tTT\tAA")
  expect_identical(f_out[4], "m3\tTT\tTT\tTA\tTA\tTT\tTT")
  
  file.remove("imaginary3")
  clear.simdata()
})

test_that("save.genotypes in transposed format without group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  save.genotypes("imaginary3", markers.as.rows=TRUE)
  f_out <- readLines("imaginary3")
  expect_identical(length(f_out), 4L)
  
  expect_identical(f_out[1], "\tG01\tG02\tG03\tG04\tG05\tG06")
  expect_identical(f_out[2], "m1\tTT\tTT\tTT\tTA\tTT\tAT")
  expect_identical(f_out[3], "m2\tAA\tAA\tAA\tAA\tTT\tAA")
  expect_identical(f_out[4], "m3\tTT\tTT\tTA\tTA\tTT\tTT")
  
  file.remove("imaginary3")
  clear.simdata()
})

test_that("save.allele.counts works with group", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  
  save.allele.counts("imaginary4", group=g, allele="T")
  #expect_warning(f_out <- readLines("imaginary4")) #warning for incomplete final line-> no longer the case, printing has been standardised
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
  
  file.remove("imaginary4")
  clear.simdata()
})

test_that("save.allele.counts works without group", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
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
  
  file.remove("imaginary4")
  clear.simdata()
})

test_that("save.pedigrees in one-step format with group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  capture_output(g2 <- make.group(c(0L,1L,2L)), print=F)
  capture_output(delete.group(g), print=F)
  capture_output(f <- make.all.unidirectional.crosses(g2, give.names = T, give.ids = T, name.prefix = "F"), print=F)
  
  save.pedigrees("imaginary5", group=f, recursive.format = FALSE)
  f_out <- readLines("imaginary5")
  expect_identical(length(f_out), 3L)
  
  expect_identical(f_out[1], "F7\tG01\tG02")
  expect_identical(f_out[2], "F8\tG01\tG03")
  expect_identical(f_out[3], "F9\tG02\tG03")
  
  file.remove("imaginary5")
  clear.simdata()
})

test_that("save.pedigrees in one-step format without group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  capture_output(g2 <- make.group(c(0L,1L,2L)), print=F)
  capture_output(delete.group(g), print=F)
  capture_output(f <- make.all.unidirectional.crosses(g2, give.names = T, give.ids = T, name.prefix = "F"), print=F)
  
  save.pedigrees("imaginary5", recursive.format=FALSE)
  f_out <- readLines("imaginary5")
  expect_identical(length(f_out), 6L)
  
  expect_identical(f_out[1], "G01\t\t")
  expect_identical(f_out[3], "G03\t\t")
  expect_identical(f_out[4], "F7\tG01\tG02")
  expect_identical(f_out[5], "F8\tG01\tG03")
  expect_identical(f_out[6], "F9\tG02\tG03")
  
  file.remove("imaginary5")
  clear.simdata()
})

test_that("save.pedigrees in recursive format with group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  capture_output(g2 <- make.group(c(0L,1L,2L)), print=F)
  capture_output(delete.group(g), print=F)
  capture_output(f <- make.all.unidirectional.crosses(g2, give.names = T, give.ids = T, name.prefix = "F"), print=F)
  capture_output(fd <- make.all.unidirectional.crosses(f, give.names = T, give.ids = T, name.prefix = "D"), print=F)
  capture_output(fe <- make.targeted.crosses(1L, 1L, give.names = T, give.ids = T, name.prefix = "E"), print=F)
  capture_output(f <- combine.groups(c(f,fd,fe)), print=F)
  
  save.pedigrees("imaginary6", group=f, recursive.format = TRUE)
  f_out <- readLines("imaginary6")
  expect_identical(length(f_out), 7L)
  
  expect_identical(f_out[1], "7\tF7=(G01,G02)")
  expect_identical(f_out[2], "8\tF8=(G01,G03)")
  expect_identical(f_out[3], "9\tF9=(G02,G03)")
  expect_identical(f_out[4], "10\tD10=(F7=(G01,G02),F8=(G01,G03))")
  expect_identical(f_out[6], "12\tD12=(F8=(G01,G03),F9=(G02,G03))")
  expect_identical(f_out[7], "13\tE13=(G02)")
  
  file.remove("imaginary6")
  clear.simdata()
})

test_that("save.pedigrees in recursive format without group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  capture_output(g2 <- make.group(c(0L,1L,2L)), print=F)
  capture_output(delete.group(g), print=F)
  capture_output(f <- make.all.unidirectional.crosses(g2, give.names = T, give.ids = T, name.prefix = "F"), print=F)
  capture_output(fd <- make.all.unidirectional.crosses(f, give.names = T, give.ids = T, name.prefix = "D"), print=F)
  capture_output(fe <- make.targeted.crosses(1L, 1L, give.names = T, give.ids = T, name.prefix = "E"), print=F)
  
  save.pedigrees("imaginary6", recursive.format=TRUE)
  f_out <- readLines("imaginary6")
  expect_identical(length(f_out), 10L)
  
  expect_identical(f_out[1], "1\tG01")
  expect_identical(f_out[4], "7\tF7=(G01,G02)")
  expect_identical(f_out[9], "12\tD12=(F8=(G01,G03),F9=(G02,G03))")
  expect_identical(f_out[10], "13\tE13=(G02)")
  
  file.remove("imaginary6")
  clear.simdata()
})

test_that("save.GEBVs with group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  g2 <- make.clones(g, inherit.names=FALSE)
  
  save.GEBVs("imaginary7", group=g)
  f_out <- readLines("imaginary7")
  expect_identical(length(f_out), 6L)
  
  expect_identical(f_out[1], "1\tG01\t1.400000")
  expect_identical(f_out[2], "2\tG02\t1.400000")
  expect_identical(f_out[3], "3\tG03\t1.600000")
  expect_identical(f_out[4], "4\tG04\t-0.100000")
  expect_identical(f_out[5], "5\tG05\t0.600000")
  expect_identical(f_out[6], "6\tG06\t-0.300000")
  
  file.remove("imaginary7")
  clear.simdata()
})

test_that("save.GEBVs without group works", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  g <- init$groupNum
  g2 <- make.clones(g, inherit.names=FALSE)
  
  save.GEBVs("imaginary8")
  f_out <- readLines("imaginary8")
  expect_identical(length(f_out), 12L)
  
  expect_identical(f_out[1], "1\tG01\t1.400000")
  expect_identical(f_out[2], "2\tG02\t1.400000")
  expect_identical(f_out[3], "3\tG03\t1.600000")
  expect_identical(f_out[4], "4\tG04\t-0.100000")
  expect_identical(f_out[5], "5\tG05\t0.600000")
  expect_identical(f_out[6], "6\tG06\t-0.300000")
  expect_identical(f_out[7], "7\t\t1.400000")
  expect_identical(f_out[8], "8\t\t1.400000")
  expect_identical(f_out[9], "9\t\t1.600000")
  expect_identical(f_out[10], "10\t\t-0.100000")
  expect_identical(f_out[11], "11\t\t0.600000")
  expect_identical(f_out[12], "12\t\t-0.300000")
  
  file.remove("imaginary8")
  clear.simdata()
})

