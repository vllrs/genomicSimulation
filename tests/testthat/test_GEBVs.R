test_that("GEBVs are correctly calculated", {})

test_that("Local GEBVs are correctly calculated and saved", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  expect_equal(save.local.GEBVs("imaginary", "helper_blocks.txt", group=g),0)
  f_out <- readLines("imaginary")
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
})