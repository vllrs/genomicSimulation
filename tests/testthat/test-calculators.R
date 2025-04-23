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

test_that("Functions to calculate local GEBVs work", {
  capture_output(g <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  capture_output(eff2 <- load.effects("helper_eff_2.txt"))
  
  b1 <- create.markerblocks.from.chrsplit(3)
  b1.dt <- see.markerblocks(b1)
  expect_equal(dim(b1.dt), c(3,2))
  expect_equal(b1.dt$marker, c("m1","m2","m3"))
  expect_equal(b1.dt$block, c(1L,3L,4L))

  b1.scores <- see.local.GEBVs(b1,1L)
  expect_equal(dim(b1.scores), c(2*6,2*3))
  expect_equal(b1.scores[1,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  expect_equal(b1.scores[2,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  expect_equal(b1.scores[3,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  expect_equal(b1.scores[4,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  expect_equal(b1.scores[5,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  expect_equal(b1.scores[6,], c(0.9,0,-0.1,0.1,0,0)) # TAA
  expect_equal(b1.scores[7,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  expect_equal(b1.scores[8,], c(-0.8,0,-0.1,0.1,0,0)) # AAA
  expect_equal(b1.scores[9,], c(0.9,0,-0.5,-0.1,0,0)) # TTT
  expect_equal(b1.scores[10,], c(0.9,0,-0.5,-0.1,0,0)) # TTT
  expect_equal(b1.scores[11,], c(-0.8,0,-0.1,-0.1,0,0)) # AAT
  expect_equal(b1.scores[12,], c(0.9,0,-0.1,-0.1,0,0)) # TAT
  
  b1.scorese2 <- see.local.GEBVs(b1,1L,eff2)
  expect_equal(dim(b1.scorese2), c(2*6,2*3))
  expect_equal(b1.scorese2[1,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  expect_equal(b1.scorese2[2,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  expect_equal(b1.scorese2[3,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  expect_equal(b1.scorese2[4,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  expect_equal(b1.scorese2[5,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  expect_equal(b1.scorese2[6,], c(2e-03,0,0.7,0.3,0,0)) # TAA
  expect_equal(b1.scorese2[7,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  expect_equal(b1.scorese2[8,], c(1.1,0,0.7,0.3,0,0)) # AAA
  expect_equal(b1.scorese2[9,], c(2e-03,0,-5e-02,-0.3,0,0)) # TTT
  expect_equal(b1.scorese2[10,], c(2e-03,0,-5e-02,-0.3,0,0)) # TTT
  expect_equal(b1.scorese2[11,], c(1.1,0,0.7,-0.3,0,0)) # AAT
  expect_equal(b1.scorese2[12,], c(2e-03,0,0.7,-0.3,0,0)) # TAT
  
  markerscores <- b1.scores[,c(1,3,4)]
  
  save.local.GEBVs("imaginarlg", b1, group=1)
  f_out <- readLines("imaginarlg")
  expect_equal(length(f_out), 2*6)
  for (row in 1:12) {
    f_out_split <- scan(text=f_out[row], what=" ", quiet=TRUE)
    expect_equal(length(f_out_split),(2*3)+1)
    expect_equal(f_out_split[1],paste0("G0", (row + 1) %/% 2, "_", (row + 1) %% 2 + 1))
    expect_equal(as.numeric(f_out_split[2:length(f_out_split)]), b1.scores[row,])
  }
  file.remove("imaginarlg")
  
  save.local.GEBVs("imaginarlg", b1, 1, eff2)
  f_out <- readLines("imaginarlg")
  expect_equal(length(f_out), 2*6)
  for (row in 1:12) {
    f_out_split <- scan(text=f_out[row], what=" ", quiet=TRUE)
    expect_equal(length(f_out_split),(2*3)+1)
    expect_equal(f_out_split[1],paste0("G0", (row + 1) %/% 2, "_", (row + 1) %% 2 + 1))
    expect_equal(as.numeric(f_out_split[2:length(f_out_split)]), b1.scorese2[row,])
  }
  file.remove("imaginarlg")
  
  # ---------------

  b1 <- create.markerblocks.from.chrsplit(1)
  b1.dt <- see.markerblocks(b1)
  expect_equal(dim(b1.dt), c(3,2))
  expect_equal(b1.dt$marker, c("m1","m2","m3"))
  expect_equal(b1.dt$block, c(1L,1L,2L))
  
  b1.scores <- see.local.GEBVs(b1,1L,1L)
  expect_equal(dim(b1.scores), c(2*6,2*1))
  for (ii in 1:12) {
    expect_equal(b1.scores[ii,], c(sum(markerscores[ii,c(1,2)]), markerscores[ii,3]))
  }
  
  save.local.GEBVs("imaginarlg", b1)
  f_out <- readLines("imaginarlg")
  expect_equal(length(f_out), 2*6)
  for (row in 1:12) {
    f_out_split <- scan(text=f_out[row], what=" ", quiet=TRUE)
    expect_equal(length(f_out_split),(2*1)+1)
    expect_equal(f_out_split[1],paste0("G0", (row + 1) %/% 2, "_", (row + 1) %% 2 + 1))
    expect_equal(as.numeric(f_out_split[2:length(f_out_split)]), b1.scores[row,])
  }
  file.remove("imaginarlg")
  
  # ---------------
  
  blocks <- data.frame("m"=c("m1","m3","m1","m2","m1213","235","m3","m1","m3","m2","m1","abc"), 
                       "b"=c(4,4,1,6,2,2,2,6,6,1,1,3))
  b2 <- create.markerblocks(blocks)
  b2.dt <- see.markerblocks(b2)
  expect_equal(dim(b2.dt),c(9,2))
  expect_equal(b2.dt$marker, c("m1","m2","m1","m3","m1","m3","m2","m1","m3"))
  expect_equal(b2.dt$block, c(1L,1L,1L,2L,3L,3L,4L,4L,4L))
  
  b2.scores <- see.local.GEBVs(b2)
  expect_equal(dim(b2.scores), c(2*6,4))
  for (ii in 1:12) {
    expect_equal(b2.scores[ii,], c(sum(markerscores[ii,c(1,2,1)]), markerscores[ii,3],
                                   sum(markerscores[ii,c(1,3)]), sum(markerscores[ii,])))
  }
  
  save.local.GEBVs("imaginarlg", b2)
  f_out <- readLines("imaginarlg")
  expect_equal(length(f_out), 2*6)
  for (row in 1:12) {
    f_out_split <- scan(text=f_out[row], what=" ", quiet=TRUE)
    expect_equal(length(f_out_split),4+1)
    expect_equal(f_out_split[1],paste0("G0", (row + 1) %/% 2, "_", (row + 1) %% 2 + 1))
    expect_equal(as.numeric(f_out_split[2:length(f_out_split)]), b2.scores[row,])
  }
  file.remove("imaginarlg")
  
  # ---------------
  
  expect_error(see.markerblocks(1L))
  expect_error(see.markerblocks(sim.data$p))
  
})

test_that("Functions to calculate local GEBVs of multiple groups work", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  trio <- make.group(c(0L,1L,4L))
  g2 <- make.doubled.haploids(trio) # make 3
  g3 <- make.targeted.crosses("G04","G06") # make 1
  
  b1 <- create.markerblocks.from.chrsplit(3)
  b1.dt <- see.markerblocks(b1)
  expect_equal(dim(b1.dt), c(3,2))
  expect_equal(b1.dt$marker, c("m1","m2","m3"))
  expect_equal(b1.dt$block, c(1L,3L,4L))
  
  g_expected_scores <- matrix(c(0.9,0,-0.1,-0.1,0,0),nrow=1) # TAT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,-0.1,0,0)) # TAT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,-0.1,0,0)) # TAT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,-0.1,0,0)) # TAT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,-0.1,0,0)) # TAT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,0.1,0,0)) # TAA
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,-0.1,0,0)) # TAT
  g_expected_scores <- rbind(g_expected_scores, c(-0.8,0,-0.1,0.1,0,0)) # AAA
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.5,-0.1,0,0)) # TTT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.5,-0.1,0,0)) # TTT
  g_expected_scores <- rbind(g_expected_scores, c(-0.8,0,-0.1,-0.1,0,0)) # AAT
  g_expected_scores <- rbind(g_expected_scores, c(0.9,0,-0.1,-0.1,0,0)) # TAT
  
  trio_expected_scores <- g_expected_scores[c(1,2,3,4,9,10),]
  
  g2_expected_scores <- g_expected_scores[c(1,1,3,3,9,9),]
  
  g3_1_possible_scores <- g_expected_scores[c(7,8,6,11),]
  g3_2_possible_scores <- g_expected_scores[c(11,12),]
  
  b1.scores <- see.local.GEBVs(b1,g2)
  expect_equal(dim(b1.scores), c(2*(3),2*3))
  expect_equal(b1.scores, g2_expected_scores)
  
  # # at the moment save.local.gebvs does not do multiple groups
  # save.local.GEBVs("imaginarlg", b1, g2)
  # f_out <- readLines("imaginarlg")
  # expect_equal(length(f_out), 2*3)
  # for (row in 1:6) {
  #   f_out_split <- scan(text=f_out[row], what=" ", quiet=TRUE)
  #   expect_equal(length(f_out_split),(2*3)+1)
  #   expect_equal(f_out_split[1],paste0("_", (row + 1) %% 2 + 1))
  #   expect_equal(as.numeric(f_out_split[2:length(f_out_split)]), b1.scores[row,])
  # }
  # file.remove("imaginarlg")
  
  # -----------------
  capture_output(eff2 <- load.effects("helper_eff_2.txt"))
  
  b1.scores <- see.local.GEBVs(b1,c(trio,g2))
  expect_equal(dim(b1.scores), c(2*(3+3),2*3))
  expect_equal(b1.scores, rbind(trio_expected_scores,g2_expected_scores))
  
  # save.local.GEBVs("imaginarlg", b1, c(trio,g2))
  # f_out <- readLines("imaginarlg")
  # expect_equal(length(f_out), 2*6)
  # for (row in 1:12) {
  #   f_out_split <- scan(text=f_out[row], what=" ", quiet=TRUE)
  #   expect_equal(length(f_out_split),(2*3)+1)
  #   expect_equal(as.numeric(f_out_split[2:length(f_out_split)]), b1.scores[row,])
  # }
  # file.remove("imaginarlg")
  
  # --------------------
  
  b1.scores <- see.local.GEBVs(b1, eff.set=init$effectID)
  expect_equal(dim(b1.scores), c(2*(6+3+1),2*3))
  expect_equal(b1.scores[1:18,], rbind(g_expected_scores,g2_expected_scores))
  expect_equal(sum(apply(g3_1_possible_scores,1,function(x) all(x == b1.scores[19,]))), 1)
  expect_equal(sum(apply(g3_2_possible_scores,1,function(x) all(x == b1.scores[20,]))), 1)
  
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