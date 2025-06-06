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


test_that("package is loading genetic map correctly", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  # what was literally read from the file
  capture_output(map.df <- read.table("helper_map.txt",header=TRUE), print=F)
  map.df$chr <- as.character(map.df$chr)
  # Chromosome's first marker's position is "0"
  chrdata <- data.frame(chr=unique(map.df$chr))
  chrdata$min <- sapply(chrdata$chr, function(chr) min(map.df$pos[map.df$chr == chr]))
  chrdata$count <- sapply(chrdata$chr, function(chr) sum(map.df$chr == chr))
  map.df$pos <- map.df$pos - chrdata$min[match(map.df$chr,chrdata$chr)] # set first pos to zero
  map.df$pos[chrdata$count[match(map.df$chr,chrdata$chr)] == 1] <- NaN # set single-marker linkage groups to not have relative positions
  # order
  map.df <- map.df[order(map.df$chr,map.df$pos),]
  rownames(map.df) <- NULL
  
  expect_identical(see.genetic.map(init$mapID), map.df)
  #and also check this manually
  expect_identical(see.genetic.map(init$mapID)$marker, c("m1","m2","m3"))
  expect_equal(see.genetic.map(init$mapID)$chr, c('1','1','3'))
  expect_equal(see.genetic.map(init$mapID)$pos, c(0,3.1,NaN))
})

test_that("package is loading allele effects (with and without centring) correctly", {
  capture_output(init <- load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt"), print=F)
  
  capture_output(eff.df <- read.table("helper_eff.txt"), print=F)
  colnames(eff.df) <- c("marker","allele","eff")
  # Order correctly (I know that in helper_eff.txt, alphabetical order of marker names is same as genome order, and T is before A)
  eff.df <- eff.df[order(eff.df$marker,eff.df$allele,decreasing=c(FALSE,TRUE),method="radix"),]
  rownames(eff.df) <- NULL
  
  expect_identical(see.marker.effects(),eff.df)
  
  # And with centering?
  ctable <- data.frame(marker=c("m1","m2","m3"),allele="(centre)",eff=c(0.7,0,0.23))
  change.eff.set.centres(data.frame(markers=c("m1","m3"),cents=c(0.7,0.23)))
  eff.df2 <- rbind(eff.df,ctable)
  expect_identical(see.marker.effects(),eff.df2)
  
  ctable$eff <- c(0.9,0.99,0.999)
  change.eff.set.centres(ctable$eff, init$effectID)
  eff.df2 <- rbind(eff.df,ctable)
  expect_identical(see.marker.effects(init$effectID),eff.df2)
  
  change.eff.set.centres.of.allele("A", data.frame(markers=c("m1","m3"),cents=c(0.6,0.13)))
  ctable$eff <- c(0.6*-0.8,0,0.13*0.1)
  eff.df2 <- rbind(eff.df,ctable)
  expect_identical(see.marker.effects(init$effectID),eff.df2)
  
  change.eff.set.centres.of.allele("T", data.frame(markers=c("m1","m2"),cents=c(0.2,0.2)),reset.centres=F)
  ctable$eff <- c(0.6*-0.8 + 0.2*0.9, 0 + 0.2*-0.5, 0.13*0.1)
  eff.df2 <- rbind(eff.df,ctable)
  expect_identical(see.marker.effects(init$effectID),eff.df2)
  
})

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
  expect_error(delete.effects(1L),"Please load.data first")
  expect_error(delete.label(1L),"Please load.data first")
  expect_error(delete.map(1L),"Please load.data first")
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
  expect_error(create.markerblocks(matrix(c("m1",2),ncol=2)),"Please load.data first")
  expect_error(create.markerblocks.from.chrsplit(1),"Please load.data first")
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
  expect_error(see.GEBVs(g0$groupNum),"Need to load at least one set of marker effects before requesting breeding values")
  b <- create.markerblocks.from.chrsplit(2)
  expect_error(see.local.GEBVs(b,g0$groupNum),"Need to load at least one set of marker effects before requesting breeding values")
  expect_error(save.GEBVs("imaginary"),"Need to load at least one set of marker effects before requesting breeding values")
  expect_error(save.local.GEBVs("imaginary",b),"Need to load at least one set of marker effects before requesting breeding values")
  expect_error(break.group.by.GEBV(g0$groupNum, number=2),"Need to load effect values before running this function")
  expect_error(see.optimal.GEBV(),"Need to load effect values before running this function")
  expect_error(see.optimal.haplotype(),"Need to load effect values before running this function")
  expect_error(see.optimal.possible.GEBV(g0$groupNum),"Need to load effect values before running this function")
  expect_error(see.optimal.possible.haplotype(g0$groupNum),"Need to load effect values before running this function")
  expect_error(see.minimal.GEBV(),"Need to load effect values before running this function")
  
  clear.simdata()
})
