## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(genomicSimulation)

## -----------------------------------------------------------------------------
# (Show the tiny example data set)
read.csv("../tests/testthat/helper_genotypes.txt", sep='\t', header=T)
read.csv("../tests/testthat/helper_map.txt", sep='\t', header=T)
read.csv("../tests/testthat/helper_eff.txt", sep='\t')

# Starting to use genomicSimulation: load the small initial data set
g0 <- load.data("../tests/testthat/helper_genotypes.txt",
                "../tests/testthat/helper_map.txt",
                "../tests/testthat/helper_eff.txt")

# Do random crosses from the progenitor lines. 
f1 <- cross.randomly(g0, n.crosses=20, give.names=TRUE, name.prefix="F1.")

#(Let's see what crossovers occured. For this few markers we don't get much information.)
save.pedigrees("a.txt", f1, type="P")
find.plot.crossovers("a.txt", "b.txt")

# Find the 25% with the top breeding value/GEBV
f1.selected <- select.by.gebv(f1, percentage=25)

# Delete groups we are not currently using, to free up some memory.
delete.group(f1)

# Make 40 random crosses from those selected 25%
f2 <- cross.randomly(f1.selected, n.crosses=5)
delete.group(f1.selected)

# Complete 4 rounds of selfing with save-as-you-go genotype saving
f6 <- self.n.times(f2, 4, file.prefix="af6", save.genotype=TRUE)

# Show current state of groups.
see.existing.groups()

# (These are the groups the above command should show)
print(c(g0, f2, f6))

# Show the save-as-you-go output file
read.csv("af6-genome", sep='\t')


# (cleanup)
file.remove("a.txt")
file.remove("b.txt")
file.remove("af6-genome")

## -----------------------------------------------------------------------------
get.top.10.phenotypes <- function(group.info, H2) {
  # calculate Ve
  Vg <- var(group.info$GEBV)
  Ve <- Vg/H2 - Vg
  
  # Simulate phenotypes
  group.info$pheno <- group.info$GEBV + rnorm(length(group.info$GEBV), mean=0, sd = sqrt(Ve))
  
  # Pick the best 10 and return their indices
  sorted <- group.info[order(group.info$pheno, decreasing=TRUE),]
  return(sorted$i[1:10])
}

## -----------------------------------------------------------------------------
g0 <- load.data("../tests/testthat/helper_genotypes.txt",
                "../tests/testthat/helper_map.txt",
                "../tests/testthat/helper_eff.txt")

# Simulate crosses
f1 <- cross.randomly(g0, n.crosses=15, offspring=3)
# Apply custom selection method
f1.info <- see.group.gebvs(f1)
f1.selected <- make.group(get.top.10.phenotypes(f1.info, 0.3))
# (delete the non-selected genotypes)
delete.group(f1)

# ... repeat for further generations
f2 <- cross.randomly(f1.selected, n.crosses=100)
f2.selected <- make.group(get.top.10.phenotypes(see.group.gebvs(f2), 0.5))

# ...
delete.group(f1.selected)
delete.group(f2)
delete.group(f2.selected)

