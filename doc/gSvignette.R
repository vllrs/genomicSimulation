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

# load the small initial data set
g0 <- load.data("../tests/testthat/helper_genotypes.txt",
                "../tests/testthat/helper_map.txt",
                "../tests/testthat/helper_eff.txt")

# Do many random crosses from the data set. 
f1 <- cross.randomly(g0, n.crosses=20, give.names=TRUE, name.prefix="F1.")

# Let's see what crossovers we got. For this few markers we don't get much information.
save.pedigrees("a.txt", f1, type="P")
find.plot.crossovers("a.txt", "b.txt")

# Make 40 random crosses from the 25% with the top GEBV
f1.selected <- select.by.gebv(f1, percentage=25)

# Delete groups we are not currently using, to free up some memory.
delete.group(f1)

f2 <- cross.randomly(f1.selected, n.crosses=5)
delete.group(f1.selected)

# Complete 4 rounds of selfing with save-as-you-go genotype saving
f6 <- self.n.times(f2, 4, file.prefix="af6", save.genotype=TRUE)

# Show current state of groups
print(c(g0, f2, f6))
see.existing.groups()

# Show the save-as-you-go output file
read.csv("af6-genome", sep='\t')


# (cleanup)
file.remove("a.txt")
file.remove("b.txt")
file.remove("af6-genome")

## -----------------------------------------------------------------------------
f <- cross.randomly(g0, n.crosses=5, offspring=3)

f.info <- see.group.gebvs(f)

# Create the masked values. In this case, a normally distributed term with sd corresponding
# to environmental variance is added to each GEBV
H2 <- 0.5  #estimate the heritability of the GEBV's trait at this step
Vg <- var(f.info$GEBV)
Ve <- Vg/H2 - Vg  # using the formula H2 = Vg/(Vg + Ve) to find the environmental variation
f.info$Pheno <- f.info$GEBV + rnorm(length(f.info$GEBV), mean=0, sd = sqrt(Ve))
  
# Sort the list by phenotypes
f.info <- f.info[order(f.info$Pheno, decreasing=TRUE),]
# Send the top 10 phenotypes to a new group
fselected <- make.group(f.info$i[1:10])

# (delete the non-selected genotypes)
delete.group(f)

# Then generate further generations from the selected genotypes
f2 <- cross.randomly(fselected, n.crosses=100)
# ...
delete.group(f2)

