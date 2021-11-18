- adding new group splitter functions. Added tests for them.

# genomicSimulation 0.2

- Introduces new version numbering scheme. This update would be numbered 0.1-5 under the old numbering scheme instead of 0.2, but it is not a more significant update than v0.1-3 -> v0.1-4. The new numbering scheme is as suggested by Hadley & Bryan in *R Packages*. It uses the second digit rather than the third for minor releases.

## Improvements

- Edited vignette for clarity. Improved vignette diagrams.
- New function `see.minimum.GEBV` returns minimum possible breeding value for the loaded set of allele effects. It is the counterpart to `see.optimal.GEBV`.
- Implement breeding value calculator speedups. Breeding value calculations are now up to 3x faster.
- Add new data types as possible outputs to `see.group.data`: "B" gives breeding values, "P1" and "P2" give the names (or IDs if none) of parents 1 and 2, "PED" gives the full known pedigree as a string. Furthermore, the preferred way to get indexes is now "X", not "I", for reducing confusion between getting ID and index data.

## Breaking changes

- `cross` is deprecated because it causes too many namespace collision issues with other packages, and its functionality is already covered by `cross.combinations`.
- `see.group.gebvs` is deprecated for no longer being anything more than a sugar wrapper around `data.frame("i"=see.group.data(...,"X"),"GEBV"=see.group.data(...,"bv"))`, now that `see.group.data` can deal with calculated datatypes like breeding values. The data.frame syntax is encouraged because it is more flexible/extensible for users who need to incorporate more types of data.


# genomicSimulation 0.1-4

## Improvements

- Add ability to load space-separated SNP matrix in addition to tab-separated
- Add local GEBV and block calculator based on chromosome lengths
- Optimal GEBV and optimal haplotype calculator functions added

## Bug fixes

- Inconsistencies in pedigree saving formats fixed
- Increased test coverage

## Breaking changes

- `cross.combinations` split into `cross.combinations.file` and `cross.combinations`
- `cross.dc.combinations` renamed to `cross.dc.combinations.file` for consistency
- `save.local.GEBVs` split into to `save.local.GEBVs.by.chr` and `save.local.GEBVs.by.file` 


# genomicSimulation 0.1-3

## Improvements

- Added short automated tests using testthat
- Removed compiler warnings
- Removed most R CMD check warnings and notes
- C code brought inline with sister repo genomicSimulationC, which ensures bugs in C code will be found faster.

## Bug fixes 

- Local GEBV calculator now only requests reasonable amounts of memory, so can be run outside of HPC environment
- Ability to load files with more than 1000 genotypes has been added.
- No additional genotypes disappear on deleting group number 1


# genomicSimulation 0.1-2

- Fixed bugs in loader functions and pedigree saver functions that caused R to crash.
- Add a function to calculate local GEBVs/block effect values.
- Updated vignette and documentation.


# genomicSimulation 0.1-1

- Initial release