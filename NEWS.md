# genomicSimulation 0.2.4

## New Features

- Add ability to load multiple effect files. Multiple sets of marker effects can now be available at one time. Loading an effect file will return a marker effect set identifier that represents that particular set of marker effects. Sets of marker effects can be removed from memory using the new function delete.effect.set.

## Improvements

- BREAKING CHANGE: When load.data is called with an effect file, the return value is now a list containing an element $groupNum (identical to the old return value) and an element $effectID (representing the identifier of the marker effect set loaded). 
- The save-as-you-go setting for saving breeding values in GenOptions now takes the identifier for a set of marker effects, if you wish to save the breeding values from that set of marker effects, or 0, to not use the save-as-you-go functionality.
 
## Bug Fixes
- Standardised file-output functions so that genotype names are consistently substituted with their PedigreeIDs, if the genotype does not have a name. 
- Stopped save-as-you-go genotype saving repeating the header row (of genetic marker names) every 1000 rows. The header row now only appears in the first row.
- Semicolon separators appear correctly in the output of save_marker_blocks
- Pedigrees are printed consistently between group-specific and whole-simulation variants of the same function. Parents are no longer skipped in favour of starting with grandparents in group-specific functions.


# genomicSimulation 0.2.3.002

## Bug Fixes

- genomicSimulation had a chance of crashing on load.data for certain marker effect files. The chance was higher for effect files listing few markers or listing many alleles. This release is a quick-fix for this bug.

# genomicSimulation 0.2.3 

## Improvements

- Added new functions to calculate the optimal haplotype and breeding value possible using only the alleles present in a particular group. 
	- `see.optimal.possible.GEBV`
	- `see.optimal.possible.haplotype`
- Increased name reading buffer from 30 to 45 as a temporary fix to allow longer marker names to be loaded.
- Added custom labels. Once a label is created, every genotype in the simulation has some (integer) value for that label. Users can create, set, and update the values of as many labels as desired. These can be used to track age or other custom subcategorisations.
	- `make.label`
	- `delete.label`
	- `change.label.to.values`
	- `change.label.to.this`
	- `change.label.by.amount`
	- `change.label.default`
	- Custom labels offer an easy way to have subcategories in a group. Those subcategories can be split from the group for individual examination with the following convenient functions, then merged back in:
		- `make.group.from.label`
		- `make.group.from.label.range`
- Added a name setter function, for users who want control over the names of their simulated genotypes.
	- `change.names.to.values`
- Updated internals of `see.group.data` and `cross.combinations` to match changes to C version.

# genomicSimulation 0.2.2

## Bug Fixes
- Removed a "negative ID" warning when trying to access the names of parents of genotypes that have no known parents
- Fixed a bug in many of the `break.` family functions which incorrectly identified some existing groups as empty.
- Fixed a bug in the crossing functions that resulted in crashes when `give.ids = FALSE`
- Fixed a bug in `make.clones` that resulted in crashes when trying to inherit the name (`inherit.names = TRUE`) of a parent with no name.

# genomicSimulation 0.2.1

## Improvements 

- `delete.group` can now be passed a vector of group ids to delete from memory. Previously it could only delete a single group at a time. 
- Added new group splitting functions:
	- `break.group.into.halfsib.families`
	- `break.group.evenly`
	- `break.group.into.buckets`
	- `break.group.randomly`
	- `break.group.with.probabilities`
- Added function `cross.randomly.between` to perform crosses where one parent is picked randomly from one group and the other from another group.
- Added option to have a cap on the number of uses of each group member as a parent of a cross in `cross.randomly` and `cross.randomly.between`. 
- Added function `make.clones` to clone or duplicate members of a group.

## Bug Fixes

- Fix segfault in select.by.gebv when trying to select more individuals than exist in the group (e.g. asking for the best 5 members of a group of 2). Now, it just moves all group members to the new selected group, and doesn't worry about the missing requested remainder.
- Fix see.group.data so it stops rather than requesting a 0-length block of heap space when it is asked to investigate a nonexistent group.
- Other minor bug fixes in the underlying C library.


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