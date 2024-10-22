
# genomicSimulation (development version)

## New Features

- `see.group.data` can now be called with multiple groups. The results from each group will be concatenated. This allows you to replace lines like `c(see.group.data(group1,"BV"),see.group.data(group2,"BV"))` with `see.group.data(c(group1,group2),"BV")`.
- `see.group.gene.data` can now be called with multiple groups. The columns from each group will be concatenated.
- Added ability to observe the pedigree IDs of parents using `see.group.data`. (Previously, `see.group.data` could be used to observe the names of parents, which would fall back to the pedigree IDs of parents for parents without names. Some use cases may appreciate the ability to access pedigree IDs consistently.)
- Added ability to override the automatically detected file layout details in `load.data` and `load.genotypes` using the new function `define.matrix.format.details`, so that files whose format is incorrectly detected are not prevented from being loaded. 
- Make missing alleles visible in the matrix returned by `see.group.gene.data` (previously, they would cause issues with sometimes hiding non-missing alleles, due to being represented internally by the null character '\0').
- The orientation of the output matrix of `save.allele.counts` can now be chosen, as has previously been possible only with `save.genotypes`. 

## Improvements

- No longer `Depends` on package `fs`. `fs` is only used if file paths are used that contain tildes `~` needing to be expanded to home directories.
- Fix a potential (untested) issue when using very long R vectors, by updating vector length calculation in R wrapper functions from `length()` to `xlength()`.
- `see.group.data` and `see.group.gene.data` now raise a warning, instead of an error, when there are no group members allocated to the group(s) they access.
- `see.group.data`, `see.group.gene.data`, `delete.group`, and `combine.groups` no longer raise an error when the group numbers passed are of numeric type instead of integer type. Instead, as long as the numbers are whole numbers, they perform the type conversion. Likewise, `make.targeted.crosses` no longer raises an error when it is passed parent indexes that are of numeric type instead of integer type, as long as the numbers are whole numbers. Likewise, `make.random.crosses`, `make.all.unidirectional.crosses`, `self.n.times` and `make.doubled.haploids` no longer raise an error when the group numbers or map IDs passed are of numeric type instead of integer type, as long as the numbers are whole numbers. 
- Column names and row names of the matrix returned by `see.group.gene.data` now reflect genotype names and genetic marker names respectively.
- Vignette and function docs edited for better clarity around compatible input file formats and function name changes made in v0.2.5
- Under-the-hood improvements and simplifications of the genotype matrix input file layout detection systems.
- Under-the-hood improvements and simplifications of the file saving output functions.
- File saving output functions with multiple possible output formats (`save.genotypes` & `save.pedigrees`) now prefer a logical/boolean parameter to select the output format, the same way as the underlying C library does, instead of needing to memorise format strings. The old string parameters will still be accepted, so no need to modify old code.
- `break.group.into.buckets` now accepts integer bucket sizes even if they are not in R integer vector format.
- Functions which accept vectors of integers (eg. `break.group.into.buckets`, `make.group`, `change.label.to.values`) provide more meaningful error messages when passed a character vector instead.
- Use of superseded progression function names (eg. `cross.combinations` for `make.targeted.crosses`) now produces warnings. 
- Superseded progression function names in tests have been replaced with their current names. 
- Clarified in documentation which progression functions accept a vector of genetic map IDs, and which only accept a single genetic map.
- More informative error message when calling `break.group.by.GEBV` on a nonexistent group.

## Bug Fixes

- In `load.data` or `load.genotypes`, the number of markers per genotype that were successfully loaded is now accurately reflected in the printed log messages. Previously, these functions incorrectly printed out the total number of markers in the stored genetic map while claiming it was the number of markers in the genotype file that had been accurately matched to that map.
- Patched a couple of memory leaks in the underlying C library.
- Fixed a crash in `make.clones` with `inherit.names = TRUE` when genotypes being cloned had no names.
- Fix an infinite loop in `make.random.crosses.between` when both groups had breeding usage caps and one of the {caps x group sizes} was exactly the number of requested offspring.

# genomicSimulation 0.2.5

## New Features

- Some culling of (unused) dependencies means we can now release genomicSimulation (R version) under the same MIT license as the C version.
- Add ability to observe the values of custom labels using `see.group.data`.
- Add ability to load multiple recombination maps into simulation at a time. The associated new functions are `load.map` and `delete.recombination.map`.
- New function `change.allele.symbol` changes the internal representation of a particular allele. It does not modify corresponding marker effects. This may be useful for tidying printed output.

## Bug Fixes

- Fixed a bug where some genotypes lost their group numbers during internal data-structure tidying in large simulations.
- Input/output file paths can now include the tilde `~` character to represent the home directory, thanks to R package `fs`.

## Improvements

- Smarter and more robust file loaders. They can now automatically detect headers and allele encodings. See R version vignette (Section: Input Files) or Templates page of C version documentation for more information.
- Increased function naming consistency between R and C versions of the package. However, this means some R functions have changed names. Old function names still work (they directly call the new function name), but sometime in the future the old names may be removed. (Old name) -> (new recommended name) pairs for the R package are as follows:
	- see.minimum.GEBV ->            see.minimal.GEBV
	- select.by.gebv ->              break.group.by.GEBV
	- make.group.from.label ->       break.group.by.label.value
	- make.group.from.label.range -> break.group.by.label.range
	- make.label ->                  create.new.label
	- load.more.genotypes ->         load.genotypes
	- load.different.effects ->      load.effects
	- cross.randomly ->              make.random.crosses
	- cross.randomly.between ->      make.random.crosses.between
	- cross.combinations ->          make.targeted.crosses
	- cross.combinations.file ->     make.crosses.from.file
	- cross.dc.combinations.file ->  make.double.crosses.from.file
	- cross.all.pairs ->             make.all.unidirectional.crosses
	- save.local.GEBVs.by.file ->    save.local.GEBVs.blocks.from.file
	- save.local.GEBVs.by.chr ->     save.local.GEBVs.blocks.from.chrsplit
- genomicSimulation now calls R_Calloc/R_Free for memory allocation rather than stdlib's malloc/free.
- Under-the-hood improvements to crossing functions and to the internal tidying function.
- Same script and same random seed will produce different genotypes post-v0.2.4.003, because gametes are now generated successively (first one, then the other) rather than simultaneously.

# genomicSimulation 0.2.4.003

## Bug Fixes

- Group modification functions no longer crash R if there are more than 10 groups in existence.

## Improvements

- Parameter `maxgroups` in `get.existing.groups` removed. A call to `get.existing.groups` will now return all of the existing groups in the simulation.
- Under-the-hood improvements to group manipulation functions.

# genomicSimulation 0.2.4.002

## Bug Fixes
- A quick-fix to group modification functions crashing R if there are more than 10 groups in existence. Full fix coming soon.


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