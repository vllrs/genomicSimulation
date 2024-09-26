# package can load files

    Code
      load.data("helper_genotypes_long.txt", "helper_map.txt", "helper_eff_2.txt")
    Output
      (Loading helper_map.txt) Format: map file with header
      (Loading helper_map.txt) 3 marker(s) with map positions were loaded. Failed to parse 0 line(s).
      (Loading helper_genotypes_long.txt) Format axis: genetic markers are -rows-, founder lines are |columns|
      (Loading helper_genotypes_long.txt) Allele format: phased allele pairs
      (Loading helper_genotypes_long.txt) Format: genotype matrix with header row
      (Loading helper_genotypes_long.txt) 1003 genotype(s) of 3 marker(s) were loaded.
      (Loading helper_eff_2.txt) Format: effect file without header
      (Loading helper_eff_2.txt) 6 effect value(s) spanning 2 allele(s) were loaded. Failed to parse 0 line(s).
      $groupNum
      [1] 1
      
      $mapID
      [1] 1
      
      $effectID
      [1] 1
      

---

    Code
      load.data("helper_genotypes.txt", "helper_map.txt")
    Output
      (Loading helper_map.txt) Format: map file with header
      (Loading helper_map.txt) 3 marker(s) with map positions were loaded. Failed to parse 0 line(s).
      (Loading helper_genotypes.txt) Format axis: genetic markers are -rows-, founder lines are |columns|
      (Loading helper_genotypes.txt) Allele format: phased allele pairs
      (Loading helper_genotypes.txt) Format: genotype matrix with header row
      (Loading helper_genotypes.txt) 6 genotype(s) of 3 marker(s) were loaded.
      $groupNum
      [1] 1
      
      $mapID
      [1] 1
      
      $effectID
      [1] NA
      

---

    Code
      load.data("helper_genotypes.txt", "helper_map.txt", "helper_eff.txt")
    Output
      (Loading helper_map.txt) Format: map file with header
      (Loading helper_map.txt) 3 marker(s) with map positions were loaded. Failed to parse 0 line(s).
      (Loading helper_genotypes.txt) Format axis: genetic markers are -rows-, founder lines are |columns|
      (Loading helper_genotypes.txt) Allele format: phased allele pairs
      (Loading helper_genotypes.txt) Format: genotype matrix with header row
      (Loading helper_genotypes.txt) 6 genotype(s) of 3 marker(s) were loaded.
      (Loading helper_eff.txt) Format: effect file without header
      (Loading helper_eff.txt) 6 effect value(s) spanning 2 allele(s) were loaded. Failed to parse 0 line(s).
      $groupNum
      [1] 1
      
      $mapID
      [1] 1
      
      $effectID
      [1] 1
      

---

    Code
      load.genotypes("helper_genotypes.txt")
    Output
      (Loading helper_genotypes.txt) Format axis: genetic markers are -rows-, founder lines are |columns|
      (Loading helper_genotypes.txt) Allele format: phased allele pairs
      (Loading helper_genotypes.txt) Format: genotype matrix with header row
      (Loading helper_genotypes.txt) 6 genotype(s) of 3 marker(s) were loaded.
      [1] 2

---

    Code
      load.genotypes("helper_genotypes_long.txt")
    Output
      (Loading helper_genotypes_long.txt) Format axis: genetic markers are -rows-, founder lines are |columns|
      (Loading helper_genotypes_long.txt) Allele format: phased allele pairs
      (Loading helper_genotypes_long.txt) Format: genotype matrix with header row
      (Loading helper_genotypes_long.txt) 1003 genotype(s) of 3 marker(s) were loaded.
      [1] 3

---

    Code
      load.effects("helper_eff_2.txt")
    Output
      (Loading helper_eff_2.txt) Format: effect file without header
      (Loading helper_eff_2.txt) 6 effect value(s) spanning 2 allele(s) were loaded. Failed to parse 0 line(s).
      [1] 2

# package can load files with manually-specified formats

    Code
      load.data("helper_genotypes_long.txt", "helper_map.txt", format = define.matrix.format.details(
        has.header = TRUE, markers.as.rows = TRUE, cell.style = "pairs"))
    Output
      (Loading helper_map.txt) Format: map file with header
      (Loading helper_map.txt) 3 marker(s) with map positions were loaded. Failed to parse 0 line(s).
      (Loading helper_genotypes_long.txt) 1003 genotype(s) of 3 marker(s) were loaded.
      $groupNum
      [1] 1
      
      $mapID
      [1] 1
      
      $effectID
      [1] NA
      

---

    Code
      load.data("helper_genotypes_inverted_counts.txt", "helper_map.txt", format = define.matrix.format.details(
        has.header = TRUE, markers.as.rows = FALSE))
    Output
      (Loading helper_map.txt) Format: map file with header
      (Loading helper_map.txt) 3 marker(s) with map positions were loaded. Failed to parse 0 line(s).
      (Loading helper_genotypes_inverted_counts.txt) Allele format: reference allele counts (phase will be randomised)
      (Loading helper_genotypes_inverted_counts.txt) 6 genotype(s) of 3 marker(s) were loaded.
      $groupNum
      [1] 1
      
      $mapID
      [1] 1
      
      $effectID
      [1] NA
      

---

    Code
      load.data("helper_genotypes_noheader_nullable_counts.txt", "helper_map.txt",
        format = define.matrix.format.details(has.header = FALSE, markers.as.rows = TRUE,
          cell.style = "C"))
    Output
      (Loading helper_map.txt) Format: map file with header
      (Loading helper_map.txt) 3 marker(s) with map positions were loaded. Failed to parse 0 line(s).
      (Loading helper_genotypes_noheader_nullable_counts.txt) 6 genotype(s) of 3 marker(s) were loaded.
      $groupNum
      [1] 1
      
      $mapID
      [1] 1
      
      $effectID
      [1] NA
      

---

    Code
      load.genotypes("helper_genotypes_slash_nocorner.txt",
        define.matrix.format.details(cell.style = "/", markers.as.rows = TRUE))
    Output
      (Loading helper_genotypes_slash_nocorner.txt) Format: genotype matrix with header row
      (Loading helper_genotypes_slash_nocorner.txt) 6 genotype(s) of 3 marker(s) were loaded.
      [1] 2

