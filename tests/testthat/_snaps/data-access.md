# see.group.data works with multiple groups

    Code
      ntmp <- see.group.data(c(g, 50), "Names")
    Output
      NOTE! Group 50 (entry 2 in the `group` vector) is an empty group

---

    Code
      ntmp <- see.group.data(c(g3, g, g2), "D")
    Output
      NOTE! Group 1 (entry 2 in the `group` vector) is an empty group

# See.group.gene.data works

    Code
      tmp <- see.group.gene.data(g2, unknown.allele = "\n")
    Output
      NOTE! Defaulting to printing '\0' alleles as a dash '-'

