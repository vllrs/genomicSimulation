# make.targeted.crosses works

    Code
      g4 <- make.targeted.crosses(first.parents = c(0, "abc"), second.parents = c(3,
        0))
    Warning <simpleWarning>
      Didn't find the name 0
      Didn't find the name abc
      Targeted crossing failed for 2 out of the 2 requested pairings due to one or both genotype indexes being invalid

---

    Code
      g5 <- make.targeted.crosses(first.parents = c(0, 10, 3), second.parents = c(30,
        100, 4))
    Warning <simpleWarning>
      Targeted crossing failed for 2 out of the 3 requested pairings due to one or both genotype indexes being invalid

