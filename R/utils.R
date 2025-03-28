#' @noRd
my.expand.path <- function(name) {
  if (name != "" && is.character(name)) { #nonempty name
    
    if (grepl("~",name,fixed=TRUE)) {
      if (requireNamespace("fs",quietly=TRUE)) { # more portable
        return(fs::path_expand(name)) 
      } else { # slightly less portable
        return(base::path.expand(name))
      }
      
    } else { 
      return(name)
    }
    
  } else { #empty name
    return(NULL)
  }
}

#' @noRd
chr.names.to.genomicSimulation.numbering <- function(chr.names) {
  chr.uniqnames <- unique(chr.names)
  data.frame("chr.names"=chr.uniqnames,
             "genomicSimulation.chr.number"=rank(strtoi(chr.uniqnames,36))-1)
}