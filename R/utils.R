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
chr.ix.to.genomicSimulation.numbering <- function(chr.names) {
  chr.uniqnames <- unique(chr.names)
  data.frame("chr.names"=chr.uniqnames,
             "genomicSimulation.chr.number"=rank(strtoi(chr.uniqnames,36))-1)
}


#' @noRd
fallback.param.names <- function(emptydefault, params.in.order) {
  for (param in params.in.order) {
    if (param != emptydefault) {
      return(param)
    }
  }
  return(emptydefault)
}

#' @noRd
convert.to.integer <- function(values, param.name, allow.na=F) {
  if (!is.integer(values)) {
    ivalues <- as.integer(values)
    if (!isTRUE(all(ivalues==values, na.rm=allow.na))) { stop(param.name, " must be integers.") }
    values <- ivalues
  }
  return(values)
}
