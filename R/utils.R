#' @noRd
expand.path <- function(name) {
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
