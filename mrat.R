mRatio <- function(g, by.strata = TRUE, rpt.size = 8:1) {
  # function takes a matrix of allele frequencies (can be more than one column)
  # rownames represent allele sizes
  calc.mratio<- function(freqs) {
    # extract first column
    freqs <- freqs[, 1]
    
    if(length(freqs) == 1) {
      warning("only one allele")
      return(NA)
    }
    
    # check if rownames are numerics
    #if(!all.is.numeric(names(freqs))) {
    #  warning("allele names are non-numeric")
    #  return(NA)
    #}
    
    if(all(freqs == 0)) { 
      warning("all frequencies are 0")
      NA
    } else {
      # sort alleles in numerical order
      freqs <- freqs[order(as.numeric(names(freqs)))]
      # convert names to numbers to get sizes
      sizes <- as.numeric(names(freqs))
      # find repeat sizes
      size.diff <- diff(sizes)
      rpt.found <- FALSE
      for(r in sort(rpt.size, decreasing = TRUE)) {
        if(all(size.diff %% r == 0)) {
          rpt.found <- TRUE
          break
        }
      }
      if(!rpt.found) {
        warning("valid repeat length not found")
        return(NA)
      }
      # find smallest and largest alleles that are present
      smallest <- min(sizes[freqs > 0])
      largest <- max(sizes[freqs > 0])
      # compute number of alleles between smallest and largest
      n <- (largest - smallest) / r
      # calculate metric
      sum(freqs > 0) / (n + 1)
    }
  }
  
  if(nStrata(g) == 1 & by.strata) {
    by.strata <- FALSE
    g <- g[, , strataNames(g)]
  }
  if(by.strata) {
    freqs <- alleleFreqs(g, by.strata = TRUE)
    do.call(rbind, lapply(freqs, function(loc) apply(loc, 3, calc.mratio)))
  } else {
    freqs <- alleleFreqs(g, by.strata = FALSE)
    sapply(freqs, calc.mratio)
  }
}