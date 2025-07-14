seal_cod.pref <-
function(x) { # x is length in mm
    pref <- (maximum / (1 + exp(-slope1 * (x - midPoint)))) - 
      (maximum / (1 + exp(-slope2 * (x - (midPoint + midPointDistance)))))
    pref <- pref/max(pref)
    return(pref)
  }
