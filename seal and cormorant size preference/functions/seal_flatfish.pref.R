seal_flatfish.pref <-
function(x) { # x is length in mm
    pref <-   maximum / (1 + exp(-slope * (x - midPoint))) # set maximum to the highest remaining observation
    pref <- pref/max(pref)
    return(pref)
  }
