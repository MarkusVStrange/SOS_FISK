corm_herring.pref <-
function(x) { # x is length in mm
    pref <- A * exp(- (x - mu)^2 / (2 * sigma^2))
    pref <- pref/max(pref)
    return(pref)
  }
