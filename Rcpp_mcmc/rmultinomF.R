rmultinomF=
  function(p) {
    return(sum(runif(1) > cumsum(p))+1)
  }