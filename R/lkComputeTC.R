lkComputeTC = function(pg, fig){

  lig = fig %*% diag(pg)
  li = log(rowSums(lig))
  lk = sum(li)

  return(out = list(li = li, lig = lig, lk = lk))
}
