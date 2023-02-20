postComputeTC = function(lig, fig, pg, G, Ti, order.time){

  n = nrow(fig)
  wig = lig / matrix(fig %*% pg, n, G)
  Wig = apply(wig, 2, function(x){rep(x, Ti)}) # estend Ti times each
  Wig = Wig[order.time,] # reorder wtr to time

  return(out = list(Wig=Wig, wig=wig))
}
