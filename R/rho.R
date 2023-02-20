rho = function(x, qtl){
  x * (qtl - ifelse(x <0, 1, 0))
}
