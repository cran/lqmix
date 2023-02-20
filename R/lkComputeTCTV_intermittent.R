lkComputeTCTV = function(delta, Gamma, pg, fithg, m, G, sbj.obs, time.obs, n, T, observed){

  time = rep(1:T, each = n)
  sbj = rep(1:n, T)

  fithg.temp = array(1, c(n*T, m, G))
  fithg.temp[observed,,] = fithg


  # FORWARD VARIABLES
  A = array(1, c(n*T, m, G))

  A[time == 1,,] = array(apply(fithg.temp[time == 1,,], 3, function(xx){xx %*% diag(delta)}), c(n,m,G))

  for (t in 2:T){
    ## forward at time t-1
    temp.A = array(A[time == (t-1),,], c(n, m, G))
    A[time == t,,] = array(apply(temp.A, 3, function(xx){xx %*% Gamma}), c(n, m, G)) * fithg.temp[time == t,,]
  }

  li = rowSums(apply(A[time == T, ,], c(1,2), function(xx){xx %*% pg}))
  lk = sum(log(li))

  return(out = list(li = li, lk = lk, A = A))

}
