lkComputeTV_intermittent = function(delta, Gamma, fith, m, sbj.obs, time.obs, n, T, observed){

  time = rep(1:T, each = n)
  sbj = rep(1:n, T)

  fith.temp = matrix(1, n*T, m)
  fith.temp[observed,] = fith

  # Forward variables
  A = matrix(1, n*T, m)
  A[time == 1,] = (fith.temp[time == 1,] %*% diag(delta))

  for (t in 2:T){
    ## forward at time t-1
    temp.A = matrix(A[time == (t-1),], n, m)
    A[time == t,] = (as.matrix(temp.A %*% Gamma) * fith.temp[time == t,])
  }

  li = rowSums(A[time == T, ])
  lk = sum(log(li))

  return(out = list(li = li, lk = lk, A = A))

}
