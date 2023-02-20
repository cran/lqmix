postComputeTV_intermittent = function(A, li, delta, Gamma, fith, m, sbj.obs, time.obs, observed){

  n = max(sbj.obs)
  T = max(time.obs)
  Ti = table(sbj.obs)

  time = rep(1:T, each = n)
  sbj = rep(1:n, T )
  fith.temp = matrix(1, n*T, m)
  fith.temp[observed,] = fith

  # compute posterior
  B = array(1, dim = c(n*T, m))
  uSingle = matrix(, n*T, m)
  uCouple = array(, c(T, m, m))
  den = matrix(1/li, n*T, m)
  #[observed,]


  wh.pre = time == T-1
  wh.now = (time == T)

  uSingle[wh.now,] = den[wh.now,] * A[wh.now,]
  uCouple[T,, ] = (t(A[wh.pre, ]) %*% (den[wh.now,] * fith.temp[wh.now,] * B[wh.now,])) * Gamma

  for(t in (T-1):1){

    wh.now = (time == t)
    wh.post = (time == t+1)
    wh.pre = time == t-1

    B[wh.now, ] = (fith.temp[wh.post,] * B[wh.post, ]) %*% t(Gamma)

    uSingle[wh.now, ] = den[wh.now,] * A[wh.now, ] * B[wh.now, ]
    if(t>1) uCouple[t,,] = (t(A[wh.pre, ]) %*% (den[wh.now,] * fith.temp[wh.now,] * B[wh.now,]))*Gamma

  }
  uSingle = uSingle[observed,]


  return(out = list(uSingle = uSingle, uCouple = uCouple))
}
