postComputeTCTV = function(A, li, delta, Gamma, pg, fithg, m, G, sbj.obs, time.obs, n, T, observed){

  time = rep(1:T, each = n)
  sbj = rep(1:n, T )
  fithg.temp = array(1, c(n*T, m, G))
  fithg.temp[observed,,] = fithg


  # compute posterior
  B = array(1, dim = c(n*T, m, G))
  uSingle = matrix(, n*T, m)
  uCouple = array(, c(T, m, m))

  wh.pre = time == T-1
  wh.now = time == T

  uSingle[wh.now,] = 1/li * apply(array(A[wh.now,,], c(n,m,G)), 2, function(xx){matrix(xx, n,G) %*% pg})

  tmp = 0
  for(g in 1:G)  tmp = tmp + (t(A[wh.pre, ,g] ) %*% (1/li * fithg.temp[wh.now,,g] * B[wh.now,,g]))*pg[g]
  uCouple[T,, ] = tmp * Gamma


  for(t in (T-1):1){
    wh.now = (time == t)
    wh.post = (time == t+1)
    wh.pre = time == t-1

    B[wh.now,,] = apply(fithg.temp[wh.post,,] * B[wh.post,,], 3, function(xx){xx %*% t(Gamma)})

    uSingle[wh.now, ] = 1/li * apply(A[wh.now,,]*B[wh.now,,], 2, function(xx){xx %*% pg})

    if(t>1){
      tmp = 0
      for(g in 1:G) tmp = tmp + (t(A[wh.pre, ,g]) %*% (1/li * fithg.temp[wh.now,,g] * B[wh.now,,g]))*pg[g]
      uCouple[t,, ] = tmp * Gamma
    }

  }

  uSingle = uSingle[observed,]

  # posterior for mixture proportions
  etag = 1/li * t(apply(apply(A[time == T,,], c(1,2), function(xx){xx * pg}), c(1,2), sum))

  # posterior for the longitudinal process

  tmp = aperm(apply(A*B, c(1,2), function(xx){xx * pg}), c(2,3,1))

  den = (array(1/li, c(n*T, m, G)))
  wgt =  tmp * den
  wgt = wgt[observed,,]

  return(out = list(uSingle = uSingle, uCouple = uCouple, etag=etag, wgt=wgt))
}
