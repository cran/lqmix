xyBuildHOM = function(formula, data){


  mff = model.frame(formula, data)

  # fixed covariates
  x.fixed = model.matrix(formula, mff)

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  namesFix = colnames(x.fixed)


  # response
  y.obs = model.response(mff)
  namesY = formula[[2]]

  return(list(y.obs = y.obs, namesY = namesY, x.fixed = x.fixed,
              nameFix = namesFix))
  }
