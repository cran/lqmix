xyBuildTV = function(formula, randomTV, data){


  # fixed model frame
  mff = model.frame(formula, data)
  # random model frame
  mfr = model.frame(randomTV, data)
  # random model matrix
  mmr = model.matrix(randomTV, mfr)

  # intercept derived from the formula
  fixInt = attr(terms(mff), "intercept") == 1
  ranInt = attr(terms(mfr), "intercept") == 1

  # if random intercept, remove it from fixed formula
  if (ranInt && fixInt) fixInt = FALSE

  # identify the type of intercept: 0 -> fixed, 1 -> random TC, 999 -> no intercept in the model
  ranInt = if (ranInt && !fixInt) 1 else if (!ranInt && fixInt) 0 else 999

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  termsRan = attr(terms(randomTV),"term.labels")

  # if slopes in common between random and fixed model formula, remove it from fixed formula
  termsFix = setdiff(termsFix, termsRan)

  # any fixed slope? + reformulate the fixed model formula
  fixed = FALSE
  if(length(termsFix)>0){ # if any fixed term
    formula = reformulate(termsFix, response=as.character(formula[2]))
    fixed = TRUE
    mmf = model.matrix(formula, mff)
    if(!fixInt) mmf = mmf[, -1, drop = FALSE]
  }else if(fixInt){ # if no fixed terms -- if only intercept
    formula = reformulate("1", response = as.character(formula[2]))
    mmf = model.matrix(formula, mff)
    fixed = TRUE
  } else formula = mmf = NULL
  namesFix = colnames(mmf)

  # any random slope?
  ranSlope = length(termsRan) > 0
  namesRan = colnames(mmr)


  # prepare covariates
  # *******************
  # random covariates
  mfr = model.frame(randomTV, data)
  x.random = model.matrix(randomTV, mfr)

# fixed covariates
  if(!is.null(formula)){
    mff = model.frame(formula, data)
    x.fixed = model.matrix(formula, mff)
    if(!fixInt) x.fixed = x.fixed[, -1, drop = FALSE]
  }else x.fixed = NULL

  # response
  y.obs = model.response(mff)
  namesY = formula[[2]]

  return(list(y.obs = y.obs, namesY = namesY, x.fixed = x.fixed, x.random = x.random,
              nameFix = namesFix, namesRan = namesRan,
              fixed = fixed, fixInt = fixInt, ranInt = ranInt, ranSlope = ranSlope))
  }
