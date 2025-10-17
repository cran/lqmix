xyBuildTCTV = function(formula, randomTC, randomTV, data){

  # fixed model frame
  mff = model.frame(formula, data)
  # random model frame TC
  mfrTC = model.frame(randomTC, data)
  # random model matrix TC
  mmrTC = model.matrix(randomTC, mfrTC)
  # random model frame TV
  mfrTV = model.frame(randomTV, data)
  # random model matrix TV
  mmrTV = model.matrix(randomTV, mfrTV)

  # intercept derived from the formula
  fixInt = attr(terms(mff), "intercept") == 1
  ranIntTC = attr(terms(mfrTC), "intercept") == 1
  ranIntTV = attr(terms(mfrTV), "intercept") == 1


  # terms specified as both TC and TV -- error
  randomTermsTC = c(ranIntTC*1, attr(terms(randomTC),"term.labels"))
  randomTermsTV = c(ranIntTV*1, attr(terms(randomTV),"term.labels"))

  wh = intersect(randomTermsTC, randomTermsTV)
  if(length(wh)>0) stop(paste("One or more terms are specified as both TC and TV random coefficients:", paste(wh, collapse=", ")))

  # if random intercept, remove it from fixed formula
  if((ranIntTC | ranIntTV) & fixInt) fixInt = FALSE

  # identify the type of intercept: 0 -> fixed, 1 -> random TC, 2 -> TV, 999 -> no intercept in the model
  ranInt = 999  # Default value
  if (ranIntTC == 1 & fixInt == 0) ranInt = 1 else if (ranIntTV == 1 & fixInt == 0) ranInt = 2 else if (ranIntTC == 0 & ranIntTV == 0 & fixInt == 1) ranInt = 0

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  termsRanTC = attr(terms(randomTC),"term.labels")
  termsRanTV = attr(terms(randomTV),"term.labels")

  # if slopes in common between random and fixed model formula, remove it from fixed formula
  termsFix = setdiff(termsFix, termsRanTC)
  termsFix = setdiff(termsFix, termsRanTV)

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
  ranSlopeTC = length(termsRanTC) > 0
  ranSlopeTV = length(termsRanTV) > 0

  namesRanTC = colnames(mmrTC)
  namesRanTV = colnames(mmrTV)


  # prepare covariates
  # *******************
  # random covariates
  mfrTC = model.frame(randomTC, data)
  x.randomTC = model.matrix(randomTC, mfrTC)

  mfrTV = model.frame(randomTV, data)
  x.randomTV = model.matrix(randomTV, mfrTV)

  # fixed covariates
  if(!is.null(formula)){
    mff = model.frame(formula, data)
    x.fixed = model.matrix(formula, mff)
    if(!fixInt) x.fixed = x.fixed[, -1, drop = FALSE]
  }else x.fixed = NULL

  # response variable
  y.obs = model.response(mff)
  namesY = formula[[2]]

  return(list(y.obs = y.obs, namesY = namesY, x.fixed = x.fixed, x.randomTC = x.randomTC, x.randomTV = x.randomTV,
              nameFix = namesFix, namesRanTC = namesRanTC, namesRanTV = namesRanTV,
              fixed = fixed, fixInt = fixInt, ranInt = ranInt, ranSlopeTC = ranSlopeTC, ranSlopeTV = ranSlopeTV))
  }
