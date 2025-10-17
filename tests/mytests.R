library(testthat)
library(lqmix)

# build a small dataset ----
set.seed(123)
n_subj = 4
times = 1:3
df = expand.grid(id = 1:n_subj, time = times)
df$meas = rnorm(nrow(df))
df$trt = sample(c(0,1), nrow(df), replace = TRUE)



# test errors in the lqmix function -----
test_that("errors for missing or wrong data arg", {
  expect_error(lqmix(formula = meas ~ trt, randomTC = ~1, group="id", time="time", G=2, data = NULL),
               "No input dataset has been given\\.")
  expect_error(lqmix(formula = meas ~ trt, randomTC = ~1, group="id", time="time", G=2, data = list()),
               "`data' must be a data frame\\.")
})

test_that("error when start=2 but parInit empty", {
  expect_error(lqmix(formula = meas ~ trt, randomTC = ~1, group="id", time="time", G=2, data = df, start=2, parInit = list()),
               "No input parameters have been given with start = 2\\.")
})

test_that("quantile out of range triggers error", {
  expect_error(lqmix(formula = meas ~ trt, randomTC = ~1, group="id", time="time", G=2, data = df, qtl = 0),
               "Quantile level out of range")
  expect_error(lqmix(formula = meas ~ trt, randomTC = ~1, group="id", time="time", G=2, data = df, qtl = 1),
               "Quantile level out of range")
})


test_that("TC branch: errors for missing arguments", {
  # randomTC specified but G missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTC = ~1,
          group = "id", time = "time", data = df),
    "Argument 'G' must be specified to fit a linear quantile mixture with TC random coefficients\\."
  )

  # G provided but randomTC missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTV = ~1, m = 2, G = 2,
          group = "id", time = "time", data = df),
    "Argument 'G' is provided, but no 'randomTC' formula is specified\\. To include TC random ciefficients, define 'randomTC' as well\\."
  )

  # m provided but randomTV missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTC = ~1, G = 2, m = 2,
          group = "id", time = "time", data = df),
    "Argument 'm' is provided, but no 'randomTV' formula is specified\\. To include TV random coefficients, define 'randomTV' as well\\."
  )
})


test_that("TV branch: errors for missing arguments", {
  # randomTV specified but m missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTV = ~1,
          group = "id", time = "time", data = df),
    "Argument 'm' must be specified to fit a linear quantile mixture with TV random coefficients\\."
  )

  # m provided but randomTV missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTC = ~1, G = 2, m = 2,
          group = "id", time = "time", data = df),
    "Argument 'm' is provided, but no 'randomTV' formula is specified\\. To include TV random coefficients, define 'randomTV' as well\\."
  )

  # G provided but randomTC missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTV = ~1, m = 2, G = 2,
          group = "id", time = "time", data = df),
    "Argument 'G' is provided, but no 'randomTC' formula is specified\\. To include TC random ciefficients, define 'randomTC' as well\\."
  )
})



test_that("TCTV branch: errors for missing arguments", {
  # --- TCTV branch errors ---
  # both randomTC and randomTV provided, but G or m missing
  expect_error(
    lqmix(formula = meas ~ trt, randomTC = ~1, randomTV = ~1, m = NULL,
          group = "id", time = "time", data = df),
    "Both 'G' and 'm' must be specified to fit a linear quantile mixture with TC and TV random coefficients\\."
  )
})

test_that("No random-coefficients model", {
  expect_error(
    lqmix(formula = meas ~ trt,
          group = "id", time = "time", data = df),
    "No random-coefficients formulas \\('randomTC' or 'randomTV'\\) have been specified\\. The model corresponds to a standard linear quantile regression without random coefficients\\. Please use the function lqr\\(\\)\\."
  )
})


# ----- test that lqmix retuns proper objects and estimates ----

test_that("TC branch calls lqmixTC() and returns an lqmix object.", {

  res = lqmix(formula = meas ~ trt, randomTC = ~1, group="id", time="time", G=2, data = df, qtl=0.5, se=FALSE)
  expect_s3_class(res, "lqmix")
  expect_named(res, c("betaf","betarTC", "pg", "sigma.e",
                      "npar", "AIC", "BIC", "qtl", "G", "nsbjs", "nobs", "postTC", "miss",
                      "model", "mmf", "mmrTC", "y",
                      "lk","scale","call","formula","randomTC","group","time"), ignore.order = TRUE, ignore.case = TRUE)
  expect_equal(res$formula, meas ~ trt)
  expect_equal(res$randomTC, ~1)
  expect_equal(c(round(res$betarTC,3)), c(-0.490, 0.663))
  expect_null(res$randomTV)
  expect_equal(res$group, "id")
  expect_equal(res$time, "time")
})


test_that("TV branch dispatches to lqmixTV and returns lqmix object", {

  res = lqmix(formula = meas ~ trt, randomTV = ~1, group="id", time="time", m=2, data = df, qtl=0.5, se=FALSE)
  expect_s3_class(res, "lqmix")
  expect_named(res, c("betaf","betarTV", "delta","Gamma", "sigma.e",
                      "npar", "AIC", "BIC", "qtl", "m", "nsbjs", "nobs", "postTV", "miss",
                      "model", "mmf", "mmrTV", "y",
                      "lk","scale","call","formula","randomTV","group","time"), ignore.order = TRUE)
  expect_equal(res$randomTV, ~1)
  expect_null(res$randomTC)
  expect_equal(c(round(res$betarTV,3)), c(-0.490, 0.663))
  expect_equal(res$group, "id")
  expect_equal(res$time, "time")
})


test_that("TCTV branch dispatches to lqmixTCTV and returns lqmix object", {

  res = lqmix(formula = meas ~ trt, randomTC = ~trt, randomTV = ~1, group="id", time="time", G=2, m=2, data = df, qtl=0.5, se=FALSE)
  expect_named(res, c("betarTV","betarTC", "pg","delta","Gamma", "sigma.e",
                      "npar", "AIC", "BIC", "qtl", "m", "G", "nsbjs", "nobs", "postTC", "postTV", "miss",
                      "model", "mmrTC", "mmrTV", "y",
                      "lk","scale","call","formula", "randomTC", "randomTV","group","time"), ignore.order = TRUE)
  expect_null(res$beta)
  expect_null(res$mmf)
  expect_s3_class(res, "lqmix")
  expect_equal(c(round(res$betarTC,3)), c(-0.289, 1.199))
  expect_equal(res$randomTC, ~trt)
  expect_equal(res$randomTV, ~1)
})



# ----- test that search_lqmix retuns proper objects and estimates ----

test_that("search_lqmix: errors for missing Gv and mv", {
  expect_error(search_lqmix(formula = y ~ x, data = df),
               "No values for both Gv and mv are provided\\.")
})

test_that("search_lqmix: errors for unsupported method", {
  expect_error(search_lqmix(formula = y ~ x, data = df, Gv = 1, method = "unknown"),
               "The method specified for selecting the optimal model is not supported\\.")
})


test_that("search_lqmix: TC branch calls lqmixTC and returns search_lqmix object", {

  res = search_lqmix(formula = meas ~ trt, randomTC = ~1, group = "id", time = "time",
                      data = df, Gv = 1:2, method = "bic", se = FALSE, nran = 0)

  expect_s3_class(res, "search_lqmix")
  expect_named(res, c("optimal","allmodels","lkv","aicv","bicv",
                      "qtl","Gv","mv","method","call","formula","randomTC","group","time"),
               ignore.order = TRUE, ignore.case = TRUE)
  expect_equal(res$formula, meas ~ trt)
  expect_equal(round(c(res$optimal$betarTC),3), c(-0.490, 0.663))
  expect_equal(res$randomTC, ~1)
  expect_null(res$randomTV)
})


test_that("search_lqmix: TV branch calls lqmixTV and returns search_lqmix object", {

  res <- search_lqmix(formula = meas ~ trt, randomTV = ~1, group = "id", time = "time",
                      data = df, mv = 1:2, method = "bic", se = FALSE, nran = 0)

  expect_s3_class(res, "search_lqmix")
  expect_named(res, c("optimal","allmodels","lkv","aicv","bicv",
                      "qtl","Gv","mv","method","call","formula","randomTV","group","time"),
               ignore.order = TRUE, ignore.case = TRUE)
  expect_equal(round(c(res$optimal$betarTV),3), c(-0.490, 0.663))
  expect_equal(res$formula, meas ~ trt)
  expect_equal(res$randomTV, ~1)
  expect_null(res$randomTC)
})


test_that("search_lqmix: TCTV branch calls lqmixTCTV and returns search_lqmix object", {

  res <- search_lqmix(formula = meas ~ trt, randomTC = ~trt, randomTV = ~1,
                      group = "id", time = "time", data = df, Gv = 1:2, mv = 1:2,
                      method = "bic", se = FALSE, nran = 0)

  expect_s3_class(res, "search_lqmix")
  expect_named(res, c("optimal","allmodels","lkv","aicv","bicv",
                      "qtl","Gv","mv","method","call","formula","randomTC", "randomTV","group","time"),
               ignore.order = TRUE, ignore.case = TRUE)
  expect_equal(round(c(res$optimal$betarTC),3), c(0.138, 1.292))
  expect_equal(res$randomTC, ~trt)
  expect_equal(res$randomTV, ~1)
})





