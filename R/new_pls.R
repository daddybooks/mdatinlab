#' Partial Least Squares regression
#'
#' @description
#' \code{pls} is used to calibrate, validate and use of partial least squares (PLS)
#' regression model.
#'
#' @param x
#' matrix with predictors.
#' @param y
#' matrix with responses.
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, center or not predictors and response values.
#' @param scale
#' logical, scale (standardize) or not predictors and response values.
#' @param cv
#' cross-validation settings (see details).
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test
#' matrix with predictors for test set.
#' @param y.test
#' matrix with responses for test set.
#' @param method
#' algorithm for computing PLS model (only 'simpls' is supported so far)
#' @param lim.type
#' which method to use for calculation of critical limits for residual distances (see details)
#' @param alpha
#' significance level for extreme limits for T2 and Q disances.
#' @param gamma
#' significance level for outlier limits for T2 and Q distances.
#' @param info
#' short text with information about the model.
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components (\code{'min'} for
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#'
#' @return
#' Returns an object of \code{pls} class with following fields:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{xcenter }{vector with values used to center the predictors (x).}
#' \item{ycenter }{vector with values used to center the responses (y).}
#' \item{xscale }{vector with values used to scale the predictors (x).}
#' \item{yscale }{vector with values used to scale the responses (y).}
#' \item{xloadings }{matrix with loading values for x decomposition.}
#' \item{yloadings }{matrix with loading values for y decomposition.}
#' \item{xeigenvals }{vector with eigenvalues of components (variance of x-scores).}
#' \item{yeigenvals }{vector with eigenvalues of components (variance of y-scores).}
#' \item{weights }{matrix with PLS weights.}
#' \item{coeffs }{object of class \code{\link{regcoeffs}} with regression coefficients calculated
#' for each component.}
#' \item{info }{information about the model, provided by user when build the model.}
#' \item{cv }{information cross-validation method used (if any).}
#' \item{res }{a list with result objects (e.g. calibration, cv, etc.)}
#'
#' @details
#' So far only SIMPLS method [1] is available. Implementation works both with one
#' and multiple response variables.
#'
#' Like in \code{\link{pca}}, \code{pls} uses number of components (\code{ncomp}) as a minimum of
#' number of objects - 1, number of x variables and the default or provided value. Regression
#' coefficients, predictions and other results are calculated for each set of components from 1
#' to \code{ncomp}: 1, 1:2, 1:3, etc. The optimal number of components, (\code{ncomp.selected}),
#' is found using first local minumum, but can be also forced to user defined value using function
#' (\code{\link{selectCompNum.pls}}). The selected optimal number of components is used for all
#' default operations - predictions, plots, etc.
#'
#' Cross-validation settings, \code{cv}, can be a number or a list. If \code{cv} is a number, it
#' will be used as a number of segments for random cross-validation (if \code{cv = 1}, full
#' cross-validation will be preformed). If it is a list, the following syntax can be used:
#' \code{cv = list("rand", nseg, nrep)} for random repeated cross-validation with \code{nseg}
#' segments and \code{nrep} repetitions or \code{cv = list("ven", nseg)} for systematic splits
#' to \code{nseg} segments ('venetian blinds').
#'
#' Calculation of confidence intervals and p-values for regression coefficients can by done
#' based on Jack-Knifing resampling. This is done automatically if cross-validation is used.
#' However it is recommended to use at least 10 segments for stable JK result. See help for
#' \code{\link{regcoeffs}} objects for more details.
#'
#' @references
#' 1. S. de Jong, Chemometrics and Intelligent Laboratory Systems 18 (1993) 251-263.
#' 2. Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), 35-48.
#' 3. Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), 103-112.
#'
#' @seealso
#' Main methods for \code{pls} objects:
#' \tabular{ll}{
#'  \code{print} \tab prints information about a \code{pls} object.\cr
#'  \code{\link{summary.pls}} \tab shows performance statistics for the model.\cr
#'  \code{\link{plot.pls}} \tab shows plot overview of the model.\cr
#'  \code{\link{pls.simpls}} \tab implementation of SIMPLS algorithm.\cr
#'  \code{\link{predict.pls}} \tab applies PLS model to a new data.\cr
#'  \code{\link{selectCompNum.pls}} \tab set number of optimal components in the model.\cr
#'  \code{\link{setDistanceLimits.pls}} \tab allows to change parameters for critical limits.\cr
#'  \code{\link{categorize.pls}} \tab categorize data rows similar to
#'    \code{\link{categorize.pca}}.\cr
#'  \code{\link{selratio}} \tab computes matrix with selectivity ratio values.\cr
#'  \code{\link{vipscores}} \tab computes matrix with VIP scores values.\cr
#' }
#'
#' Plotting methods for \code{pls} objects:
#' \tabular{ll}{
#'  \code{\link{plotXScores.pls}} \tab shows scores plot for x decomposition.\cr
#'  \code{\link{plotXYScores.pls}} \tab shows scores plot for x and y decomposition.\cr
#'  \code{\link{plotXLoadings.pls}} \tab shows loadings plot for x decomposition.\cr
#'  \code{\link{plotXYLoadings.pls}} \tab shows loadings plot for x and y decomposition.\cr
#'  \code{\link{plotXVariance.pls}} \tab shows explained variance plot for x decomposition.\cr
#'  \code{\link{plotYVariance.pls}} \tab shows explained variance plot for y decomposition.\cr
#'  \code{\link{plotXCumVariance.pls}} \tab shows cumulative explained variance plot for y
#'  decomposition.\cr
#'  \code{\link{plotYCumVariance.pls}} \tab shows cumulative explained variance plot for y
#'  decomposition.\cr
#'  \code{\link{plotXResiduals.pls}} \tab shows distance/residuals plot for x decomposition.\cr
#'  \code{\link{plotXYResiduals.pls}} \tab shows joint distance plot for x and y decomposition.\cr
#'  \code{\link{plotWeights.pls}} \tab shows plot with weights.\cr
#'  \code{\link{plotSelectivityRatio.pls}} \tab shows plot with selectivity ratio values.\cr
#'  \code{\link{plotVIPScores.pls}} \tab shows plot with VIP scores values.\cr
#' }
#'
#' Methods inherited from \code{regmodel} object (parent class for \code{pls}):
#' \tabular{ll}{
#'  \code{\link{plotPredictions.regmodel}} \tab shows predicted vs. measured plot.\cr
#'  \code{\link{plotRMSE.regmodel}} \tab shows RMSE plot.\cr
#'  \code{\link{plotYResiduals.regmodel}} \tab shows residuals plot for y values.\cr
#'  \code{\link{getRegcoeffs.regmodel}} \tab returns matrix with regression coefficients.\cr
#' }
#'
#' Most of the methods for plotting data (except loadings and regression coefficients) are also
#' available for PLS results (\code{\link{plsres}}) objects. There is also a randomization test
#' for PLS-regression (\code{\link{randtest}}) and implementation of interval PLS algorithm
#' for variable selection (\code{\link{ipls}})
#'

#'

#'
#' @export
pls <- function(x, y, ncomp = min(nrow(x) - 1, ncol(x), 20), center = TRUE, scale = FALSE,
  cv = NULL, exclcols = NULL, exclrows = NULL, x.test = NULL, y.test = NULL, method = "simpls",
  info = "", ncomp.selcrit = "min", lim.type = "ddmoments", alpha = 0.05, gamma = 0.01) {

  # if y is a vector, convert it to matrix
  if (is.null(dim(y))) {
    dim(y) <- c(length(y), 1)
  }

  # check calibration data and process excluded rows and columns
  x <- prepCalData(x, exclrows = exclrows, exclcols = exclcols, min.nrows = 2, min.ncols = 1)
  y <- prepCalData(y, exclrows = exclrows, exclcols = NULL, min.nrows = 2, min.ncols = 1)

  # build a model and apply to calibration set
  model <- pls.cal(x, y, ncomp, center = center, scale = scale, method = method, cv = cv)
  model$info <- info
  model$call <- match.call()

  # get calibration results
  model$res <- list()
  model$res[["cal"]] <- predict.pls(model, x, y)
  model$res[["cal"]]$info <- "calibration results"
  model$calres <- model$res[["cal"]]

  # compute critical limit parameters
  model$limParams <- list(
    "Q" = ldecomp.getLimParams(model$res[["cal"]]$xdecomp$Q),
    "T2" = ldecomp.getLimParams(model$res[["cal"]]$xdecomp$T2),
    "Z" = ldecomp.getLimParams(model$res[["cal"]]$ydecomp$Q)
  )

  # do cross-validation if needed
  if (!is.null(cv)) {
    cvres <- crossval.regmodel(model, x, y, cv, cal.fun = pls.cal)
    model$res[["cv"]] <- plsres(cvres$y.pred, cvres$y.ref, ncomp.selected = model$ncomp)
    model$res[["cv"]]$info <- "cross-validation results"
    model$cvres <- model$res[["cv"]]
    model$coeffs <- regcoeffs(model$coeffs$values, cvres$jk.coeffs)
  }

  # do test set validation if provided
  if (!is.null(x.test) && !is.null(y.test)) {
    model$res[["test"]] <- predict.pls(model, x.test, y.test)
    model$res[["test"]]$info <- "test set validation results"
    model$testres <- model$res[["test"]]
  }

  # select optimal number of components
  model$cv <- cv
  model$ncomp.selcrit <- ncomp.selcrit
  model <- selectCompNum(model, selcrit = ncomp.selcrit)

  # set distance limits
  model <- setDistanceLimits(model, lim.type = lim.type, alpha = alpha, gamma = gamma)

  return(model)
}

#' Select optimal number of components for PLS model
#'
#' @description
#' Allows user to select optimal number of components for PLS model
#'
#' @param obj
#' PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to select
#' @param selcrit
#' criterion for selecting optimal number of components (\code{'min'} for
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#' @param ...
#' other parameters if any
#'
#' @return
#' the same model with selected number of components
#'
#' @details
#'
#' The method sets \code{ncomp.selected} parameter for the model and return it back. The parameter
#' points out to the optimal number of components in the model. You can either specify it manually,
#' as argument \code{ncomp}, or use one of the algorithms for automatic selection.
#'
#' Automatic selection by default based on cross-validation statistics. If no cross-validation
#' results are found in the model, the method will use test set validation results. If they are
#' not available as well, the model will use calibration results and give a warning as in this case
#' the selected number of components will lead to overfitted model.
#'
#' There are two algorithms for automatic selection you can chose between: either first local
#' minimum of RMSE (`selcrit="min"`) or Wold's rule (`selcrit="wold"`).
#'
#' The first local minimum criterion finds at which component, A, error of prediction starts
#' raising and selects (A - 1) as the optimal number. The Wold's criterion finds which component A
#' does not make error smaller at least by 5% comparing to the previous value and selects (A - 1)
#' as the optimal number.
#'
#' If model is PLS2 model (has several response variables) the method computes optimal number of
#' components for each response and returns the smallest value. For example, if for the first
#' response 2 components give the smallest error and for the second response this number is 3,
#' A = 2 will be selected as a final result.
#'
#' It is not recommended to use automatic selection for real applications, always investigate
#' your model (via RMSE, Y-variance plot, regression coefficients) to make correct decision.
#'
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
selectCompNum.pls <- function(obj, ncomp = NULL, selcrit = obj$ncomp.selcrit, ...) {

  # returns index based on Wold's R criterion
  # finds first number of components which does not make error smaller by 5%
  # comparing to the previous value
  fWold <- function(press, ncomp) {
    r <- press[, 2:ncomp, drop = FALSE] / press[, 1:(ncomp - 1), drop = FALSE]
    return(which(r > 0.95, arr.ind = TRUE))
  }

  # returns index based on first local minimum
  fMin <- function(press, ncomp) {
    r <- press[, 2:ncomp, drop = FALSE] - press[, 1:(ncomp - 1), drop = FALSE]
    return(which(r > 0, arr.ind = TRUE))
  }

  # returns number of components based on PRESS and index values
  f <- function(res, indFun) {
    press <- res$rmse^2 * nrow(res$y.pred)
    ncomp <- ncol(press)

    # if number of components is 2 - for every row find which component
    # gives smallest error and take the smallest number
    if (ncomp < 3) return(min(apply(press, 1, which.min)))

    # otherwise use dedicated function
    ind <- indFun(press, ncomp)

    return((if (length(ind) == 0) ncomp else min(ind[, 2])))
  }

  # if only one component in the model do nothing
  if (obj$ncomp == 1) return(obj)

  # if user provided ncomp use it
  if (!is.null(ncomp)) selcrit <- ""

  # get name of first result available in the sequence
  name <- intersect(c("cv", "test", "cal"), names(obj$res))[1]
  res <- obj$res[[name]]

  if (name == "cal") {
    warning("No validation results were found.")
  }

  # estimate number of optimal components
  if (is.null(selcrit)) selcrit <- ""
  ncomp <- switch(selcrit,
    "wold" = f(res, fWold),
    "min" = f(res, fMin),
    ncomp
  )

  # if NULL - somthing went wrong
  if (is.null(ncomp)) {
    stop("Can not estimate correct number of PLS components.")
  }

  # if not, check that value is meaningful
  if (ncomp > obj$ncomp || ncomp < 0) {
    stop("Wrong number of selected components.")
  }

  # correct number of model and calibration results
  obj$ncomp.selected <- ncomp
  obj$res[["cal"]]$ncomp.selected <- ncomp
  obj$calres <- obj$res[["cal"]]

  # correct number of components for cross-validation results
  if (!is.null(obj$res[["cv"]])) {
    obj$res[["cv"]]$ncomp.selected <- ncomp
    obj$cvres <- obj$res[["cv"]]
  }

  # correct number of components for test set results
  if (!is.null(obj$res[["test"]])) {
    obj$res[["test"]]$ncomp.selected <- ncomp
    obj$testres <- obj$res[["test"]]
  }

  obj$call <- match.call()
  return(obj)
}

#' Compute and set statistical limits for residual distances.
#'
#' @description
#' Computes statisticsl limits for orthogonal and score distances (x-decomposition) and
#' orthogonal distance (y-decomposition) based on calibration set and assign the calculated
#' values as model properties.
#'
#' @param obj
#' object with PLS model
#' @param lim.type
#' type of limits ("jm", "chisq", "ddmoments", "ddrobust")
#' @param alpha
#' significance level for detection of extreme objects
#' @param gamma
#' significance level for detection of outliers (for data driven approach)
#' @param ...
#' other arguments
#'
#' @details
#'
#' The limits can be accessed as fields of model objects: \code{$Qlim}, \code{$T2lim}, and
#' \code{$Zlim}. Each is a matrix with four rows and \code{ncomp} columns. In case of limits
#' for x-decomposition, first row contains critical limits for extremes, second row - for outliers,
#' third row contains mean value for corresponding distances (or its robust estimate in case of
#' \code{lim.type = "ddrobust"}) and last row contains the degrees of freedom.
#'
#' @return
#' Object models with the three fields updated.
#'
#' @export
setDistanceLimits.pls <- function(obj, lim.type = obj$lim.type, alpha = obj$alpha,
  gamma = obj$gamma, ...) {

  obj$T2lim <- ldecomp.getT2Limits(lim.type, alpha, gamma, obj$limParams)
  obj$Qlim <- ldecomp.getQLimits(lim.type, alpha, gamma, obj$limParams,
    obj$res[["cal"]]$xdecomp$residuals, obj$xeigenvals)
  obj$Zlim <- pls.getZLimits(lim.type, alpha, gamma, obj$limParams)

  obj$alpha <- alpha
  obj$gamma <- gamma
  obj$lim.type <- lim.type

  attr(obj$res[["cal"]]$xdecomp$Q, "u0") <- obj$Qlim[3, ]
  attr(obj$res[["cal"]]$xdecomp$Q, "Nu") <- obj$Qlim[4, ]

  attr(obj$res[["cal"]]$xdecomp$T2, "u0") <- obj$T2lim[3, ]
  attr(obj$res[["cal"]]$xdecomp$T2, "Nu") <- obj$T2lim[4, ]

  attr(obj$res[["cal"]]$ydecomp$Q, "u0") <- obj$Zlim[3, ]
  attr(obj$res[["cal"]]$ydecomp$Q, "Nu") <- obj$Zlim[4, ]

  obj$calres <- obj$res[["cal"]]

  if (!is.null(obj$res$test)) {
    attr(obj$res[["test"]]$xdecomp$Q, "u0") <- obj$Qlim[3, ]
    attr(obj$res[["test"]]$xdecomp$Q, "Nu") <- obj$Qlim[4, ]

    attr(obj$res[["test"]]$xdecomp$T2, "u0") <- obj$T2lim[3, ]
    attr(obj$res[["test"]]$xdecomp$T2, "Nu") <- obj$T2lim[4, ]

    attr(obj$res[["test"]]$ydecomp$Q, "u0") <- obj$Zlim[3, ]
    attr(obj$res[["test"]]$ydecomp$Q, "Nu") <- obj$Zlim[4, ]

    obj$testres <- obj$res[["test"]]
  }

  return(obj)
}

#' PLS predictions
#'
#' @description
#' Applies PLS model to a new data set
#'
#' @param object
#' a PLS model (object of class \code{pls})
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with reference y values (responses)
#' @param cv
#' logical, shall predictions be made for cross-validation procedure or not
#' @param ...
#' other arguments
#'
#' @return
#' PLS results (an object of class \code{\link{plsres}})
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
predict.pls <- function(object, x, y = NULL, cv = FALSE, ...) {

  # get names
  prednames <- rownames(object$xloadings)
  respnames <- rownames(object$yloadings)
  objnames <- rownames(x)

  # preprocess x and calculate scores, total and full variance
  x.attrs <- mda.getattr(x)
  y.attrs <- mda.getattr(y)

  # set names for y-axis (rows if it is empty)
  if (is.null(x.attrs$yaxis.name)) {
    x.attrs$yaxis.name <- "Objects"
  }

  # check datasets and convert to matrix if needed
  x <- prepCalData(x, min.nrows = 1, min.ncols = nrow(object$xloadings) - length(x.attrs$exclcols))

  # get dimensions
  nresp <- dim(object$coeffs$values)[3]
  ncomp <- dim(object$coeffs$values)[2]

  # check dimensions of predictors
  if (ncol(x) != dim(object$coeffs$values)[1]) {
    stop("Wrong number of columns in matrix with predictors (x).")
  }

  # autoscale x
  x <- prep.autoscale(x, center = object$xcenter, scale = object$xscale)

  # compute x scores and residuals
  xscores <- x %*% (object$weights %*% solve(crossprod(object$xloadings, object$weights)))
  xresiduals <- x - tcrossprod(xscores, object$xloadings)

  # set attributes
  xscores <- mda.setattr(xscores, x.attrs, "row")
  xresiduals <- mda.setattr(xresiduals, x.attrs)
  attr(xscores, "name") <- "X-scores"
  attr(xscores, "xaxis.name") <- "Components"
  attr(xresiduals, "name") <- "Residuals"

  # set names
  rownames(xscores) <- rownames(xresiduals) <- objnames
  colnames(xscores) <- colnames(object$xloadings)
  colnames(xresiduals) <- prednames

  # make predictions
  yp <- apply(object$coeffs$values, 3, function(x, y)(y %*% x), x)
  dim(yp) <- c(nrow(x), ncomp, dim(object$coeffs$values)[3])

  # if reference values are provided calculate and set up names for ydecomp
  y.ref <- NULL
  ydecomp <- NULL
  if (!is.null(y)) {

    if (is.null(dim(y))) dim(y) <- c(length(y), 1)

    if (nrow(x) != nrow(y)) {
      stop("Matrices with predictors (x) and response (y) should have the same number of rows.")
    }

    if (ncol(y) != nresp) {
      stop("Wrong number of columns in matrix with response values (y).")
    }

    # keep the original y values as reference
    y.ref <- y

    # autoscale y-values
    y <- prep.autoscale(y, center = object$ycenter, scale = object$yscale)

    # compute and orthogonalize y-scores
    yscores <- as.matrix(y) %*% object$yloadings
    for (a in seq_len(ncomp)) {
      for (n in 1:2) {
        for (j in seq_len(a - 1)) {
          yscores[, a] <- yscores[, a] -
            tcrossprod(xscores[, j], yscores[, a]) %*% xscores[, j]
        }
      }
    }

    # compute y-residuals
    yresiduals <- y - yp[, ncomp, ]

    # set names
    rownames(yscores) <- rownames(yresiduals) <- objnames
    colnames(yscores) <- colnames(object$yloadings)
    colnames(yresiduals) <- respnames

    # set attributes
    yscores <- mda.setattr(yscores, x.attrs, "row")
    yresiduals <- mda.setattr(yresiduals, y.attrs)
    attr(yscores, "exclrows") <- attr(yresiduals, "exclrows") <- x.attrs$exclrows
    attr(yscores, "name") <- "Y-scores"
    attr(yscores, "xaxis.name") <- "Components"
    attr(yresiduals, "name") <- "Residuals"

    # create ydecomp object (we use xscores as residuals for different components are computed
    # as xscores %*% t(yloadings)), but then we assign correct residuals
    ydecomp <- ldecomp(scores = xscores, loadings = object$yloadings, residuals = yresiduals,
      eigenvals = object$yeigenvals, ncomp.selected = object$ncomp.selected)
    ydecomp$scores <- yscores
    attr(ydecomp$Q, "u0") <- object$Zlim[3, ]
    attr(ydecomp$Q, "Nu") <- object$Zlim[4, ]
  }

  # unscale predicted y values
  yp <- if (is.numeric(object$yscale)) sweep(yp, 3, object$yscale, "*") else yp

  # uncenter predicted y values
  yp <- if (is.numeric(object$ycenter)) sweep(yp, 3, object$ycenter, "+") else yp

  # if predictions for cross-validation - return
  if (cv) {
    return(list(y.pred = yp))
  }

  # set up all attributes and names
  yp <- mda.setattr(yp, x.attrs, "row")
  attr(yp, "exclrows") <- x.attrs$exclrows
  attr(yp, "name") <- "Response values, predicted"
  dimnames(yp) <- c(list(rownames(x)), dimnames(object$coeffs$values)[2:3])

  # create xdecomp object
  xdecomp <- ldecomp(scores = xscores, residuals = xresiduals, loadings = object$xloadings,
    eigenvals = object$xeigenvals, ncomp.selected = object$ncomp.selected)

  # add u0 and Nu parameters as arguments, so the residuals can be normalized
  attr(xdecomp$Q, "u0") <- object$Qlim[3, ]
  attr(xdecomp$Q, "Nu") <- object$Qlim[4, ]

  attr(xdecomp$T2, "u0") <- object$T2lim[3, ]
  attr(xdecomp$T2, "Nu") <- object$T2lim[4, ]

  return(
    plsres(yp, y.ref = y.ref, ncomp.selected = object$ncomp.selected,
      xdecomp = xdecomp, ydecomp = ydecomp)
  )
}

#' Categorize data rows based on PLS results and critical limits for total distance.
#'
#' @description
#' The method uses full distance for decomposition of X-data and squared Y-residuals of PLS results
#' from \code{res} with critical limits computed for the PLS model and categorizes the
#' corresponding objects as "regular", "extreme" or "outlier".
#'
#' @param obj
#' object with PCA model
#' @param res
#' object with PCA results
#' @param ncomp
#' number of components to use for the categorization
#' @param ...
#' other parameters
#'
#' @details
#' The method does not categorize hidden values if any. It is based on the approach described in
#' [1] and works only if data driven approach is used for computing critical limits.
#'
#' @return
#' vector (factor) with results of categorization.
#'
#' @references
#' 1. Rodionova O. Ye., Pomerantsev A. L. Detection of Outliers in Projection-Based Modeling.
#' Analytical Chemistry (2020, in publish). doi: 10.1021/acs.analchem.9b04611
#'
#' @export
categorize.pls <- function(obj, res = obj$res$cal, ncomp = obj$ncomp.selected, ...) {

  create_categories <- function(extremes_ind, outliers_ind) {
    categories <- rep(1, length(extremes_ind))
    categories[extremes_ind] <- 2
    categories[outliers_ind] <- 3
    return(factor(categories, levels = 1:3, labels = c("regular", "extreme", "outlier")))
  }

  # if not data driven - quit
  if (!(obj$lim.type %in% c("ddmoments", "ddrobust"))) {
    stop("categorize.pls() works only with data driven limit types ('ddoments' or 'ddrobust').")
  }

  # get distance values for selected number of components
  h <- res$xdecomp$T2[, ncomp]
  q <- res$xdecomp$Q[, ncomp]
  z <- res$ydecomp$Q[, ncomp]

  # remove excluded values if any
  rows_excluded <- attr(res$xdecomp$Q, "exclrows")
  if (length(rows_excluded) > 0) {
    h <- h[-rows_excluded]
    q <- q[-rows_excluded]
    z <- z[-rows_excluded]
  }

  # get DoF
  Nh <- obj$T2lim[4, ncomp]
  Nq <- obj$Qlim[4, ncomp]
  Nz <- obj$Zlim[4, ncomp]
  Nf <- Nq + Nh

  # get scale factor
  h0 <- obj$T2lim[3, ncomp]
  q0 <- obj$Qlim[3, ncomp]
  z0 <- obj$Zlim[3, ncomp]

  # process degrees of freedom for (Z)
  Nz <- round(Nz)
  Nz[Nz < 1] <- 1
  Nz[Nz > 250] <- 250

  # process degrees of freedom for (F)
  Nf <- round(Nf)
  Nf[Nf < 1] <- 1
  Nf[Nf > 250] <- 250

  # compute total distance and DoF for it
  g <- Nh * h / h0 + Nq * q / q0 + Nz * z / z0
  Ng <- Nh + Nq + Nz
  nobj <- nrow(obj$res$cal$xdecomp$scores)

  # compute limits for total distance
  ext_lim <- qchisq(1 - obj$alpha, Ng)
  out_lim <- qchisq((1 - obj$gamma) ^ (1 / nobj), Ng)

  outliers_ind <- g > out_lim
  extremes_ind <- g > ext_lim & g < out_lim

  return(create_categories(extremes_ind, outliers_ind))
}

#' Summary method for PLS model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to count.
#' @param ny
#' which y variables to show the summary for (can be a vector)
#' @param ...
#' other arguments
#'
#' @export
summary.pls <- function(object, ncomp = object$ncomp.selected,
  ny = seq_len(nrow(object$yloadings)), ...) {

  if (length(ncomp) != 1 || ncomp < 0 || ncomp > object$ncomp) {
    stop("Wrong value for the 'ncomp' parameter.")
  }

  cat("\nPLS model (class pls) summary\n")
  cat("-------------------------------\n")
  fprintf("Info: %s\n", object$info)
  fprintf("Number of selected components: %d\n", ncomp)
  fprintf("Cross-validation: %s\n", crossval.str(object$cv))

  cat("\n")
  for (y in ny) {
    fprintf("Response variable: %s\n", rownames(object$yloadings)[y])
    out <- do.call(rbind, lapply(object$res, as.matrix, ncomp = ncomp, ny = y))
    rownames(out) <- capitalize(names(object$res))

    if (!any(is.na(out[, 1:4]))) out[, 1:4] <- round(out[, 1:4], 3)
    out[, 5] <- round(out[, 5], 3)
    out[, 6] <- mdaplot.formatValues(out[, 6], round.only = T)
    out[, 7] <- round(out[, 7], 3)
    out[, 8] <- round(out[, 8], 4)
    out[, 9] <- round(out[, 9], 2)

    print(out[, -c(1, 3), drop = FALSE])
    cat("\n")
  }
}

#' Print method for PLS model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ...
#' other arguments
#'
#' @export
print.pls <- function(x, ...) {
  cat("\nPLS model (class pls)\n")
  cat("\nCall:\n")
  print(x$call)

  cat("\nMajor fields:\n")
  cat("$ncomp - number of calculated components\n")
  cat("$ncomp.selected - number of selected components\n")
  cat("$coeffs - object (regcoeffs) with regression coefficients\n")
  cat("$xloadings - vector with x loadings\n")
  cat("$yloadings - vector with y loadings\n")
  cat("$weights - vector with weights\n")
  cat("$res - list with results (calibration, cv, etc)\n")

  cat("\nTry summary(model) and plot(model) to see the model performance.\n")
}


################################
#  Plotting methods            #
################################


#' Explained X variance plot for PLS
#'
#' @description
#' Shows plot with explained X variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param main
#' title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXVariance.pls <- function(obj, type = "b", main = "Variance (X)", ...) {
  plotVariance(obj, decomp = "xdecomp", type = type, main = main, ...)
}

#' Explained Y variance plot for PLS
#'
#' @description
#' Shows plot with explained Y variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param main
#' title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotYVariance.pls <- function(obj, type = "b", main = "Variance (Y)", ...) {
  plotVariance(obj, decomp = "ydecomp", type = type, main = main, ...)
}

#' Cumulative explained X variance plot for PLS
#'
#' @description
#' Shows plot with cumulative explained X variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param main
#' title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXCumVariance.pls <- function(obj, type = "b", main = "Cumulative variance (X)", ...) {
  plotVariance(obj, decomp = "xdecomp", variance = "cumexpvar", type = type, main = main, ...)
}

#' Cumulative explained Y variance plot for PLS
#'
#' @description
#' Shows plot with cumulative explained Y variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param main
#' title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotYCumVariance.pls <- function(obj, type = "b", main = "Cumulative variance (Y)", ...) {
  plotVariance(obj, decomp = "ydecomp", variance = "cumexpvar", type = type, main = main, ...)
}

#' Variance plot for PLS
#'
#' @description
#' Shows plot with variance values vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param decomp
#' which decomposition to use ("xdecomp" for x or "ydecomp" for y)
#' @param variance
#' which variance to use ("expvar", "cumexpvar")
#' @param type
#' type of the plot("b", "l" or "h")
#' @param labels
#' what to show as labels for plot objects.
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ylab
#' label for y-axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotVariance.pls <- function(obj, decomp = "xdecomp", variance = "expvar", type = "b",
  labels = "values", res = obj$res, ylab = "Explained variance, %", ...) {

  plot_data <- lapply(res, plotVariance, decomp = decomp, variance = variance, show.plot = FALSE)
  mdaplotg(plot_data, labels = labels, type = type, ylab = ylab, ...)
}

#' X scores plot for PLS
#'
#' @description
#' Shows plot with X scores values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param main
#' main plot title
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXScores.pls <- function(obj, comp = c(1, 2), show.axes = T,  main = "Scores (X)",
  res = obj$res, ...) {


  # set up values for showing axes lines
  show.lines <- FALSE
  if (show.axes) {
    show.lines <- if (length(comp) == 2) c(0, 0) else c(NA, 0)
  }

  plot_data <- lapply(res, plotXScores, comp = comp, type = "p", show.plot = FALSE)
  mdaplotg(plot_data, show.lines = show.lines, type = "p", main = main, ...)
}

#' XY scores plot for PLS
#'
#' @description
#' Shows plot with X vs. Y scores values for selected component.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' which component to show the plot for
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXYScores.pls <- function(obj, ncomp = 1, show.axes = T,  res = obj$res, ...) {

  show.lines <- if (show.axes) c(0, 0) else FALSE
  plot_data <- lapply(res, plotXYScores, ncomp = ncomp, type = "p", show.plot = FALSE)
  mdaplotg(plot_data, show.lines = show.lines, type = "p", ...)
}

#' Residual distance plot for decomposition of X data
#'
#' @description
#' Shows a plot with orthogonal distance vs score distance for PLS decomposition of X data.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (by default optimal value selected for the model will be used)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param cgroup
#' color grouping of plot points (works only if one result object is available)
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
#' @param show.legend
#' logical, show or not a legend on the plot (needed if several result objects are available)
#' @param show.limits
#' vector with two logical values defining if limits for extreme and/or outliers must be shown
#' @param lim.col
#' vector with two values - line color for extreme and outlier limits
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier limits
#' @param lim.lty
#' vector with two values - line type for extreme and outlier limits
#' @param main
#' title for the plot
#' @param legend.position
#' position of legend (if shown)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The function is almost identical to \code{\link{plotResiduals.pca}}.
#'
#' @export
plotXResiduals.pls <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
  main = sprintf("X-distances (ncomp = %d)", ncomp), cgroup = NULL, xlim = NULL, ylim = NULL,
  show.limits = c(TRUE, TRUE), lim.col = c("darkgray", "darkgray"), lim.lwd = c(1, 1),
  lim.lty = c(2, 3), show.legend = TRUE, legend.position = "topright", res = obj$res, ...) {

  # get xdecomp from list with result objects
  res <- lapply(res, function(x) if ("ldecomp" %in% class(x$xdecomp)) x$xdecomp)
  res <- res[!sapply(res, is.null)]

  ldecomp.plotResiduals(res, obj$Qlim, obj$T2lim, ncomp = ncomp, log = log, norm = norm,
    cgroup = cgroup, xlim = xlim, ylim = ylim, show.limits = show.limits, lim.col = lim.col,
    lim.lwd = lim.lwd, show.legend = show.legend, main = main, ...)
}

#' Residual XY-distance plot
#'
#' @description
#' Shows a plot with full X-distance (f) vs. orthogonal Y-distance (z) for PLS model results.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (by default optimal value selected for the model will be used)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param cgroup
#' color grouping of plot points (works only if one result object is available)
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
#' @param show.legend
#' logical, show or not a legend on the plot (needed if several result objects are available)
#' @param show.limits
#' vector with two logical values defining if limits for extreme and/or outliers must be shown
#' @param lim.col
#' vector with two values - line color for extreme and outlier limits
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier limits
#' @param lim.lty
#' vector with two values - line type for extreme and outlier limits
#' @param main
#' title for the plot
#' @param legend.position
#' position of legend (if shown)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The function presents a way to identify extreme objects and outliers based on both full distance
#' for X-decomposition (known as f) and squared residual distance for Y-decomposition (z). The
#' approach has been proposed in [1].
#'
#' The plot is available only if data driven methods (classic or robust) have been used for
#' computing of critical limits.
#'
#' @references
#' 1. Rodionova O. Ye., Pomerantsev A. L. Detection of Outliers in Projection-Based Modeling.
#' Analytical Chemistry (2020, in publish). doi: 10.1021/acs.analchem.9b04611
#'
#' @export
plotXYResiduals.pls <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
  main = sprintf("XY-distances (ncomp = %d)", ncomp), cgroup = NULL, xlim = NULL, ylim = NULL,
  show.limits = c(TRUE, TRUE), lim.col = c("darkgray", "darkgray"), lim.lwd = c(1, 1),
  lim.lty = c(2, 3), show.legend = TRUE, legend.position = "topright", res = obj$res, ...) {

  if (!(obj$lim.type %in% c("ddmoments", "ddrobust"))) {
    stop("plotXYResiduals() works only with data driven limit types ('ddoments' or 'ddrobust')")
  }

  # generate values for cgroup if categories should be used
  if (length(cgroup) == 1 && cgroup == "categories") {
    cgroup <- categorize(obj, res[[1]], ncomp = ncomp)
  }

  # get xdecomp from list with result objects
  res <- lapply(res, function(x) if ("ldecomp" %in% class(x$xdecomp)) x)
  res <- res[!sapply(res, is.null)]

  # function to compute plot limits
  getPlotLim <- function(lim, pd, ld, dim, show.limits) {
    if (!is.null(lim) || all(!show.limits)) return(lim)
    limits <- if (show.limits[[2]]) ld$outliers else ld$extremes
    return(c(0, max(sapply(pd, function(x) max(x[, dim])), limits[, dim])) * 1.05)
  }

  # check that show.limits is logical
  if (!all(is.logical(show.limits))) {
    stop("Parameter 'show.limits' must have logical value(s).")
  }

  # if show.limits has only one value - duplicate it
  if (length(show.limits) == 1) {
    show.limits <- rep(show.limits, 2)
  }

  # compute plot data for each result object
  plot_data <- lapply(res, plotXYResiduals.plsres, ncomp = ncomp, norm = norm, log = log,
    show.plot = FALSE)

  # get coordinates for critical limits
  lim_data <- pls.getLimitsCoordinates(obj$Qlim, obj$T2lim, obj$Zlim,
    ncomp = ncomp, nobj = obj$limParams$Q$moments$nobj, norm = norm, log = log)

  xlim <- getPlotLim(xlim, plot_data, lim_data, 1, show.limits)
  ylim <- getPlotLim(ylim, plot_data, lim_data, 2, show.limits)

  # make plot
  if (length(plot_data) == 1) {
    mdaplot(plot_data[[1]], type = "p", xlim = xlim, ylim = ylim, cgroup = cgroup, main = main, ...)
  } else {
    mdaplotg(plot_data, type = "p", xlim = xlim, ylim = ylim, show.legend = show.legend, main = main,
      legend.position = legend.position, ...)
  }

  # show critical limits
  if (show.limits[[1]]) {
    lines(lim_data$extremes[, 1], lim_data$extremes[, 2],
      col = lim.col[1], lty = lim.lty[1], lwd = lim.lwd[1])
  }

  if (show.limits[[2]]) {
    lines(lim_data$outliers[, 1], lim_data$outliers[, 2],
      col = lim.col[2], lty = lim.lty[2], lwd = lim.lwd[2])
  }
}

#' X loadings plot for PLS
#'
#' @description
#' Shows plot with X loading values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param type
#' type of the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.legend
#' logical, show or not legend on the plot (when it is available)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXLoadings.pls <- function(obj, comp = c(1, 2), type = "p", show.axes = TRUE,
  show.legend = TRUE, ...) {

  plot_data <- mda.subset(obj$xloadings, select = comp)
  colnames(plot_data) <- sprintf("Comp %d (%.2f%%)", comp, obj$res[["cal"]]$xdecomp$expvar[comp])
  attr(plot_data, "name") <- "Loadings (X)"

  # set up values for showing axes lines
  show.lines <- FALSE
  if (show.axes) {
    show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
  }

  if (type == "p") {
    return(mdaplot(plot_data, type = type, show.lines = show.lines, ...))
  }

  plot_data <- mda.t(plot_data)
  attr(plot_data, "yaxis.name") <- "Loading"
  mdaplotg(plot_data, show.legend = show.legend, type = type, show.lines = show.lines, ...)
}

#' X loadings plot for PLS
#'
#' @description
#' Shows plot with X loading values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param type
#' type of the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.legend
#' logical, show or not a legend
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotWeights.pls <- function(obj, comp = 1, type = (if (nrow(obj$weights) < 20) "h" else "l"),
  show.axes = TRUE, show.legend = TRUE, ...) {

  plot_data <- mda.subset(obj$weights, select = comp)
  colnames(plot_data) <- sprintf("Comp %d (%.2f%%)", comp, obj$res[["cal"]]$xdecomp$expvar[comp])
  attr(plot_data, "name") <- "Weights"

  # set up values for showing axes lines
  show.lines <- FALSE
  if (show.axes) {
    show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
  }

  if (type == "p") {
    return(mdaplot(plot_data, type = type, show.lines = show.lines, ...))
  }

  plot_data <- mda.t(plot_data)
  attr(plot_data, "yaxis.name") <- "Weight"
  mdaplotg(plot_data, show.legend = show.legend, type = type, show.lines = show.lines, ...)
}

#' XY loadings plot for PLS
#'
#' @description
#' Shows plot with X and Y loading values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXYLoadings.pls <- function(obj, comp = c(1, 2), show.axes = TRUE, ...) {

  if (length(comp) != 2) {
    stop("This plot can be made for only two components.")
  }

  plot_data <- list(
    "X" = mda.subset(obj$xloadings, select = comp),
    "Y" = mda.subset(obj$yloadings, select = comp)
  )
  colnames(plot_data[[1]]) <- sprintf("Comp %d (%.2f%%)", comp,
    obj$res[["cal"]]$xdecomp$expvar[comp])

  attr(plot_data, "name") <- "Loadings (XY)"
  show.lines <- if (show.axes) c(0, 0) else FALSE
  mdaplotg(plot_data, type = "p", show.lines = show.lines, ...)
}

#' VIP scores plot for PLS model
#'
#' @description
#' Shows a plot with VIP scores values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to plot the values for (if y is multivariate), can be a vector.
#' @param ncomp
#' number of components to count
#' @param type
#' type of the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See \code{\link{vipscores}} for more details.
#'
#' @export
plotVIPScores.pls <- function(obj, ny = 1, ncomp = obj$ncomp.selected,
  type = "l", ...) {

  vipscores <- vipscores(obj, ncomp = ncomp)
  mdaplotg(mda.t(mda.subset(vipscores, select = ny)), type = type, ...)
}

#' Selectivity ratio plot for PLS model
#'
#' @description
#' Computes and shows a plot for Selectivity ratio values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to plot the values for (if y is multivariate), can be a vector.
#' @param ncomp
#' number of components to count
#' @param type
#' type of the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See \code{\link{vipscores}} for more details.
#'
#' @export
plotSelectivityRatio.pls <- function(obj, ny = 1,
  ncomp = obj$ncomp.selected, type = "l", ...) {

  selratio <- selratio(obj, ncomp = ncomp)
  mdaplotg(mda.t(mda.subset(selratio, select = ny)), type = type, ...)
}

#' Model overview plot for PLS
#'
#' @description
#' Shows a set of plots (x residuals, regression coefficients, RMSE and predictions) for PLS model.
#'
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plot.pls <- function(x, ncomp = x$ncomp.selected, ny = 1, show.legend = TRUE, ...) {

  if (!is.null(ncomp) && (ncomp <= 0 || ncomp > x$ncomp)) {
    stop("Wrong value for number of components.")
  }

  par(mfrow = c(2, 2))
  plotXResiduals(x, ncomp = ncomp, show.legend = show.legend)
  plotRegcoeffs(x, ncomp = ncomp, ny = ny)
  plotRMSE(x, ny = ny, show.legend = show.legend)
  plotPredictions(x, ncomp = ncomp, ny = ny, show.legend = show.legend)
  par(mfrow = c(1, 1))
}


################################
#  Static methods              #
################################


#' Runs selected PLS algorithm
#'
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param y
#' a matrix with y values (responses from calibration set)
#' @param ncomp
#' how many components to compute
#' @param method
#' algorithm for computing PLS model
#' @param cv
#' logical, is this for CV or not
#'
#' @export
pls.run <- function(x, y, ncomp = min(nrow(x) - 1, ncol(x)), method = "simpls", cv = FALSE) {

  if (ncomp < 1 || ncomp > min(nrow(x) - 1, ncol(x))) {
    stop("Wrong value for 'ncomp' parameter.")
  }

  methods <- list("simpls" = pls.simpls)

  if (!(method %in% names(methods))) {
    stop("Method with this name is not supported.")
  }

  return(methods[[method]](x, y, ncomp, cv = cv))
}

#' SIMPLS algorithm
#'
#' @description
#' SIMPLS algorithm for calibration of PLS model
#'
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param cv
#' logical, is model calibrated during cross-validation or not
#'
#' @return
#' a list with computed regression coefficients, loadings and scores for x and y matrices,
#' and weights.
#'
#' @references
#' [1]. S. de Jong. SIMPLS: An Alternative approach to partial least squares regression.
#' Chemometrics and Intelligent Laboratory Systems, 18, 1993 (251-263).
#'
pls.simpls <- function(x, y, ncomp, cv = FALSE) {

  x <- as.matrix(x)
  y <- as.matrix(y)

  npred <- ncol(x)
  nresp <- ncol(y)

  # initial estimation
  A <- crossprod(x, y)
  M <- crossprod(x, x)
  C <- diag(npred)

  # prepare space for results
  B <- array(0, dim = c(npred, ncomp, nresp))
  W <- matrix(0, nrow = npred, ncol = ncomp)
  P <- matrix(0, nrow = npred, ncol = ncomp)
  Q <- matrix(0, nrow = nresp, ncol = ncomp)

  # loop for each components
  for (n in seq_len(ncomp)) {
    # get the dominate eigenvector of A'A
    e <- eigen(crossprod(A))
    q <- e$vectors[seq_len(nresp)]

    # calculate and store weights
    w <- A %*% q
    c <- as.numeric(crossprod(w, (M %*% w)))

    # stop cycle since c-value is very small and can result in singular matrix
    if (c < .Machine$double.eps) {
      n <- n - 1
      warning(paste0(
        "PLS can not compute more than ", n, " components (eigenvalues are too small). "
      ), call. = FALSE)
      break
    }

    w <- w / sqrt(c)
    W[, n] <- w

    # calculate and store x loadings
    p <- M %*% w
    P[, n] <- p

    # calculate and store y loadings
    q <- crossprod(A, w)
    Q[, n] <- q

    v <- C %*% p
    v <- v / sqrt(as.numeric(crossprod(v)))

    # compute coefficients for current component
    B[, n, ] <- tcrossprod(W[, seq_len(n), drop = FALSE], Q[, seq_len(n), drop = FALSE])

    # recalculate matrices for the next compnonent
    C <- C - tcrossprod(v)
    M <- M - tcrossprod(p)
    A <- C %*% A
  }

  # truncate results if n is smaller than ncomp
  W <- W[, seq_len(n), drop = FALSE]
  P <- P[, seq_len(n), drop = FALSE]
  Q <- Q[, seq_len(n), drop = FALSE]
  B <- B[, seq_len(n), , drop = FALSE]

  return(list(coeffs = B, weights = W, xloadings = P, yloadings = Q, ncomp = n))
}

#' PLS model calibration
#'
#' @description
#' Calibrates (builds) a PLS model for given data and parameters
#'
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' algorithm for computing PLS model (only 'simpls' is supported so far)
#' @param cv
#' logical, is model calibrated during cross-validation or not (or cv settings for calibration)
#'
#' @return model
#' an object with calibrated PLS model
#'
pls.cal <- function(x, y, ncomp, center, scale, method = "simpls", cv = FALSE) {

  # prepare empty list for model object and assign
  # several properties, which do not depend on calculations below
  model <- list()
  model$center <- center
  model$scale <- scale
  model$method <- method
  class(model) <- c("pls", "regmodel")

  # get attributes
  x.attrs <- mda.getattr(x)
  y.attrs <- mda.getattr(y)

  # get names of variables
  prednames <- colnames(x)
  respnames <- colnames(y)

  # if y is a vector convert it to a matrix
  if (is.null(dim(y))) dim(y) <- c(length(y), 1)

  # check dimensions
  if (nrow(x) != nrow(y)) {
    stop("Number of rows for predictors and responses should be the same.")
  }

  # convert data to a matrix
  x <- mda.df2mat(x)
  y <- mda.df2mat(y)

  # get dimension of original data
  x.ncols <- ncol(x)
  y.ncols <- ncol(y)

  # check if data has missing values
  if (any(is.na(x))) {
    stop("Predictors have missing values, try to fix this using pca.mvreplace.")
  }

  if (any(is.na(y))) {
    stop("Responses have missing values, try to fix this using pca.mvreplace.")
  }

  # set column names for predictors if missing
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("y", seq_len(ncol(y)))
  }

  # correct x-axis name
  if (is.null(x.attrs$xaxis.name)) {
    x.attrs$xaxis.name <- "Predictors"
  }

  if (is.null(y.attrs$xaxis.name)) {
    y.attrs$xaxis.name <- "Responses"
  }

  # remove excluded rows
  if (length(x.attrs$exclrows) > 0) {
    x <- x[-x.attrs$exclrows, , drop = FALSE]
    y <- y[-x.attrs$exclrows, , drop = FALSE]
  }

  # autoscale and save the mean and std values for predictors
  x <- prep.autoscale(x, center = center, scale = scale)
  model$xcenter <- attr(x, "prep:center")
  model$xscale <- attr(x, "prep:scale")

  # autoscale and save the mean and std values for responses
  y <- prep.autoscale(y, center = center, scale = scale)
  model$ycenter <- attr(y, "prep:center")
  model$yscale <- attr(y, "prep:scale")

  # remove excluded columns for predictors
  if (length(x.attrs$exclcols) > 0) {
    x <- x[, -x.attrs$exclcols, drop = FALSE]
  }

  # remove excluded columns for responses
  if (length(y.attrs$exclcols) > 0) {
    y <- y[, -y.attrs$exclcols, drop = FALSE]
  }

  # get dimensions of reduced datasets
  xc.nrows <- nrow(x)
  xc.ncols <- ncol(x)
  yc.ncols <- ncol(y)

  # find maximum number of objects in a segment
  nobj.cv <- 0
  if (!is.logical(cv) && !is.null(cv)) {
    nseg <- if (is.numeric(cv)) cv else cv[[2]]
    nobj.cv <- if (nseg == 1) 1 else ceiling(xc.nrows / nseg)

    # we set cv to FALSE so fitting knows that it is not a part of cross-validation
    cv <- FALSE
  }

  # set cv to FALSE also if it was null (needed for correct call of pls.run() method)
  if (is.null(cv)) cv <- FALSE

  # correct maximum number of components
  ncomp <- min(xc.ncols, xc.nrows - 1 - nobj.cv, ncomp)

  # fit model
  fit <- pls.run(x, y, method = method, ncomp = ncomp, cv = cv)
  model$ncomp <- ncomp <- fit$ncomp

  # if it is for cross-validation return the results as is
  if (is.logical(cv) && cv) {
    model$coeffs <- regcoeffs(fit$coeffs)
    model$xloadings <- fit$xloadings
    model$weights <- fit$weights
    return(model)
  }

  # compute x-scores and residuals
  xscores <- x %*% (fit$weights %*% solve(crossprod(fit$xloadings, fit$weights)))
  yscores <- as.matrix(y) %*% fit$yloadings

  # compute eigenvalues
  xeigenvals <- colSums(xscores^2) / (xc.nrows - 1)
  attr(xeigenvals, "DoF") <- (xc.nrows - 1)
  yeigenvals <- colSums(yscores^2) / (xc.nrows - 1)
  attr(yeigenvals, "DoF") <- (xc.nrows - 1)

  # correct results related to predictors for missing columns in x
  # corresponding rows will be set to 0 and excluded
  xloadings <- matrix(0, nrow = x.ncols, ncol = ncomp)
  yloadings <- matrix(0, nrow = y.ncols, ncol = ncomp)
  weights <- matrix(0, nrow = x.ncols, ncol = ncomp)
  coeffs <- array(0, dim = c(x.ncols, ncomp, yc.ncols))

  pred_ind <- seq_len(x.ncols)
  if (length(x.attrs$exclcols) > 0) pred_ind <- pred_ind[-x.attrs$exclcols]

  resp_ind <- seq_len(y.ncols)
  if (length(y.attrs$exclcols) > 0) resp_ind <- resp_ind[-y.attrs$exclcols]

  # x-loadings
  xloadings[pred_ind, ] <- fit$xloadings
  xloadings <- mda.exclrows(xloadings, x.attrs$exclcols)

  # y-loadings
  yloadings[resp_ind, ] <- fit$yloadings
  yloadings <- mda.exclrows(yloadings, y.attrs$exclcols)

  # weights
  weights[pred_ind, ] <- fit$weights
  weights <- mda.exclrows(weights, x.attrs$exclcols)

  # coeffs
  coeffs[pred_ind, , ] <- fit$coeffs
  coeffs <- mda.exclrows(coeffs, x.attrs$exclcols)

  # set names and attributes
  compnames <- paste("Comp", seq_len(ncomp))

  # x-loadings and weights
  rownames(xloadings) <- rownames(weights) <- prednames
  colnames(xloadings) <- colnames(weights) <- compnames
  attr(xloadings, "name") <- "X loadings"
  attr(weights, "name") <- "Weights"
  attr(xloadings, "xaxis.name") <- attr(weights, "xaxis.name") <- "Components"
  attr(xloadings, "yaxis.name") <- attr(weights, "yaxis.name") <- x.attrs$xaxis.name
  attr(xloadings, "yaxis.values") <- attr(weights, "yaxis.values") <- x.attrs$xaxis.values

  # coefficients
  dimnames(coeffs) <- list(prednames, compnames, colnames(y))
  attr(coeffs, "yaxis.name") <- x.attrs$xaxis.name
  attr(coeffs, "yaxis.values") <- x.attrs$xaxis.values

  # y-loadings
  rownames(yloadings) <- respnames
  colnames(yloadings) <- colnames(xloadings)
  attr(yloadings, "name") <- "Y loadings"
  attr(yloadings, "xaxis.name") <- "Components"
  attr(yloadings, "yaxis.name") <- y.attrs$xaxis.name
  attr(yloadings, "yaxis.values") <- y.attrs$xaxis.values

  # set up and return model parameters
  model$xloadings <- xloadings
  model$yloadings <- yloadings
  model$weights <- weights
  model$coeffs <- regcoeffs(coeffs)

  model$ncomp.selected <- ncomp
  model$exclrows <- x.attrs$exclrows
  model$exclcols <- x.attrs$exclcols
  model$xeigenvals <- xeigenvals
  model$yeigenvals <- yeigenvals

  return(model)
}

#' VIP scores for PLS model
#'
#' @description
#' Calculates VIP (Variable Importance in Projection) scores for predictors for given number
#' of components and response variable.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to count
#'
#' @return
#' matrix \code{nvar x ny} with VIP score values (columns correspond to responses).
#'
#' @details
#' May take some time in case of large number of predictors Returns results as a column-vector,
#' with all necessary attributes inherited (e.g. xaxis.values, excluded variables, etc.). If you
#' want to make a plot use for example: \code{mdaplot(mda.t(v), type = "l")}, where \code{v} is
#' a vector with computed VIP scores. Or just try \code{\link{plotVIPScores.pls}}.
#'
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#'
#' @export
vipscores <- function(obj, ncomp = obj$ncomp.selected) {

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
    stop("Wrong value for the 'ncomp' parameter.")
  }

  # subset needed model parameters
  comp <- seq_len(ncomp)
  weights <- obj$weights[, comp, drop = FALSE]
  yloads <- obj$yloadings[, comp, drop = FALSE];

  # get eigenvalues and multiply them to degrees of freedom
  xeigenvals <- obj$xeigenvals[comp]

  # get number and indices of variables and adjust dimension for regcoeffs
  nvar <- nrow(weights)
  var_ind <- seq_len(nvar)

  # remove hidden variables
  if (length(obj$exclcols) > 0) {
    weights <- weights[-obj$exclcols, , drop = FALSE]
    var_ind <- var_ind[-obj$exclcols]
  }

  # prepare matrix for vipscores
  vipscores <- matrix(0, nrow = nvar, ncol = nrow(yloads))

  # normalize scores
  wnorm <- weights %*% diag(1 / sqrt(colSums(weights^2)), nrow = ncomp, ncol = ncomp)

  # compute sum of squares for explained y variance and normalize it
  ssq <- yloads^2 %*% diag(xeigenvals, nrow = ncomp, ncol = ncomp)
  ssq <- ssq %*% diag(1 / rowSums(ssq), nrow = ncomp, ncol = ncomp)

  # compute VIP scores
  vipscores[var_ind, ] <- sqrt(nvar * wnorm^2 %*% t(ssq))

  rownames(vipscores) <- rownames(obj$xloadings)
  colnames(vipscores) <- rownames(obj$yloadings)

  attr(vipscores, "exclrows") <- obj$exclcols
  attr(vipscores, "yaxis.values") <- attr(obj$xloadings, "yaxis.values")
  attr(vipscores, "yaxis.name") <- attr(obj$xloadings, "yaxis.name")
  attr(vipscores, "xaxis.name") <- ""
  attr(vipscores, "name") <- sprintf("VIP scores (ncomp = %d)", ncomp)

  return(vipscores)
}

#' VIP scores for PLS model
#'
#' @description
#' Returns vector with VIP scores values. This function is a proxy for \code{\link{vipscores}}
#' and will be removed in future releases.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to count
#' @param ...
#' other parameters
#'
#' @return
#' matrix \code{nvar x 1} with VIP score values
#'
#' @export
getVIPScores.pls <- function(obj, ncomp = obj$ncomp.selected, ...) {

  warning("This function is deprecated and will be removed in future. Use 'vipscores()' insted.")
  return(vipscores(obj, ncomp = ncomp))
}

#' Selectivity ratio calculation
#'
#' @description
#' Calculates selectivity ratio for each component and response variable in
#' the PLS model
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to count
#'
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#'
#' @return
#' array \code{nvar x ncomp x ny} with selectivity ratio values
#'
#' @export
selratio <- function(obj, ncomp = obj$ncomp.selected) {

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
    stop("Wrong value for the 'ncomp' parameter.")
  }

  # get number and indices of variables and adjust dimension for regcoeffs
  nvar <- nrow(obj$weights)
  nresp <- nrow(obj$yloadings)
  var_ind <- seq_len(nvar)

  # reproduce x values
  xresiduals <- obj$res[["cal"]]$xdecomp$residuals
  xscores <- obj$res[["cal"]]$xdecomp$scores
  x <- xresiduals + tcrossprod(xscores, obj$xloadings)

  # remove excluded rows
  if (length(obj$exclrows) > 0) {
    x <- x[-obj$exclrows, , drop = FALSE]
  }

  # subset needed model parameters
  coeffs <- obj$coeffs$values[, ncomp, , drop = FALSE]

  # correct dimension for coefficients
  dim(coeffs) <- c(nvar, nresp)

  # remove hidden variables
  if (length(obj$exclcols) > 0) {
    x <- x[, -obj$exclcols, drop = FALSE]
    coeffs <- coeffs[-obj$exclcols, , drop = FALSE]
    var_ind <- var_ind[-obj$exclcols]
  }

  # prepare matrix for vipscores
  selratio <- matrix(0, nrow = nvar, ncol = nresp)

  # get norm value for regression coefficients
  bnorm <- sqrt(colSums(coeffs^2))

  # compute target projections
  ttp <-  x %*% (coeffs %*% diag(1 / bnorm, nrow = nresp, ncol = nresp))
  ptp <- t(crossprod(x, ttp) %*% diag(1 / colSums(ttp^2), nrow = nresp, ncol = nresp))

  # compute selectivity ratio
  for (y in seq_len(nresp)) {
    expvar <- ttp[, y, drop = FALSE] %*% ptp[y, , drop = FALSE]
    selratio[var_ind, y] <- colSums(expvar^2) / colSums((x - expvar)^2)
  }

  rownames(selratio) <- rownames(obj$xloadings)
  colnames(selratio) <- rownames(obj$yloadings)

  attr(selratio, "exclrows") <- obj$exclcols
  attr(selratio, "yaxis.values") <- attr(obj$xloadings, "yaxis.values")
  attr(selratio, "yaxis.name") <- attr(obj$xloadings, "yaxis.name")
  attr(selratio, "xaxis.name") <- ""
  attr(selratio, "name") <- sprintf("Selectivity ratio (ncomp = %d)", ncomp)

  return(selratio)
}

#' Selectivity ratio for PLS model
#'
#' @description
#' Returns vector with Selectivity ratio values. This function is a proxy for \code{\link{selratio}}
#' and will be removed in future releases.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to get the values for (if NULL user selected as optimal will be used)
#' @param ...
#' other parameters
#'
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#'
#' @return
#' vector with selectivity ratio values
#'
#' @export
getSelectivityRatio.pls <- function(obj, ncomp = obj$ncomp.selected, ...) {
  warning("This function is deprecated and will be removed in future. Use 'selratio()' insted.")
  return(selratio(obj, ncomp = ncomp))
}

#' Compute critical limits for orthogonal distances (Q)
#'
#' @param lim.type
#' which method to use for calculation of critical limits for residuals
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#' @param params
#' distribution parameters returned by ldecomp.getLimParams
#'
#' @export
pls.getZLimits <- function(lim.type, alpha, gamma, params) {

  if (!(lim.type %in% c("ddmoments", "ddrobust"))) {
    return(NULL)
  }

  pZ <- if (regexpr("robust", lim.type) > 0) params$Z$robust else params$Z$moments
  DoF <- round(pZ$Nu)
  DoF[DoF < 1] <- 1
  DoF[DoF > 250] <- 250

  ncomp <- length(pZ$u0)
  lim <- rbind(0, 0, pZ$u0, DoF)

  colnames(lim) <- paste("Comp", seq_len(ncomp))
  rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
  attr(lim, "name") <- "Critical limits for orthogonal distance (Z)"
  attr(lim, "alpha") <- alpha
  attr(lim, "gamma") <- gamma
  attr(lim, "lim.type") <- lim.type

  return(lim)
}

#' Compute coordinates of lines or curves with critical limits
#'
#' @param Qlim
#' matrix with critical limits for orthogonal distances (X)
#' @param T2lim
#' matrix with critical limits for score distances (X)
#' @param Zlim
#' matrix with critical limits for orthogonal distances (Y)
#' @param nobj
#' number of objects to compute the limits for
#' @param ncomp
#' number of components for computing the coordinates
#' @param norm
#' logical, shall distance values be normalized or not
#' @param log
#' logical, shall log transformation be applied or not
#'
#' @return
#' list with two matrices (x and y coordinates of corresponding limits)
#'
#' @export
pls.getLimitsCoordinates <- function(Qlim, T2lim, Zlim, nobj, ncomp, norm, log) {

  # get DoF
  Nh <- T2lim[4, ncomp]
  Nq <- Qlim[4, ncomp]
  Nz <- Zlim[4, ncomp]
  Nf <- Nq + Nh

  # get scaling factor
  z0 <- Zlim[3, ncomp]
  f0 <- Nf

  # process degrees of freedom for (Z)
  Nz <- round(Nz)
  Nz[Nz < 1] <- 1
  Nz[Nz > 250] <- 250

  # process degrees of freedom for (F)
  Nf <- round(Nf)
  Nf[Nf < 1] <- 1
  Nf[Nf > 250] <- 250

  # get limit parameters
  alpha <- attr(Qlim, "alpha")
  gamma <- attr(Qlim, "gamma")

  ## slope and intercepts
  eB <- qchisq(1 - alpha, Nf + Nz) / Nz * z0
  oB <- qchisq((1 - gamma) ^ (1 / nobj), Nf + Nz) / Nz * z0
  eA <- oA <- -1 * (z0 / f0) * (Nf / Nz)

  fE <- seq(-0.95, -eB / eA, length.out = 100)
  fO <- seq(-0.95, -oB / oA, length.out = 100)
  zE <- eA * fE + eB
  zO <- oA * fO + oB

  if (norm) {
    fE <- fE / f0
    zE <- zE / z0
    fO <- fO / f0
    zO <- zO / z0
  }

  if (log) {
    fE <- log(1 + fE)
    zE <- log(1 + zE)
    fO <- log(1 + fO)
    zO <- log(1 + zO)
  }

  return(list(
    extremes = cbind(fE, zE),
    outliers = cbind(fO, zO)
  ))
}

#' Convert image to data matrix
#'
#' @param img
#' an image (3-way array)
#'
#' @export
mda.im2data <- function(img) {
  width <- dim(img)[2]
  height <- dim(img)[1]
  nchannels <- dim(img)[3]

  npixels <- width * height
  dim(img) <- c(npixels, nchannels)
  attr(img, "width") <- width
  attr(img, "height") <- height

  return(img)
}

#' Convert data matrix to an image
#'
#' @param data
#' data matrix
#'
#' @export
mda.data2im <- function(data) {
  width <- attr(data, "width", exact = TRUE)
  height <- attr(data, "height", exact = TRUE)
  bgpixels <- attr(data, "bgpixels", exact = TRUE)

  if (length(bgpixels) > 0) {
    img <- matrix(NA, nrow = nrow(data) + length(bgpixels), ncol = ncol(data))
    img[-bgpixels, ] <- data
  } else {
    img <- data
  }

  dim(img) <- c(height, width, ncol(data))
  return(img)
}

#' Remove background pixels from image data
#'
#' @param data
#' a matrix with image data
#' @param bgpixels
#' vector with indices or logical values corresponding to background pixels
#'
#' @export
mda.setimbg <- function(data, bgpixels) {
  attrs <- mda.getattr(data)

  if (length(attrs$exclrows) > 0) {
    stop("You can not set background pixels if some of them have been already excluded.")
  }

  # unfold bgpixels to a vector
  dim(bgpixels) <- NULL

  # get indices instead of logical values
  if (is.logical(bgpixels)) {
    bgpixels <- which(bgpixels)
  }

  # correct indices of bgpixels if some of the pixels were already removed
  if (length(attrs$bgpixels) > 0) {
    npixels <- attrs$width * attrs$height
    row.ind <- seq_len(npixels)
    row.ind <- row.ind[-attrs$bgpixels]
    bgpixels <- row.ind[bgpixels]
  }

  # remove corresponding rows and correct attributes
  data <- data[-bgpixels, , drop = FALSE]
  attrs$bgpixels <- unique(c(attrs$bgpixels, bgpixels))

  data <- mda.setattr(data, attrs)
  return(data)
}

#' show image data as an image
#'
#' @param data
#' data with image
#' @param channels
#' indices for one or three columns to show as image channels
#' @param show.excluded
#' logical, if TRUE the method also shows the excluded (hidden) pixels
#' @param main
#' main title for the image
#' @param colmap
#' colormap using to show the intensity levels
#'
#' @export
imshow <- function(data, channels = 1, show.excluded = FALSE,
  main = paste0(" ", colnames(data)[channels]), colmap = "jet") {

  attrs <- mda.getattr(data)

  data <- mda.subset(data, select = channels)
  data <- (data - min(data)) / (max(data) - min(data))
  data <- mda.data2im(data)

  bg <- is.na(data)

  nrows <- dim(data)[1]
  ncols <- dim(data)[2]

  if (is.character(colmap) && length(colmap) == 1) {
    colmap <- if (colmap == "gray") colorRampPalette(c("#000000", "#ffffff"), space = "Lab")(256)
    else mdaplot.getColors(256, NULL, colmap)
  }

  if (length(channels) == 1) {
    nrows <- nrow(data)
    image(t(data[seq(nrows, 1, -1), , 1]), xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1),
      main = main, useRaster = TRUE, col = colmap, axes = FALSE)

    if (any(bg)) {
      bgimg <- matrix(NA, nrows, ncols)
      bgimg[bg[, , 1]] <- 0
      rasterImage(bgimg, 0, 0, 1, 1)
    }
  } else {
    if (any(bg)) data[bg] <- 0
    plot(0, main = main, type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "",
      xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    rasterImage(data, 0, 0, 1, 1)
  }

  # hide excluded pixels with dark gray color
  if (show.excluded == FALSE && length(attrs$exclrows) > 0) {
    npixels <- nrows * ncols
    ind <- seq_len(npixels)
    if (length(attrs$bgpixels) > 0) ind <- ind[-attrs$bgpixels]
    eximage <- rep(NA, npixels)
    eximage[ind[attrs$exclrows]] <- 0.25
    dim(eximage) <- c(nrows, ncols)
    rasterImage(eximage, 0, 0, 1, 1)
  }

}

#' Wrapper for show() method
#'
#' @param x
#' data set
#' @param n
#' number of rows to show
#'
#' @export
mda.show <- function(x, n = 50) {
  exclrows <- attr(x, "exclrows", exact = TRUE)
  exclcols <- attr(x, "exclcols", exact = TRUE)

  name <- attr(x, "name", exact = TRUE)

  if (!is.null(name) && nchar(name) > 0) {
    fprintf("%s\n%s\n", name, paste(rep("-", nchar(name)), collapse = ""))
  }

  if (!is.null(exclrows)) {
    x <- x[-exclrows, , drop = FALSE]
  }

  if (!is.null(exclcols)) {
    x <- x[, -exclcols, drop = FALSE]
  }

  if (n > nrow(x)) n <- nrow(x)
  show(x[seq_len(n), , drop = FALSE])
}

#' A wrapper for subset() method with proper set of attributed
#'
#' @param x
#' dataset (data frame or matrix)
#' @param subset
#' which rows to keep (indices, names or logical values)
#' @param select
#' which columns to select (indices, names or logical values)
#'
#' @return
#' a data with the subset
#'
#' @details
#' The method works similar to the standard \code{subset()} method, with minor differences. First
#' of all it keeps (and correct, if necessary) all important attributes. If only columns are
#' selected, it keeps all excluded rows as excluded. If only rows are selected, it keeps all
#' excluded columns. If both rows and columns are selected it removed all excluded elements first
#' and then makes the subset.
#'
#' The parameters \code{subset} and \code{select} may each be a vector with numbers or nanes
#' without excluded elements, or a logical expression.
#'
#' @export
mda.subset <- function(x, subset = NULL, select = NULL) {

  if (is.null(x)) return(NULL)

  attrs <- mda.getattr(x)

  if (!is.null(subset)) {
    if (is.logical(subset) & !is.null(attrs$exclrows))
      subset <- subset[-attrs$exclrows]

    # remove excluded rows first
    if (!is.null(attrs$exclrows))
      x <- x[-attrs$exclrows, , drop = F]

    # get numeric indices for the rows and subset them
    subset <- mda.getexclind(subset, rownames(x), nrow(x))
    x <- x[subset, , drop = FALSE]

    # correct attributes
    if (!is.null(attrs$yaxis.values)) {
      if (!is.null(attrs$exclrows)) attrs$yaxis.values <- attrs$yaxis.values[-attrs$exclrows]
      attrs$yaxis.values <- attrs$yaxis.values[subset]
    }

    attrs$exclrows <- NULL
  }

  if (!is.null(select)) {
    if (is.logical(select) && !is.null(attrs$exclcols))
      select <- select[-attrs$exclcols]

    # remove excluded rows first
    if (!is.null(attrs$exclcols))
      x <- x[, -attrs$exclcols, drop = F]

    # get numeric indices for the rows and subset them
    select <- mda.getexclind(select, colnames(x), ncol(x))
    x <- x[, select, drop = F]

    # correct attributes
    if (!is.null(attrs$xaxis.values)) {
      if (!is.null(attrs$exclcols)) attrs$xaxis.values <- attrs$xaxis.values[-attrs$exclcols]
      attrs$xaxis.values <- attrs$xaxis.values[select]
    }

    attrs$exclcols <- NULL
  }

  x <- mda.setattr(x, attrs)
  return(x)
}

#' A wrapper for rbind() method with proper set of attributes
#'
#' @param ...
#' datasets (data frames or matrices) to bind
#'
#' @return
#' the merged datasets
#'
#' @export
mda.rbind <- function(...) {
  objects <- list(...)
  nobj <- length(objects)

  attrs <- mda.getattr(objects[[1]])
  out.exclrows <- attrs$exclrows
  out.yaxis.values <- attrs$yaxis.values

  out.x <- objects[[1]]
  for (i in 2:nobj) {
    x <- objects[[i]]
    exclrows <- attr(x, "exclrows", exact = TRUE)
    yaxis.values <- attr(x, "yaxis.values")
    if (!is.null(exclrows)) out.exclrows <- c(out.exclrows, exclrows + nrow(out.x))
    if (is.null(out.yaxis.values) || is.null(yaxis.values))
      out.yaxis.values <- NULL
    else
      out.yaxis.values <- c(out.yaxis.values, yaxis.values)
    out.x <- rbind(out.x, x)
  }

  out.x <- mda.setattr(out.x, attrs)
  attr(out.x, "exclrows") <- out.exclrows
  attr(out.x, "yaxis.values") <- out.yaxis.values

  return(out.x)
}

#' A wrapper for cbind() method with proper set of attributes
#'
#' @param ...
#' datasets (data frames or matrices) to bind
#'
#' @return
#' the merged datasets
#'
#' @export
mda.cbind <- function(...) {
  objects <- list(...)
  nobj <- length(objects)

  attrs <- mda.getattr(objects[[1]])
  out.exclcols <- attrs$exclcols
  out.xaxis.values <- attrs$xaxis.values
  out.x <- objects[[1]]

  for (i in 2:nobj) {
    x <- objects[[i]]
    exclcols <- attr(x, "exclcols")
    xaxis.values <- attr(x, "xaxis.values")
    if (!is.null(exclcols))
      out.exclcols <- c(out.exclcols, exclcols + ncol(out.x))
    if (is.null(out.xaxis.values) || is.null(xaxis.values))
      out.xaxis.values <- NULL
    else
      out.xaxis.values <- c(out.xaxis.values, xaxis.values)
    out.x <- cbind(out.x, x)
  }

  out.x <- mda.setattr(out.x, attrs)
  attr(out.x, "exclcols") <- out.exclcols
  attr(out.x, "xaxis.values") <- out.xaxis.values

  return(out.x)
}

#' A wrapper for t() method with proper set of attributes
#'
#' @param x
#' dataset (data frames or matrices) to transpose
#'
#' @return
#' the transposed dataset
#'
#' @export
mda.t <- function(x) {
  attrs <- mda.getattr(x)
  out.attrs <- attrs
  out.attrs$exclrows <- attrs$exclcols
  out.attrs$exclcols <- attrs$exclrows
  out.attrs$xaxis.name <- attrs$yaxis.name
  out.attrs$yaxis.name <- attrs$xaxis.name
  out.attrs$xaxis.values <- attrs$yaxis.values
  out.attrs$yaxis.values <- attrs$xaxis.values

  x <- t(x)
  x <- mda.setattr(x, out.attrs)
}

#' Exclude/hide rows in a dataset
#'
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' indices of rows to exclude (numbers, names or logical values)
#'
#' @return
#' dataset with excluded rows
#'
#' @details
#' The method assign attribute 'exclrows', which contains number of rows, which should be
#' excluded/hidden from calculations and plots (without removing them physically). The
#' argument \code{ind} should contain rows numbers (excluding already hidden), names or logical
#' values.
#'
#' @export
mda.exclrows <- function(x, ind) {

  if (length(ind) < 1) return(x)

  excl.rows <- attr(x, "exclrows", exact = TRUE)
  nrows.tot <- nrow(x)
  nrows.excl <- length(excl.rows)

  if (nrows.excl == 0) {
    # no objects are excluded yet
    attr(x, "exclrows") <- mda.getexclind(ind, rownames(x), nrows.tot)
  } else {
    # some objects were excluded before
    if (is.logical(ind))
      ind <- ind[-excl.rows]
    ind <- mda.getexclind(ind, rownames(x)[-excl.rows], nrows.tot - nrows.excl)
    ind.tot <- seq_len(nrows.tot)
    ind.tot <- ind.tot[-excl.rows]
    attr(x, "exclrows") <- sort(unique(c(ind.tot[ind], excl.rows)))
  }

  # check that number of rows is still sufficient
  if (is.null(nrow(x)) || nrow(x) == 0) {
    stop("No rows left when excluded hidden values.")
  }

  return(x)
}

#' include/unhide the excluded rows
#'
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' number of excluded rows to include
#'
#' @return
#' dataset with included rows
#'
#' @description
#' include rows specified by user (earlier excluded using mda.exclrows)
#'
#' @export
mda.inclrows <- function(x, ind) {
  excl.rows <- attr(x, "exclrows", exact = TRUE)
  ind.log <- excl.rows %in% ind
  attr(x, "exclrows") <- excl.rows[!ind.log]

  return(x)
}

#' Exclude/hide columns in a dataset
#'
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' indices of columns to exclude (numbers, names or logical values)
#'
#' @return
#' dataset with excluded columns
#'
#' @details
#' The method assign attribute 'exclcols', which contains number of columns, which should be
#' excluded/hidden from calculations and plots (without removing them physically). The argument
#' \code{ind} should contain column numbers (excluding already hidden), names or logical values.
#'
#' @export
mda.exclcols <- function(x, ind) {
  if (length(ind) < 1) return(x)

  excl.cols <- attr(x, "exclcols", exact = TRUE)
  ncols.tot <- ncol(x)
  ncols.excl <- length(excl.cols)

  if (ncols.excl == 0) {
    # no objects are excluded yet
    attr(x, "exclcols") <- mda.getexclind(ind, colnames(x), ncols.tot)
    return(x)
  }

  # some objects were excluded before
  if (is.logical(ind)) ind <- ind[-excl.cols]
  ind <- mda.getexclind(ind, colnames(x)[-excl.cols], ncols.tot - ncols.excl)
  ind.tot <- seq_len(ncols.tot)
  ind.tot <- ind.tot[-excl.cols]
  attr(x, "exclcols") <- sort(unique(c(ind.tot[ind], excl.cols)))
  return(x)
}

#' Include/unhide the excluded columns
#'
#' @param x
#' dataset (data frame or matrix).
#' @param ind
#' number of excluded columns to include
#'
#' @return
#' dataset with included columns.
#'
#' @description
#' include colmns specified by user (earlier excluded using mda.exclcols)
#'
#' @export
mda.inclcols <- function(x, ind) {
  excl.cols <- attr(x, "exclcols", exact = TRUE)
  ind.log <- excl.cols %in% ind
  attr(x, "exclcols") <- excl.cols[!ind.log]

  return(x)
}

#' Set data attributes
#'
#' @description
#' Set most important data attributes (name, xvalues, excluded rows and columns, etc.) to a dataset
#'
#' @param x
#' a dataset
#' @param attrs
#' list with attributes
#' @param type
#' a text variable telling which attributes to set ('all', 'row', 'col')
#'
#' @export
mda.setattr <- function(x, attrs, type = "all") {

  attr(x, "name") <- attrs$name
  attr(x, "width") <- attrs$width
  attr(x, "height") <- attrs$height
  attr(x, "bgpixels") <- attrs$bgpixels

  if (type != "col") {
    attr(x, "yaxis.name") <- attrs$yaxis.name
    attr(x, "yaxis.values") <- attrs$yaxis.values
    attr(x, "exclrows") <- attrs$exclrows
  }

  if (type != "row") {
    attr(x, "xaxis.name") <- attrs$xaxis.name
    attr(x, "xaxis.values") <- attrs$xaxis.values
    attr(x, "exclcols") <- attrs$exclcols
  }

  return(x)
}

#'
#' Get data attributes
#'
#' @description
#' Returns a list with important data attributes (name, xvalues, excluded rows and columns, etc.)
#'
#' @param x
#' a dataset
#'
#' @export
mda.getattr <- function(x) {
  attrs <- list()

  attrs$name <- attr(x, "name", exact = TRUE)
  attrs$exclrows <- attr(x, "exclrows", exact = TRUE)
  attrs$exclcols <- attr(x, "exclcols", exact = TRUE)
  attrs$xaxis.values <- attr(x, "xaxis.values", exact = TRUE)
  attrs$yaxis.values <- attr(x, "yaxis.values", exact = TRUE)
  attrs$xaxis.name <- attr(x, "xaxis.name", exact = TRUE)
  attrs$yaxis.name <- attr(x, "yaxis.name", exact = TRUE)
  attrs$width <- attr(x, "width", exact = TRUE)
  attrs$height <- attr(x, "height", exact = TRUE)
  attrs$bgpixels <- attr(x, "bgpixels", exact = TRUE)

  return(attrs)
}

#' Get indices of excluded rows or columns
#'
#' @param excl
#' vector with excluded values (logical, text or numbers)
#' @param names
#' vector with names for rows or columns
#' @param n
#' number of rows or columns
#'
#' @export
mda.getexclind <- function(excl, names, n) {
  nitems <- if (is.logical(excl)) sum(excl) else length(excl)

  if (is.character(excl))
    excl <- which(names %in% excl)
  if (is.logical(excl))
    excl <- which(excl)

  if (length(excl) < nitems)
    stop("At least one index or name is incorrect.")

  if (length(excl) > 0 && (!is.numeric(excl) || min(excl) < 1 || max(excl) > n))
    stop("At least one index or name is incorrect.")

  return(excl)
}

#' Convert data frame to a matrix
#'
#' @description
#' The function converts data frame to a numeric matrix.
#'
#' @param x
#' a data frame
#' @param full
#' logical, if TRUE number of dummy variables for a factor will be the same as number of levels,
#' otherwise by one smaller
#'
#' @details
#' If one or several columns of the data frame are factors they will be converted to a set of dummy
#' variables. If any columns/rows were hidden in the data frame they will remain hidden in the
#' matrix. If there are factors among the hidden columns, the corresponding dummy variables will be
#' hidden as well.
#'
#' All other attributes (names, axis names, etc.) will be inherited.
#'
#' @return
#' a numeric matrix
#'
#' @export
mda.df2mat <- function(x, full = FALSE) {
  attrs <- mda.getattr(x)
  if (is.null(x) || is.matrix(x) || is.vector(x)) return(x)

  if (is.factor(x)) {
    x <- data.frame(x)
  }

  if (any(sapply(x, is.character))) {
    stop("At least one column in the provided data frame has text values.", call. = FALSE)
  }

  # get indices of factor and numeric columns
  col.fac <- unlist(lapply(x, is.factor))
  col.num <- which(!col.fac)
  col.fac <- which(col.fac)

  dummy <- function(i, x, full = FALSE) {
    x <- x[, i]
    n <- if (full) nlevels(x) else nlevels(x) - 1
    y <- matrix(seq_len(n), nrow = length(x), ncol = n, byrow = TRUE)
    d <- y == as.numeric(x)
    colnames(d) <- levels(x)[seq_len(n)]
    attr(d, "cols.info") <- c(i, n)
    return(d)
  }

  if (is.null(col.fac) || length(col.fac) == 0) {
    # no factors among columns - easy job
    x <- as.matrix(x)
    x <- mda.setattr(x, attrs)
    return(x)
  }

  exclcols.fac.ind <- NULL
  exclcols.num.ind <- NULL
  if (!is.null(attrs$exclcols)) {
    if (is.character(attrs$exclcols)) attrs$exclcols <- which(colnames(x) %in% attrs$exclcols)
    if (is.logical(attrs$exclcols)) attrs$exclcols <- which(attrs$exclcols)

    exclcols.fac.ind <- which(col.fac %in% attrs$exclcols) # hidden factors
    exclcols.num.ind <- which(col.num %in% attrs$exclcols) # hidden numeric columns
  }

  # split data to numeric columns and factors
  num.data <- if (length(col.fac) < ncol(x)) as.matrix(x[, -col.fac, drop = FALSE])
  fac.data <- x[, col.fac, drop = FALSE]

  fac.data.hidden <- NULL
  if (length(exclcols.fac.ind) > 0) {
    fac.data.hidden <- fac.data[, exclcols.fac.ind, drop = FALSE]
    fac.data <- fac.data[, -exclcols.fac.ind, drop = FALSE]
  }

  # convert all non-excluded factors to dummy variables
  fac.data <- lapply(seq_len(ncol(fac.data)), dummy, x = fac.data, full = full)
  fac.data <- do.call(cbind, fac.data)

  # convert all excluded factors to numeric values
  exclcols.fac.ind <- NULL
  if (!is.null(fac.data.hidden)) {
    fac.data.hidden <- as.matrix(as.data.frame(lapply(fac.data.hidden, as.numeric)))
    n.incl.col <- ncol(num.data) + ncol(fac.data)
    exclcols.fac.ind <- (n.incl.col + 1):(n.incl.col + ncol(fac.data.hidden))
  }

  # combine the data values and set attributes
  x <- cbind(num.data, fac.data, fac.data.hidden)

  # correct and set arguments
  attrs$exclcols <- c(exclcols.num.ind, exclcols.fac.ind)
  x <- mda.setattr(x, attrs)
  return(x)
}

#' Removes excluded (hidden) rows from data
#'
#' @param data
#' data frame or matrix with data
#'
#' @export
mda.purgeRows <- function(data) {
  attrs <- mda.getattr(data)
  if (length(attrs$exclrows) == 0) return(data)
  new_data <- data[-attrs$exclrows, , drop = FALSE]
  attrs$yaxis.values <- if (!is.null(attrs$yaxis.values)) attrs$yaxis.values[-attrs$exclrows]
  attrs$exclrows <- NULL
  new_data <- mda.setattr(new_data, attrs)
  return(new_data)
}

#' Removes excluded (hidden) colmns from data
#'
#' @param data
#' data frame or matrix with data
#'
#' @export
mda.purgeCols <- function(data) {
  attrs <- mda.getattr(data)
  if (length(attrs$exclcols) == 0) return(data)

  new_data <- data[, -attrs$exclcols, drop = FALSE]
  attrs$xaxis.values <- if (!is.null(attrs$xaxis.values)) attrs$xaxis.values[-attrs$exclcols]
  attrs$exclcols <- NULL
  new_data <- mda.setattr(new_data, attrs)
  return(new_data)
}

#' Removes excluded (hidden) rows and colmns from data
#'
#' @param data
#' data frame or matrix with data
#'
#' @export
mda.purge <- function(data) {
  return(mda.purgeCols(mda.purgeRows(data)))
}

#' Get selected components
#'
#' @description
#' returns number of components depending on a user choice
#'
#' @param obj
#' an MDA model or result object (e.g. \code{pca}, \code{pls}, \code{simca}, etc)
#' @param ncomp
#' number of components to select, provided by user
#'
#' @details
#' Depedning on a user choice it returns optimal number of component for the model (if
#' use did not provide any value) or check the user choice for correctness and returns
#' it back
#'
getSelectedComponents <- function(obj, ncomp = NULL) {
  if (!is.null(ncomp)) return(ncomp)
  return(if (is.null(obj$ncomp.selected)) 1 else obj$ncomp.selected)
}

#' Get main title
#'
#' @description
#' returns main title for a plot depending on a user choice
#'
#' @param main
#' main title of a plot, provided by user
#' @param ncomp
#' number of components to select, provided by user
#' @param default
#' default title for the plot
#'
#' @details
#' Depedning on a user choice it returns main title for a plot
#'
getMainTitle <- function(main, ncomp, default) {
  if (!is.null(main)) return(main)
  return(if (is.null(ncomp)) default else sprintf("%s (ncomp = %d)", default, ncomp))
}

#' Imitation of fprinf() function
#'
#' @param ...
#' arguments for sprintf function
#'
#' @export
fprintf <- function(...) {
  cat(sprintf(...))
}

#' Return list with valid results
#'
#' @param res
#' list with results
#' @param classname
#' name of class (for result object) to look for
#'
#' @export
getRes <- function(res, classname = "ldecomp") {

  if (!is.list(res)) {
    stop("Parameter 'res' should be a list with result objects.")
  }

  res <- res[sapply(res, function(x) classname %in% class(x))]

  if (length(res) == 0) {
    stop("No valid results provided.")
  }

  return(res)
}

#' Capitalize text or vector with text values
#'
#' @param str
#' text of vector with text values
#'
#' @export
capitalize <- function(str) {
  return(sapply(str,  function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))))
}

#' Replicate matric x
#'
#' @param x
#' original matrix
#' @param nrows
#' number of times replicate matrix row wise
#' @param ncols
#' number of times replicate matrix columns wise
#'
#' @export
repmat <- function(x, nrows, ncols = nrows) {
  x <- as.matrix(x)
  return(matrix(1, nrows, ncols) %x% x)
}


#' Prepares calibration data
#'
#' @param x
#' matrix or data frame with values (calibration set)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param min.nrows
#' smallest number of rows which must be in the dataset
#' @param min.ncols
#' smallest number of columns which must be in the dataset
#'
#' @export
prepCalData <- function(x, exclrows = NULL, exclcols = NULL, min.nrows = 1, min.ncols = 2) {

  # check that x has a dimension
  stopifnot("Data values must be provided in form of a matrix or a data frame." = !is.null(dim(x)))

  if (is.data.frame(x) && any(sapply(x, is.character))) {
    stop("At least one column in the provided data frame has text values.", call. = FALSE)
  }

  # exclude columns if "exclcols" is provided
  if (length(exclcols) > 0) {
    x <- mda.exclcols(x, exclcols)
  }

  # exclude rows if "exclrows" is provided
  if (length(exclrows) > 0) {
    x <- mda.exclrows(x, exclrows)
  }

  # check number of rows
  if (nrow(x) - length(attr(x, "exclrows")) < min.nrows) {
    stop(sprintf("Dataset should contain at least %d measurements (rows).", min.nrows))
  }

  # check number of columns
  if (ncol(x) - length(attr(x, "exclcols")) < min.ncols) {
    stop(sprintf("Dataset should contain at least %d variables (columns).", min.ncols))
  }

  return(x)
}

#' Define parameters based on 'cv' value
#'
#' @param cv
#' settings for cross-validation provided by user
#' @param nobj
#' number of objects in calibration set
#'
crossval.getParams <- function(cv, nobj) {

  nrep <- 1

  # random
  if (is.numeric(cv)) {
    return(
      list(
        type = "rand",
        nrep = 1,
        nseg = if (cv == 1) nobj else cv
      )
    )
  }

  # leave one out
  type <- cv[[1]]
  if (type == "loo") {
    return(
      list(
        type = "rand",
        nrep = nrep,
        nseg = nobj
      )
    )
  }

  # venetian blinds
  nseg <- cv[[2]]
  if (type == "ven") {
    return(
      list(
        type = "ven",
        nrep = nrep,
        nseg = nseg
      )
    )
  }

  nrep <- if (length(cv) == 3) cv[[3]] else 1
  return(
    list(
      type = type,
      nrep = nrep,
      nseg = nseg
    )
  )
}

#' Generate sequence of indices for cross-validation
#'
#' @description
#' Generates and returns sequence of object indices for each segment in random segmented
#' cross-validation
#'
#' @param nobj
#' number of objects in a dataset
#' @param cv
#' cross-validation settings, can be a number or a list. If cv is a number, it will be
#' used as a number of segments for random cross-validation (if cv = 1, full cross-validation
#' will be preformed), if it is a list, the following syntax can be used:
#' cv = list('rand', nseg, nrep) for random repeated cross-validation with nseg segments and nrep
#' repetitions or cv = list('ven', nseg) for systematic splits to nseg segments ('venetian blinds').
#' @param resp
#' vector with response values to use in case of venetian blinds
#'
#' @return
#' matrix with object indices for each segment
#'
#' @export
crossval <- function(cv = 1, nobj = NULL, resp = NULL) {

  # get cross-validation parameters
  if (is.null(nobj)) nobj <- length(resp)

  p <- crossval.getParams(cv = cv, nobj = nobj)
  if (!(p$type %in% c("rand", "ven", "loo"))) {
    stop("Wrong name for cross-validation method.")
  }

  # check number of repetitions
  if (p$nrep < 1 || p$nrep > 100) {
    stop("Wrong value for cv repetitions (should be between 1 and 100).")
  }

  # check number of segments
  if (p$nseg < 2 || p$nseg > nobj) {
    stop("Wrong value for number of segments (should be between 2 and number of objects).")
  }

  seglen <- ceiling(nobj / p$nseg)
  fulllen <- seglen * p$nseg
  ind <- array(0, dim = c(p$nseg, seglen, p$nrep))

  if (p$type == "rand") {
    for (i in seq_len(p$nrep)) {
      v <- c(sample(nobj), rep(NA, fulllen - nobj))
      ind[, , i] <- matrix(v, nrow = p$nseg, byrow = TRUE)
    }
    return(ind)
  }

  if (p$type == "ven") {
    v <- c(order(resp), rep(NA, fulllen - nobj))
    ind[, , 1] <- matrix(v, nrow = p$nseg, byrow = FALSE)
    return(ind)
  }

  stop("Something went wrong.")
}

#' String with description of cross-validation method
#'
#' @param cv
#' a list with cross-validation settings
#'
#' @return
#' a string with the description text
#'
crossval.str <- function(cv) {

  if (length(cv) == 0) return("none")

  if (is.numeric(cv)) {
    return(
      if (cv == 1) "full (leave one out)"
      else sprintf("random with %.0f segments", cv)
    )
  }

  type <- cv[[1]]
  if (type == "loo") {
    return("full (leave one out)")
  }

  if (type == "ven") {
    return(sprintf("venetian blinds with %.0f segments", cv[[2]]))
  }

  return(
    sprintf("random with %.0f segments%s",
      cv[[2]], if (length(cv) == 3) paste(" and", cv[[3]], "repetitions") else "")
  )
}


#' Plot purity spectra
#' @param obj
#' object with mcr pure case
#' @param ...
#' other parameters
#'
#' @export
plotPuritySpectra <- function(obj, ...) {
  UseMethod("plotPuritySpectra")
}

#' Plot purity values
#' @param obj
#' object with mcr pure case
#' @param ...
#' other parameters
#'
#' @export
plotPurity <- function(obj, ...) {
  UseMethod("plotPurity")
}

#' Plot resolved spectra
#' @param obj
#' object with mcr case
#' @param ...
#' other parameters
#'
#' @export
plotSpectra <- function(obj, ...) {
  UseMethod("plotSpectra")
}

#' Plot resolved contributions
#' @param obj
#' object with mcr case
#' @param ...
#' other parameters
#'
#' @export
plotContributions <- function(obj, ...) {
  UseMethod("plotContributions")
}

#' Categorize PCA results
#' @param obj
#' object with PCA model
#' @param ...
#' other parameters
#'
#' @export
categorize <- function(obj, ...) {
  UseMethod("categorize")
}

#' Confusion matrix for classification results
#'
#' @details
#' Returns confusion matrix for classification results represented by the object.
#'
#' @param obj
#' classification results (object of class \code{simcares}, \code{simcamres}, etc)
#' @param ...
#' other parameters.
#'
#' @export
getConfusionMatrix <- function(obj, ...) {
  UseMethod("getConfusionMatrix")
}

#' Plot for class belonging probability
#'
#' @description
#' Makes a plot with class belonging probabilities for each object of the classification results.
#' Works only with classification methods, which compute this probability (e.g. SIMCA).
#'
#' @param obj
#' an object with classification results (e.g. SIMCA)
#' @param ...
#' other parameters
#'
#' @export
plotProbabilities <- function(obj, ...) {
  UseMethod("plotProbabilities")
}

#' Get class belonging probability
#'
#' @description
#' Compute class belonging probabilities for classification results.
#'
#' @param obj
#' an object with classification results (e.g. SIMCA)
#' @param ...
#' other parameters
#'
#' @export
getProbabilities <- function(obj, ...) {
  UseMethod("getProbabilities")
}

#' Set residual distance limits
#'
#' @description
#' Calculates and set critical limits for residuals of PCA model
#'
#' @param obj
#' a model object
#' @param ...
#' other parameters
#'
#' @export
setDistanceLimits <- function(obj, ...) {
  UseMethod("setDistanceLimits")
}

#' Show residual distance limits
#'
#' @description
#' Calculates and set critical limits for residuals of PCA model
#'
#' @param obj
#' a model object
#' @param ...
#' other parameters
#'
#' @export
showDistanceLimits <- function(obj, ...) {
  UseMethod("showDistanceLimits")
}

#' Shows extreme plot for SIMCA model
#'
#' @description
#' Generic function for creating extreme plot for SIMCA model
#'
#' @param obj
#' a SIMCA model
#' @param ...
#' other parameters
#'
#' @export
plotExtreme <- function(obj, ...) {
  UseMethod("plotExtreme")
}

#' Get regression coefficients
#'
#' @description
#' Generic function for getting regression coefficients from PLS model
#'
#' @param obj
#' a PLS model
#' @param ...
#' other parameters
#'
#' @export
getRegcoeffs <- function(obj, ...) {
  UseMethod("getRegcoeffs")
}

#' VIP scores plot
#'
#' @description
#' Generic function for plotting VIP scores values for regression model (PCR, PLS, etc)
#'
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#'
#' @export
plotVIPScores <- function(obj, ...) {
  UseMethod("plotVIPScores")
}

#' VIP scores
#'
#' @description
#' Generic function for returning VIP scores values for regression model (PCR, PLS, etc)
#'
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#'
#' @export
getVIPScores <- function(obj, ...) {
  UseMethod("getVIPScores")
}

#' Selectivity ratio plot
#'
#' @description
#' Generic function for plotting selectivity ratio values for regression model (PCR, PLS, etc)
#'
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#'
#' @export
plotSelectivityRatio <- function(obj, ...) {
  UseMethod("plotSelectivityRatio")
}

#' Selectivity ratio
#'
#' @description
#' Generic function for returning selectivity ratio values for regression model (PCR, PLS, etc)
#'
#' @param obj
#' a regression model
#' @param ...
#' other parameters
#'
#' @export
getSelectivityRatio <- function(obj, ...) {
  UseMethod("getSelectivityRatio")
}

#' Select optimal number of components for a model
#'
#' @description
#' Generic function for selecting number of components for multivariate models (e.g. PCA, PLS, ...)
#'
#' @param obj
#' a model object
#' @param ncomp
#' number of components to select
#' @param ...
#' other arguments
#'
#' @export
selectCompNum <- function(obj, ncomp = NULL, ...) {
  UseMethod("selectCompNum")
}

#' Cooman's plot
#'
#' @details
#' Generic function for Cooman's plot
#'
#' @param obj
#' classification model or result object
#' @param ...
#' other arguments
#'
#' @export
plotCooman <- function(obj, ...) {
  UseMethod("plotCooman")
}

#' Model distance plot
#'
#' @details
#' Generic function for plotting distance from object to a multivariate model
#'
#' @param obj
#' a model object
#' @param ...
#' other arguments
#'
#' @export
plotModelDistance <- function(obj, ...) {
  UseMethod("plotModelDistance")
}

#' Discrimination power plot
#'
#' @details
#' Generic function for plotting discrimination power values for classification model
#'
#' @param obj
#' a model object
#' @param ...
#' other arguments
#'
#' @export
plotDiscriminationPower <- function(obj, ...) {
  UseMethod("plotDiscriminationPower")
}

#' Calibration data
#'
#' @details
#' Generic function getting calibration data from a linear decomposition model (e.g. PCA)
#'
#' @param obj
#' a model object
#' @param ...
#' other arguments
#'
#' @export
getCalibrationData <- function(obj, ...) {
  UseMethod("getCalibrationData")
}

#' Modelling power plot
#'
#' @details
#' Generic function for plotting modelling power values for classification model
#'
#' @param obj
#' a model object
#' @param ...
#' other arguments
#'
#' @export
plotModellingPower <- function(obj, ...) {
  UseMethod("plotModellingPower")
}

#' Misclassification ratio plot
#'
#' @details
#' Generic function for plotting missclassification values for classification model or results
#'
#' @param obj
#' a model or a result object
#' @param ...
#' other arguments
#'
#' @export
plotMisclassified <- function(obj, ...) {
  UseMethod("plotMisclassified")
}

#' Specificity plot
#'
#' @details
#' Generic function for plotting specificity values for classification model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotSpecificity <- function(obj, ...) {
  UseMethod("plotSpecificity")
}

#' Sensitivity plot
#'
#' @details
#' Generic function for plotting sensitivity values for classification model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotSensitivity <- function(obj, ...) {
  UseMethod("plotSensitivity")
}

#' Classification performance plot
#'
#' @details
#' Generic function for plotting classification performance for model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotPerformance <- function(obj, ...) {
  UseMethod("plotPerformance")
}

#' Predictions
#'
#' @details
#' Generic function for showing predicted values for classification or regression model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
showPredictions <- function(obj, ...) {
  UseMethod("showPredictions")
}

#' X residuals plot
#'
#' @details
#' Generic function for plotting x residuals for classification or regression model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXResiduals <- function(obj, ...) {
  UseMethod("plotXResiduals")
}

#' Y residuals plot
#'
#' @details
#' Generic function for plotting y residuals for classification or regression model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotYResiduals <- function(obj, ...) {
  UseMethod("plotYResiduals")
}

#' X variance plot
#'
#' @details
#' Generic function for plotting explained variance for decomposition of x data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXVariance <- function(obj, ...) {
  UseMethod("plotXVariance")
}

#' Y variance plot
#'
#' @details
#' Generic function for plotting explained variance for decomposition of y data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotYVariance <- function(obj, ...) {
  UseMethod("plotYVariance")
}

#' Biplot
#'
#' @details
#' Generic function for biplot
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotBiplot <- function(obj, ...) {
  UseMethod("plotBiplot")
}

#' Scores plot
#'
#' @details
#' Generic function for scores values for data decomposition
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotScores <- function(obj, ...) {
  UseMethod("plotScores")
}

#' X scores plot
#'
#' @details
#' Generic function for plotting scores values for decomposition of x data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXScores <- function(obj, ...) {
  UseMethod("plotXScores")
}

#' XY scores plot
#'
#' @details
#' Generic function for plotting scores values for decomposition of x and y data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXYScores <- function(obj, ...) {
  UseMethod("plotXYScores")
}

#' Selected intervals plot
#'
#' @details
#' Generic function for plotting selected intervals or variables
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotSelection <- function(obj, ...) {
  UseMethod("plotSelection")
}

#' RMSE plot
#'
#' @details
#' Generic function for plotting RMSE values vs. complexity of a regression model
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotRMSE <- function(obj, ...) {
  UseMethod("plotRMSE")
}

#' Variance plot
#'
#' @details
#' Generic function for plotting explained variance for data decomposition
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotCumVariance <- function(obj, ...) {
  UseMethod("plotCumVariance")
}

#' X cumulative variance plot
#'
#' @details
#' Generic function for plotting cumulative explained variance for decomposition of x data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXCumVariance <- function(obj, ...) {
  UseMethod("plotXCumVariance")
}

#' Y cumulative variance plot
#'
#' @details
#' Generic function for plotting cumulative explained variance for decomposition of y data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotYCumVariance <- function(obj, ...) {
  UseMethod("plotYCumVariance")
}

#' Loadings plot
#'
#' @details
#' Generic function for plotting loadings values for data decomposition
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotLoadings <- function(obj, ...) {
  UseMethod("plotLoadings")
}

#' Predictions plot
#'
#' @details
#' Generic function for plotting predicted values for classification or regression model or results
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotPredictions <- function(obj, ...) {
  UseMethod("plotPredictions")
}

#' Regression coefficients plot
#'
#' @details
#' Generic function for plotting regression coefficients values for a regression model
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotRegcoeffs <- function(obj, ...) {
  UseMethod("plotRegcoeffs")
}

#' Residuals plot
#'
#' @details
#' Generic function for plotting residual values for data decomposition
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotResiduals <- function(obj, ...) {
  UseMethod("plotResiduals")
}

#' Variance plot
#'
#' @details
#' Generic function for plotting explained variance for data decomposition
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotVariance <- function(obj, ...) {
  UseMethod("plotVariance")
}

#' X loadings plot
#'
#' @details
#' Generic function for plotting loadings values for decomposition of x data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXLoadings <- function(obj, ...) {
  UseMethod("plotXLoadings")
}

#' X loadings plot
#'
#' @details
#' Generic function for plotting loadings values for decomposition of x and y data
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXYLoadings <- function(obj, ...) {
  UseMethod("plotXYLoadings")
}

#' Statistic histogram
#'
#' @details
#' Generic function for plotting statistic histogram plot
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotHist <- function(obj, ...) {
  UseMethod("plotHist")
}

#' Correlation plot
#'
#' @details
#' Generic function for correlation plot
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotCorr <- function(obj, ...) {
  UseMethod("plotCorr")
}

#' Plot for PLS weights
#'
#' @details
#' Generic function for weight plot
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotWeights <- function(obj, ...) {
  UseMethod("plotWeights")
}

#' Plot for XY-residuals
#'
#' @details
#' Generic function for XY-residuals plot
#'
#' @param obj
#' a model or result object
#' @param ...
#' other arguments
#'
#' @export
plotXYResiduals <- function(obj, ...) {
  UseMethod("plotXYResiduals")
}

#' Class for storing and visualising linear decomposition of dataset (X = TP' + E)
#'
#' @description
#' Creates an object of ldecomp class.
#'
#' @param scores
#' matrix with score values (I x A).
#' @param loadings
#' matrix with loading values (J x A).
#' @param residuals
#' matrix with data residuals (I x J)
#' @param eigenvals
#' vector with eigenvalues for the loadings
#' @param ncomp.selected
#' number of selected components
#'
#' @return
#' Returns an object (list) of \code{ldecomp} class with following fields:
#' \item{scores }{matrix with score values (I x A).}
#' \item{residuals }{matrix with data residuals (I x J).}
#' \item{T2 }{matrix with score distances (I x A).}
#' \item{Q }{matrix with orthogonal distances (I x A).}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#'
#' @details
#' \code{ldecomp} is a general class for storing results of decomposition of dataset in
#' form X = TP' + E. Here, X is a data matrix, T - matrix with scores, P - matrix with
#' loadings and E - matrix with residuals. It is used, for example, for PCA results
#' (\code{\link{pcares}}), in PLS and other methods. The class also includes methods for
#' calculation of residual distances and explained variance.
#'
#' There is no need to use the \code{ldecomp} manually. For example, when build PCA model
#' with \code{\link{pca}} or apply it to a new data, the results will automatically inherit
#' all methods of \code{ldecomp}.
#'
#' @importFrom methods show
#' @importFrom stats convolve cor lm na.exclude predict pt qf qnorm qt sd var
#'
#' @export
ldecomp <- function(scores, loadings, residuals, eigenvals, ncomp.selected = ncol(scores)) {

  ncomp <- ncol(scores)

  obj <- list(
    scores = scores,
    residuals = residuals,
    ncomp = ncomp,
    ncomp.selected = ncomp.selected,
    categories = NULL
  )

  # get distances and add them to the object
  dist <- ldecomp.getDistances(scores, loadings, residuals, eigenvals)
  obj <- c(obj, dist)

  # get variance and add it to the object
  var <- ldecomp.getVariances(scores, loadings, residuals, dist$Q)
  obj <- c(obj, var)

  obj$call <- match.call()
  class(obj) <- "ldecomp"

  return(obj)
}

#' Cumulative explained variance plot
#'
#' @description
#' Shows a plot with cumulative explained variance vs. number of components.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param labels
#' what to show as labels for plot objects
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotCumVariance.ldecomp <- function(obj, type = "b", labels = "values", show.plot = TRUE, ...) {

  return(
    plotVariance(obj, variance = "cumexpvar", type = type, labels = labels,
      show.plot = show.plot, ...)
  )
}

#' Explained variance plot
#'
#' @description
#' Shows a plot with explained variance vs. number of components.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param variance
#' string, which variance to make the plot for ("expvar", "cumexpvar")
#' @param labels
#' what to show as labels for plot objects.
#' @param xticks
#' vector with ticks for x-axis
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ylab
#' label for y-axis
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotVariance.ldecomp <- function(obj, type = "b", variance = "expvar", labels = "values",
  xticks = seq_len(obj$ncomp), show.plot = TRUE, ylab = "Explained variance, %", ...) {

  if (!show.plot) return(obj[[variance]])

  return(
    mdaplot(obj[[variance]], xticks = xticks, labels = labels, type = type, ylab = ylab, ...)
  )
}

#' Scores plot
#'
#' @description
#' Shows a plot with scores values for data objects.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param comp
#' which components to show the plot for (can be one value or vector with two values).
#' @param type
#' type of the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotScores.ldecomp <- function(obj, comp = c(1, 2), type = "p", show.axes = TRUE,
  show.plot = TRUE, ...) {

  # get scores for given components and generate column names with explained variance
  plot_data <- mda.subset(obj$scores, select = comp)
  colnames(plot_data) <- paste0("Comp ", comp, " (", round(obj$expvar[comp], 2), "%)")
  attr(plot_data, "name") <- "Scores"

  # if no plot required - return plot series object
  if (!show.plot) {
    return(plot_data)
  }

  # set up values for showing axes lines
  show.lines <- FALSE
  if (show.axes) {
    show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
  }

  # scatter plot
  if (type == "p") {
    p <- mdaplot(plot_data, type = type, show.lines = show.lines, ...)
    return(invisible(p))
  }

  # line or bar plot
  plot_data <- mda.t(plot_data)
  attr(plot_data, "yaxis.name") <- "Score"
  return(mdaplotg(plot_data, type = type, show.lines = show.lines, ...))
}

#' Residual distance plot
#'
#' @description
#' Shows a plot with orthogonal (Q, q) vs. score (T2, h) distances for data objects.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' number of components to show the plot for (if NULL, selected by model value will be used).
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels if necessary
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotResiduals.ldecomp <- function(obj, ncomp = obj$ncomp.selected, norm = FALSE, log = FALSE,
  show.labels = FALSE, labels = "names", show.plot = TRUE, ...) {

  attrs <- mda.getattr(obj$Q)

  # function for transforming distances
  transform <- function(u, u0, norm, log) {
    if (norm) u <- u / u0
    if (log) u <- log(1 + u)
    return(u)
  }

  # function for creating labels depending on transformation
  get_label <- function(lab, norm, log) {
    if (norm) lab <- paste0(lab, "/", lab, "0")
    if (log) lab <- paste0("log(1 + ", lab, ")")
    return(lab)
  }

  # get scale factors
  h0 <- if (!is.null(attr(obj$T2, "u0"))) attr(obj$T2, "u0")[[ncomp]]
  q0 <- if (!is.null(attr(obj$Q, "u0"))) attr(obj$Q, "u0")[[ncomp]]

  # check that scaling values exist
  if (norm && (is.null(h0) || is.null(q0))) {
    warning("Can not normalize distances as scaling values are absent.")
    norm <- FALSE
  }

  # prepare plot data
  h <- transform(obj$T2[, ncomp], h0, norm, log)
  q <- transform(obj$Q[, ncomp], q0, norm, log)

  # default values for local labels
  lxlab <- get_label("h", norm, log)
  lylab <- get_label("q", norm, log)

  # combine everything to dataset and assign attributes
  plot_data <- mda.cbind(h, q)
  plot_data <- mda.setattr(plot_data, attrs, "row")
  rownames(plot_data) <- rownames(obj$Q)
  colnames(plot_data) <- c(
    paste0("Score distance, ", lxlab),
    paste0("Orthogonal distance, ", lylab)
  )

  attr(plot_data, "name") <- sprintf("Distances (ncomp = %d)", ncomp)

  # if no plot required - return plot series object
  if (!show.plot) return(plot_data)

  # show plot
  return(mdaplot(plot_data, ...))
}

#' Print method for linear decomposition
#'
#' @description
#' Generic \code{print} function for linear decomposition. Prints information about
#' the \code{ldecomp} object.
#'
#' @param x
#' object of class \code{ldecomp}
#' @param str
#' user specified text to show as a description of the object
#' @param ...
#' other arguments
#'
#' @export
print.ldecomp <- function(x, str = NULL, ...) {
  if (is.null(str)) {
    str <- "Results of data decomposition (class ldecomp)."
  }

  if (nchar(str) > 0) {
    fprintf("\n%s\n", str)
  }

  cat("\nMajor fields:\n")
  cat("$scores - matrix with score values\n")
  cat("$T2 - matrix with T2 distances\n")
  cat("$Q - matrix with Q residuals\n")
  cat("$ncomp.selected - selected number of components\n")
  cat("$expvar - explained variance for each component\n")
  cat("$cumexpvar - cumulative explained variance\n")
}

#' as.matrix method for ldecomp object
#'
#' @description
#' Generic \code{as.matrix} function for linear decomposition. Returns a matrix with information
#' about the decomposition.
#'
#' @param x
#' object of class \code{ldecomp}
#' @param ncomp
#' number of components to get the result for (if NULL will return for each available)
#' @param ...
#' other arguments
#'
#' @export
as.matrix.ldecomp <- function(x, ncomp = NULL, ...) {

  out <- cbind(x$expvar, x$cumexpvar)
  rownames(out) <- colnames(x$Q)
  colnames(out) <- c("Expvar", "Cumexpvar")

  if (!is.null(ncomp)) {
    out <- out[ncomp, , drop = FALSE]
  }

  return(out)
}

#' Summary statistics for linear decomposition
#'
#' @description
#' Generic \code{summary} function for linear decomposition. Prints statistic about
#' the decomposition.
#'
#' @param object
#' object of class \code{ldecomp}
#' @param str
#' user specified text to show as a description of the object
#' @param ...
#' other arguments
#'
#' @export
summary.ldecomp <- function(object, str = NULL, ...) {
  if (is.null(str)) {
    str <- "Summary for data decomposition (class ldecomp)."
  }

  fprintf("\n%s\n", str)
  fprintf("\nSelected components: %d\n\n", object$ncomp.selected)

  print(round(as.matrix(object), 2))
}


##########################
# Static methods         #
##########################


#' Compute explained variance
#'
#' @description
#' Computes explained variance and cumulative explained variance for data decomposition.
#'
#' @param scores
#' matrix with scores (T).
#' @param loadings
#' matrix with loadings (P).
#' @param residuals
#' matrix with residuals (E).
#' @param Q
#' matrix with squared orthogonal distances.
#'
#' @return
#' Returns a list with two vectors.
#'
ldecomp.getVariances <- function(scores, loadings, residuals, Q) {

  # get names and attributes
  rows_excluded <- attr(scores, "exclrows")
  cols_excluded <- attr(scores, "exclcols")


  # remove excluded columns from loadings and residuals
  if (length(cols_excluded) > 0) {
    loadings <- loadings[-cols_excluded, , drop = FALSE]
    residuals <- residuals[-cols_excluded, , drop = FALSE]
  }

  # remove excluded rows from scores, residuals and Q
  if (length(rows_excluded) > 0) {
    scores <- scores[-rows_excluded, , drop = FALSE]
    residuals <- residuals[-rows_excluded, , drop = FALSE]
    Q <- Q[-rows_excluded, , drop = FALSE]
  }

  # compute total variance
  totvar <- sum(tcrossprod(scores, loadings)^2) + sum(residuals^2)

  # compute explained variance
  cumexpvar <- 100 * (1 - colSums(Q) / totvar)
  expvar <- c(cumexpvar[1], diff(cumexpvar))

  names(cumexpvar) <- names(expvar) <- colnames(Q)
  attr(expvar, "name") <- "Variance"
  attr(cumexpvar, "name") <- "Cumulative variance"
  attr(expvar, "xaxis.name") <- attr(cumexpvar, "xaxis.name") <- "Components"

  return(list(expvar = expvar, cumexpvar = cumexpvar))
}

#' Compute score and residual distances
#'
#' @description
#' Compute orthogonal Euclidean distance from object to PC space (Q, q) and Mahalanobis
#' squared distance between projection of the object to the space and its origin (T2, h).
#'
#' @param scores
#' matrix with scores (T).
#' @param loadings
#' matrix with loadings (P).
#' @param residuals
#' matrix with residuals (E).
#' @param eigenvals
#' vector with eigenvalues for the components
#'
#' @details
#' The distances are calculated for every 1:n components, where n goes from 1 to ncomp
#' (number of columns in scores and loadings).
#'
#' @return
#' Returns a list with Q, T2 and tnorm values for each component.
#'
ldecomp.getDistances <- function(scores, loadings, residuals, eigenvals) {

  # get names and attributes
  rows_excluded <- attr(scores, "exclrows")
  cols_excluded <- attr(loadings, "exclrows")

  # get sizes
  ncomp <- ncol(scores)
  nobj <- nrow(scores)

  # remove excluded variables from loadings and residuals
  if (length(cols_excluded) > 0) {
    loadings <- loadings[-cols_excluded, , drop = FALSE]
    residuals <- residuals[, -cols_excluded, drop = FALSE]
  }

  # get rid of hidden scores and residuals (needed for some calculations)
  scores_visible <- scores
  residuals_visible <- residuals
  if (length(rows_excluded) > 0) {
    scores_visible <- scores_visible[-rows_excluded, , drop = FALSE]
    residuals_visible <- residuals_visible[-rows_excluded, , drop = FALSE]
  }

  # normalize the scores
  scoresn <- scale(scores, center = FALSE, scale = sqrt(eigenvals))

  # prepare zero matrices for the and model power
  T2 <- matrix(0, nrow = nobj, ncol = ncomp)
  Q <- matrix(0, nrow = nobj, ncol = ncomp)

  # calculate distances and model power for each possible number of components in model
  for (i in seq_len(ncomp)) {
    res <- residuals
    if (i < ncomp) {
      res <- res +
        tcrossprod(
          scores[, (i + 1):ncomp, drop = F],
          loadings[, (i + 1):ncomp, drop = F]
        )
    }

    Q[, i] <- rowSums(res^2)
    T2[, i] <- rowSums(scoresn[, seq_len(i), drop = F]^2)
  }

  # set attributes for Q
  Q <- mda.setattr(Q, mda.getattr(scores), type = "row")
  attr(Q, "name") <- "Squared residual distance (q)"
  attr(Q, "xaxis.name") <- "Components"

  # set attributes for T2
  T2 <- mda.setattr(T2, mda.getattr(Q))
  attr(T2, "name") <- "Score distance (h)"

  colnames(Q) <- colnames(T2) <- colnames(loadings)
  rownames(Q) <- rownames(T2) <- rownames(scores)

  # return the results
  return(list(Q = Q, T2 = T2))
}

###############################
# Methods for critical limits #
###############################

#' Calculate critical limits for distance values using Jackson-Mudholkar approach
#'
#' @param residuals
#' matrix with PCA residuals
#' @param eigenvals
#' vector with eigenvalues for PCA components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @return
#' vector with four values: critical limits for given alpha and gamma, mean distance and DoF.
#'
#' @export
jm.crit <- function(residuals, eigenvals, alpha = 0.05, gamma = 0.01) {

  # if not all eigenvalues available - ise residuals to compute the rest
  ncomp <- length(eigenvals)
  nobj <- nrow(residuals)
  max_ncomp <- min(nrow(residuals) - 1, ncol(residuals))
  if (length(eigenvals) < max_ncomp) {
    eigenvals <- c(eigenvals, svd(residuals)$d[seq_len(max_ncomp - ncomp)]^2 / (nobj - 1))
  }

  # since it is residuals we do not need eigenvalue for PC1
  eigenvals <- eigenvals[-1]
  t1 <- rev(cumsum(rev(eigenvals)))[seq_len(ncomp)]
  t2 <- rev(cumsum(rev(eigenvals)^2))[seq_len(ncomp)]
  t3 <- rev(cumsum(rev(eigenvals)^3))[seq_len(ncomp)]

  h0 <- 1 - 2 * t1 * t3 / 3 / (t2^2);
  ifelse(h0 < 0.001, h0 <- 0.001, h0)

  # inverse error function
  erfinv <- function(x) qnorm((1 + x) / 2) / sqrt(2)
  gcl <- 1 - (1 - gamma) ^ (1 / nobj)
  ca <- sqrt(2) * erfinv(c(1 - 2 * alpha, (1 - 2 * gcl)))

  # compute h1 for alpha and gamma
  h1a <- ca[1] * sqrt(2 * t2 * h0^2) / t1
  h1g <- ca[2] * sqrt(2 * t2 * h0^2) / t1
  h2 <- t2 * h0 * (h0 - 1) / (t1 ^ 2)

  out <- rbind(
    t1 * (1 + h1a + h2) ^ (1 / h0),
    t1 * (1 + h1g + h2) ^ (1 / h0)
  )

  if (ncomp == max_ncomp) {
    out[, max_ncomp] <- 0
  }

  attr(out, "eigenvals") <- eigenvals
  return(out)
}

#' Calculate probabilities for distance values and given parameters using Hotelling T2 distribution
#'
#' @param u
#' vector with distances
#' @param eigenvals
#' vector with eigenvalues for PCA components
#' @param ncomp
#' number of components
#'
#' @export
jm.prob <- function(u, eigenvals, ncomp) {

  erf <- function(x) 1 - pnorm(-x * sqrt(2)) * 2

  t1 <- rev(cumsum(rev(eigenvals)))[ncomp]
  t2 <- rev(cumsum(rev(eigenvals)^2))[ncomp]
  t3 <- rev(cumsum(rev(eigenvals)^3))[ncomp]

  h0 <- 1 - 2 * t1 * t3 / 3 / (t2^2);
  ifelse(h0 < 0.001, h0 <- 0.001, h0)

  h1 <- (u / t1)^h0
  h2 <- t2 * h0 * (h0 - 1) / t1^2
  d <- t1 * (h1 - 1 - h2) / (sqrt(2 * t2) * h0)

  return(0.5 * (1 + erf(d / sqrt(2))))
}

#' Calculate critical limits for distance values using Hotelling T2 distribution
#'
#' @param nobj
#' number of objects in calibration set
#' @param ncomp
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @return
#' vector with four values: critical limits for given alpha and gamma, mean distance and DoF.
#'
#' @export
hotelling.crit <- function(nobj, ncomp, alpha = 0.05, gamma = 0.01) {
  return(
    rbind(
      (ncomp * (nobj - 1) / (nobj - ncomp)) * qf(1 - alpha, ncomp, (nobj - ncomp)),
      (ncomp * (nobj - 1) / (nobj - ncomp)) * qf((1 - gamma) ^ (1 / nobj), ncomp, (nobj - ncomp))
    )
  )
}

#' Calculate probabilities for distance values and given parameters using Hotelling T2 distribution
#'
#' @param u
#' vector with distances
#' @param ncomp
#' number of components
#' @param nobj
#' number of objects in calibration set
#'
#' @export
hotelling.prob <- function(u, ncomp, nobj) {
  return(pf(u * (nobj - ncomp) / (ncomp * (nobj - 1)), ncomp, (nobj - ncomp)))
}

#' Calculates critical limits for distance values using Chi-square distribution
#'
#' @description
#' The method is based on Chi-squared distribution with DF = 2 * (m(u)/s(u)^2
#'
#' @param param
#' matrix with distribution parameters
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
chisq.crit <- function(param, alpha = 0.05, gamma = 0.01) {

  u0 <- param$u0
  nobj <- param$nobj
  Nu <- param$Nu

  DoF <- floor(Nu)
  DoF <- ifelse(DoF < 1, 1, DoF)

  return(
    rbind(
      qchisq(1 - alpha, DoF) * u0 / Nu,
      qchisq((1 - gamma) ^ (1 / nobj), DoF) * u0 / Nu
    )
  )
}

#' Calculate probabilities for distance values using Chi-square distribution
#'
#' @param u
#' vector with distances
#' @param param
#' vector with distribution parameters
#'
#' @export
chisq.prob <- function(u, param) {
  u0 <- param[1]
  Nu <- param[2]

  DoF <- floor(Nu)
  DoF[DoF == 0] <- 1
  return(pchisq(Nu * u / u0, DoF))
}

#' Calculates critical limits for distance values using Data Driven moments approach
#'
#' @param paramQ
#' matrix with parameters for distribution of Q distances
#' @param paramT2
#' matrix with parameters for distribution of T2 distances
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
dd.crit <- function(paramQ, paramT2, alpha = 0.05, gamma = 0.01) {

  nobj <- paramQ$nobj
  Nq <- round(paramQ$Nu)
  Nq[Nq < 1] <- 1
  Nq[Nq > 250] <- 250

  Nh <- round(paramT2$Nu)
  Nh[Nh < 1] <- 1
  Nh[Nh > 250] <- 250

  return(
    rbind(
      qchisq(1 - alpha, Nq + Nh),
      qchisq((1 - gamma) ^ (1 / nobj), Nq + Nh)
    )
  )
}

#' Calculates critical limits for distance values using Data Driven moments approach
#'
#' @param U
#' matrix or vector with distance values
#'
#' @export
ddmoments.param <- function(U) {

  if (is.null(dim(U))) dim(U) <- c(length(U), 1)

  u0 <- apply(U, 2, mean)
  su <- apply(U, 2, sd)
  Nu <- 2 * (u0 / su)^2

  return(list(u0 = u0, Nu = Nu, nobj = nrow(U)))
}

#' Calculates critical limits for distance values using Data Driven robust approach
#'
#' @param U
#' matrix or vector with distance values
#' @param ncomp
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
ddrobust.param <- function(U, ncomp, alpha, gamma) {

  if (is.null(dim(U))) dim(U) <- c(length(U), 1)

  Mu <- apply(U, 2, median)
  Su <- apply(U, 2, IQR)

  RM <- Su / Mu
  Nu <- round(exp((1.380948 * log(2.68631 / RM)) ^ 1.185785))
  Nu[RM > 2.685592117] <- 1
  Nu[RM < 0.194565995] <- 100

  u0 <- 0.5 * Nu * (Mu / qchisq(0.50, Nu) + Su / (qchisq(0.75, Nu) - qchisq(0.25, Nu)))
  return(list(u0 = u0, Nu = Nu, nobj = nrow(U)))
}

#' Compute parameters for critical limits based on calibration results
#'
#' @param U
#' matrix with residual distances
#'
#' @export
ldecomp.getLimParams <- function(U) {

  U <- mda.purgeRows(U)

  return(
    list(
      "moments" = ddmoments.param(U),
      "robust" = ddrobust.param(U),
      "nobj" = nrow(U)
    )
  )
}

#' Compute critical limits for orthogonal distances (Q)
#'
#' @param lim.type
#' which method to use for calculation of critical limits for residuals
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#' @param params
#' distribution parameters returned by ldecomp.getLimParams
#' @param residuals
#' matrix with residuals (E)
#' @param eigenvals
#' egenvalues for the components used to decompose the data
#'
#' @export
ldecomp.getQLimits <- function(lim.type, alpha, gamma, params, residuals, eigenvals) {

  pQ <- if (regexpr("robust", lim.type) > 0) params$Q$robust else params$Q$moments
  ncomp <- length(pQ$u0)

  if (lim.type == "jm") {
    # methods based on Jackson-Mudholkar approach
    residuals <- mda.purge(residuals)
    lim <- jm.crit(residuals, eigenvals, alpha, gamma)
    eigenvals <- attr(lim, "eigenvals")
    lim <- rbind(lim, pQ$u0, nrow(residuals))
    attr(lim, "eigenvals") <- eigenvals

  } else {
    # methods based on chi-square distribution
    pT2 <- if (regexpr("robust", lim.type) > 0) params$T2$robust else params$T2$moments
    DoF <- round(pQ$Nu)
    DoF[DoF < 1] <- 1
    DoF[DoF > 250] <- 250

    lim <- switch(lim.type,
      "chisq" = chisq.crit(pQ, alpha, gamma),
      "ddmoments" = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pQ$u0),
      "ddrobust"  = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pQ$u0),
      stop("Wrong value for 'lim.type' parameter.")
    )

    lim <- rbind(lim, pQ$u0, DoF)
  }

  colnames(lim) <- paste("Comp", seq_len(ncomp))
  rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
  attr(lim, "name") <- "Critical limits for orthogonal distances (Q)"
  attr(lim, "alpha") <- alpha
  attr(lim, "gamma") <- gamma
  attr(lim, "lim.type") <- lim.type

  return(lim)
}

#' Compute critical limits for score distances (T2)
#'
#' @param lim.type
#' which method to use for calculation ("chisq", "ddmoments", "ddrobust")
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#' @param params
#' distribution parameters returned by ldecomp.getLimParams
#'
#' @export
ldecomp.getT2Limits <- function(lim.type, alpha, gamma, params) {

  pQ <- if (regexpr("robust", lim.type) > 0) params$Q$robust else params$Q$moments
  pT2 <- if (regexpr("robust", lim.type) > 0) params$T2$robust else params$T2$moments
  ncomp <- length(pT2$u0)

  DoF <- round(pT2$Nu)
  DoF[DoF < 1] <- 1
  DoF[DoF > 250] <- 250

  if (lim.type %in% c("jm", "chisq")) DoF <- pT2$nobj

  lim <- switch(lim.type,
    "jm" = hotelling.crit(pT2$nobj, seq_len(ncomp), alpha, gamma),
    "chisq" = hotelling.crit(pT2$nobj, seq_len(ncomp), alpha, gamma),
    "ddmoments" = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pT2$u0),
    "ddrobust"  = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pT2$u0),
    stop("Wrong value for 'lim.type' parameter.")
  )

  lim <- rbind(lim, pT2$u0, DoF)
  colnames(lim) <- paste("Comp", seq_len(ncomp))
  rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
  attr(lim, "name") <- "Critical limits for score distances (T2)"
  attr(lim, "alpha") <- alpha
  attr(lim, "gamma") <- gamma
  attr(lim, "lim.type") <- lim.type

  return(lim)
}

#' Compute coordinates of lines or curves with critical limits
#'
#' @param Qlim
#' matrix with critical limits for orthogonal distances
#' @param T2lim
#' matrix with critical limits for score distances
#' @param ncomp
#' number of components for computing the coordinates
#' @param norm
#' logical, shall distance values be normalized or not
#' @param log
#' logical, shall log transformation be applied or not
#' @param show.limits
#' vector with two logical values defining if limits for extreme and/or outliers must be shown
#'
#' @return
#' list with two matrices (x and y coordinates of corresponding limits)
#'
#' @export
ldecomp.getLimitsCoordinates <- function(Qlim, T2lim, ncomp, norm, log,
  show.limits = c(TRUE, TRUE)) {

  # get parameters
  h0 <- T2lim[3, ncomp]
  q0 <- Qlim[3, ncomp]
  Nh <- T2lim[4, ncomp]
  Nq <- Qlim[4, ncomp]
  lim.type <- attr(Qlim, "lim.type")

  # check that show.limits is logical
  if (!all(is.logical(show.limits))) {
    stop("Parameter 'show.limits' must have logical value(s).")
  }

  # if show.limits has only one value - duplicate it
  if (length(show.limits) == 1) {
    show.limits <- rep(show.limits, 2)
  }

  # compute the limits
  if (lim.type %in% c("jm", "chisq")) {

    # quadratic limits
    hE <- c(0, T2lim[1, ncomp], T2lim[1, ncomp])
    hO <- c(0, T2lim[2, ncomp], T2lim[2, ncomp])

    qE <- c(Qlim[1, ncomp], Qlim[1, ncomp], 0)
    qO <- c(Qlim[2, ncomp], Qlim[2, ncomp], 0)

  } else {

    ## slope and intercepts
    eB <- Qlim[1, ncomp]
    oB <- Qlim[2, ncomp]
    eA <- oA <- -1 * (q0 / h0) * (Nh / Nq)

    hE <- seq(-0.95, -eB / eA, length.out = 100)
    hO <- seq(-0.95, -oB / oA, length.out = 100)
    qE <- eA * hE + eB
    qO <- oA * hO + oB
  }

  if (norm) {
    hE <- hE / h0
    qE <- qE / q0
    hO <- hO / h0
    qO <- qO / q0
  }

  if (log) {
    hE <- log(1 + hE)
    qE <- log(1 + qE)
    hO <- log(1 + hO)
    qO <- log(1 + qO)
  }

  return(list(
    extremes = if (show.limits[[1]]) cbind(hE, qE),
    outliers = if (show.limits[[2]]) cbind(hO, qO)
  ))
}

#' Residuals distance plot for a set of ldecomp objects
#'
#' @description
#' Shows a plot with score (T2, h) vs orthogonal (Q, q) distances and corresponding critical
#' limits for given number of components.
#'
#' @param res
#' list with result objects to show the plot for
#' @param Qlim
#' matrix with critical limits for orthogonal distance
#' @param T2lim
#' matrix with critical limits for score distance
#' @param ncomp
#' how many components to use (by default optimal value selected for the model will be used)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param cgroup
#' color grouping of plot points (works only if one result object is available)
#' @param xlim
#' limits for x-axis (if NULL will be computed automatically)
#' @param ylim
#' limits for y-axis (if NULL will be computed automatically)
#' @param show.legend
#' logical, show or not a legend on the plot (needed if several result objects are available)
#' @param show.limits
#' vector with two logical values defining if limits for extreme and/or outliers must be shown
#' @param lim.col
#' vector with two values - line color for extreme and outlier limits
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier limits
#' @param lim.lty
#' vector with two values - line type for extreme and outlier limits
#' @param show.legend
#' logical, show or not legend on the plot (if more than one result object)
#' @param legend.position
#' if legend must be shown, where it should be
#' @param show.excluded
#' logical, show or hide rows marked as excluded (attribute `exclrows`).
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The function is a bit more advanced version of \code{\link{plotResiduals.ldecomp}}. It allows to
#' show distance values for several result objects (e.g. calibration and test set or calibration
#' and new prediction set) as well as display the correspondng critical limits in form of lines
#' or curves.
#'
#' Depending on how many result objects your model has or how many you specified manually,
#' using the \code{res} parameter, the plot behaves in a bit different way.
#'
#' If only one result object is provided, then it allows to colorise the points using \code{cgroup}
#' parameter. If two or more result objects are provided, then the function show
#' distances in groups, and adds corresponding legend.
#'
#' The function can show distance values normalised (h/h0 and q/q0) as well as with log
#' transformation (log(1 + h/h0), log(1 + q/q0)). The latter is useful if distribution of the
#' points is skewed and most of them are densely located around bottom left corner.
#'
#' @export
ldecomp.plotResiduals <- function(res, Qlim, T2lim, ncomp, log = FALSE, norm = FALSE,
  cgroup = NULL, xlim = NULL, ylim = NULL, show.limits = c(TRUE, TRUE),
  lim.col = c("darkgray", "darkgray"), lim.lwd = c(1, 1), lim.lty = c(2, 3),
  show.legend = TRUE, legend.position = "topright", show.excluded = FALSE, ...) {

  # return column with values either with or without excluded outliers
  getValues <- function(x, dim) {
    return(if (show.excluded) x[, dim] else mda.purgeRows(x)[, dim])
  }

  # compute limits fo axis depending on values and position of critical limits
  getPlotLim <- function(lim, pd, ld, dim) {
    if (!is.null(lim) || all(!show.limits)) return(lim)
    limits <- if (show.limits[[2]]) max(ld$outliers[, dim]) else max(ld$extremes[, dim])
    return(c(0, max(sapply(pd, function(x) max(c(getValues(x, dim), limits)) * 1.05))))
  }

  # check that show.limits is logical
  if (!all(is.logical(show.limits))) {
    stop("Parameter 'show.limits' must have logical value(s).")
  }

  # if show.limits has only one value - duplicate it
  if (length(show.limits) == 1) {
    show.limits <- rep(show.limits, 2)
  }

  # keep only ojects of class "ldecomp" in result list
  res <- getRes(res, "ldecomp")

  # compute plot data for each result object
  plot_data <- lapply(res, plotResiduals.ldecomp, ncomp = ncomp, norm = norm, log = log,
    show.plot = FALSE)

  # get coordinates for critical limits
  lim_data <- ldecomp.getLimitsCoordinates(Qlim, T2lim, ncomp = ncomp, norm = norm, log = log)
  xlim <- getPlotLim(xlim, plot_data, lim_data, 1)
  ylim <- getPlotLim(ylim, plot_data, lim_data, 2)

  # make plot
  if (length(plot_data) == 1) {
    mdaplot(plot_data[[1]], type = "p", xlim = xlim, ylim = ylim, cgroup = cgroup,
      show.excluded = show.excluded, ...)
  } else {
    mdaplotg(plot_data, type = "p", xlim = xlim, ylim = ylim, show.legend = show.legend,
      show.excluded = show.excluded, legend.position = legend.position, ...)
  }

  # show critical limits
  if (show.limits[[1]]) {
    lines(lim_data$extremes[, 1], lim_data$extremes[, 2],
      col = lim.col[1], lty = lim.lty[1], lwd = lim.lwd[1])
  }

  if (show.limits[[2]]) {
    lines(lim_data$outliers[, 1], lim_data$outliers[, 2],
      col = lim.col[2], lty = lim.lty[2], lwd = lim.lwd[2])
  }
}

#' PLS results
#'
#' @description
#' \code{plsres} is used to store and visualize results of applying a PLS model to a new data.
#'
#' @param y.pred
#' predicted y values.
#' @param y.ref
#' reference (measured) y values.
#' @param ncomp.selected
#' selected (optimal) number of components.
#' @param xdecomp
#' PLS decomposition of X data (object of class \code{ldecomp}).
#' @param ydecomp
#' PLS decomposition of Y data (object of class \code{ldecomp}).
#' @param info
#' information about the object.
#'
#' @details
#' Do not use \code{plsres} manually, the object is created automatically when one applies a PLS
#' model to a new data set, e.g. when calibrate and validate a PLS model (all calibration and
#' validation results in PLS model are stored as objects of \code{plsres} class) or use function
#' \code{\link{predict.pls}}.
#'
#' The object gives access to all PLS results as well as to the plotting methods for visualisation
#' of the results. The \code{plsres} class also inherits all properties and methods of \code{regres}
#'  - general class for regression results.
#'
#' If no reference values provided, regression statistics will not be calculated and most of the
#' plots not available. The class is also used for cross-validation results, in this case some of
#' the values and methods are not available (e.g. scores and scores plot, etc.).
#'
#' All plots are based on \code{\link{mdaplot}} function, so most of its options can be used (e.g.
#' color grouping, etc.).
#'
#' RPD is ratio of standard deviation of response values to standard error of prediction (SDy/SEP).
#'
#' @return
#' Returns an object of \code{plsres} class with following fields:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{y.ref }{a matrix with reference values for responses.}
#' \item{y.pred }{a matrix with predicted values for responses.}
#' \item{rmse }{a matrix with root mean squared error values for each response and component.}
#' \item{slope }{a matrix with slope values for each response and component.}
#' \item{r2 }{a matrix with determination coefficients for each response and component.}
#' \item{bias }{a matrix with bias values for each response and component.}
#' \item{sep }{a matrix with standard error values for each response and component.}
#' \item{rpd }{a matrix with RPD values for each response and component.}
#' \item{xdecomp }{decomposition of predictors (object of class \code{ldecomp}).}
#' \item{ydecomp }{decomposition of responses (object of class \code{ldecomp}).}
#' \item{info }{information about the object.}
#'
#' @seealso
#'
#' Methods for \code{plsres} objects:
#' \tabular{ll}{
#'    \code{print} \tab prints information about a \code{plsres} object.\cr
#'    \code{\link{summary.plsres}} \tab shows performance statistics for the results.\cr
#'    \code{\link{plot.plsres}} \tab shows plot overview of the results.\cr
#'    \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
#'    \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
#'    \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
#'    \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
#'    \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y
#'    decomposition.\cr
#'    \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y
#'    decomposition.\cr
#'    \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'    \code{\link{plotYResiduals.plsres}} \tab shows residuals plot for y values.\cr
#' }
#'
#' Methods inherited from \code{regres} class (parent class for \code{plsres}):
#' \tabular{ll}{
#'    \code{\link{plotPredictions.regres}} \tab shows predicted vs. measured plot.\cr
#'    \code{\link{plotRMSE.regres}} \tab shows RMSE plot.\cr
#' }
#'
#' See also \code{\link{pls}} - a class for PLS models.
#'
#'
#' @export
plsres <- function(y.pred, y.ref = NULL, ncomp.selected = dim(y.pred)[2], xdecomp = NULL,
  ydecomp = NULL, info = "") {

  obj <- regres(y.pred, y.ref = y.ref, ncomp.selected = ncomp.selected)
  obj$ncomp <- dim(y.pred)[2]
  obj$xdecomp <- xdecomp
  obj$ydecomp <- ydecomp
  obj$info <- info
  obj$ncomp.selected <- ncomp.selected

  obj$call <- match.call()
  class(obj) <- c("plsres", "regres")

  return(obj)
}

#' as.matrix method for PLS results
#'
#' @description
#' Returns a matrix with model performance statistics for PLS results
#'
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' number of components to calculate the statistics for
#' @param ny
#' for which response variable calculate the statistics for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.plsres <- function(x, ncomp = NULL, ny = 1, ...) {

  xdecomp_res <- if (!is.null(x$xdecomp)) as.matrix(x$xdecomp) else matrix(NA, x$ncomp, 2)
  ydecomp_res <- if (!is.null(x$ydecomp)) as.matrix(x$ydecomp) else matrix(NA, x$ncomp, 2)

  out <- cbind(
    xdecomp_res,
    ydecomp_res,
    as.matrix.regres(x, ny = ny)
  )

  rownames(out) <- paste("Comp", seq_len(x$ncomp))
  colnames(out)[1:4] <- c("X expvar", "X cumexpvar", "Y expvar", "Y cumexpvar")

  if (!is.null(ncomp)) {
    out <- out[ncomp, , drop = FALSE]
  }

  return(out)
}

#' summary method for PLS results object
#'
#' @description
#' Shows performance statistics for the results.
#'
#' @param object
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' for which response variable show the summary for
#' @param ...
#' other arguments
#'
#' @export
summary.plsres <- function(object, ny = seq_len(object$nresp), ncomp = NULL, ...) {

  cat("\nPLS regression results (class plsres) summary\n")
  fprintf("Info: %s\n", object$info)
  fprintf("Number of selected components: %d\n", object$ncomp.selected)

  if (is.null(object$y.ref)) {
    cat("No reference data provided to calculate prediction performance.")
    return()
  }

  if (length(ncomp) == 1) {
    fprintf("\nNumber of selected components: %d\n", ncomp)
  }

  for (y in ny) {
    fprintf("\nResponse variable %s:\n", dimnames(object$y.pred)[[3]][y])
    out <- as.matrix.plsres(object, ny = y, ncomp = ncomp)
    if (!any(is.na(out[, 1:4]))) out[, 1:4] <- round(out[, 1:4], 3)
    out[, 5] <- round(out[, 5], 3)
    out[, 6] <- mdaplot.formatValues(out[, 6], round.only = T)
    out[, 7] <- round(out[, 7], 3)
    out[, 8] <- round(out[, 8], 4)
    out[, 9] <- round(out[, 9], 2)
    print(out)
  }
}

#' print method for PLS results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ...
#' other arguments
#'
#' @export
print.plsres <- function(x, ...) {
  cat("\nPLS results (class plsres)\n")
  cat("\nCall:\n")
  print(x$call)

  cat("\nMajor fields:\n")
  cat("$ncomp.selected - number of selected components\n")
  cat("$y.pred - array with predicted y values\n")

  if (!is.null(x$y.ref)) {
    cat("$y.ref - matrix with reference y values\n")
    cat("$rmse - root mean squared error\n")
    cat("$r2 - coefficient of determination\n")
    cat("$slope - slope for predicted vs. measured values\n")
    cat("$bias - bias for prediction vs. measured values\n")
    cat("$ydecomp - decomposition of y values (ldecomp object)\n")
  }

  cat("$xdecomp - decomposition of x values (ldecomp object)\n")
}

################################
#  Plotting methods            #
################################

#' Explained X variance plot for PLS results
#'
#' @description
#' Shows plot with explained X variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param decomp
#' which dcomposition to use ("xdecomp" or "ydecomp")
#' @param variance
#' which variance to use ("expvar", "cumexpvar")
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotVariance.plsres <- function(obj, decomp = "xdecomp", variance = "expvar", ...) {
  if (is.null(obj[[decomp]])) return(NULL)
  return(plotVariance.ldecomp(obj[[decomp]], variance = variance, ...))
}

#' Explained X variance plot for PLS results
#'
#' @description
#' Shows plot with explained X variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXVariance.plsres <- function(obj, main = "Variance (X)", ...) {
  return(plotVariance.plsres(obj, decomp = "xdecomp", main = main, ...))
}

#' Explained Y variance plot for PLS results
#'
#' @description
#' Shows plot with explained Y variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotYVariance.plsres <- function(obj, main = "Variance (Y)", ...) {
  return(plotVariance.plsres(obj, decomp = "ydecomp", main = main, ...))
}

#' Explained cumulative X variance plot for PLS results
#'
#' @description
#' Shows plot with cumulative explained X variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXCumVariance.plsres <- function(obj, main = "Cumulative variance (X)", ...) {
  return(plotVariance.plsres(obj, decomp = "xdecomp", variance = "cumexpvar", main = main, ...))
}

#' Explained cumulative Y variance plot for PLS results
#'
#' @description
#' Shows plot with cumulative explained Y variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotYCumVariance.plsres <- function(obj, main = "Cumulative variance (Y)", ...) {
  return(plotVariance.plsres(obj, decomp = "ydecomp", variance = "cumexpvar", main = main, ...))
}

#' X scores plot for PLS results
#'
#' @description
#' Shows plot with scores values for PLS decomposition of x data.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXScores.plsres <- function(obj, comp = c(1, 2), main = "Scores (X)", ...) {
  if (is.null(obj$xdecomp)) return(invisible(NULL))
  return(plotScores.ldecomp(obj$xdecomp, comp = comp, main = main, ...))
}

#' XY scores plot for PLS results
#'
#' @description
#' Shows plot with X vs. Y scores values for PLS results.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' which component to show the plot for
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXYScores.plsres <- function(obj, ncomp = 1, show.plot = TRUE, ...) {

  if (is.null(obj$xdecomp) || is.null(obj$ydecomp)) return(invisible(NULL))

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
    stop("Wrong value for ncomp argument.")
  }

  plot_data <- cbind(
    obj$xdecomp$scores[, ncomp, drop = FALSE],
    obj$ydecomp$scores[, ncomp, drop = FALSE]
  )

  plot_data <- mda.setattr(plot_data, mda.getattr(obj$xdecomp$scores))
  rownames(plot_data) <- rownames(obj$xdecomp$scores)
  colnames(plot_data) <- c(
    sprintf("X-scores (Comp %d, %.2f%%)", ncomp, obj$xdecomp$expvar[ncomp]),
    sprintf("Y-scores (Comp %d, %.2f%%)", ncomp, obj$ydecomp$expvar[ncomp])
  )

  attr(plot_data, "name") <- "Scores (XY)"

  if (!show.plot) {
    return(plot_data)
  }

  return(mdaplot(plot_data, type = "p", ...))
}

#' X residuals plot for PLS results
#'
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for PLS decomposition of x data.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param main
#' main title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXResiduals.plsres <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
  main = sprintf("X-distances (ncomp = %d)", ncomp), ...) {

  if (is.null(obj$xdecomp)) return(invisible(NULL))

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
    stop("Wrong value for ncomp argument.")
  }

  return(plotResiduals.ldecomp(obj$xdecomp, ncomp = ncomp, main = main,
    norm = norm, log = log, ...))
}

#' Y residuals plot for PLS results
#'
#' @description
#' Shows a plot with Y residuals vs reference Y values for selected component.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' Proxy for \code{\link{plotResiduals.regres}} function.
#'
#' @export
plotYResiduals.plsres <- function(obj, ncomp = obj$ncomp.selected, ...) {

  if (is.null(obj$y.ref)) return(invisible(NULL))
  return(plotResiduals.regres(obj, ncomp = ncomp, ...))
}


#' Overview plot for PLS results
#'
#' @description
#' Shows a set of plots for PLS results.
#'
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plot.plsres <- function(x, ncomp = x$ncomp.selected, ny = 1, show.labels = FALSE, ...) {

  if (is.null(x$y.ref)) {
    par(mfrow = c(1, 2))
    plotXResiduals(x, ...)
    plotPredictions.regres(x, ncomp = ncomp, ny = ny, ...)
    par(mfrow = c(1, 1))
    return()
  }

  par(mfrow = c(2, 2))
  plotXResiduals(x, ncomp = ncomp, ...)
  plotYVariance(x, ...)
  plotRMSE(x, ny = ny, ...)
  plotPredictions.regres(x, ncomp = ncomp, ny = ny, ...)
  par(mfrow = c(1, 1))
}

#' Residual distance plot
#'
#' @description
#' Shows a plot with orthogonal (Q, q) vs. score (T2, h) distances for data objects.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' number of components to show the plot for (if NULL, selected by model value will be used).
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels if necessary
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotXYResiduals.plsres <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
  show.labels = FALSE, labels = "names", show.plot = TRUE, ...) {

  attrs <- mda.getattr(obj$xdecomp$Q)

  # function for transforming distances
  transform <- function(u, u0, norm, log) {
    if (norm) u <- u / u0
    if (log) u <- log(1 + u)
    return(u)
  }

  # function for creating labels depending on transformation
  get_label <- function(lab, norm, log) {
    if (norm) lab <- paste0(lab, "/", lab, "0")
    if (log) lab <- paste0("log(1 + ", lab, ")")
    return(lab)
  }

  # get scale factors
  h0 <- if (!is.null(attr(obj$xdecomp$T2, "u0"))) attr(obj$xdecomp$T2, "u0")[[ncomp]]
  q0 <- if (!is.null(attr(obj$xdecomp$Q, "u0"))) attr(obj$xdecomp$Q, "u0")[[ncomp]]
  z0 <- if (!is.null(attr(obj$ydecomp$Q, "u0"))) attr(obj$ydecomp$Q, "u0")[[ncomp]]

  # get DoF factors
  Nh <- if (!is.null(attr(obj$xdecomp$T2, "Nu"))) attr(obj$xdecomp$T2, "Nu")[[ncomp]]
  Nq <- if (!is.null(attr(obj$xdecomp$Q, "Nu"))) attr(obj$xdecomp$Q, "Nu")[[ncomp]]

  # get distances
  h <- obj$xdecomp$T2[, ncomp]
  q <- obj$xdecomp$Q[, ncomp]
  z <- obj$ydecomp$Q[, ncomp]

  # compute full distance for X
  f <- Nh * h / h0 + Nq * q / q0
  f0 <- Nh + Nq

  # prepare plot data
  f <- transform(f, f0, norm, log)
  z <- transform(z, z0, norm, log)

  # default values for local labels
  lxlab <- get_label("f", norm, log)
  lylab <- get_label("z", norm, log)

  # combine everything to dataset and assign attributes
  plot_data <- mda.cbind(f, z)
  plot_data <- mda.setattr(plot_data, attrs, "row")

  rownames(plot_data) <- rownames(obj$xdecomp$Q)
  colnames(plot_data) <- c(
    paste0("Full X-distance, ", lxlab),
    paste0("Y-distance, ", lylab)
  )

  attr(plot_data, "name") <- sprintf("XY-distances (ncomp = %d)", ncomp)

  # if no plot required - return plot series object
  if (!show.plot) return(plot_data)

  # show plot
  return(mdaplot(plot_data, ...))
}

#' Regression coefficients
#'
#' @description
#' class for storing and visualisation of regression coefficients
#' for regression models
#'
#' @param coeffs
#' array (npred x ncomp x nresp) with regression coefficients
#' @param ci.coeffs
#' array (npred x ncomp x nresp x cv) with regression coefficients for
#' computing confidence intervals (e.g. from cross-validation) using Jack-Knifing method
#' @param use.mean
#' logical, tells how to compute standard error for regression coefficients. If \code{TRUE}
#' mean values for ci.coeffs is computed first. If \code{FALSE}, \code{values} (coefficients
#' computed for global model) are used as mean.
#'
#' @return
#' a list (object of \code{regcoeffs} class) with fields, including:
#' \tabular{ll}{
#'    \code{values} \tab an array (nvar x ncomp x ny) with regression coefficients \cr
#'    \code{se} \tab an array (nvar x ncomp x ny) with standard errors for the coefficients \cr
#'    \code{t.values} \tab an array (nvar x ncomp x ny) with t-values for the coefficients \cr
#'    \code{p.values} \tab an array (nvar x ncomp x ny) with p-values for coefficients \cr
#' }
#'
#' last three fields are available if parameter \code{ci.coeffs} was provided.
#'
#' Check also \code{\link{confint.regcoeffs}}, \code{\link{summary.regcoeffs}} and
#' \code{\link{plot.regcoeffs}}.
#'
#' @export
regcoeffs <- function(coeffs, ci.coeffs = NULL, use.mean = TRUE) {

  if (is.null(dim(coeffs)) || length(dim(coeffs)) != 3) {
    stop("Coefficients must be provided as 3-way array.")
  }

  obj <- list()
  obj$values <- coeffs
  obj$nvar <- dim(coeffs)[1]

  # assign response number and names
  obj$nresp <- dim(coeffs)[3]
  obj$respnames <- dimnames(coeffs)[[3]]
  if (is.null(obj$respnames)) obj$respnames <- paste0("y", seq_len(obj$nresp))

  # add statistics and class name
  obj <- c(obj, regcoeffs.getStats(coeffs, ci.coeffs, use.mean))
  obj$call <- match.call()
  class(obj) <- "regcoeffs"

  return(obj)
}

#' Confidence intervals for regression coefficients
#'
#' @description
#' returns matrix with confidence intervals for regression coeffocoents
#' for given response number and number of components.
#'
#' @param object
#' regression coefficients object (class \code{regcoeffs})
#' @param parm
#' not used, needed for compatiility with general method
#' @param level
#' confidence level
#' @param ncomp
#' number of components (one value)
#' @param ny
#' index of response variable (one value)
#' @param ...
#' other arguments
#'
#' @export
confint.regcoeffs <- function(object, parm = NULL, level = 0.95, ncomp = 1, ny = 1, ...) {

  if (length(ncomp) != 1) {
    stop("Parameter 'ncomp' should be just one value.")
  }

  if (ncomp < 1 || ncomp > dim(object$values)[2]) {
    stop("Wrong value for parameter 'ncomp'.")
  }

  if (length(ny) != 1) {
    stop("Parameter 'ny' should be just one value.")
  }

  if (ny < 1 || ny > dim(object$values)[3]) {
    stop("Wrong value for parameter 'ny'.")
  }

  alpha <- 1 - level
  t <- -qt(alpha / 2, object$DoF)
  ci <- repmat(object$se[, ncomp, ny], 1, 2) %*% diag(c(-t, t)) +
    repmat(object$values[, ncomp, ny], 1, 2)

  if (length(attr(object$values, "exclrows"))) {
    ci[attr(object$values, "exclrows"), ] <- 0
  }

  ci <- mda.setattr(ci, mda.getattr(object$values))
  rownames(ci) <- dimnames(object$values)[[1]]
  colnames(ci) <- paste0(round(c(alpha / 2, alpha / 2 + level) * 100, 1), "%")
  attr(ci, "name") <- "Confidence interval"

  return(ci)
}

#' as.matrix method for regression coefficients class
#'
#' @description
#' returns matrix with regression coeffocoents for given response number and amount of components
#'
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' number of response variable to return the coefficients for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.regcoeffs <- function(x, ncomp = 1, ny = 1, ...) {
  return(x$values[, ncomp, ny])
}

#' Summary method for regcoeffs object
#'
#' @description
#' Shows estimated coefficients and statistics (if available).
#'
#' @param object
#' object of class \code{regcoeffs}
#' @param ncomp
#' how many components to use
#' @param ny
#' which y variable to show the summary for
#' @param alpha
#' significance level for confidence interval (if statistics available)
#' @param ...
#' other arguments
#'
#' @details
#' Statistcs are shown if Jack-Knifing was used when model is calibrated.
#'
#' @export
summary.regcoeffs <- function(object, ncomp = 1, ny = 1, alpha = 0.05, ...) {

  if (length(ncomp) != 1) {
    stop("Parameter 'ncomp' should be just one value.")
  }

  if (ncomp < 1 || ncomp > dim(object$values)[2]) {
    stop("Wrong value for parameter 'ncomp'.")
  }

  if (length(ny) != 1) {
    stop("Parameter 'ny' should be just one value.")
  }

  if (ny < 1 || ny > dim(object$values)[3]) {
    stop("Wrong value for parameter 'ny'.")
  }

  attrs <- mda.getattr(object$values)
  coeffs <- object$values[, ncomp, ny, drop = F]
  dim(coeffs) <- c(dim(object$values)[1], 1)
  colnames(coeffs)[1] <- "Coeffs"
  if (!is.null(object$se)) {
    coeffs <- cbind(
      coeffs,
      object$se[, ncomp, ny],
      round(object$t.values[, ncomp, ny], 2),
      round(object$p.values[, ncomp, ny], 3),
      confint(object, ncomp = ncomp, ny = ny, level = 1 - alpha)
    )
    colnames(coeffs)[1:4] <- c("Coeffs", "Std. err.", "t-value", "p-value")
  }

  rownames(coeffs) <- dimnames(object$values)[[1]]
  attr(coeffs, "exclrows") <- attrs$exclrows
  attr(coeffs, "name") <- paste0(
    "Regression coefficients for ", object$respnames[ny], " (ncomp = ", ncomp, ")"
  )

  cat("\n")
  mda.show(coeffs)
  if (ncol(coeffs) > 1) {
    cat(sprintf("\nDegrees of freedom (Jack-Knifing): %d\n", object$DoF))
  }
  cat("\n")
}

#' print method for regression coefficients class
#'
#' @description
#' prints regression coeffocoent values for given response number and amount of components
#'
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ...
#' other arguments
#'
#' @export
print.regcoeffs <- function(x, ...) {
  cat("\nRegression coefficients (class regcoeffs)\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nMajor fields:\n")
  cat("$values - array with regression coefficients\n")
  cat("$se - array with standard errors\n")
  cat("$t.values - array with t-values\n")
  cat("$p.values - array with p-values\n")
  cat("$DoF - degrees of freedom for Jack-Knifing\n")
  cat("\nThe last four fields available only if Jack-Knifing was used.\n\n")
}


################################
#  Static methods              #
################################


#' Distribution statistics for regression coeffificents
#'
#' @description
#' calculates standard error, t-values and p-values for
#' regression coefficients based on Jack-Knifing method.
#'
#' @param coeffs
#' array (npred x ncomp x nresp) with regression coefficients
#' @param ci.coeffs
#' array (npred x ncomp x nresp x cv) with regression coefficients for
#' computing confidence intervals (e.g. from cross-validation) using Jack-Knifing method
#' @param use.mean
#' logical, tells how to compute standard error for regression coefficients. If \code{TRUE}
#' mean values for ci.coeffs is computed first. If \code{FALSE}, \code{values} (coefficients
#' computed for global model) are used as mean.
#'
#' @return
#' a list with statistics three arrays: srandard error, t-values and p-values computed for
#' each regression coefficient.
#'
#' @export
regcoeffs.getStats <- function(coeffs, ci.coeffs = NULL, use.mean = TRUE) {

  if (is.null(ci.coeffs)) return()

  if (is.null(dim(ci.coeffs)) || length(dim(ci.coeffs)) != 4) {
    stop("Coefficients for distribution statistics must be provided as 4-way array.")
  }

  # get attributes and prepare arrays
  nseg <- dim(ci.coeffs)[4]
  dim_names <- dimnames(coeffs)
  DoF <- dim(ci.coeffs)[4] - 1
  attrs <- mda.getattr(coeffs)
  se <- p.values <- t.values <- array(0, dim = dim(coeffs))

  # to make sure p-values for excluded predictors are 1 (not important)
  p.values <- p.values + 1

  # prepare correct indices for predictors taking into accound excluded variables
  row_ind <- seq_len(dim(coeffs)[1])
  if (length(attrs$exclrows) > 0) {
    row_ind <- row_ind[-attrs$exclrows]
    coeffs <- coeffs[-attrs$exclrows, , , drop = FALSE]
  }

  # check the dimension
  if (any(dim(coeffs) != dim(ci.coeffs)[1:3])) {
    stop("Dimension of coefficients for distribution statistics does not much the 'coeffs'.")
  }

  # compute mean if needed
  if (use.mean) {
    coeffs <- apply(ci.coeffs, 1:3, mean)
  }

  # compute main statistics
  err <- sweep(ci.coeffs, 1:3, coeffs, "-")
  ssq <- apply(err^2, 1:3, sum)
  se[row_ind, , ] <- sqrt((nseg - 1) / nseg * ssq)
  t.values[row_ind, , ] <- coeffs / se[row_ind, , , drop = FALSE]
  p.values[row_ind, , ] <- 2 * (1 - pt(abs(t.values[row_ind, , , drop = FALSE]), nseg - 1))

  # add names and attributes
  dimnames(se) <- dimnames(t.values) <- dimnames(p.values) <- dim_names
  t.values <- mda.setattr(t.values, attrs)
  p.values <- mda.setattr(p.values, attrs)
  se <- mda.setattr(se, attrs)
  attr(se, "name") <- "Standard error (Jack-knifing)"
  attr(t.values, "name") <- "t-values (Jack-knifing)"
  attr(p.values, "name") <- "p-values (Jack-knifing)"

  return(
    list(
      se = se,
      t.values = t.values,
      p.values = p.values,
      DoF = DoF
    )
  )
}


################################
#  Plotting methods            #
################################


#' Regression coefficients plot
#'
#' @description
#' Shows plot with regression coefficient values for every predictor variable (x)
#'
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to use for creating the plot
#' @param ny
#' index of response variable to make the plot for
#' @param type
#' type of the plot
#' @param col
#' vector with two colors for the plot (one is used to show real coefficient and another one to
#' show confidence intervals)
#' @param show.lines
#' allows to show horizontal line at c(NA, 0)
#' @param show.ci
#' logical, show or not confidence intervals if they are available
#' @param alpha
#' significance level for confidence intervals (a number between 0 and 1, e.g.
#' for 95\% alpha = 0.05)
#' @param ylab
#' label for y-axis
#' @param ...
#' other arguments for plotting methods (e.g. main, xlab, etc)
#'
#' @export
plot.regcoeffs <- function(x, ncomp = 1, ny = 1, type = (if (x$nvar > 30) "l" else "h"),
  col = c(mdaplot.getColors(1), "lightgray"), show.lines = c(NA, 0), show.ci = FALSE,
  alpha = 0.05, ylab = paste0("Coefficients (", x$respnames[ny], ")"), ...) {

  if (length(ncomp) != 1) {
    stop("Parameter 'ncomp' should be just one value.")
  }

  if (ncomp < 1 || ncomp > dim(x$values)[2]) {
    stop("Wrong value for parameter 'ncomp'.")
  }

  if (length(ny) != 1) {
    stop("Parameter 'ny' should be just one value.")
  }

  if (ny < 1 || ny > dim(x$values)[3]) {
    stop("Wrong value for parameter 'ny'.")
  }

  plot_data <- matrix(x$values[, ncomp, ny], nrow = 1)
  attr(plot_data, "exclcols") <- attr(x$values, "exclrows")
  attr(plot_data, "xaxis.name") <- attr(x$values, "yaxis.name")
  attr(plot_data, "xaxis.values") <- attr(x$values, "yaxis.values")

  colnames(plot_data) <- rownames(x$values)
  attr(plot_data, "name") <- paste0("Regression coefficients (ncomp = ", ncomp, ")")

  if (!(show.ci && !is.null(x$se))) {
    return(mdaplot(plot_data, col = col[1], type = type, show.lines = show.lines, ylab = ylab, ...))
  }

  ci <- confint(x, ncomp = ncomp, ny = ny, level = 1 - alpha)
  if (type == "l") {
    return(
      mdaplot(mda.rbind(plot_data, mda.t(ci)), type = "l", col = col[c(2, 1, 1)],
        show.lines = show.lines, ylab = ylab, ...)
    )
  }

  err <- (ci[, 2] - ci[, 1]) / 2
  mdaplotg(list(plot_data, mda.rbind(plot_data, err)), type = c(type, "e"), show.legend = FALSE,
    col = col[c(2, 1)], show.grid = T, show.lines = show.lines, ylab = ylab, ...)
}

regmodel <- function(...) {
}

#' Cross-validation of a regression model
#'
#' @description
#' Does cross-validation of a regression model
#'
#' @param obj
#' a regression model (object of class \code{regmodel})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param y
#' a matrix with y values (responses from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param cal.fun
#' reference to function for model calibration
#'
#' @return
#' object of class \code{plsres} with results of cross-validation
#'
#' @export
crossval.regmodel <- function(obj, x, y, cv, cal.fun) {

  # get attributes
  x.attrs <- attributes(x)
  y.attrs <- attributes(y)

  # remove excluded rows
  if (length(x.attrs$exclrows) > 0) {
    x <- x[-x.attrs$exclrows, , drop = FALSE]
    y <- y[-x.attrs$exclrows, , drop = FALSE]
    attr(x, "exclrows") <- NULL
    attr(y, "exclrows") <- NULL
  }

  # remove excluded columns
  if (length(x.attrs$exclcols) > 0) {
    x <- x[, -x.attrs$exclcols, drop = FALSE]
    attr(x, "exclcols") <- NULL
  }

  if (length(y.attrs$exclcols) > 0) {
    y <- y[, -y.attrs$exclcols, drop = FALSE]
    attr(y, "exclcols") <- NULL
  }

  y.ref <- y

  # get main data parameters
  nvar <- ncol(x)
  nobj <- nrow(x)
  nresp <- ncol(y)
  ncomp <- obj$ncomp

  # get matrix with indices for cv segments
  cv_ind <- crossval(cv, nobj = nobj, resp = y[, 1])
  nseg <- nrow(cv_ind);
  nrep <- dim(cv_ind)[3]

  # prepare arrays for results
  yp.cv <- array(0, dim = c(nobj, ncomp, nresp))
  jk.coeffs <- array(0, dim = c(nvar, ncomp, nresp, nseg))

  # loop over segments and repetitions
  for (ir in seq_len(nrep)) {
    for (is in seq_len(nseg)) {
      ind <- na.exclude(cv_ind[is, , ir])
      if (length(ind) == 0) next

      xc <- x[-ind, , drop = FALSE]
      yc <- y[-ind, , drop = FALSE]
      xt <- x[ind, , drop = FALSE]

      # create a model
      m.loc <- cal.fun(xc, yc, ncomp, method = obj$method, center = obj$center,
        scale = obj$scale, cv = TRUE)

      if (m.loc$ncomp < ncomp) {
        stop(
          "Local model inside cross-validation can not be computed with the same number of\n",
          "components as used for calibration. Limit the number by using parameter 'ncomp'\n",
          "and run the code again.\n", call. = FALSE
        )
      }

      r.loc <- predict(m.loc, xt, cv = TRUE)

      # if any have NA values quit
      if (any(is.na(r.loc$y.pred))) {
        stop("NA results produced during cross-validation.")
      }

      # save results
      yp.cv[ind, , ] <- yp.cv[ind, , , drop = FALSE] + r.loc$y.pred
      jk.coeffs[, , , is] <- jk.coeffs[, , , is, drop = FALSE] +
        array(m.loc$coeffs$values, dim = c(dim(m.loc$coeffs$values), 1))
    }
  }

  # average results over repetitions
  yp.cv <- yp.cv / nrep
  jk.coeffs <- jk.coeffs / nrep

  # set up names
  dimnames(jk.coeffs) <- list(
    colnames(x),
    colnames(obj$coeffs$values),
    colnames(y),
    seq_len(nseg)
  )

  dimnames(yp.cv) <- list(
    rownames(x),
    colnames(obj$coeffs$values),
    colnames(y)
  )

  # make pls results and return
  return(list(y.pred = yp.cv, y.ref = y.ref, jk.coeffs = jk.coeffs))
}

#' Regression coefficients for PLS model'
#'
#' @description
#' Returns a matrix with regression coefficients for
#' the PLS model which can be applied to a data directly
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' if y is multivariate which variables you want to see the coefficients for
#' @param full
#' if TRUE the method also shows p-values and t-values as well as confidence intervals for the
#' coefficients (if available)
#' @param alpha
#' significance level for confidence intervals (a number between 0 and 1, e.g. 0.05)
#' @param ...
#' other parameters
#'
#' @details
#' The method recalculates the regression coefficients found by the PLS algorithm
#' taking into account centering and scaling of predictors and responses, so the
#' matrix with coefficients can be applied directly to original data (yp = Xb).
#'
#' If number of components is not specified, the optimal number, selected by user
#' or identified by a model will be used.
#'
#' If Jack-knifing method was used to get statistics for the coefficient the method
#' returns all statistics as well (p-value, t-value, confidence interval). In this case user
#' has to specified a number of y-variable (if there are many) to get the statistics and
#' the coefficients for. The confidence interval is computed for unstandardized coefficients.
#'
#' @return
#' A matrix  with regression coefficients and (optinally) statistics.
#'
#' @export
getRegcoeffs.regmodel <- function(obj, ncomp = obj$ncomp.selected, ny = 1, full = FALSE,
  alpha = 0.05, ...) {

  if (length(ncomp) != 1 || ncomp <= 0 || ncomp > obj$ncomp) {
    stop("Wrong value for number of components.")
  }

  attrs <- mda.getattr(obj$coeffs$values)
  out <- obj$coeffs$values[, ncomp, ny]

  # get center values and scale factors
  sx <- if (is.logical(obj$xscale)) rep(1, nrow(out)) else obj$xscale
  mx <- if (is.logical(obj$xcenter)) rep(0, nrow(out)) else obj$xcenter
  sy <- if (is.logical(obj$yscale)) rep(1, ncol(out)) else obj$yscale[ny]
  my <- if (is.logical(obj$ycenter)) rep(0, ncol(out)) else obj$ycenter[ny]


  # rescale coefficients and find intercept
  s <- sy / sx
  out <- matrix(c(my - sum(s * out * mx), s * out), ncol = 1)
  colnames(out) <- "Estimated"
  rownames(out) <- c("Intercept", dimnames(obj$coeffs$values)[[1]])


  if (full && !is.null(obj$coeffs$se)) {
    ci <- rbind(c(NA, NA), confint(obj$coeffs, level = 1 - alpha))
    se <- c(NA, obj$coeffs$se[, ncomp, ny] * s)
    t  <- c(NA, obj$coeffs$t.values[, ncomp, ny])
    p  <- c(NA, obj$coeffs$p.values[, ncomp, ny])
    out <- cbind(out, se, t, p, ci)
    colnames(out)[2:6] <- c("Std. err.", "t-value", "p-value",
      paste0(round(c(alpha / 2, 1 - alpha / 2) * 100, 2), "%"))
  }

  attr(out, "exclrows") <- attrs$exclrows
  attr(out, "name") <- paste("Regression coefficients for ", obj$coeffs$respnames[ny])

  return(out)
}

#' Summary method for regression model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a regression model (object of class \code{regmodel})
#' @param ncomp
#' number of components to show summary for
#' @param ny
#' which y variables to show the summary for (can be a vector)
#' @param res
#' list of results to show summary for
#' @param ...
#' other arguments
#'
#' @export
summary.regmodel <- function(object, ncomp = object$ncomp.selected,
  ny = seq_len(object$res$cal$nresp), res = object$res, ...) {

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > dim(object$res[["cal"]]$y.pred)[2]) {
    stop("Wrong value for 'ncomp' parameter.")
  }

  cat("\nRegression model (class regmodel) summary\n")
  cat("\nPerformance and validation:\n")
  fprintf("Number of selected components: %d\n", ncomp)

  for (y in ny) {
    sum_data <- do.call(rbind, lapply(res, as.matrix, ny = y, ncomp = ncomp))
    rownames(sum_data) <- capitalize(names(res))

    sum_data[, "R2"] <- round(sum_data[, "R2"], 3)
    sum_data[, "RMSE"] <- mdaplot.formatValues(sum_data[, "RMSE"], round.only = T)
    sum_data[, "Slope"] <- round(sum_data[, "Slope"], 3)
    sum_data[, "Bias"] <- round(sum_data[, "Bias"], 4)
    sum_data[, "RPD"] <- round(sum_data[, "RPD"], 1)

    attr(sum_data, "name") <- sprintf("\nResponse variable #%d (%s)", y, res[[1]]$respnames[y])
    mda.show(sum_data)
  }
  cat("\n")
}

#' Print method for PLS model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a regression model (object of class \code{regmodel})
#' @param ...
#' other arguments
#'
#' @export
print.regmodel <- function(x, ...) {
  cat("\nRegression model (class regmodel)\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nMajor fields:\n")
  cat("$ncomp - number of calculated components\n")
  cat("$ncomp.selected - number of selected components\n")
  cat("$coeffs - object (regcoeffs) with regression coefficients\n")
  cat("$res - list with result objects\n")
  cat("\nTry summary(model) and plot(model) to see the model performance.\n")
}


################################
#  Plotting methods            #
################################


#' RMSE plot for regression model
#'
#' @description
#' Shows plot with root mean squared error values vs. number of components for PLS model.
#'
#' @param obj
#' a regression model (object of class \code{regmodel})
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param type
#' type of the plot("b", "l" or "h")
#' @param labels
#' what to show as labels (vector or name, e.g. "names", "values", "indices")
#' @param xticks
#' vector with ticks for x-axis values
#' @param res
#' list with result objects
#' @param ylab
#' label for y-axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @export
plotRMSE.regmodel <- function(obj, ny = 1, type = "b", labels = "values",
  xticks = seq_len(obj$ncomp), res = obj$res, ylab = paste0("RMSE (", obj$res$cal$respnames[ny], ")"), ...) {

  plot_data <- lapply(res, plotRMSE, ny = ny, show.plot = FALSE)
  mdaplotg(plot_data, type = type, xticks = xticks, labels = labels, ylab = ylab, ...)
}

#' Predictions plot for regression model
#'
#' @description
#' Shows plot with predicted vs. reference (measured) y values for selected components.
#'
#' @param obj
#' a regression model (object of class \code{regmodel})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param legend.position
#' position of legend on the plot (if shown)
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param res
#' list with result objects
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @export
plotPredictions.regmodel <- function(obj, ncomp = obj$ncomp.selected, ny = 1,
  legend.position = "topleft", show.line = TRUE, res = obj$res, ...) {

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > dim(obj$res[["cal"]]$y.pred)[2]) {
    stop("Wrong value for 'ncomp' parameter.")
  }

  plot_data <- lapply(res, plotPredictions.regres, ny = ny, ncomp = ncomp, show.plot = FALSE)
  attr(plot_data[[1]], "name") <- sprintf("Predictions (ncomp = %d)", ncomp)
  plots <- mdaplotg(plot_data, type = "p", legend.position = legend.position, ...)

  if (show.line) {
    for (p in plots) {
      plotRegressionLine(p)
    }
  }
}

#' Y residuals plot for regression model
#'
#' @description
#' Shows plot with y residuals (predicted vs. reference values) for selected components.
#'
#' @param obj
#' a regression model (object of class \code{regmodel})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param show.lines
#' allows to show the horizonta line at 0 level
#' @param res
#' list with result objects
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @export
plotYResiduals.regmodel <- function(obj, ncomp = obj$ncomp.selected, ny = 1, show.lines = c(NA, 0),
  res = obj$res, ...) {

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > dim(obj$res[["cal"]]$y.pred)[2]) {
    stop("Wrong value for 'ncomp' parameter.")
  }

  plot_data <- lapply(res, plotResiduals, ny = ny, ncomp = ncomp, show.plot = FALSE)
  attr(plot_data[[1]], "name") <- sprintf("Y-residuals (ncomp = %d)", ncomp)
  mdaplotg(plot_data, show.lines = show.lines, ...)
}

#' Regression coefficient plot for regression model
#'
#' @description
#' Shows plot with regression coefficient values. Is a proxy for \code{link{plot.regcoeffs}} method.
#'
#' @param obj
#' a regression model (object of class \code{regmodel})
#' @param ncomp
#' number of components to show the plot for
#' @param ...
#' other plot parameters (see \code{link{plot.regcoeffs}} for details)
#'
#' @export
plotRegcoeffs.regmodel <- function(obj, ncomp = obj$ncomp.selected, ...) {
  plot(obj$coeffs, ncomp = ncomp, ...)
}

#' Regression results
#'
#' @description
#' Class for storing and visualisation of regression predictions
#'
#' @param y.pred
#' vector or matrix with y predicted values
#' @param y.ref
#' vector with reference (measured) y values
#' @param ncomp.selected
#' if y.pred calculated for different components, which to use as default
#'
#' @return
#' a list (object of \code{regres} class) with fields, including:
#' \tabular{ll}{
#'    \code{y.pred} \tab a matrix with predicted values \cr
#'    \code{y.ref} \tab a vector with reference (measured) values \cr
#'    \code{ncomp.selected} \tab selected column/number of components for predictions \cr
#'    \code{rmse} \tab root mean squared error for predicted vs measured values \cr
#'    \code{slope} \tab slope for predicted vs measured values \cr
#'    \code{r2} \tab coefficient of determination for predicted vs measured values \cr
#'    \code{bias} \tab bias for predicted vs measured values \cr
#'    \code{rpd} \tab RPD values \cr
#' }
#'
#' @export
regres <- function(y.pred, y.ref = NULL, ncomp.selected = 1) {

  if (is.null(y.pred) || length(dim(y.pred)) != 3) {
    stop("Parameter 'y.pred' should be a 3-way array.")
  }

  if (ncomp.selected > dim(y.pred)[2]) {
    stop("Wrong value for 'ncomp.selected' parameter.")
  }

  if (!is.null(y.ref)) y.ref <- as.matrix(y.ref)

  obj <- list()
  obj$y.pred <- y.pred
  obj$y.ref <- y.ref
  obj$ncomp <- dim(y.pred)[2]
  obj$ncomp.selected <- ncomp.selected
  obj$nresp <- dim(y.pred)[3]
  obj$respnames <- dimnames(y.pred)[[3]]
  if (is.null(obj$respnames)) obj$respnames <- paste0("y", seq_len(obj$nresp))

  obj <- c(obj, regres.getPerformanceStats(y.pred, y.ref))
  obj$call <- match.call()
  class(obj) <- "regres"

  return(obj)
}

#' as.matrix method for regression results
#'
#' @description
#' Returns a matrix with model performance statistics for regression results
#'
#' @param x
#' regression results (object of class \code{regres})
#' @param ncomp
#' model complexity (number of components) to calculate the statistics for (can be a vector)
#' @param ny
#' for which response variable calculate the statistics for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.regres <- function(x, ncomp = NULL, ny = 1, ...) {
  if (is.null(x$y.ref)) return()

  out <- cbind(x$r2[ny, ], x$rmse[ny, ], x$slope[ny, ],
    x$bias[ny, ], x$rpd[ny, ])

  colnames(out) <- c("R2", "RMSE", "Slope", "Bias", "RPD")
  rownames(out) <- dimnames(x$y.pred)[[2]]

  if (!is.null(ncomp)) {
    out <- out[ncomp, , drop = FALSE]
  }

  return(out)
}

#' summary method for regression results object
#'
#' @description
#' Shows performance statistics for the regression results.
#'
#' @param object
#' regression results (object of class \code{regres})
#' @param ncomp
#' model complexity to show the summary for (if NULL - shows for all available values)
#' @param ny
#' for which response variable show the summary for (one value or a vector)
#' @param ...
#' other arguments
#'
#' @export
summary.regres <- function(object, ncomp = object$ncomp.selected, ny = seq_len(object$nresp), ...) {

  cat("\nRegression results (class regres) summary\n")
  if (is.null(object$y.ref)) {
    cat("No reference data provided to calculate prediction performance.")
    return()
  }

  if (!is.null(ncomp)) {
    fprintf("\nNumber of selected components: %d\n\n", ncomp)
  }

  for (i in ny) {
    fprintf("\nResponse variable %s:\n", object$respnames[i])
    print(as.matrix.regres(object, ny = i, ncomp = ncomp))
  }
}

#' print method for regression results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' regression results (object of class \code{regres})
#' @param ...
#' other arguments
#'
#' @export
print.regres <- function(x, ...) {
  cat("\nRegression results (class regres)\n")
  cat("\nCall:\n")
  print(x$call)

  cat("\nMajor fields:\n")
  cat("$y.pred - matrix or vector with predicted y values\n")
  if (!is.null(x$y.ref)) {
    cat("$y.ref - vector with reference y values\n")
    cat("$rmse - root mean squared error\n")
    cat("$r2 - coefficient of determination\n")
    cat("$slope - slope for predicted vs. measured values\n")
    cat("$bias - bias for prediction vs. measured values\n")
  }

  if (ncol(x$y.pred) > 1) {
    cat("$ncomp.selected - number of selected components\n")
  }
}


################################
#  Static methods              #
################################


regres.getPerformanceStats <- function(y.pred, y.ref) {
  if (is.null(y.ref)) return(NULL)

  attrs <- mda.getattr(y.pred)

  # remove excluded rows so they are not counted
  # when calculating statistics
  if (length(attrs$exclrows) > 0) {
    y.pred <- y.pred[-attrs$exclrows, , , drop = F]
    y.ref <- y.ref[-attrs$exclrows, , drop = F]
  }

  # residuals (errors) based statistics
  err <- regres.err(y.pred, y.ref)
  ytot <- colSums(scale(y.ref, center = TRUE, scale = FALSE)^2)

  stats <- list(
    "r2" = regres.r2(err, ytot),
    "bias" = regres.bias(err),
    "rmse" = regres.rmse(err)
  )

  stats$slope <- regress.addattrs(regres.slope(y.pred, y.ref), attributes(err), "Slope")
  stats$sep <- regress.addattrs(sqrt(stats$rmse^2 - stats$bias^2), attributes(err), "SEP")
  stats$rpd <- regress.addattrs(apply(y.ref, 2, sd) / stats$sep, attributes(err), "RPD")

  return(stats)
}

#' Add names and attributes to matrix with statistics
#'
#' @param stat
#' matrix with statistics
#' @param attrs
#' attributes from error matrix
#' @param name
#' name of statistic
#'
regress.addattrs <- function(stat, attrs, name) {

  attr(stat, "name") <- name
  dimnames(stat) <- attrs$dimnames[c(3, 2)]
  attr(stat, "xaxis.name") <- attrs$xaxis.name
  attr(stat, "yaxis.name") <- attrs$yaxis.name

  return(stat)
}

#' Error of prediction
#'
#' @description
#' Calculates array of differences between predicted and reference values.
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.err <- function(y.pred, y.ref) {
  ncomp <- dim(y.pred)[2]
  err <- array(repmat(y.ref, ncomp, 1), dim = dim(y.pred)) - y.pred

  attr(err, "name") <- "Error of prediction"
  attr(err, "xaxis.name") <- "Components"
  attr(err, "yaxis.name") <- "Predictors"

  return(err)
}

#' Determination coefficient
#'
#' @description
#' Calculates matrix with coeffient of determination for every response and components
#'
#' @param err
#' vector with difference between reference and predicted y-values
#' @param ytot
#' total variance for y-values
#'
#' @export
regres.r2 <- function(err, ytot) {
  r2 <- t(1 - scale(colSums(err^2), center = F, scale = ytot))
  return(regress.addattrs(r2, attributes(err), "Coefficient of determination"))
}

#' Prediction bias
#'
#' @description
#' Calculates matrix with bias (average prediction error) for every response and components
#'
#' @param err
#' vector with difference between reference and predicted y-values
#'
regres.bias <- function(err) {
  bias <- t(colSums(err) / nrow(err))
  return(regress.addattrs(bias, attributes(err), "Bias"))
}

#' RMSE
#'
#' @description
#' Calculates matrix with root mean squared error of prediction for every response and components.
#'
#' @param err
#' vector with difference between reference and predicted y-values
#'
regres.rmse <- function(err) {
  rmse <- t(sqrt(colSums(err^2) / nrow(err)))
  return(regress.addattrs(rmse, attributes(err), "RMSE"))
}

#' Slope
#'
#' @description
#' Calculates matrix with slope of predicted and measured values for every response and components.
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.slope <- function(y.pred, y.ref) {
  nresp <- ncol(y.ref)
  ncomp <- ncol(y.pred)
  slope <- matrix(0, nrow = nresp, ncol = ncomp)
  for (i in seq_len(nresp)) {
    slope[i, ] <- matrix(coefficients(lm(y.pred[, , i] ~ y.ref[, i])), nrow = 2)[2, ]
  }

  return(slope)
}


################################
#  Plotting methods            #
################################


#' Predictions plot for regression results
#'
#' @description
#' Shows plot with predicted y values.
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param ncomp
#' complexity of model (e.g. number of components) to show the plot for
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param show.stat
#' logical, show or not legend with statistics on the plot
#' @param stat.col
#' color of text in legend with statistics
#' @param stat.cex
#' size of text in legend with statistics
#' @param xlim
#' limits for x-axis (if NULL will be computed automatically)
#' @param ylim
#' limits for y-axis (if NULL will be computed automatically)
#' @param axes.equal
#' logical, make limits for x and y axes equal or not
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' If reference values are available, the function shows a scatter plot with predicted vs.
#' reference values, otherwise predicted values are shown vs. object numbers.
#'
#' @export
plotPredictions.regres <- function(obj, ny = 1, ncomp = obj$ncomp.selected, show.line = TRUE,
  show.stat = FALSE, stat.col = "#606060", stat.cex = 0.85, xlim = NULL, ylim = NULL,
  axes.equal = TRUE, show.plot = TRUE, ...) {

  if (length(ny) != 1) {
    stop("You can show prediction plot only for one selected response variable.")
  }

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
    stop("Wrong value for ncomp argument.")
  }

  if (is.null(obj$y.ref)) {
    plot_data <- matrix(obj$y.pred[, ncomp, ny], ncol = 1)
    attr(plot_data, "xaxis.name") <- colnames(plot_data) <-
      paste0(obj$respnames[ny], ", predicted")
    attr(plot_data, "yaxis.name") <- attr(obj$y.pred, "yaxis.name")
    attr(plot_data, "yaxis.values") <- attr(obj$y.pred, "yaxis.values")
  } else {
    plot_data <- cbind(obj$y.ref[, ny], obj$y.pred[, ncomp, ny])
    colnames(plot_data) <- c(
      paste0(obj$respnames[ny], ", reference"),
      paste0(obj$respnames[ny], ", predicted")
    )
  }

  plot_data <- mda.setattr(plot_data, mda.getattr(obj$y.pred))
  rownames(plot_data) <- rownames(obj$y.pred)
  attr(plot_data, "name") <- paste0("Predictions (ncomp = ", ncomp, ")")

  if (!show.plot) {
    return(plot_data)
  }

  if (axes.equal && !is.null(obj$yref)) {
    xlim <- ylim <- range(plot_data)
  }

  p <- mdaplot(plot_data, type = "p", xlim = xlim, ylim = ylim, ...)

  if (is.null(obj$y.ref)) {
    return(invisible(p))
  }

  if (show.stat) {
    stat.text <- sprintf("nLV = %d\nRMSE = %.3f\nR2 = %.3f\n",
      ncomp, obj$rmse[ny, ncomp], obj$r2[ny, ncomp])
    text(p$xlim[1], p$ylim[2], stat.text, adj = c(0, 1), col = stat.col, cex = stat.cex)
  }

  if (show.line) {
    plotRegressionLine(p)
  }

  return(invisible(p))
}

#' Residuals plot for regression results
#'
#' @description
#' Shows plot with Y residuals (difference between predicted and reference values) for selected
#' response variable and complexity (number of components).
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param ncomp
#' complexity of model (e.g. number of components) to show the plot for
#' @param show.lines
#' allows to show the horisontal line at y = 0
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @export
plotResiduals.regres <- function(obj, ny = 1, ncomp = obj$ncomp.selected,
  show.lines = c(NA, 0), show.plot = TRUE, ...) {

  if (is.null(obj$y.ref)) {
    stop("Y-residuals can not be plotted without reference values.")
  }

  if (length(ny) != 1) {
    stop("You can make residuals plot only for one selected response variable.")
  }

  if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
    stop("Wrong value for ncomp argument.")
  }

  plot_data <- cbind(obj$y.ref[, ny], obj$y.ref[, ny] - obj$y.pred[, ncomp, ny])
  plot_data <- mda.setattr(plot_data, mda.getattr(obj$y.pred))
  colnames(plot_data) <- c(
    sprintf("%s, reference", obj$respnames[ny]),
    sprintf("%s, residuals", obj$respnames[ny])
  )
  attr(plot_data, "name") <- sprintf("Y-residuals (ncomp = %d)", ncomp)

  if (!show.plot) {
    return(plot_data)
  }

  return(mdaplot(plot_data, type = "p", show.lines = show.lines, ...))
}

#' RMSE plot for regression results
#'
#' @description
#' Shows plot with RMSE values vs. model complexity (e.g. number of components).
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param type
#' type of the plot
#' @param xticks
#' vector with ticks for x-axis
#' @param labels
#' what to use as labels ("names", "values" or "indices")
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ylab
#' label for y-axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @export
plotRMSE.regres <- function(obj, ny = 1, type = "b", xticks = seq_len(obj$ncomp),
  labels = "values", show.plot = TRUE, ylab = paste0("RMSE (", obj$respnames[ny], ")"), ...) {

  if (is.null(obj$rmse)) {
    stop("RMSE values are not available.")
  }

  if (length(ny) != 1) {
    stop("You can make residuals plot only for one selected response variable.")
  }

  plot_data <- mda.subset(obj$rmse, ny)
  attr(plot_data, "name") <- "RMSE"

  if (!show.plot) {
    return(plot_data)
  }

  return(mdaplot(plot_data, type = type, xticks = xticks, labels = labels, ylab = ylab, ...))
}

#' Plot method for regression results
#'
#' @details
#' This is a shortcut for \code{\link{plotPredictions.regres}}
#'
#' @param x
#' regression results (object of class \code{regres})
#' @param ...
#' other arguments
#'
#' @export
plot.regres <- function(x, ...) {
  plotPredictions.regres(x, ...)
}
