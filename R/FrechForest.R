#' Frechet Random Forest
#'
#' This function builds Frechet random Forest introduced by Capitaine et.al, this includes the OOB predictions, OOB errors and variable importance computations.
#'
#'
#' @param Curve [list]: A list that contains the different input curves. It must
#'   contain the following elements (no choice): \code{X} the matrix of the
#'   different curves, each column code for a different curve variable;
#'   \code{id} is the vector of the identifiers for the different trajectories
#'   contained in \code{X}; \code{time} is the vector of the measurement times
#'   associated with the trajectories contained in \code{X}.
#' @param Scalar [list]: A list that contains the different input scalars. It
#'   must contain the following elements (no choice):  \code{X} the matrix of
#'   the scalars, each column code for a different variable; \code{id} is the
#'   vector of the identifiers for each individual.
#' @param Factor [list]: A list that contains the different input factors. It
#'   must contain the following elements (no choice):  \code{X} the matrix of
#'   the factors, each column code for a different variable; \code{id} is the
#'   vector of the identifiers for each individual.
#' @param Shape [list]: A list that contains the different input shapes. It must
#'   contain the following elements (no choice):  \code{X} the array of the
#'   shapes of dimension \code{n}x2x\code{l}x\code{p} where \code{n} is the
#'   number of points for composing each shape, \code{l} is the number of shapes
#'   and \code{p} is the number of shapes variables, \code{id} is the vector of
#'   the identifiers for each individual.
#' @param Image [list]: A list that contains the different input images. It must
#'   contain the following elements (no choice):  \code{X} the array of the
#'   images of dimension \code{n}x\code{m}x\code{l}x\code{p} where
#'   \code{n}*\code{m} is the size of each image, \code{l} is the number of
#'   images and \code{p} is the number of shapes variables; \code{id} is the
#'   vector of the identifiers for each individual.
#' @param Y [list]: A list that contains the output, It must contain the
#'   following elements (no choice): \code{type} defines the nature of the
#'   output, can be "\code{curve}", "\code{sclalar}", "\code{factor}",
#'   "\code{shape}", "\code{image}"; \code{Y} is the output variable; \code{id}
#'   is the vector of the identifiers for each individuals, they should be the
#'   same as the identifiers of the inputs.
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at
#'   each split. The default value \code{p/3}
#' @param ntree [numeric]: Number of trees to grow. This should not be set to
#'   too small a number, to ensure that every input row gets predicted at least
#'   a few times.
#' @param ncores [numeric]: Number of cores used to build Frechet randomized
#'   trees in parallel, defaulting to number of cores of the computer minus 1.
#' @param ERT [logical]: If \code{TRUE} uses Extremly Randomized Frechet Trees
#'   to build the Frechet forest.
#' @param ntry [numeric]: Only with \code{ERT=TRUE}, allows to manage with
#'   randomness of the trees.
#' @param timeScale [numeric]: Allow to modify the time scale, increasing or
#'   decreasing the cost of the horizontal shift. If timeScale is very big, then
#'   the Frechet mean tends to the Euclidean distance. If timeScale is very
#'   small, then it tends to the Dynamic Time Warping. Only used when there are
#'   trajectories either in input or output.
#' @param importance [logical]: TRUE to compute the variables importance FALSE
#'   otherwise (default \code{importance=}TRUE)
#' @param nodesize [numeric]: minimal number of observations in a node.
#' @param d_out [numeric]: Time scale for the output curves (\code{d_out=0.1} by
#'   default).
#' @param ... : optional parameters to be passed to the low level function
#'
#'
#' @return A Frechet random forest which is a list of the following elements:
#' \itemize{
#' \item \code{rf:} a list of the \code{ntree} randomized Frechet trees that
#' compose the forest.
#'
#' \item \code{oobTrees :} a vector containing the OOB error of each tree
#' of the Frechet random forest.
#'
#' \item \code{OOB.err: } a vector containing the OOB prediction error of each
#' individual in the learning sample.
#'
#' \item \code{OOB.pred: } a list of the OOB prediction for each individual in
#' the learning set.
#'
#' \item \code{Importance: } A vector containing the variables importance.
#'
#' \item \code{pseudoR2: } “pseudo R-squared”: Percentage of variance explained.
#' }
#'
#' @import emdist
#' @import Evomorph
#' @import doParallel
#' @import foreach
#' @import geomorph
#' @import kmlShape
#' @import parallel
#' @import pbapply
#' @import stringr
#' @import stats
#'
#' @export
#'
FrechForest <- function(Curve = NULL, Scalar = NULL, Factor = NULL,
  Shape = NULL, Image = NULL, Y, mtry = NULL, ntree = 100, ncores = NULL,
  ERT = FALSE, timeScale = 0.1, ntry = 3, nodesize = 1, importance = TRUE,
  d_out = 0.1, ...){

  startAll <- Sys.time()

  ### On va regarder les differentes entrees:
  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  if (Y$type=="shape"){
    Y$Y <- gpagen(Y$Y,print.progress = FALSE)$coords
  }

  # On recupère le nombre de variables au total :
  nvar <- 0
  for (k in Inputs){
    nvar <- nvar + dim(get(k)$X)[length(dim(get(k)$X))]
  }

  if (is.null(mtry)==TRUE || mtry> nvar){
    mtry <- floor(nvar/3)*(floor(nvar/3)>=1) + 1*(floor(nvar/3)<1)
  }

  if (is.null(Shape)!=TRUE || is.null(Image)!=TRUE) ERT <- TRUE

  size <- NULL
  if (Y$type=="shape" || Y$type=="image"){
    size <- c(dim(Y$Y)[1],dim(Y$Y)[2])
  }

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()-1
  }

  print("Building the maximal Frechet trees")

  startFRF <- Sys.time()

  rf <-  rf_shape_para(
    Curve = Curve, Scalar = Scalar, Factor = Factor, Shape = Shape,
    Image = Image, Y = Y, mtry = mtry, ntree = ntree, ERT = ERT, ntry = ntry,
    timeScale = timeScale, nodesize = nodesize, ncores = ncores, ...
  )

  if (Y$type=="shape" || Y$type=="image"){
    rf <- list(type=Y$type, rf=rf, size = dim(Y$Y))
  }
  else {
    rf <- list(type=Y$type, rf=rf, levels=levels(Y$Y))
  }

  print("Computing the OOB error of the forest")

  oob.err <- OOB.rfshape(rf, Curve = Curve, Scalar =Scalar, Factor = Factor,
    Shape = Shape, Image = Image, Y = Y, timeScale = timeScale, d_out = d_out,
    ncores = ncores, ...)

  if (Y$type == "image"){
    varTot <- apply(Y$Y, 2, "var")
    pseudoR2 <- 1 - apply(oob.err$err, 2, "mean") / varTot
  } else {
    varTot <- impurity(Y = Y, timeScale = d_out, ...)
    pseudoR2 <- 1 - mean(oob.err$err) / varTot
  }

  computTimeFRF <- Sys.time() - startFRF

  if (importance == TRUE) {

    print("Computing permutation-based variable importance")

    imp <- Importance(rf = rf, Curve = Curve, Scalar = Scalar,
      Factor = Factor, Shape = Shape, Image = Image, Y = Y, ncores = ncores,
      timeScale = timeScale, d_out = d_out, ...)
    varImp <- imp$varImp
    oobTrees <- imp$oobTrees
    computTimeImp <- imp$computTimeImp
  } else {
    varImp <- NULL
    oobTrees <- NULL
    computTimeImp <- NULL
  }

  computTimeTot <- Sys.time() - startAll

  FRF <- list(
    rf = rf$rf,
    type = rf$type,
    levels = rf$levels,
    oob.err = oob.err$err,
    oob.pred = oob.err$oob.pred,
    pseudoR2 = pseudoR2,
    varImp = varImp,
    oobTrees = oobTrees,
    computTimeFRF = computTimeFRF,
    computTimeImp = computTimeImp,
    computTimeTot = computTimeTot,
    size = size
  )
  class(FRF) <- c("FrechForest")
  return(FRF)
}
