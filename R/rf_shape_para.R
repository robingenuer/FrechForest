#' @inheritParams FrechForest
#'
#' @keywords internal
rf_shape_para <- function(Curve = NULL, Scalar = NULL, Factor = NULL,
  Shape = NULL, Image = NULL, Y , mtry, ntree, ncores, ERT = FALSE, ntry = 3,
  nodesize = 1, timeScale = 0.1, ...){

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  trees <- pbsapply(1:ntree, FUN=function(i){
    Rtmax(Curve = Curve, Scalar = Scalar, Factor = Factor, Shape = Shape,
          Image = Image, Y = Y, mtry, ERT = ERT, ntry = ntry,
          nodesize = nodesize, timeScale = timeScale, ...)
  }, cl = cl)

  parallel::stopCluster(cl)

  return(trees)
}
