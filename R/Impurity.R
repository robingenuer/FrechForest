#' Impurity
#'
#' Compute the impurity of a given vector
#'
#' @inheritParams FrechForest
#'
#' @keywords internal
impurity <- function(Y, d_out = 0.1, FrechetSumOrMax = "max", ...){

  if (Y$type=="curve"){
    trajLong <- data.frame(id = Y$id, time = Y$time, traj = Y$Y)
    meanF <- meanFrechet(trajLong = trajLong, timeScale = d_out,
                         FrechetSumOrMax = FrechetSumOrMax, ...)
    allDistF <- sapply(unique(Y$id), FUN = function(i) {
      distF <- distFrechet(
        meanF$times, meanF$traj, Y$time[which(Y$id == i)], Y$Y[which(Y$id == i)],
        timeScale = d_out, FrechetSumOrMax = FrechetSumOrMax)^2
      return(distF)
    })
    imp <- mean(allDistF)
  }

  if (Y$type == "image"){
    if (length(Y$id) == 1){
      imp <- 0
    }
    imp <- mean(apply(Y$Y,2,stats::var))
  }

  if (Y$type == "scalar"){
    if (length(Y$Y) == 1){
      imp <- 0
    }
    imp <- stats::var(Y$Y)
  }

  if (Y$type == "factor"){
    p <- table(Y$Y) / length(Y$Y)
    imp <- -sum(p * log2(p))
  }

  if (Y$type == "shape"){
    ms <- mshape(Y$Y[ , , , drop=FALSE])
    imp <- mean(ShapeDist(Y$Y,ms)^2)
  }

  return(imp)
}
