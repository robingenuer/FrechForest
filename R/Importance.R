#' @inheritParams FrechForest
#' @inheritParams OOB.rfshape
#' @keywords internal
Importance <- function(rf, Curve = NULL, Scalar = NULL, Factor = NULL,
                       Shape = NULL, Image = NULL, Y, ncores = NULL,
                       timeScale = 0.1, d_out = 0.1, ...) {

  startImp <- Sys.time()

  ntree <- ncol(rf$rf)
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  print("Computing the OOB error of each tree")

  oobTrees <- pbsapply(1:ntree, FUN = function(numTree){
    OOB.tree(rf$rf[, numTree], Curve = Curve, Scalar = Scalar,
             Factor = Factor, Shape = Shape, Image = Image, Y = Y,
             timeScale = timeScale, d_out = d_out)
  }, cl = cl)

  Curve.perm <- Curve
  Scalar.perm <- Scalar
  Factor.perm <- Factor
  Shape.perm <- Shape
  Image.perm <- Image

  Importance.Curve <- NULL
  Importance.Scalar <- NULL
  Importance.Factor <- NULL
  Importance.Shape <- NULL
  Importance.Image <- NULL

  if (is.element("curve", inputs) == TRUE) {
    print('Computing the importance for curves input variables')
    Curve.err <- matrix(NA, ntree, dim(Curve$X)[2])

    Importance.Curve <- pbsapply(1:dim(Curve$X)[2], FUN = function(p) {
      for (k in 1:ntree) {
        BOOT <- rf$rf[, k]$boot
        nboot <- length(unique(Y$id)) - length(BOOT)

        id_boot_Curve <- NULL
        for (i in 1:length(BOOT)) {
          id_boot_Curve <- c(id_boot_Curve, which(Curve$id == BOOT[i]))
        }

        Curve.perm$X[-id_boot_Curve, p] <- permutation_courbes(
          Curve$X[-id_boot_Curve, p], Curve$id[-id_boot_Curve])

        Curve.err[k, p] <- OOB.tree(
          rf$rf[, k], Curve = Curve.perm, Scalar = Scalar, Factor = Factor,
          Shape = Shape, Image = Image, Y = Y, timeScale = timeScale,
          d_out = d_out, ...)

        Curve.perm$X[, p] <- Curve$X[, p]
      }
      return(mean(Curve.err[, p] - oobTrees))
    }, cl = cl)
  }

  if (is.element("scalar",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of scalars')
    Scalar.err <- matrix(NA, ntree, dim(Scalar$X)[2])

    Importance.Scalar <- foreach::foreach(
      p = 1:dim(Scalar$X)[2], .packages = "kmlShape" ,
      .combine = "c") %dopar% {

        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)

          id_boot_Scalar <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Scalar <- c(id_boot_Scalar, which(Scalar$id==BOOT[i]))
          }


          Scalar.perm$X[-id_boot_Scalar,p] <- sample(
            Scalar.perm$X[-id_boot_Scalar,p])

          Scalar.err[k,p] <- OOB.tree(rf$rf[,k], Curve = Curve,
                                      Scalar = Scalar.perm, Factor = Factor, Shape = Shape,
                                      Image = Image, Y, timeScale=timeScale, d_out=d_out, ...)

        }
        Scalar.perm$X[,p] <- Scalar$X[,p]
        res <- mean(Scalar.err[,p]- oobTrees)
      }
  }

  if (is.element("factor",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of factors')
    Factor.err <- matrix(NA, ntree, dim(Factor$X)[2])

    Importance.Factor <- foreach::foreach(
      p = 1:dim(Factor$X)[2], .packages = "kmlShape" ,
      .combine = "c") %dopar% {

        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)

          id_boot_Factor <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Factor <- c(id_boot_Factor, which(Factor$id==BOOT[i]))
          }

          # Il faut maintenant faire la permutation :

          Factor.perm$X[-id_boot_Factor,p] <- sample(
            Factor.perm$X[-id_boot_Factor,p])

          Factor.err[k,p] <- OOB.tree(
            rf$rf[,k], Curve = Curve, Scalar = Scalar, Factor = Factor.perm,
            Shape = Shape, Image = Image, Y = Y, timeScale = timeScale,
            d_out = d_out, ...)

        }
        ##on remet la variable en place :::
        Factor.perm$X[,p] <- Factor$X[,p]
        res <- mean(Factor.err[,p]- oobTrees)
      }
  }

  if (is.element("shape",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of shapes')
    Shape.err <- matrix(NA, ntree, dim(Shape$X)[length(dim(Shape$X))])

    Importance.Shape <- foreach::foreach(
      p = 1:dim(Shape$X)[length(dim(Shape$X))], .packages = "kmlShape",
      .combine = "c") %dopar% {

        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)

          id_boot_Shape <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Shape <- c(id_boot_Shape, which(Shape$id==BOOT[i]))
          }

          # Il faut maintenant faire la permutation :

          Shape.perm$X[,,-id_boot_Shape,p] <- permutation_shapes(
            Shape.perm$X[,,-id_boot_Shape, p], Shape.perm$id[-id_boot_Shape])

          Shape.err[k,p] <- OOB.tree(
            rf$rf[,k], Curve = Curve, Scalar = Scalar, Factor = Factor,
            Shape = Shape.perm, Image = Image, Y = Y, timeScale = timeScale,
            d_out = d_out, ...)

        }
        ##on remet la variable en place :::
        Shape.perm$X[,,,p] <- Shape$X[,,,p]
        res <- mean(Shape.err[,p]- oobTrees)
      }
  }

  if (is.element("image",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of images')
    Image.err <- matrix(NA, ntree, dim(Image$X)[3])

    Importance.Image <- foreach::foreach(
      p = 1:dim(Image$X)[3], .packages = "kmlShape", .combine = "c") %dopar% {

        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)

          id_boot_Image <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Image <- c(id_boot_Image, which(Image$id==BOOT[i]))
          }

          # Il faut maintenant faire la permutation :

          Image.perm$X[-id_boot_Image,,p] <-
            Image.perm$X[-id_boot_Image,,p][sample(nboot),]

          Image.err[k,p] <- OOB.tree(
            rf$rf[,k], Curve = Curve, Scalar = Scalar, Factor = Factor,
            Shape = Shape, Image = Image.perm, Y = Y, timeScale = timeScale,
            d_out = d_out, ...)

        }
        ##on remet la variable en place :::
        Image.perm$X[,,p] <- Image$X[,,p]
        res <- mean(Image.err[,p]- oobTrees)
      }
  }

  parallel::stopCluster(cl)

  varImp <- list(
    Curve = as.vector(Importance.Curve),
    Scalar = as.vector(Importance.Scalar),
    Factor = as.vector(Importance.Factor),
    Shape = as.vector(Importance.Shape),
    Image = as.vector(Importance.Image))

  computTimeImp <- Sys.time() - startImp

  return(list(varImp = varImp,
              oobTrees = oobTrees,
              computTimeImp = computTimeImp))
}
