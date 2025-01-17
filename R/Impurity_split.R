#' Impurity Split
#'
#' @param split
#'
#' @inheritParams FrechForest
#'
#' @keywords internal
impurity_split <- function(Y, split, d_out = 0.1, ...){
  impur <- 0
  imp <- list()
  for (i in 1:2){
    fils <- unique(Y$id)[which(split==i)]
    prop <- length(fils)/length(unique(Y$id))
    if (Y$type=="curve"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="curve",Y=Y$Y[w],id=Y$id[w],time=Y$time[w]), ...)
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type=="image"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w,which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="image", Y=Y$Y[w,], id=Y$id[w]), ...)
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type=="shape"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      if (length(w)>1){imp[[i]] <- impurity(list(type=Y$type,Y=Y$Y[,,w],id=Y$id[w]), ...)
      impur <- impur + imp[[i]]*prop}
      else {imp[[i]] <- 0}
    }

    if (Y$type=="scalar" || Y$type=="factor") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type=Y$type,Y=Y$Y[w],id=Y$id[w]), ...)
      impur <- impur + imp[[i]]*prop
    }

  }
  return(list(impur=impur, imp_list=imp))
}

