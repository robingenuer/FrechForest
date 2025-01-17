#' Extremely randomized split
#'
#' @param X Inputs data
#' @inheritParams FrechForest
#' @inheritParams pred.FT
#'
#' @keywords internal
ERvar_split <- function(X, Y, ntry = 3, timeScale = 0.1, d_out = 0.1,
  FrechetSumOrMax = "max", ...){

  impur <- rep(0,dim(X$X)[length(dim(X$X))])
  toutes_imp <- list()
  impur_list = list()
  split <- list()
  Pure <- FALSE

  Imp_shape <- Inf
  var_shape <- Inf

  for (i in 1:dim(X$X)[length(dim(X$X))]){

    if (X$type=="factor"){

      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()

        # On tire une partition au hasard
        tirage <- sample(1:length(L), 1)
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::

        split[[i]] <- rep(2,length(X$id))
        for (l in L[[tirage]]){
          split[[i]][which(X$id==l)] <- 1
        }
        # Il faut maintenant regarder la qualite du decoupage ::
        impurete <- impurity_split(Y,split[[i]], ...)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }

    if( X$type=="curve"){

      # Il faut commencer par tirer les multiples centres ::

      id_centers <- matrix(NA,ntry,2)
      for (l in 1:ntry){
        id_centers[l,] <- sample(unique(X$id),2)
      }

      ### Il faut ensuite boucler sur le ntry
      split_prime <- matrix(2,ntry,length(unique(X$id)))
      u <- 0
      impurete2 <- list()
      qui <- NULL
      imp <- NULL

      for (c in 1:ntry){
        w_gauche <- which(X$id==id_centers[c,1])
        w_droit <- which(X$id==id_centers[c,2])

        for (l in 1:length(unique(X$id))){
          w <- which(X$id==unique(X$id)[l])
          dg <- distFrechet(
            X$time[w_gauche],X$X[w_gauche,i],X$time[w],X$X[w,i],
            timeScale = timeScale, FrechetSumOrMax = FrechetSumOrMax)
          dd <- distFrechet(
            X$time[w_droit],X$X[w_droit,i],X$time[w],X$X[w,i],
            timeScale = timeScale, FrechetSumOrMax = FrechetSumOrMax)
          if (dg<=dd) split_prime[c,l] <- 1
        }

        if (length(unique(split_prime[c,]))>1){
          u <- u + 1
          qui <- c(qui, c)
          impurete2[[c]] <- impurity_split(
            Y, split_prime[c,], d_out = d_out, ...)
          imp <- c(imp,impurete2[[c]]$impur)
        }

      }

      if (u>0){
        gagnant <- qui[which.min(imp)]
        split[[i]] <- split_prime[gagnant,]
        impurete <- impurete2[[gagnant]]
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }

      else{
        impur[i] <- Inf
        split[[i]] <- Inf}
    }


    if (X$type=="shape"){
      n_elem = dim(X$X)[3]
      if (n_elem>2){

        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(X$id,2)
        }

        split_prime <- matrix(2,ntry,length(X$id))

        dd = rep(NA,n_elem)
        dg = rep(NA,n_elem)

        for (c in 1:ntry){

          for (k in 1:n_elem){
            dg[k] <- emd2d(X$X[,,k,i],X$X[,,which(X$id==id_centers[c,1]),i])
            dd[k] <- emd2d(X$X[,,k,i],X$X[,,which(X$id==id_centers[c,2]),i])
          }

          for (l in 1:length(unique(X$id))){

            if (is.nan(dg[l]) || is.nan(dd[l])) split_prime[c,l] <- sample(c(1,2),1)
            else if (dg[l]<=dd[l]) split_prime[c,l] <- 1
          }
          if (length(split_prime[c,])>1){
            impurete2 <- impurity_split(Y, split_prime[c,], d_out = d_out, ...)

            if (impurete2$impur <Imp_shape && is.na(impurete2$impur)==FALSE){
              Imp_shape <- impurete2$impur
              var_shape <- i
              gauche = id_centers[c,1]
              droite = id_centers[c,2]

              impur_list = impurete2$imp_list

              split = split_prime[c,]
              Pure = FALSE
            }
          }
        }

      }
    }


    if (X$type=="image"){
      if (nrow(X$X)>2){
        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(X$id,2)
        }

        split_prime <- matrix(2,ntry,length(X$id))


        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL

        for (c in 1:ntry){

          w_g <- which(X$id==id_centers[c,1])
          w_d <- which(X$id==id_centers[c,2])
          ### Il nous faut calculer la distance :
          dg = apply(apply(X$X[,,i],1,"-",X$X[w_g,,i])^2,2,"mean")
          dd = apply(apply(X$X[,,i],1,"-",X$X[w_d,,i])^2,2,"mean")

          split_prime[c,which((dg<=dd)==TRUE)]=1
          if (length(unique(split_prime[c,]))>1){
            u <-u+1
            qui <- c(qui,c)
            impurete2[[c]] <- impurity_split(
              Y, split_prime[c,], d_out = d_out, ...)
            imp <- c(imp,impurete2[[c]]$impur)
          }
        }

        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }

        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }

      }

      else{
        split[[i]] <- c(1,2)
        impurete <- impurity_split(Y, split[[i]], d_out = d_out, ...)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
    }

    if(X$type=="scalar"){
      if (length(unique(X$X[,i]))>2){

        ### On doit tier les centres
        #centers <- sample(X$X[,i],2)

        centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          centers[l,] <- sample(X$X[,i],2)
        }

        #split[[i]] <- rep(2,length(X$X[,i]))
        split_prime <- matrix(2,ntry,length(X$X[,i]))

        for (l in 1:length(X$X[,i])){
          for (k in 1:ntry){
            if (abs(centers[k,1]-X$X[l,i])<= abs(centers[k,2]-X$X[l,i])) {
              split_prime[k,l] <- 1
            }
          }
        }

        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL
        for (k in 1:ntry){
          if (length(unique(split_prime[k,]))>1){
            u <- u+1
            qui <- c(qui,k)
            impurete2[[k]] <- c(
              impurete2, impurity_split(Y, split_prime[k,], d_out = d_out, ...))
            imp <- c(imp, impurete2[[k]]$impur)
          }
        }

        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }

        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }
      }

      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
  }

  if (Imp_shape<Inf){
    return(list(split = split, impurete = Imp_shape, gauche = gauche, droite= droite , variable = var_shape,Pure = Pure, impur_list = impur_list))
  }

  if (length(unique(impur))==1 & is.element(Inf,impur)==TRUE){
    return(list(Pure=TRUE))
  }
  true_split <- which.min(impur)
  split <- split[[true_split]]
  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur), Pure=Pure))
}

