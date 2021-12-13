#' Randomized Frechet tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param mtry [integer]:
#' @param timeScale [numeric]:
#' @param ntry [numeric]:
#' @param ... : option
#'
#' @import kmlShape
#' @import stringr
#' @import stats
#'
#' @export
Rtmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Image=NULL,Y,mtry,ntry=3, timeScale, ...){

  V <- NULL
  if (is.null(Curve)==FALSE)  V <- c(V, rep("curve",length(Curve)))
  if (is.null(Scalar)==FALSE) V <- c(V,rep("scalar",ncol(Scalar$X)))
  if (is.null(Factor)==FALSE) V <- c(V,rep("factor",ncol(Factor$X)))
  if (is.null(Image)==FALSE) V <- c(V,rep("image",dim(Image$X)[3]))


  impurity_feuilles <- NULL
  V_split <- NULL
  hist_nodes <- list()
  id_boot <- unique(sample(unique(Y$id), length(unique(Y$id)), replace=TRUE))
  boot <- id_boot
  decoupe <- 1

  wY = which(Y$id %in% id_boot)

  Y_pred <- list()
  Y_pred_surv  <- list()


  id_feuille  <- rbind(unique(Y$id),rep(NA,length(unique(Y$id))))
  id_feuille[2,which(id_feuille[1,] %in% id_boot)] <- 1
  qui <- id_feuille #### localisation des feuilles de l'arbre
  id_feuille_prime = id_feuille

  for (p in 1:(length(id_boot)/2-1)){
    count_split <- 0
    for (i in unique(na.omit(id_feuille[2,]))){

      variables <- sample(V,mtry) # Maintenant on sait combien on doit en tirer dans chaque espace
      # On ne va regarder que les espaces tirés :
      split.spaces <- unique(variables)

      w <- which(id_feuille[2,]==i)
      qui[2,w]=1
      qui[2,-w]=NA

      if (length(unique(w))>1){

        if (is.element("curve",split.spaces)==TRUE) tiragecurve <- sample(1:length(Curve),length(which(variables=="curve")))

        if (is.element("scalar",split.spaces)==TRUE) tiragescalar <- sample(1:ncol(Scalar$X),length(which(variables=="scalar")))

        if (is.element("factor",split.spaces)==TRUE) tiragefactor <- sample(1:ncol(Factor$X),length(which(variables=="factor")))

        if (is.element("image",split.spaces)==TRUE) tirageimage <- sample(1:dim(Image$X)[3],length(which(variables=="image")))

        F_SPLIT <- data.frame(vari = as.character(split.spaces), hetero= rep(Inf,length(split.spaces)))

        if (is.element("factor",split.spaces)==TRUE){

          feuille_split_factor = list(Pure = TRUE)
          tryCatch({
            feuille_split_factor <- ERvar_split(X=Factor,Y=Y,V=tiragefactor,qui=qui, type="factor",timeScale=timeScale,ntry = ntry)
          }, error = function(sp){feuille_split_factor = list(Pure = TRUE)})

          if (feuille_split_factor$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="factor"),2] <- feuille_split_factor$impurete}
        }

        if (is.element("curve",split.spaces)==TRUE){
          feuille_split_curve = list(Pure = TRUE)
          tryCatch({
            feuille_split_curve <- ERvar_split(X=Curve,Y,V=tiragecurve,qui=qui,type="curve",timeScale=timeScale, ntry=ntry)
            }, error = function(sp){feuille_split_curve = list(Pure = TRUE)})

          if (feuille_split_curve$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="curve"),2] <- feuille_split_curve$impurete
          }

        }

        if (is.element("scalar",split.spaces)==TRUE){

          feuille_split_scalar = list(Pure = TRUE)
          tryCatch({
            feuille_split_scalar <- ERvar_split(X=Scalar,Y=Y,V=tiragescalar,qui=qui,type="scalar",timeScale=timeScale, ntry=ntry)
          }, error = function(sp){feuille_split_scalar = list(Pure = TRUE)})

          if (feuille_split_scalar$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="scalar"),2] <- feuille_split_scalar$impurete
          }
        }


        if (is.element("image",split.spaces)==TRUE){


          feuille_split_image = list(Pure = TRUE)
          tryCatch({
            feuille_split_image <- ERvar_split(X=Image,Y=Y,V=tirageimage,qui=qui,type="image",timeScale=timeScale, ntry=ntry)
          }, error = function(sp){feuille_split_image = list(Pure = TRUE)})

          if (feuille_split_image$Pure==FALSE){
            F_SPLIT[which(F_SPLIT[,1]=="image"),2] <- feuille_split_image$impurete
          }
        }


        if (min(F_SPLIT[,2])<Inf){

          TYPE <- F_SPLIT[which.min(F_SPLIT[,2]),1]

          feuille_split <- get(paste("feuille_split_",TYPE, sep=""))

          Vari = get(paste(str_to_upper(str_sub(TYPE,1,1)),str_sub(TYPE,2), sep=""))

          vsplit_space <- feuille_split$variable

          gauche_id <- feuille_split$split[1,which(feuille_split$split[2,]==1)]
          droit_id <- feuille_split$split[1,which(feuille_split$split[2,]==2)]


          V_split <- rbind(V_split,c(as.character(TYPE),i,vsplit_space,feuille_split$centers[1],feuille_split$centers[2]))


          id_feuille_prime[2,which(feuille_split$split[2,]==1)] <- 2*i
          id_feuille_prime[2,which(feuille_split$split[2,]==2)] <- 2*i+1

          #print(paste("Split on the variable", vsplit_space, "on the space of ", paste(TYPE,"s",sep="")))


          if (TYPE=="factor"){
            w_gauche <- which(Vari$id %in% gauche_id)
            w_droit <- which(Vari$id %in% droit_id)
            hist_nodes[[2*i]] <- unique(Vari$X[w_gauche, vsplit_space])
            hist_nodes[[2*i+1]] <- unique(Vari$X[w_droit,vsplit_space])
          }

          count_split <- count_split+1
        }
      }
    }

    id_feuille <- id_feuille_prime

    if (count_split ==0 ){

      V_split <- data.frame(V_split)
      names(V_split) <- c("type","num_noeud", "var_split","idG","idD")
      for (q in unique(na.omit(id_feuille[2,]))){

        wq = id_feuille[1,which(id_feuille[2,]==q)]
        w <- which(Y$id %in% wq)

        if (Y$type=="curve"){
          if (length(q)>1){
            Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y$id[w], Y$time[w], Y$Y[w]), timeScale = Y$timeScale)
          }
          else Y_pred[[q]] <- data.frame(time=Y$time[w],traj=Y$Y[w])
        }

        if(Y$type=="scalar"){
          Y_pred[[q]]<- mean(Y$Y[w])
        }

        if(Y$type=="factor"){
          Table <- which.max(table(Y$Y[w]))
          Y_pred[[q]] <-  as.factor(attributes(Table)$names)
        }

        if (Y$type=="image"){
          Y_pred[[q]] <- apply(Y$Y[w,,drop=FALSE], 2, "mean")
        }

      }
      if (Y$type=="factor"){
        Ylevels <- unique(Y$Y)
        return(list(feuilles = id_feuille,Ytype=Y$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, boot=boot, Ylevels=Ylevels, timeScale=timeScale))
      }
      return(list(feuilles = id_feuille, Ytype=Y$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, boot=boot, timeScale=timeScale))
    }
  }


  V_split <- data.frame(V_split)
  names(V_split) <- c("type","num_noeud", "var_split", "idG", "idD")
  for (q in unique(na.omit(id_feuille[2,]))){

    wq = id_feuille[1,which(id_feuille[2,]==q)]
    w <- which(Y$id %in% wq)

    if (Y$type=="curve"){
      if (length(q)>1){
        Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y$id[w], Y$time[w], Y$Y[w]), timeScale = Y$timeScale)
      }
      else Y_pred[[q]] <- data.frame(time=Y$time[w],traj=Y$Y[w])
    }

    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y$Y[w])
    }

    if(Y$type=="factor"){
      Table <- which.max(table(Y$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }

    if (Y$type=="image"){
      Y_pred[[q]] <- apply(Y$Y[w,,drop=FALSE], 2, "mean")
    }

  }

  if (Y$type=="factor"){
    Ylevels <- unique(Y$Y)
    return(list(feuilles = id_feuille,Ytype=Y$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, boot=boot, Ylevels=Ylevels, timeScale=timeScale))
  }
  return(list(feuilles = id_feuille, Ytype=Y$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, boot=boot, timeScale=timeScale))
}
