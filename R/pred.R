#' Predict with Extra Frechet random forests
#'
#' @param object : Frechet random forest
#' @param Curve [list]: bla bla
#' @param Scalar [list]: bla bla
#' @param Factor [list]: bla bla bla
#' @param Image [list]: bla bla
#' @param Curve.app [list]: bla bla
#' @param Factor.app [list]: bla bla
#' @param Image.app [list]: bla bla
#' @param ... : bla bla
#'
#'
#' @import kmlShape
#' @import stringr
#' @import stats
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#'
#' @export
#'
predict.ExtraFrech <- function(object, Curve=NULL,Scalar=NULL,Factor=NULL, Image=NULL,
                              Curve.app=NULL, Factor.app= NULL, Image.app=NULL,...){
  # La première étape est de toujours lire les prédicteurs ::
  if (is.null(Curve)==FALSE) Id.pred <- unique(Curve[[1]]$id)
  else if (is.null(Scalar)==FALSE) Id.pred <- unique(Scalar$id)
  else if (is.null(Factor)==FALSE) Id.pred <- unique(Factor$id)
  else  Id.pred <- unique(Image$id)

  ncores <- detectCores()-1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)


  pred.feuille <- foreach(t=1:ncol(object$rf), .packages = "kmlShape" , .combine = "rbind") %dopar% {
    pred.FT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor,Image=Image, timeScale=object$timeScale,
                           Curve.app=Curve.app, Factor.app=Factor.app, Image.app=Image.app)$pred}

  parallel::stopCluster(cl)


  if (object$type=="scalar"){
    pred <- apply(pred.feuille, 2, "mean")
    return(pred)
  }

  if (object$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=object$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/ncol(object$rf)
    }
    prediction <- data.frame(pred=val, prob=proba)
    return(prediction)
  }

  if (object$type=="image"){
    pred <- matrix(0, length(Id.pred), ncol(object$Y$Y))
    for (l in 1:length(Id.pred)){
      pred_courant <- matrix(0,ncol(object$rf),ncol(object$Y$Y))
      for(k in 1:nrow(pred.feuille)){
        pred_courant[k,] <- object$Y$Y[which(object$Y$id== pred.feuille[k,l]),]
      }
      pred[l,] <-  apply(pred_courant, 2, "mean")
    }
  }

  if (object$type=="curve"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[,k]$Y_pred[[pred.feuille[k,l]]])[1]),object$rf[,k]$Y_pred[[pred.feuille[k,l]]]))
      }
      predouille <- kmlShape::meanFrechet(pred_courant, timeScale = object$Y$timeScale)
      predouille <- cbind(predouille, rep(Id.pred[l],dim(predouille)[1]))
      pred <- rbind(pred, predouille)
    }
    names(pred) <- c("times", "traj", "ID")
  }

  return(pred)
}
