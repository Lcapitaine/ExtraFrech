#' OOB tree
#'
#' @param tree
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Image
#' @param Y
#' @param timeScale
#'
#' @import kmlShape
#' @import stringr
#'
#' @keywords internal
OOB.tree <- function(tree, Curve=NULL, Scalar=NULL, Factor=NULL, Image=NULL ,Y, timeScale){

  BOOT <- tree$boot
  OOB <- setdiff(unique(Y$id), BOOT)
  xerror <- rep(NA,length(OOB))
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Curve_courant <- NULL
  Image_courant <- NULL

  if (Y$type=="curve"){

    if (is.null(Image)==FALSE){
      id_wXImage <- which(Image$id %in% OOB)
      Image_courant <- list(X=Image$X[id_wXImage,,,drop=FALSE], id=Image$id[id_wXImage])
    }

    if (is.null(Factor)==FALSE){
      id_wXFactor <- which(Factor$id %in% OOB)
      Factor_courant <- list(X=Factor$X[id_wXFactor,,drop=FALSE], id=Factor$id[id_wXFactor])
    }

    if (is.null(Scalar)==FALSE){
      id_wXScalar <- which(Scalar$id %in% OOB)
      Scalar_courant <- list(X=Scalar$X[id_wXScalar,,drop=FALSE], id=Scalar$id[id_wXScalar])
    }



    if (is.null(Curve)==FALSE){
      Curve_courant = Curve
      for (k in 1:length(Curve_courant)){
        id_wXCurve <- which(Curve[[k]]$id %in% OOB)
        Curve_courant[[k]] <- Curve[[k]][id_wXCurve,]
      }
    }

    pred_courant <- pred.FT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Image=Image_courant,
                            Curve.app = Curve, Factor.app = Factor, Image.app = Image ,timeScale=timeScale)

    for (k in 1:nrow(pred_courant)){
      xerror[k] <- kmlShape::distFrechet(tree$Y_pred[[pred_courant$pred[k]]]$time,
                                                     tree$Y_pred[[pred_courant$pred[k]]]$traj,
                                                     Y$time[which(Y$id %in% pred_courant$id[k])], Y$Y[which(Y$id %in% pred_courant$id[k])])^2
    }
  }
  else {

    if (is.null(Image)==FALSE){
      id_wXImage <- which(Image$id %in% OOB)
      Image_courant <- list(X=Image$X[id_wXImage,,,drop=FALSE], id=Image$id[id_wXImage])
    }

    if (is.null(Factor)==FALSE){
      id_wXFactor <- which(Factor$id %in% OOB)
      Factor_courant <- list(X=Factor$X[id_wXFactor,,drop=FALSE], id=Factor$id[id_wXFactor])
    }

    if (is.null(Scalar)==FALSE){
      id_wXScalar <- which(Scalar$id %in% OOB)
      Scalar_courant <- list(X=Scalar$X[id_wXScalar,,drop=FALSE], id=Scalar$id[id_wXScalar])
    }



    if (is.null(Curve)==FALSE){
      Curve_courant = Curve
      for (k in 1:length(Curve_courant)){
        id_wXCurve <- which(Curve_courant[[k]]$id %in% OOB)
        Curve_courant[[k]] <- Curve_courant[[k]][id_wXCurve,]
      }
    }

    pred <- pred.FT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Image=Image_courant,
                            Curve.app = Curve, Factor.app = Factor, Image.app = Image ,timeScale=timeScale)

    ### Maintenant on calcule l'erreur de prédiction :::


    if (Y$type=="scalar"){xerror <- (Y$Y[which(Y$id %in% OOB)][order(Y$id[which(Y$id %in% OOB)])]-pred$pred[order(pred$id)])^2}
    if (Y$type=="factor"){xerror <- 1*(Y$Y[which(Y$id %in% OOB)][order(Y$id[which(Y$id %in% OOB)])] != pred$pred[order(pred$id)])}

    if (Y$type=="image"){
      xerror <- apply((Y$Y[which(Y$id %in% OOB),][order(Y$id[which(Y$id %in% OOB)]),]-pred$pred[order(pred$id)])^2,1,"mean")
    }
  }
  return(mean(xerror))
}



