#' Title
#'
#' @param time.init
#' @param traj.init
#' @param time.new
#'
#'
#' @keywords internal
Curve.reduc.times <- function(time.init , traj.init, time.new){
  new.Curve <- matrix(NA,length(time.new),2)
  for (j in 1:length(time.new)){
    w.time <- which.min(abs(time.new[j]-time.init))
    if (round(time.init[w.time]-time.new[j],5)==0){
      new.Curve[j,] <- c(time.new[j], traj.init[w.time])
    }
    else {
      t_g <- (time.new[j]>time.init[w.time])*(time.init[w.time]) + (time.new[j]<time.init[w.time])*(time.init[w.time-1])
      t_d <- (time.new[j]<time.init[w.time])*(time.init[w.time]) + (time.new[j]>time.init[w.time])*(time.init[w.time+1])
      Y_g <- (time.new[j]>time.init[w.time])*(traj.init[w.time]) + (time.new[j]<time.init[w.time])*(traj.init[w.time-1])
      Y_d <- (time.new[j]<time.init[w.time])*(traj.init[w.time]) + (time.new[j]>time.init[w.time])*(traj.init[w.time+1])
      d1 <- time.new[j]-t_g
      d2 <- t_d - time.new[j]
      new.Curve[j,] <- c(time.new[j], (1 - (d1/(d1+d2)))*Y_g + (1 - (d2/(d1+d2)))*Y_d)
    }
  }
  return(new.Curve)
}

#' OOB predictions and error for Extra Frechet rf
#'
#' @param rf
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Image
#' @param timeScale
#' @param Y
#'
#' @import kmlShape
#'
#' @keywords internal
OOB.rfshape <- function(rf, Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL, Y, timeScale){

  err <- rep(NA,length(unique(Y$id)))

  Curve_courant <- NULL
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Image_courant <- NULL


  if (Y$type=="curve"){
    oob.pred <- list()

    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf)){

        BOOT <- rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)

        if (is.element(indiv, oob)== TRUE){

          if (is.null(Curve)==FALSE){

            Curve_courant = Curve
            for (k in 1:length(Curve_courant)){
              id_wXCurve <- which(Curve[[k]]$id == indiv)
              Curve_courant[[k]] <- Curve[[k]][id_wXCurve,]
            }

          }

          if (is.null(Scalar)==FALSE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list( X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.null(Factor)==FALSE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.null(Image)==FALSE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Image=Image_courant,
                          Curve.app = Curve, Factor.app = Factor, Image.app = Image,timeScale = timeScale)

          courbe <- rf[,t]$Y_pred[[pred$pred]]

          pred_courant <- rbind(cbind(rep(t,nrow(courbe)[1]),courbe),pred_courant)
        }
      }
      mean_pred <- meanFrechet(pred_courant, timeScale = Y$timeScale)
      dp <- as.data.frame(Curve.reduc.times(mean_pred$times, mean_pred$traj, Y$time[w_y]))
      names(dp) <- c("x","y")
      oob.pred[[i]] <- dp
      names(oob.pred)[i]= indiv
      err[i] <- distFrechet(dp$x, dp$y, Y$time[w_y], Y$Y[w_y], timeScale = Y$timeScale)^2
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="factor"){
    oob.pred <- rep(NA, length(unique(Y$id)))

    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- rep(NA,ncol(rf))
      for (t in 1:ncol(rf)){

        BOOT <- rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)

        if (is.element(indiv, oob)== TRUE){

          if (is.null(Curve)==FALSE){

            Curve_courant = Curve
            for (k in 1:length(Curve_courant)){
              id_wXCurve <- which(Curve[[k]]$id == indiv)
              Curve_courant[[k]] <- Curve[[k]][id_wXCurve,]
            }

          }

          if (is.null(Scalar)==FALSE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list( X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.null(Factor)==FALSE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.null(Image)==FALSE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Image=Image_courant,
                          Curve.app = Curve, Factor.app = Factor, Image.app = Image,timeScale = timeScale)

          pred_courant[t] <- pred$pred
        }
      }
      oob.pred[i] <- as.factor(attributes(which.max(table(na.omit(pred_courant)))))
    }
    return(list(err=1*(oob.pred!=Y$Y),oob.pred=oob.pred))
  }

  if (Y$type=="scalar"){
    oob.pred <- rep(NA, length(unique(Y$id)))

    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- rep(NA,ncol(rf))
      for (t in 1:ncol(rf)){

        BOOT <- rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)

        if (is.element(indiv, oob)== TRUE){

          if (is.null(Curve)==FALSE){

            Curve_courant = Curve
            for (k in 1:length(Curve_courant)){
              id_wXCurve <- which(Curve[[k]]$id == indiv)
              Curve_courant[[k]] <- Curve[[k]][id_wXCurve,]
            }

          }

          if (is.null(Scalar)==FALSE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list( X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.null(Factor)==FALSE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          if (is.null(Image)==FALSE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }

          pred <- pred.FT(rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Image=Image_courant,
                          Curve.app = Curve, Factor.app = Factor, Image.app = Image,timeScale = timeScale)

          pred_courant[t] <- pred$pred
        }
      }
      oob.pred[i] <- mean(na.omit(pred_courant))
    }
    return(list(err=(oob.pred[i]-Y$Y[w_y])^2,oob.pred=oob.pred))
  }

}

