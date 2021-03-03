#' Impurity Split
#'
#' @param Y
#' @param split
#' @param timeScale
#'
#' @import DescTools
#' @keywords internal
impurity_split <- function(Y,split){

  if (Y$type=="curve"){
    gauche = split[1,which(split[2,]==1)]
    droit = split[1,which(split[2,]==2)]
    w_gauche = which(Y$id %in% gauche)
    w_droit = which(Y$id %in% droit)

    if (length(droit)>0 & length(gauche)>0){
      trajLong_gauche <- data.frame(id=Y$id[w_gauche],time=Y$time[w_gauche],traj=Y$Y[w_gauche])
      trajLong_droit <- data.frame(id=Y$id[w_droit],time=Y$time[w_droit],traj=Y$Y[w_droit])

      meanF_gauche = meanFrechet(trajLong = trajLong_gauche, timeScale = Y$timeScale)
      meanF_droit = meanFrechet(trajLong = trajLong_droit, timeScale = Y$timeScale)

      imp_gauche=0
      imp_droit=0

      for (i in 1:(max(length(droit),length(gauche)))){

        w_gauche2 = which(Y$id==gauche[i])
        w_droit2 = which(Y$id==droit[i])

        if (length(w_gauche2)>0){
          imp_gauche <- imp_gauche + distFrechet(meanF_gauche$times, meanF_gauche$traj,Y$time[w_gauche2], Y$Y[w_gauche2], timeScale = Y$timeScale)^2
        }

        if (length(w_droit2)>0){
          imp_droit <- imp_droit + distFrechet(meanF_droit$times, meanF_droit$traj,Y$time[w_droit2], Y$Y[w_droit2], timeScale = Y$timeScale)^2
        }
      }

      return(imp_gauche*(length(gauche)/ncol(split)) + imp_droit*(length(droit)/ncol(split)))
    }
    else return(Inf)
  }


  if (Y$type=="image"){

    gauche = which(Y$id %in% split[1,which(split[2,]==1)])
    droit = which(Y$id %in% split[1,which(split[2,]==2)])

    if (length(gauche)>0 & length(droit)>0){
      return( mean(apply(Y$Y[gauche,],2,"var"))*(length(gauche)/length(gauche)) +
                mean(apply(Y$Y[droit,],2,"var"))*(length(droit)/length(droit)) )
    }

    else return(Inf)
  }


  if (Y$type=="scalar") {

    gauche = which(Y$id %in% split[1,which(split[2,]==1)])
    droit = which(Y$id %in% split[1,which(split[2,]==2)])

    if (length(gauche)>0 & length(droit)>0){
      return( var(Y$Y[gauche])*(length(gauche)/length(gauche)) + var(Y$Y[droit])*(length(droit)/length(droit)))
    }
    else return(Inf)
  }

  if (Y$type=="factor"){

    if (length(unique(Y$Y[which(Y$id %in% split[1,is.na(split[2,])==FALSE])]))==1) return(Inf)
    gauche = which(Y$id %in% split[1,which(split[2,]==1)])
    droit = which(Y$id %in% split[1,which(split[2,]==2)])

    if (length(gauche)>0 & length(droit)>0){
      return( Entropy(table(Y$Y[gauche]))*(length(gauche)/length(gauche)) + Entropy(table(Y$Y[droit]))*(length(droit)/length(droit)))
    }
    else return(Inf)
  }
}
