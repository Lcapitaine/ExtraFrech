
#' Predict Frechet tree
#'
#' @param tree : Frechet tree.
#' @param Curve [list]: A list that contains the input curves.
#' @param Scalar [list]: A list that contains the input scalars.
#' @param Factor [list]: A list that contains the input factors.
#' @param Image [list]: A list that contains the input images.
#' @param timeScale [numeric]: Time scale for the input and output curves (\code{timeScale=0.1} by default)
#' @param Curve.app [list]: A list that contains the learning input curves.
#' @param Factor.app [list]: A list that contains the learning input factors.
#' @param Image.app [list]: A list that contains the learning input images.
#'
#' @import stringr
#' @import kmlShape
#'
#' @export
#'
pred.FT <- function(tree, Curve=NULL,Scalar=NULL,Factor=NULL,Image=NULL,
                    Curve.app=NULL,Factor.app=NULL,Image.app=NULL,timeScale){



  if (is.null(Curve)==FALSE) id.pred <- unique(Curve[[1]]$id)
  else if (is.null(Scalar)==FALSE) id.pred <- unique(Scalar$id)
  else if (is.null(Factor)==FALSE) id.pred <- unique(Factor$id)
  else  id.pred <- unique(Image$id)

  pred <- rep(NA,length(id.pred))

  if (tree$Y$type=="factor"){
    pred <- factor(rep(NA,length(id.pred)), levels=levels(tree$Y$Y))

  }

  for (i in 1:length(id.pred)){

    noeud_courant <- 1

    while (is.element(noeud_courant, tree$feuilles[2,])==FALSE){

      type <- as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1])
      var.split <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),3]))
      centers = tree$V_split[which(tree$V_split[,2]==noeud_courant),4:5]

      # Maintenant il nous faut regarder la différence entre la moyenne à gauche et a droite et conclure :

      if (type=="curve"){
        meanG = Curve.app[[var.split]][which(Curve.app[[var.split]]$id==centers$idG),c(2,1)]
        meanD = Curve.app[[var.split]][which(Curve.app[[var.split]]$id==centers$idD),c(2,1)]

        w = which(Curve[[var.split]]$id == id.pred[i])
        distG <- distFrechet(meanG[,1], meanG[,2],
                             Curve[[var.split]]$time[w], Curve[[var.split]]$traj[w], timeScale = timeScale[var.split])
        distD <- distFrechet(meanD[,1], meanD[,2],
                             Curve[[var.split]]$time[w], Curve[[var.split]]$traj[w], timeScale = timeScale[var.split])
      }

      if (type=="scalar"){

        w= which(Scalar$id==id.pred[i])
        distG <- abs(as.numeric(as.character(centers$idG)) - Scalar$X[w,var.split])
        distD <- abs(as.numeric(as.character(centers$idD)) -Scalar$X[w,var.split])

      }

      if (type=="image"){

        meanG = Image.app$X[which(Image.app$id==centers$idG),,var.split]
        meanD = Image.app$X[which(Image.app$id==centers$idD),,var.split]

        w= which(Image$id==id.pred[i])
        distG <- mean((Image$X[w,,var.split]-meanG)^2)
        distD <- mean((Image$X[w,,var.split]-meanD)^2)
      }

      if (type=="factor"){

        meanG <- tree$hist_nodes[[2*noeud_courant]]
        meanD <- tree$hist_nodes[[2*noeud_courant+1]]

        w= which(Factor$id==id.pred[i])
        distG <- -1*(is.element(Factor$X[w,var.split],meanG))
        distD <- -1*(is.element(Factor$X[w,var.split],meanD))
      }

      if (distG <= distD) { noeud_courant <- 2*noeud_courant}
      if (distD < distG) {noeud_courant <- 2*noeud_courant +1}

    }

    if(tree$Y$type=="curve" || tree$Y$type=="image"){
      pred[i] <- noeud_courant
    }

    else{
      pred[i] <- tree$Y_pred[[noeud_courant]]
    }
  }
  return(data.frame(id=id.pred,pred=pred))
}


