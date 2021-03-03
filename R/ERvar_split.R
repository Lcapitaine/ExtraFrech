#' Extremely randomized split
#'
#' @param X
#' @param Y
#' @param timeScale
#' @param ntry
#'
#' @import kmlShape
#' @import stats
#'
#' @keywords internal
ERvar_split <- function(X,Y,V,qui,type,ntry=3, timeScale){

  p=length(V)
  impur = rep(Inf,p)

  split <- matrix(NA,(p+1),ncol(qui)) ## On construit la matrice des splits:
  split_prime <- matrix(NA,(ntry+1),ncol(qui))

  split_prime[1,]=split[1,]=qui[1,]
  Centers = matrix(NA, p,2)
  Pure <- FALSE
  indiv = qui[1,which(qui[2,]==1)]


  for (i in 1:p){

    if (type=="factor"){

      if (length(unique(X$X[which(X$id %in% indiv),V[i]]))>1){
        L <- Fact.partitions(X$X[which(X$id %in% indiv),V[i]],X$id[which(X$id %in% indiv)])
        tirage <- sample(1:length(L), 1)
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::

        split[(i+1),which(qui[2,]==1)] <- rep(2,length(which(qui[2,]==1)))
        split[(i+1),which(qui[1,] %in% L[[tirage]])] = 1

        # Il faut maintenant regarder la qualité du découpage ::
        impur[i] <- impurity_split(Y,split[c(1,(i+1)),])
      }
      else {impur[i] <- Inf}
    }

    if(type=="curve"){


      id_centers <- matrix(NA,ntry,2)
      for (l in 1:ntry){
        id_centers[l,] <- sample(qui[1,which(qui[2,]==1)],2)
      }

      ### Il faut ensuite boucler sur le ntry
      impurete2 <- rep(Inf,ntry)

      for (k in 1:ntry){

        w_gauche <- which(X[[V[i]]]$id==id_centers[k,1])
        w_droit <- which(X[[V[i]]]$id==id_centers[k,2])

        for (l in qui[1,which(qui[2,]==1)]){

          w <- which(X[[V[i]]]$id==l)
          dg <- distFrechet(X[[V[i]]]$time[w_gauche],X[[V[i]]]$traj[w_gauche],X[[V[i]]]$time[w],X[[V[i]]]$traj[w], timeScale = timeScale[p])
          dd <- distFrechet(X[[V[i]]]$time[w_droit],X[[V[i]]]$traj[w_droit],X[[V[i]]]$time[w],X[[V[i]]]$traj[w], timeScale = timeScale[p])
          if (dg<=dd) split_prime[(k+1),which(qui[1,]==l)] <- 1
          else split_prime[(k+1),which(qui[1,]==l)] <- 2

        }

        if (length(na.omit(unique(split_prime[(k+1),])))>1){
          impurete2[k] <- impurity_split(Y,split_prime[c(1,(k+1)),])
        }
      }

      gagnant <- which.min(impurete2)
      split[(i+1),] <- split_prime[(gagnant+1),]
      impur[i] <- impurete2[gagnant]
      Centers[i,]=id_centers[gagnant,]
      }




    if (type=="image"){

      indiv = which(X$id %in% qui[1,which(qui[2,]==1)])
      if (nrow(X$X)>2){
        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(qui[1,which(qui[2,]==1)],2)
        }

        impurete2 = rep(Inf,ntry)

        for (k in 1:ntry){

          w_g <- which(X$id==id_centers[k,1])
          w_d <- which(X$id==id_centers[k,2])

          ### Il nous faut calculer la distance :
          dg = apply(apply(X$X[indiv,,V[i]],1,"-",X$X[w_g,,V[i]])^2,2,"mean")
          dd = apply(apply(X$X[indiv,,V[i]],1,"-",X$X[w_d,,V[i]])^2,2,"mean")

          split_prime[(k+1),indiv]=2
          split_prime[(k+1),indiv[which((dg<=dd)==TRUE)]]=1

          if (length(unique(split_prime[(k+1),]))>1){
            impurete2[k] <- impurity_split(Y,split_prime[c(1,(k+1)),])
          }
        }

        gagnant <- which.min(impurete2)
        split[(i+1),] <- split_prime[(gagnant+1),]
        impur[i] <- impurete2[gagnant]
        Centers[i,]= id_centers[gagnant, ]
      }

      else{impur[i] <- Inf}
    }

    if(type=="scalar"){

      indiv = which(X$id %in% qui[1,which(qui[2,]==1)])
      impurete2= rep(NA,ntry)

      if (length(indiv)>2){

        ### On doit tier les centres
        #centers <- sample(X$X[,i],2)

        centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          centers[l,] <- sample(X$X[indiv,V[i]],2)
        }


        for (k in 1:ntry){
          split_prime[(k+1),indiv]=2
          split_prime[(k+1),indiv[abs(centers[k,1]-X$X[indiv,V[i]])<= abs(centers[k,2]-X$X[indiv,V[i]])]] <- 1
          impurete2[k] = impurity_split(Y,split_prime[c(1,k+1),])
        }

        gagnant <- which.min(impurete2)
        split[(i+1),] <- split_prime[(gagnant+1),]
        impur[i] <- impurete2[gagnant]
        Centers[i,] = centers[gagnant,]
      }

      else{
        impur[i] <- Inf
        split[(i+1),indiv] = c(1,2)
      }
    }

  }

  if (min(impur)==Inf){
    Pure =TRUE
  }


  true_split <- which.min(impur)
  split <- split[c(1,(true_split+1)),]
  return(list(split=split, impurete=min(impur), variable=V[true_split], centers = Centers[true_split,],Pure=Pure))
}

