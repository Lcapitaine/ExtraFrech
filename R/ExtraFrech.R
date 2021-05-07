#' Parallelized Frechet random Forest
#'
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Image
#' @param Y
#' @param mtry
#' @param ntree
#' @param ncores
#' @param timeScale
#' @param ntry
#'
#' @import foreach
#' @import kmlShape
#' @import doParallel
#' @import pbapply
#'
#' @keywords internal
rf_para <- function(Curve=NULL, Scalar=NULL, Factor=NULL,Image=NULL,Y,mtry,ntree, ncores,ntry=3,timeScale){

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  Curve=Curve
  Scalar = Scalar
  Factor = Factor
  Image=Image
  Y=Y
  mtry=mtry
  ntry=ntry
  timeScale=timeScale

  trees <- pbsapply(1:ntree, FUN=function(i){
    Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor,Image=Image,Y=Y,mtry=mtry,ntry=ntry,timeScale=timeScale)
  },cl=cl)

  parallel::stopCluster(cl)
  return(trees)
}



#' Extra Frechet Random Forest
#'
#' This function builds Frechet random Forest introduced by Capitaine et.al, this includes the OOB predictions, OOB errors and variable importance computations.
#'
#'
#' @param Curve [list]: A list that contains the different input curves. It must contain the following elements (no choice): \code{X} the matrix of the different curves, each column code for a different curve variable; \code{id} is the vector of the identifiers for the different trajectories contained in \code{X}; \code{time} is the vector of the measurement times associated with the trajectories contained in \code{X}.
#' @param Scalar [list]: A list that contains the different input scalars. It must contain the following elements (no choice):  \code{X} the matrix of the scalars, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Factor [list]: A list that contains the different input factors. It must contain the following elements (no choice):  \code{X} the matrix of the factors, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Image [list]: A list that contains the different input images. It must contain the following elements (no choice):  \code{X} the array of the images of dimension \code{n}x\code{m}x\code{l}x\code{p} where \code{n}*\code{m} is the size of each image, \code{l} is the number of images and \code{p} is the number of shapes variables; \code{id} is the vector of the identifiers for each individual.
#' @param Y [list]: A list that contains the output, It must contain the following elements (no choice): \code{type} defines the nature of the output, can be "\code{curve}", "\code{sclalar}", "\code{factor}", "\code{image}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value \code{p/3}
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param ncores [numeric]: Number of cores used to build Frechet randomized trees in parallel, defaulting to number of cores of the computer minus 1.
#' @param ntry [numeric]: Only with \code{ERT=TRUE}, allows to manage with randomness of the trees.
#' @param timeScale [numeric]: Allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping. Only used when there are trajectories either in input or output.
#' @param sigma [numeric]: vector of noise
#' @param ... : parameters to be passed at low levels
#'
#' @import stringr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import pbapply
#' @import DescTools
#'
#' @return A Frechet random forest which is a list of the following elements: \itemize{
#' \item \code{rf:} a list of the \code{ntree} randomized Frechet trees that compose the forest.
#' \item \code{xerror :} a vector containing the OOB prediction error of each randomized Frechet tree composing the forest.
#' \item \code{OOB.err: } a vector containing the OOB prediction error of each individual in the learning sample.
#' \item \code{OOB.pred: } a list of the OOB prediction for each individual in the learning set.
#' \item \code{Importance: } A vector containing the variables importance.
#' \item \code{varex: } “pseudo R-squared”: Percentage of variance explained.
#' }
#'
#' @export
#'
#'
ExtraFrech <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Image=NULL ,
                       Y, mtry=NULL, ntree=100,ncores=NULL, timeScale=NULL, sigma=NULL ,ntry=3, ...){

  Inputs = NULL
  if (is.null(Image)==FALSE){
    Inputs = rep("image", dim(Image$X)[3])
    im = Image$X
    for (j in 1:dim(im)[3]){
      M=N=sqrt(dim(im)[2])
      Q= Image_transfo(M,N, sigma=sigma[j])
      im[,,j]= transfo_lin(im[,,j],Q)
    }
    Image <- list(X=Image$X,id=Image$id)
  }


  # On récupère le nombre de variables au total :

  if (is.null(Curve)==FALSE) Inputs = c(Inputs, rep("curve",length(Curve)))
  if (is.null(Scalar)==FALSE) Inputs = c(Inputs, rep("scalar",ncol(Scalar$X)))
  if (is.null(Factor)==FALSE) Inputs = c(Inputs, rep("factor",ncol(Factor$X)))


  if (is.null(mtry)==TRUE || mtry> length(Inputs)){
    mtry <- floor(length(Inputs)/3)*(floor(length(Inputs)/3)>=1) + 1*(floor(length(Inputs)/3)<1)
  }

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()-1
  }

  print("Building the maximal Frechet trees...")

  debut <- Sys.time()
  rf <-  rf_para(Curve=Curve,Scalar=Scalar, Factor=Factor, Image=Image,Y=Y, mtry=mtry, ntree=ntree,ntry = ntry,timeScale = timeScale,ncores=ncores)
  temps <- Sys.time() - debut

  print("Forest constucted !")

  print("OOB error computation")

  xerror = rep(NA, ncol(rf))
  for (i in 1:ncol(rf)){
    xerror[i] = OOB.tree(rf[,i], Curve= Curve, Scalar=Scalar, Factor=Factor, Image=Image, Y=Y,timeScale=timeScale)
  }

  print("OOB predictions...")
  OOB.rf = OOB.rfshape(rf=rf, Curve= Curve, Scalar=Scalar, Factor=Factor, Image=Image, Y=Y,timeScale=timeScale)


  ### Il nous faut calculer calmement l'impureté initiale :::

  if (Y$type=="curve"){
    trajLong <- data.frame(id=Y$id,time=Y$time,traj=Y$Y)
    meanF = meanFrechet(trajLong = trajLong, timeScale = Y$timeScale)
    dp <- as.data.frame(Curve.reduc.times(meanF$times, meanF$traj, Y$time[w]))
    imp = 0
    for (k in unique(Y$id)){
      w = which(Y$id==k)
      imp <- imp + distFrechet(dp[,1], dp[,2],Y$time[w], Y$Y[w], timeScale = Y$timeScale)^2
    }
    imp.init = imp/length(unique(Y$id))
    varex = 1- mean(OOB.rf$err)/imp.init
  }

  if (Y$type=="scalar"){
    varex = 1 - mean(OOB.rf$err)/var(Y$Y)
  }

  if (Y$type == "factor"){
    varex = 1- mean(OOB.rf$err)/Entropy(table(Y$Y))
  }

  if (Y$type == "image"){
    varex = 1-  mean(OOB.rf$err)/mean(apply(Y$Y,2,"var"))
  }


  ### Maintenant on va calculer l'importance des variables ::

  Curve.perm <- Curve
  Scalar.perm <- Scalar
  Factor.perm <- Factor
  Image.perm <- Image

  Importance.Curve <- NULL
  Importance.Scalar <- NULL
  Importance.Factor <- NULL
  Importance.Image <- NULL

  if (is.null(Curve)==FALSE){
    print('Computing the importance on the space of curves')
    Curve.err <- matrix(NA, ncol(rf), length(Curve))

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    p=1
    Importance.Curve <- foreach::foreach(p=1:length(Curve),.packages = "kmlShape" ,.combine = "c") %dopar% {

    for (k in 1:ncol(rf)){

      BOOT <- rf[,k]$boot
      OOB <- setdiff(unique(Y$id), BOOT)

      id_boot_Curve <- which(Curve[[p]]$id %in% OOB)

      perm = sample(OOB)
      new_id = Curve[[p]]$id
      new_id2 = new_id

      for (z in 1:length(OOB)){
        new_id[which(new_id2 == OOB[z])] = perm[z]
      }

      Curve.perm[[p]][,3] = new_id

      Curve.err[k,p] <- OOB.tree(rf[,k], Curve=Curve.perm, Scalar = Scalar, Factor=Factor, Image=Image, Y=Y, timeScale=timeScale)
      Curve.perm[[p]] <- Curve[[p]]
    }

    res <- mean(Curve.err[,p]- xerror)
  }

  parallel::stopCluster(cl)
  }


  if (is.null(Scalar)==FALSE){
    print('Computing the importance on the space of scalars')
    Scalar.err <- matrix(NA, ncol(rf), ncol(Scalar$X))

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    p=1

    Importance.Scalar <- foreach::foreach(p=1:ncol(Scalar$X) ,.combine = "c") %dopar% {

      for (k in 1:ncol(rf)){

        BOOT <- rf[,k]$boot
        OOB <- setdiff(unique(Y$id), BOOT)

        id_boot_Scalar <- which(Scalar$id %in% OOB)

        Scalar.perm$X[id_boot_Scalar,p] = sample(Scalar.perm$X[id_boot_Scalar,p])


        Scalar.err[k,p] <- OOB.tree(rf[,k], Curve=Curve, Scalar = Scalar.perm, Factor=Factor, Image=Image, Y=Y, timeScale=timeScale)
        Scalar.perm$X[,p] <- Scalar$X[,p]
      }

      res <- mean(Scalar.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }

  if (is.null(Factor)==FALSE){
    print('Computing the importance on the space of factors')
    Factor.err <- matrix(NA, ncol(rf), ncol(Factor$X))

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    p=1

    Importance.Factor <- foreach::foreach(p=1:ncol(Factor$X) ,.combine = "c") %dopar% {

      for (k in 1:ncol(rf)){

        BOOT <- rf[,k]$boot
        OOB <- setdiff(unique(Y$id), BOOT)

        id_boot_Factor <- which(Factor$id %in% OOB)

        Factor.perm$X[id_boot_Factor,p] = sample(Factor.perm$X[id_boot_Factor,p])


        Factor.err[k,p] <- OOB.tree(rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor.perm, Image=Image, Y=Y, timeScale=timeScale)
        Factor.perm$X[,p] <- Factor$X[,p]
      }

      res <- mean(Factor.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }

  if (is.null(Image)==FALSE){
    print('Computing the importance on the space of images')
    Image.err <- matrix(NA, ncol(rf), dim(Image$X)[3])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    p=1

    Importance.Image <- foreach::foreach(p=1:dim(Image$X)[3] ,.combine = "c") %dopar% {

      for (k in 1:ncol(rf)){

        BOOT <- rf[,k]$boot
        OOB <- setdiff(unique(Y$id), BOOT)

        id_boot_Image <- which(Image$id %in% OOB)

        Image.perm$X[id_boot_Image,,p] = sample(Image.perm$X[id_boot_Image,,p])


        Image.err[k,p] <- OOB.tree(rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor, Image=Image.perm, Y=Y, timeScale=timeScale)
        Image.perm$X[,,p] <- Image$X[,,p]
      }

      res <- mean(Image.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }


  Imp=list(Curve=Importance.Curve, Scalar=Importance.Scalar,
           Factor=Importance.Factor, Image= Importance.Image)

  frf <- list(rf=rf, type=Y$type, Y=Y, levels=levels(Y$Y), time=temps, xerror=xerror, importance = Imp, varex=varex, OOB.err= OOB.rf$err, OOB.pred=OOB.rf$oob.pred)
  class(frf) <- c("ExtraFrech")
  return(frf)
}





