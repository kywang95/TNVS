###############################################################################

## The Main function of the TNVS

###############################################################################
TNVS_main <- function(data, class_col = 1, 
                      alpha_nonzero = 1e-2, alpha_rel = 0, alpha_red = 0.05,
                      dmax = NULL, d.early = TRUE, numCores = 1, simulation, iteration, i_fold){
  stopifnot(require(stringr) & require(infotheo) & require(snowfall) & require(rlecuyer))
  sfInit(parallel=TRUE, cpus=numCores,type="SOCK" )
  core_time <- Sys.time()
  # reduce the dimension with Gram-Schmidt transformation. Train the parameters 
  # using only the training set
  y <- data.frame(class = data[, class_col])
  x = data.frame(data[,-class_col])
  name_x <- colnames(x)
  y$class <- as.numeric(y$class)
  p = ncol(x)
  colnames(x) <- paste0("V", c(1:p))

  # delete the variables without information (sd(x) < 1e-5)
  entropy_x <- apply(infotheo::discretize(x), 2, infotheo::entropy)
  non_zero <- which( entropy_x > alpha_nonzero)
  
  zero <- colnames(x[,-non_zero])
  zero_info <- entropy_x[entropy_x <= alpha_nonzero]
  zero_data <- data.frame(zero, entropy=zero_info)
  zero_data$zero <- as.numeric(stringr::str_remove(string = zero_data$zero, pattern = "V"))
  filename = paste0(simulation,"_delete_lack_info_TNVS_", iteration, "_", i_fold, ".csv")
  write.csv(zero_data, filename, row.names = FALSE)
  
  name_x <- colnames(x)[non_zero]
  x <- data.frame(x[, non_zero])
  colnames(x) <- name_x
  # print(zero)
  # print(dim(x))
  
  x <- scale(x)
  N <- nrow(x)
  p = ncol(x)
  z = NULL
  stepT = c()
  red_variables <- NULL
  k_GS <- 1
  
  if (is.null(dmax)) d.stop <- p
  else d.stop <- dmax
  
  while (p>0 && length(z)<d.stop){
    cor_yx <- array(0, dim = p)
    cor_S <- codec_S(y[,1], z)
    if (cor_S <= 0) break
    codec_yx <- function(i){
      return(codec_Q(y[,1], x[,i], z)/cor_S)
    }
    sfExport("y","x","z", "cor_S", "codec_Q", 
             ".estimateQ", ".estimateConditionalQ", ".randomNN")
    sfClusterSetupRNG()
    cor_yx = sfLapply(1:p, codec_yx) 
    cor_yx = unlist(cor_yx)
    k = order(cor_yx, decreasing = TRUE)[1]
    t = max(cor_yx)
    stepT = c(stepT, t)
    D = data.frame(x[,k])
    colnames(D) <- colnames(x)[k]
    # print(c(round(t, 3), colnames(D)))
    if (is.null(z)){
      z = data.frame(D)
      
      D_GS <- z
    }else{
      name_z <- colnames(z)
      z = cbind(z, D)
      colnames(z) <- c(name_z, colnames(D))
      
      D_GS <- cbind(D_GS, D)
      colnames(D_GS) <- colnames(z)
      for (i in 1:(ncol(D_GS)-1)){
        d_GS <- t(D) %*% D_GS[,i] / (t(D_GS[,i]) %*% D_GS[,i])
        D_GS[, ncol(D_GS)] <- D_GS[, ncol(D_GS)] - c(d_GS) * D_GS[, i]
      }
    }
    if (t < alpha_rel && d.early) break
    if(ncol(x)>2){
      x <- x[, -k]
    }else
      if (ncol(x)==2){
        xnames <- colnames(x)
        x <- data.frame(x = x[,-k])
        colnames(x) <- xnames[-k]
      }else break
    
    # delete the variables which is highly correlated with the selected x[,k]
    # using Gram-Schimidt Orthogonalization
    Gk_data <- Gram_Schmidt_delete(D_GS, x, alpha = alpha_red)
    if (length(Gk_data)>0){
      Gk <- Gk_data$red
      # print(length(Gk))
      # x[, Gk] <- list(NULL)
      x <- x[, (!colnames(x) %in% Gk), drop = FALSE]
      red_variables <- rbind(red_variables, Gk_data)
    }
    p = ncol(x)
    
    # print(ncol(z))
    #    print(dim(x))
    
    # Save the parameters of the Gram-Schimidt Orthogonalization transformation.
    Gk_data$red <- as.numeric(str_remove(string = Gk_data$red, pattern = "V"))
    k_GS <- k_GS + 1
  }
  
  # Save the parameters of the Gram-Schimidt Orthogonalization transformation.
  red_variables$red <- as.numeric(str_remove(string = red_variables$red, pattern = "V"))
  filename = paste0(simulation,"_variables_redundant_TNVS_", iteration, "_", i_fold, ".csv")
  write.csv(red_variables, filename, row.names = FALSE)
  
  col_remain <- colnames(z)
  col_remain <- as.numeric(str_remove(string = col_remain, pattern = "V"))
  output_col <- data.frame(col_remain, stepT = stepT)
  
  # print(col_remain)
  # print(delete_variables)
  if (d.early && t<alpha_rel)
    output_col <- output_col[-length(col_remain),]
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = output_col, 
                  time = core_time)
  sfStop()
  return(results)
}


# codec_S --------------------------------------------------------------------------
#' Estimate the conditional dependence coefficient (CODEC)
codec_S <- function(Y, X = NULL){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    
    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    return(.estimateS(Y))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  
  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  return(.estimateConditionalS(Y, X))
}

# codec_Q --------------------------------------------------------------------------
codec_Q <- function(Y, Z, X = NULL, na.rm = TRUE){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    if(!is.matrix(Z)) {
      Z = as.matrix(Z)
    }
    if((length(Y) != nrow(Z))) stop("Number of rows of Y and X should be equal.")
    if (na.rm == TRUE) {
      # NAs are removed here:
      ok = complete.cases(Y,Z)
      Z = as.matrix(Z[ok,])
      Y = Y[ok]
    }
    
    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    q = ncol(Z)
    return(.estimateQ(Y, Z))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  if((length(Y) != nrow(Z))) stop("Number of rows of Y and Z should be equal.")
  if((nrow(Z) != nrow(X))) stop("Number of rows of Z and X should be equal.")
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(Y,Z,X)
    Z = as.matrix(Z[ok,])
    Y = Y[ok]
    X = as.matrix(X[ok,])
  }
  
  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  q = ncol(Z)
  p = ncol(X)
  
  return(.estimateConditionalQ(Y, X, Z))
}


# .estimateConditionalQ ------------------------------------------------------------
# Estimate Q(Y, Z | X), the numinator of the measure of conditional dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Z: Matrix of predictors (n by q)
# @param Y: Vector (length n)
#
# @return estimation \eqn{Q(Y, Z|X)}
.estimateConditionalQ <- function (Y, X, Z) {
  
  id <- group <- rnn <- NULL
  
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
  
  n = length(Y)
  
  W = cbind(X, Z)
  
  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  # handling repeated data
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  
  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  # compute the nearest neighbor of W
  nn_W = RANN::nn2(W, query = W, k = 3)
  nn_index_W = nn_W$nn.idx[, 2]
  repeat_data = which(nn_W$nn.dists[, 2] == 0)
  
  df_W = data.table::data.table(id = repeat_data, group = nn_W$nn.idx[repeat_data])
  df_W[, rnn := .randomNN(id), by = "group"]
  
  nn_index_W[repeat_data] = df_W$rnn
  # nearest neighbors with ties
  ties = which(nn_W$nn.dists[, 2] == nn_W$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_W[ties] = sapply(ties, helper_ties)
  }
  
  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_W]), 2, min),
            -apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  return(Q_n)
}




# .estimateConditionalS ------------------------------------------------------------
# Estimate S(Y, X), the denominator of the measure of dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
#
# @return estimation \eqn{S(Y, X)}
.estimateConditionalS <- function (Y, X){
  
  id <- group <- rnn <- NULL
  
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  n = length(Y)
  
  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  
  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  # estimate S
  R_Y = rank(Y, ties.method = "max")
  S_n = sum(R_Y - apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  
  return(S_n)
}


# estimateConditionalT -------------------------------------------------------------
# Estimate T(Y, Z | X), the measure of dependence of Y on Z given X
#
# @param Y: Vector (length n)
# @param Z: Matrix of predictors (n by q)
# @param X: Matrix of predictors (n by p)
#
# @return estimation of \eqn{T(Y, Z|X)}.
.estimateConditionalT <- function(Y, Z, X){
  S = .estimateConditionalS(Y, X)
  
  # happens only if Y is constant
  if (S == 0) {
    return(1)
  } else {
    return(.estimateConditionalQ(Y, X, Z) / S)
  }
}




# .estimateQ -----------------------------------------------------------------------
# Estimate Q(Y, X), the numinator of the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p).
# @param Y: Vector (length n).
#
# @return estimation of \eqn{Q(Y, X)}.
.estimateQ <- function(Y, X) {
  
  id <- group <- rnn <- NULL
  
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  
  n = length(Y)
  nn_X = RANN::nn2(X, query = X, k = 3)
  # remove the first nearest neighbor for each x which is x itself in case of no repeat data
  # when there is repeated data this is wrong but for repeated data we find the nearest
  # neighbors separately.
  nn_index_X = nn_X$nn.idx[, 2]
  
  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  # for the repeated data points, choose one of their identicals at random and set its index
  # as the index of the nearest neighbor
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  nn_index_X[repeat_data] = df_X$rnn
  
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  
  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  L_Y = rank(-Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_X]), 2, min) - (L_Y^2)/n) / (n^2)
  
  return(Q_n)
}



# .estimateS -------------------------------------------------------------------------
# Estimate S(Y) , the denominator of the measure of dependence of Y on X
#
# @param Y: Vector (length n).
# @return estimation of \eqn{S(Y)}.
.estimateS <- function (Y) {
  n = length(Y)
  L_Y = rank(-Y, ties.method = "max")
  S_n = gmp::asNumeric(sum(gmp::as.bigz(L_Y) * gmp::as.bigz(n - L_Y))) / (n^3)
  return(S_n)
}



# .estimateT -------------------------------------------------------------------------
# Estimate T(Y, X), the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
# @return estimation of \eqn{T(Y, X) = Q(Y, X) / S(Y)}.
.estimateT <- function(Y, X){
  # May 15, Mona removed FOCI:::
  S = .estimateS(Y)
  # happens only if Y is a constant vector.
  if (S == 0) {
    return(1)
  } else {
    return(.estimateQ(Y, X) / S)
  }
}


# randomNN -------------------------------------------------------------------------
# Find the random nearest neighbors.
#
# Gives the indices of nearest neighbors of points that are not unique in the data set.
# For each point we know that to which cluster of points it belongs (repeated points),
# we chose one of those indices that is not equal to the same index of our given point
# at random as its nearest neighbor.
#
# @param ids: a vector of ids of groups that each index is a member of.
#
# @return a vector of indices.
.randomNN <- function(ids) {
  m <- length(ids)
  
  x <- sample(x = (m - 1), m, replace = TRUE)
  x <- x + (x >= (1:m))
  
  return(ids[x])
}


# The Gram-Schmidt used to delete the redundant predictors--------------------------
Gram_Schmidt_delete <- function(D_GS, x, alpha = 0.05){
  stopifnot(require(snowfall))
  N = nrow(x)
  p = ncol(x)
  index_x <- array(0, dim = p)
  u <- x
  t_value <- rep(0, times = p)
  
  GS_ux <- function(i){
    for (j in 1:ncol(D_GS)){
      d <- t(u[,i]) %*% D_GS[,j] / (t(D_GS[,j]) %*% D_GS[,j])
      u[,i] <- u[,i]-c(d)*D_GS[,j]
    }
    if (var(u[,i])==0 | var(x[,i])==0){
      t_value = 0
    }else{
      t_value = var(u[,i])
    }
    return(t_value)
  }
  sfExport("u","x","D_GS")
  sfClusterSetupRNG()
  t_value = sfLapply(1:p, GS_ux) 
  t_value = unlist(t_value)
  
  Gk <- colnames(u)[which(t_value < alpha)]
  # print(round(t_value, 3))
  Gk_data <- data.frame(red=Gk, alpha = t_value[t_value < alpha])
  #  print(Gk)
  
  # t_bound = qt(1-alpha, df = N-2)
  # Gk <- colnames(u[, which(abs(t_value) < t_bound)])
  # print(Gk)
  # print(t_value[abs(t_value) < t_bound])
  
  return(Gk_data)
}



