library(glasso)
library(Matrix)
Kappa <-
  function(beta, X, Y, lambda_kappa){
    n <- nrow(Y)  
    SigmaR <- 1/n * t(Y - X %*% beta) %*% (Y - X %*% beta)
    # if (any(eigen(SigmaR,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    #   stop("Residual covariance matrix is not non-negative definite")
    # }
    res <- glasso(SigmaR, lambda_kappa, penalize.diagonal = FALSE)
    return(as.matrix(forceSymmetric(res$wi)))
  }

beta_ridge_C <- function(X, Y, lambda_beta) {
  .Call('graphicalVAR_beta_ridge_C', PACKAGE = 'graphicalVAR', X, Y, lambda_beta)
}

Beta_C <- function(kappa, beta, X, Y, lambda_beta, lambda_beta_mat, convergence, maxit) {
  .Call('graphicalVAR_Beta_C', PACKAGE = 'graphicalVAR', kappa, beta, X, Y, lambda_beta, lambda_beta_mat, convergence, maxit)
}

VAR_logLik_C <- function(X, Y, kappa, beta) {
  .Call('graphicalVAR_VAR_logLik_C', PACKAGE = 'graphicalVAR', X, Y, kappa, beta)
}

LogLik_and_BIC <- function(X, Y, estimates) {
  .Call('graphicalVAR_LogLik_and_BIC', PACKAGE = 'graphicalVAR', X, Y, estimates)
}

Rothmana <-
  function(X, Y, lambda_beta, lambda_kappa, convergence = 1e-4, gamma = 0.5, maxit.in = 100, maxit.out = 100,
           penalize.diagonal, # if FALSE, penalizes the first diagonal (assumed to be auto regressions), even when ncol(X) != ncol(Y) !
           interceptColumn = 1 # Set to NULL or NA to omit
  ){
    # Algorithm 2 of Rothmana, Levinaa & Ji Zhua
    
    if(is.null(ncol(Y))){Y <- t(as.matrix(Y, ncol = length(Y)))}
    if(is.null(ncol(X))){X <- t(as.matrix(X, ncol = length(X)))}
    nY <- ncol(Y)
    nX <- ncol(X)
    
    if (missing(penalize.diagonal)){
      penalize.diagonal <- (nY != nX-1) & (nY != nX )
    }
    
    lambda_mat <- matrix(lambda_beta,nX, nY)
    if (!penalize.diagonal){
      if (nY == nX){
        add <- 0
      } else if (nY == nX - 1){
        add <- 1
      } else {
        stop("Beta is not P x P or P x P+1, cannot detect diagonal.")
      }
      for (i in 1:min(c(nY,nX))){
        lambda_mat[i+add,i] <- 0
      }
    }
    if (!is.null(interceptColumn) && !is.na(interceptColumn)){
      lambda_mat[interceptColumn,] <- 0
    }
    
    n <- nrow(X)
    beta_ridge <- beta_ridge_C(X, Y, lambda_beta)
    
    # Starting values:
    beta <- matrix(0, nX, nY)  
    
    # Algorithm:
    it <- 0
    
    repeat{
      it <- it + 1
      kappa <- Kappa(beta, X, Y, lambda_kappa)
      beta_old <- beta
      beta <- Beta_C(kappa, beta, X, Y, lambda_beta, lambda_mat, convergence, maxit.in) 
      
      if (sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){
        break
      }
      
      if (it > maxit.out){
        #warning("Model did NOT converge in outer loop")
        break
      }
    }
    
    ## Compute unconstrained kappa (codes from SparseTSCGM):
    ZeroIndex <- which(kappa==0, arr.ind=TRUE) ## Select the path of zeros
    WS <-  (t(Y)%*%Y - t(Y) %*% X  %*% beta - t(beta) %*% t(X)%*%Y + t(beta) %*% t(X)%*%X %*% beta)/(nrow(X))
    
    # if (any(eigen(WS,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    #   stop("Residual covariance matrix is not non-negative definite")
    # }
    
    if (nrow(ZeroIndex)==0){
      out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
    } else {
      out4 <- suppressWarnings(glasso(WS, rho = 0, zero = ZeroIndex,trace = FALSE))
    }
    lik1  = determinant( out4$wi)$modulus[1]
    lik2 <- sum(diag( out4$wi%*%WS))
    
    pdO = sum(sum(kappa[upper.tri(kappa,diag=FALSE)] !=0))
    pdB = sum(sum(beta[lambda_mat!=0] !=0)) # CHECK WITH LOURENS
    
    LLk <-  (n/2)*(lik1-lik2) 
    LLk0 <-  (n/2)*(-lik2)
    
    EBIC <-  -2*LLk + (log(n))*(pdO +pdB) + (pdO  + pdB)*4*gamma*log(2*nY)
    
    ### TRANSPOSE BETA!!!
    return(list(beta=t(beta), kappa=kappa, EBIC = EBIC))
  }

invGlasso <- function(x){
  if (all(eigen(x)$values > sqrt(.Machine$double.eps))){
    Xinv <- solve(x)
  } else {
    Xglas <- glasso(x,0.05,penalize.diagonal=FALSE)    
    Xinv <- Xglas$wi
  }
  Xinv
}
generate_lambdas <- function(
  X,
  Y,
  nLambda_kappa = 10,
  nLambda_beta = 10,
  lambda_min_kappa = 0.05,
  lambda_min_beta = 0.05, 
  penalize.diagonal=TRUE         
){
  N <- nrow(Y)
  P <- ncol(Y)
  
  #### Lambda sequence for Kappa:
  corY <- cov2cor(t(Y)%*%Y/nrow(Y))
  lam_K_max = max(abs(corY))
  lam_K_min = lambda_min_kappa*lam_K_max
  lam_K = exp(seq(log(lam_K_max), log(lam_K_min), length = nLambda_kappa))
  
  #### Lambda sequence for Beta
  # Initial estimate for Kappa:
  # Yinv <- invGlasso(t(Y) %*% Y / N)
  # Xinv <- invGlasso(t(X) %*% X / N)
  #   beta <- t(Y) %*% X %*% Xinv
  #   S <- 1/(nrow(Y)) * (
  #     t(Y) %*% Y -
  #       t(Y) %*% X %*% t(beta) -
  #       beta %*% t(X) %*% Y +
  #       beta %*% t(X) %*% X %*% t(beta)
  #   )
  #   S <- (S + t(S)) / 2
  #   if (any(eigen(S)$value < -sqrt(.Machine$double.eps))) stop("Residual covariances not postive definite")
  #   kappa <- invGlasso(S)
  #   kappa <- (kappa + t(kappa)) / 2
  # lam_B_max = max(abs((1/N)*t(X)%*%Y%*%Yinv))
  
  Yinv <- invGlasso(t(Y) %*% Y)
  lam_B_max = max(abs(t(X)%*%Y%*%Yinv))
  lam_B_min = lambda_min_beta*lam_B_max
  lam_B = exp(seq(log(lam_B_max), log(lam_B_min), length = nLambda_beta))
  
  return(list(lambda_kappa = lam_K, lambda_beta = lam_B))
}

summary.graphicalVAR <- function(object,...) print(object,...)

# Print method
print.graphicalVAR <- function(x, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  
  
  cat("=== graphicalVAR results ===")
  cat("\nNumber of nodes:",nrow(x[['kappa']]),
      "\nNumber of tuning parameters tested:",nrow(x[['path']]),
      "\nEBIC hyperparameter:",x[['gamma']],
      "\nOptimal EBIC:",x[['EBIC']],
      
      
      "\n\nNumber of non-zero Partial Contemporaneous Correlations (PCC):",sum(x[['PCC']][upper.tri(x[['PCC']],diag=FALSE)]==0) ,
      "\nPCC Sparsity:",mean(x[['PCC']][upper.tri(x[['PCC']],diag=FALSE)]==0) ,
      "\nNumber of PCC tuning parameters tested:",length(unique(x$path$kappa)),
      paste0("\nPCC network stored in ",name,"$PCC"),
      
      "\n\nNumber of non-zero Directed Contemporaneous Correlations (PDC):",sum(x[['PDC']][upper.tri(x[['PDC']],diag=FALSE)]==0) ,
      "\nPDC Sparsity:",mean(x[['PDC']][upper.tri(x[['PDC']],diag=FALSE)]==0) ,
      "\nNumber of PDC tuning parameters tested:",length(unique(x$path$beta)),
      paste0("\nPDC network stored in ",name,"$PDC"),
      
      paste0("\n\nUse plot(",name,") to plot the estimated networks.")
  )
}

# Plot method
plot.graphicalVAR <- plot.gVARmodel <- function(x, include = c("PCC","PDC"), repulsion = 1, horizontal = TRUE, titles = TRUE, sameLayout = TRUE, unweightedLayout = FALSE,...){
  qtitle <-  function (x) 
  {
    text(par("usr")[1] + (par("usr")[2] - par("usr")[1])/40, 
         par("usr")[4] - (par("usr")[4] - par("usr")[3])/40, x, 
         adj = c(0, 1))
  }
  
  if (length(include)>1){
    if (horizontal){
      layout(t(seq_along(include))) 
    } else {
      layout(seq_along(include))
    }
  }
  
  # Choose directed or undirected:
  if (unweightedLayout){
    wPCC <- 1*(x$PCC!=0)
    wPDC <- 1*(x$PDC!=0)
  } else {
    wPCC <- x$PCC
    wPDC <- x$PDC
  }
  
  if (sameLayout & all(c("PCC","PDC") %in% include)){
    Layout <- qgraph::averageLayout(as.matrix(wPCC), as.matrix(wPDC), repulsion=repulsion)
  }
  
  Res <- list()
  
  for (i in seq_along(include)){
    if ("PCC" == include[i]){
      if (sameLayout & all(c("PCC","PDC") %in% include)){
        
        Res[[i]] <- qgraph::qgraph(x$PCC, layout = Layout, ..., repulsion=repulsion)
      } else {
        L <- qgraph:::qgraph(wPCC,DoNotPlot=TRUE,...,repulsion=repulsion)$layout
        Res[[i]] <- qgraph::qgraph(x$PCC, layout = L,..., repulsion=repulsion)
      }
      
      if (titles){
        qtitle("Partial Contemporaneous Correlations")
      }
    }
    
    if ("PDC" == include[i]){
      if (sameLayout & all(c("PCC","PDC") %in% include)){
        Res[[i]] <- qgraph::qgraph(x$PDC, layout = Layout, ..., repulsion=repulsion, directed=TRUE)
      } else {
        L <- qgraph:::qgraph(wPDC,DoNotPlot=TRUE,...,repulsion=repulsion, directed=TRUE)$layout
        Res[[i]] <- qgraph::qgraph(x$PDC,layout=L, ..., repulsion=repulsion, directed=TRUE)
      }
      
      if (titles){
        qtitle("Partial Directed Correlations")
      }
    }
  }
  
  invisible(Res)
}

VARglm <-
  function(x,family,vars,adjacency,icfun = BIC,...)
  {
    # Returns estimated weights matrix of repeated measures data x
    ## x must be matrix, rows indicate measures and columns indicate variables
    # If adjacency is missing, full adjacency is tested
    # 'family' can be assigned family function (see ?family), list of such
    ## functions for each variable in x or character vector with names of the
    ## family functions.
    # 'vars' must be a vector indicating which variables are predicted, can be useful for parallel implementation.
    
    if (missing(x)) stop("'x' must be assigned")
    x <- as.matrix(x)
    
    Ni <- ncol(x)
    Nt <- nrow(x)
    
    # Check input:
    if (missing(vars)) vars <- 1:Ni
    No <- length(vars)
    
    if (missing(adjacency)) adjacency <- matrix(1,Ni,No)
    if (is.vector(adjacency)) adjacency <- as.matrix(adjacency)
    if (!is.matrix(adjacency) && ncol(adjacency)!=No && nrow(adjacency)!=Ni) stop("'adjacency' must be square matrix with a row for each predictor and column for each outcome variable.")
    
    if (any(apply(x,2,sd)==0))
    {
      adjacency[apply(x,2,sd)==0,] <- 0
      adjacency[,apply(x,2,sd)==0] <- 0
      warning("Adjacency matrix adjusted to not include nodes with 0 variance.")
    }
    
    if (missing(family)) 
    {
      if (identical(c(0,1),sort(unique(c(x))))) family <- rep("binomial",No) else family <- rep("",No)
    }
    if (length(family)==1)
    {
      family <- list(family)
      if (No > 1) for (i in 2:No) family[[i]] <- family[[1]]
    }
    if (length(family)!=No) stop("Length of family is not equal to number of outcome variables.")
    
    ## Output:
    Out <- list() 
    Out$graph <- matrix(0,Ni,No)
    Out$IC <- 0
    
    # Run glms:
    for (i in 1:No)
    {
      if (is.function(family[[i]])) fam <- family[[i]] else fam <- get(family[[i]])
      if (any(as.logical(adjacency[,i]))) 
      {
        tryres <- try(Res <- glm(x[-1,vars[i]] ~ x[-nrow(x),as.logical(adjacency[,i])],family=fam))
        if (is(tryres, 'try-error')) Res <- glm(x[-1,vars[i]] ~ NULL,family=fam)
      } else {
        Res <- glm(x[-1,vars[i]] ~ NULL,family=fam)
      }
      Out$graph[as.logical(adjacency[,i]),i] <- coef(Res)[-1]
      Out$IC <- Out$IC + icfun(Res,...)
    }
    Out$graph[is.na(Out$graph)] <- 0
    return(Out)
  }


computePCC <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- as.matrix(forceSymmetric(x))
  return(x)
}

computePDC <- function(beta,kappa){
  if (ncol(beta) == nrow(beta)+1){
    beta <- beta[,-1,drop=FALSE]
  }
  sigma <- solve(kappa)
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

graphicalVAR_test <- function (data, nLambda = 50, verbose = TRUE, gamma = 0.5, scale = TRUE, 
          lambda_beta, lambda_kappa, maxit.in = 100, maxit.out = 100, 
          deleteMissings = TRUE, penalize.diagonal = TRUE, lambda_min_kappa = 0.05, 
          lambda_min_beta = 0.05) 
{
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  stopifnot(is.matrix(data))
  Nvar <- ncol(data)
  Ntime <- nrow(data)
  #data <- scale(data, TRUE, scale)
  data_c <- data[-1, , drop = FALSE]
  data_l <- cbind(1, data[-nrow(data), , drop = FALSE])
  if (any(is.na(data_c)) || any(is.na(data_l))) {
    if (deleteMissings) {
      warnings("Data with missings deleted")
      missing <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 
        0
      data_c <- data_c[!missing, ]
      data_l <- data_l[!missing, ]
    }
    else {
      stop("Missing data not supported")
    }
  }
  if (missing(lambda_beta) | missing(lambda_kappa)) {
    lams <- generate_lambdas(data_l, data_c, nLambda, nLambda, 
                             lambda_min_kappa = lambda_min_kappa, lambda_min_beta = lambda_min_beta, 
                             penalize.diagonal = penalize.diagonal)
    if (missing(lambda_beta)) {
      lambda_beta <- lams$lambda_beta
    }
    if (missing(lambda_kappa)) {
      lambda_kappa <- lams$lambda_kappa
    }
  }
  Nlambda_beta <- length(lambda_beta)
  Nlambda_kappa <- length(lambda_kappa)
  lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
  Estimates <- vector("list", nrow(lambdas))
  if (verbose) {
    pb <- txtProgressBar(0, nrow(lambdas), style = 3)
  }
  for (i in seq_len(nrow(lambdas))) {
    if (lambdas$beta[i] == 0 & lambdas$kappa[i] == 0) {
      X <- data_l
      Y <- data_c
      nY <- ncol(Y)
      nX <- ncol(X)
      n <- nrow(X)
      beta <- t(Y) %*% X %*% solve(t(X) %*% X)
      S <- 1/(nrow(Y)) * (t(Y) %*% Y - t(Y) %*% X %*% 
                            t(beta) - beta %*% t(X) %*% Y + beta %*% t(X) %*% 
                            X %*% t(beta))
      S <- (S + t(S))/2
      if (any(eigen(S)$value < 0)) 
        warning("Residual covariances not postive definite")
      kappa <- solve(S)
      kappa <- (kappa + t(kappa))/2
      lik1 = determinant(kappa)$modulus[1]
      lik2 <- sum(diag(kappa %*% S))
      pdO = sum(sum(kappa[upper.tri(kappa, diag = FALSE)] != 
                      0))
      pdB = sum(sum(beta != 0))
      LLk <- (n/2) * (lik1 - lik2)
      LLk0 <- (n/2) * (-lik2)
      EBIC <- -2 * LLk + (log(n)) * (pdO + pdB) + (pdO + 
                                                     pdB) * 4 * gamma * log(2 * nY)
      Estimates[[i]] <- list(beta = beta, kappa = kappa, 
                             EBIC = EBIC)
    }
    else {
      tryres <- try(Rothmana(data_l, data_c, lambdas$beta[i], 
                             lambdas$kappa[i], gamma = gamma, maxit.in = maxit.in, 
                             maxit.out = maxit.out, penalize.diagonal = penalize.diagonal))
      if (is(tryres, "try-error")) {
        Estimates[[i]] <- list(beta = matrix(NA, Nvar, 
                                             Nvar + 1), kappa = matrix(NA, Nvar, Nvar), 
                               EBIC = Inf, error = tryres)
      }
      else {
        Estimates[[i]] <- tryres
      }
    }
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose) {
    close(pb)
  }
  lambdas$ebic <- sapply(Estimates, "[[", "EBIC")
  if (all(lambdas$ebic == Inf)) {
    stop("No model estimated without error")
  } else{}
  min <- which.min(lambdas$ebic)
  Results <- Estimates[[min]]
  Results$PCC <- computePCC(Results$kappa)
  Results$PDC <- computePDC(Results$beta, Results$kappa)
  Results$path <- lambdas
  Results$labels <- colnames(data)
  if (is.null(Results$labels)) {
    Results$labels <- paste0("V", seq_len(ncol(data)))
  }
  colnames(Results$beta) <- c("1", Results$labels)
  rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <- colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <- Results$labels
  Results$gamma <- gamma
  Results$allResults <- Estimates
  class(Results) <- "graphicalVAR"
  return(Results)
}