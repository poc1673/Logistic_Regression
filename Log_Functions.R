
# Binomial Regression -----------------------------------------------------
# Implementation of binary logistic regression:
# Inputs:
# [1]  y - a column vector of N elements containing the actual classification results.
# [2]  X - a (p+1) by N design matrix.
# [3]  beta - a p-length column vector of coeffients.

# Outputs:
# [1]  f - The log-likelihood function evaluated using the inputs.
# [2]  df - p-length column vector representing the log-likelihood's gradient.
# [3]  ddf - a p by p Hessian matrix for the log-likelihood function.

L_L <- function(y,beta,X){
  # Calculate the vector of probabilities:
  pi <- exp(X%*%beta)/(1+exp(X%*%beta))
  # print(  X%*%beta)
  f <- prod(ifelse(y==1,pi,1-pi))
  # First derivative:
  df <- t(t(y- pi)%*%X )
  # Second derivative:
  W <- diag(as.vector(pi*(1-pi)))
  ddf <- -t(X)%*%W%*%X
  # Output:
  # output_list <- list(f = f, df = df)
  output_list <- list(f = f, df = df,ddf = ddf)
  return(output_list)
}


# Simple implementation of Newton-Raphson:
# Inputs:
# [1] func - The function to be optimized by the algorithm.
# [2] initial_values - An initial guess to be passed to "func".
# [3] tol - The stopping criteria to the model.  For this implementation, the tolerance is compared to the 2-norm
#           of the gradient:

# Outputs:
# [1]  params - the output of the logistic regression.

newton <- function(func, initial_values,tol = 1e-16){
  params <- initial_values
  check <-1
  while(check > tol){
    func_eval <- func(params)
    # print(func_eval)
    params <- params -  solve(func_eval$ddf)%*%func_eval$df
    check <- sqrt(t(func_eval$df)%*%func_eval$df)
    print(check)
  }
  return(params)
}

# Simple implementation of Gradient Descent:
# Inputs:
# [1] func - The function to be optimized by the algorithm.
# [2] initial_values - An initial guess to be passed to "func".
# [3] tol - The stopping criteria to the model.  For this implementation, the tolerance is compared to the 2-norm
#           of the gradient:

# Outputs:
# [1]  params - the output of the logistic regression.

grad_desc <- function(b0,fun,tol = 1e-10,step = .1){
  cont <- 1
  
  while(cont > tol){
    f_b1 <-  fun(b0)
    b1 <- b0 + step*f_b1$df
    b0 <- b1
    cont <- as.numeric(sqrt(t(f_b1$df)%*%f_b1$df))
    # print(paste(f_b1$df[1],"  ",f_b1$df[2],"  ",f_b1$df[3],sep = ""))
    print(cont)
  }
  return(list("res" = b1,"func_val" = fun(b1)))
}



# Implementation of multi-class regression --------------------------------

# I implemented the multi-class version of the probabiliity function to produce a matrix of the class probabilities.  Let K be the number of classes.
# Inputs:
# [1]  X - A N*(K-1) by M design matrix.
# [2]  init_betas - A N*(K-1) length design matrix.
# [3]  classes - The number of classes we are predicting.

# Outputs:
# A K by M matrix containing the class probabilities calculated as:

find_pi_multi <- function(X,beta,classes){
  vars <- ncol(X)
  prob_mat <- matrix(rep(0,nrow(X)*(classes-1)),nrow = nrow(X),ncol = (classes-1))
  for(k in 1:(classes-1)){
    rel_indices <- 1:vars+(k-1)*vars
    prob_mat[,k] <- exp(X%*%beta[rel_indices])
  }
  denom <- apply(X = prob_mat,MARGIN = 1,FUN = sum)+1
  prob_mat2 <- cbind(prob_mat,as.matrix(rep(1,nrow(X))))/denom
  return(prob_mat2)
}


# I'm sure R has a better way to form a block matrix.  Pardon my MATLABese here.
# Checks out with arbitrary number of classes,  Should workdk  When doing the implementation try to use the dimensions of the matrix of
# p probabilities for the class size.
X_tilde <- function(X,classes){
  r <- nrow(X)
  N <- r*(classes -1)
  X_tilde <- matrix(rep(0,N*(classes-1)*ncol(X)),nrow = N, ncol = (classes-1)*ncol(X) )
  m <- 1
  n <- 1
  
  for(i in 1:(classes-1)){
    # block_diags[[i]] <- X
    X_tilde[n:(n+r-1),m:(m+(ncol(X)-1)  )] <- X
    n <- n + r
    m <- m + (ncol(X)) 
    n:(n+r-1)
    m:(m+(ncol(X)-1)  )
  }
  return(X_tilde)
}






# Some notes:
# Clean up the implementation of the transposition section of this matrix.  Doing the off diagonals first and then transposing makes more sense.
form_W <- function(p){
  # W <-
  # I'm not too worried about speed with this seeing as I am just forming block matrices.
  # Form the block diagonals in the W matrix:
  N <- nrow(p)*(ncol(p)  -1)
  r <- nrow(p)
  W <- matrix(rep(0,N^2),nrow = N, ncol = N )
  block_diags <- list()
  upper_diag <- list() 
  m <- 1
  n <- 1
  for(i in 1:(ncol(p)-1)){
    block_diags[[i]] <- diag((p[,i])*(1-p[,i]))
    W[n:(n+r-1),m:(m+r-1)] <- block_diags[[i]]
    n <- n +r
    m <- m + r
  }
  
  # Form the off-diagonals:
  # W_ii = -p_k*p_m
  
  k <- 1
  m <- 1
  n <- 1
  
      if(ncol(p)>2){
      for(i in 1:(ncol(p)-2)){ 
        # print(W)
         # print(paste("i is ",i))
         for(j in  c((i+1):(ncol(p)-1)))  {
           m <- r*(j-1)+1
           upper_diag[[k]] <- diag(-p[,i]*p[,j])
           W[(n:(n+r-1)),m:(m+r-1)] <- diag(-p[,i]*p[,j])
           k <- k+1
         }

        
        
      }
    W <- W + t(W)
    diag(W) <- diag(W)/2
  }
  
  # The block matrices should be the same across the diagonal so I'm just adding the transpose
  # to the original W and correcting the main diagonal.
  
  return("W" = W) 
}   




# Generate a N(K-1) length vector of indicator functions based on class.

form_y <- function(classes,y, obs){
  y0 <- as.matrix(rep(0,obs))
  
  for(i in 1:(length(classes)-1)){
    # print(i)
    if(i==1){yp <- y0
    yp[which(y==classes[i])] <- 1
    }else{y_temp <- y0
    y_temp[which(y==classes[i])] <- 1
    yp <- rbind(yp,y_temp)
    }
  }
  return(yp)
}

# Multi-class Regression -----------------------------------------------------
# Implementation of binary logistic regression:
# Inputs:
# [1]  y - a column vector of N elements containing the actual classification results.
# [2]  X - a (p+1) by N design matrix.
# [3]  beta - a p-length column vector of coeffients.

# Outputs:
# [1]  f - The log-likelihood function evaluated using the inputs.
# [2]  df - p-length column vector representing the log-likelihood's gradient.
# [3]  ddf - a p by p Hessian matrix for the log-likelihood function.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LL_Multi <- function(Y,X,beta){
  class_types <- unique(Y)[order(-unique(Y))]
  classes <- length(class_types)
  class_probs <- as.matrix(rep(0,nrow(X))) 
  p <- find_pi_multi(X,beta,classes)
  # Get the esimated probability of the observed class.
  p2 <- as.matrix(unlist(as.data.frame(p[,1:(dim(p)[2]-1)])))
  y <- form_y(classes = class_types,y = Y,obs = length(Y))
  
  for(k in 1:length(class_types)){
    class_probs[which(Y==class_types[k])] <- p[which(Y==class_types[k]),k]
  }
  XT <- X_tilde(X,ncol(p))
  ws <- form_W(p = p)
  f <- sum(class_probs)
  # Calculate the gradient:
  df <- t(XT)%*%(y-p2)
  # Calculate the Hessian:
  ddf <- -t(XT)%*%ws%*%XT
  return(list("f" = f,"df" = df, "ddf" = ddf))
  # return(list("f" = f,"df" = df))
}
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
