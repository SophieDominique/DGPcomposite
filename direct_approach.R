#-------------------------------------------------------------------------------
# DGP : direct approach with generated proxies
#-------------------------------------------------------------------------------

generatePLS <- function(N, B, blocks, w, phi = NULL, R2 = NULL)
{
  library(mvtnorm)
  
  # check param
  nbLV = length(blocks)
  nbMV = rep(0, nbLV+1)
  
  for(b in 1:nbLV)
  {
    nbMV[b+1] = length(blocks[[b]])
  }
 
  stopifnot(length(w) == sum(nbMV))
  # {
  #   print(paste("the nb of weights must be equal to the nb of indicators",err))
  #   return()
  # }
  
  #-------------------------------------------------------------------------------
  # generating exogenous proxies
  
  isExo = which(rowSums(abs(B)) == 0)
  nbExo = length(isExo)
  
  if(is.null(phi))
  {
    phi = diag(nbExo)
  }
  
  if(is.null(phi))
  {
    Zexo = apply(matrix(rnorm(N*nbExo), N, nbExo),2,scale)
  } else {
    Zexo = rmvnorm(N, mean = rep(0, nbExo), sigma = phi,
                   method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)
  }
  
  Zexo = apply(Zexo,2,scale)
  Z = Zexo
  allErr = c()
  
  #-------------------------------------------------------------------------------
  # generating endogenous proxies
  
  isEndo = which(rowSums(abs(B)) != 0)
  nbEndo = length(isEndo)
  
  for(j in 1:nbEndo)
  {
    ind = isEndo[j]
    isPred = which(B[ind,] != 0)
    beta = matrix(B[ind,isPred],ncol = 1)
    
    zPred = Z[,isPred, drop = FALSE]
    zzold = zPred%*%beta
    # print(var(zzold))
    
    if(!is.null(R2))
    {
      newbeta = beta * as.numeric(sqrt(R2[j]/var(zzold)))
      zPred = Z[,isPred, drop = FALSE]
      zzold = zPred%*%newbeta
      # print(var(zzold))
    } else {
      if(var(zzold) > 1)
      {
        message("Coefficient constraints not respected")
        return()
      }
    }
    
    if(ncol(Z) > 1)
    {
      U = svd(Z)$u
      errZ = (diag(N) - U%*%t(U))%*%rnorm(N)
      errZ = scale(errZ)
      if(length(isPred) > 1)
      {
        z_chap = zzold + as.vector(sqrt(1-(t(zzold)%*%zzold/(N-1))))*errZ
      } else {
        z_chap = zzold + as.numeric(sqrt(1-beta^2))*errZ
      }
      
    } else{
      errZ = (diag(N) - (1/(N-1))*Z%*%t(Z))%*%rnorm(N)
      errZ = scale(errZ)
      
      if(length(isPred) > 1)
      {
        z_chap = zzold + as.vector(sqrt(1-(t(zzold)%*%zzold/(N-1))))*errZ
      } else {
        z_chap = zzold + as.numeric(sqrt(1-beta^2))*errZ
      }
    }
    
    allErr = cbind(allErr, errZ)
    Z = cbind(Z,z_chap)
  }
  
  #-------------------------------------------------------------------------------
  # generating indicators
  
  indicatorx = NULL
  for(i in 1:nbExo)
  {
    indicatorx = c(indicatorx, rep(i,length(blocks[[i]])))
  }
  p = length(indicatorx)
  Wx = matrix(indicatorx,p,nbExo)
  Wx <- 1*(Wx == matrix(c(1:nbExo),p,nbExo,byrow=T))
  Wx <- Wx*w[1:p]
  
  indicatory = NULL
  for(i in 1:nbEndo)
  {
    indicatory = c(indicatory, rep(i,length(blocks[[nbExo + i]])))
  }
  q = length(indicatory)
  Wy = matrix(indicatory,q,nbEndo)
  Wy <- 1*(Wy == matrix(c(1:nbEndo),q,nbEndo,byrow=T))
  Wy <- Wy*w[(p+1):(p+q)]
  
  W = cbind(rbind(Wx,matrix(0,nrow = q, ncol = nbExo)),
            rbind(matrix(0, nrow = p, ncol = nbEndo),Wy))
  
  H = W%*%solve(t(W)%*%W)
  
  X = Z%*%t(H)
  
  #-------------------------------------------------------------------------------
  return(list(B = B, scores = Z, data = X, zeta = allErr))
}





#-------------------------------------------------------------------------------
# example
# n=100
# # variables of each block of expected data
# p1 = 3
# p2 = 4
# p3 = 5
# p4 = 3
# # p5 = 6
# # p6 = 4
# blocks=list(X1=1:p1,X2=(p1+1):(p1+p2),X3=(p1+p2+1):(p1+p2+p3),
#             X4=(p1+p2+p3+1):(p1+p2+p3+p4))
# #            X5=(p1+p2+p3+p4+1):(p1+p2+p3+p4+p5),X6=(p1+p2+p3+p4+p5+1):(p1+p2+p3+p4+p5+p6))
# 
# B = matrix(c(0,0,0,0,
# 0,0,0,0,
# 0.6,0,0,0,
# 0,0.3,0.4,0),4,4, byrow = TRUE)
# 
# allP = p1+p2+p3+p4
# weight = rep(0.9,allP)
# phi = matrix(c(1,0,
#                0,1),2,2,byrow =TRUE)
# #-------------------------------------------------------------------------------
# 
# test = generatePLS(N = n, B = B, blocks = blocks, w = weight, phi = phi,
#                    R2 = NULL, d_err = 0.1)
# 
# # test = generatePLS(N = n, B = B, blocks = blocks, w = weight, phi = phi,
# #                    R2 = c(0.5, 0.4), d_err = 0)
# #-------------------------------------------------------------------------------
# library(matrixcalc)
# matrix.rank(cov(test$data))
# 
# indicatorx <- c(rep(1,p1),rep(2,p2))
# indicatory <- c(rep(1,p3),rep(2,p4))
# 
# library(cbsem)
# 
# resPLS = plspath(test$data, B = matrix(c(0,0,0,0,
#                                          0,0,0,0,
#                                          1,0,0,0,
#                                          0,1,1,0),4,4, byrow = TRUE),
#                  indicatorx = indicatorx, indicatory = indicatory,
#                  modex = "A", modey = "A")
# #-------------------------------------------------------------------------------
# 
# as.vector(resPLS$Bhat[which(resPLS$Bhat!=0)])
# resPLS$w
# resPLS$lambdahat
# # resPLS$R2
# 
# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------





