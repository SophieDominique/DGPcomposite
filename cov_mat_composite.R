#-------------------------------------------------------------------------------
# DGP : common-factor approach revised for composite models
# with options data generation
#-------------------------------------------------------------------------------

generateCOMP <- function(N, B, blocks, w, phi, R2 = NULL, DGP = 2, d_err = 0.01)
{
  
  # paramètres
  nbExo = length(which(rowSums(B) == 0))
  nbEndo = length(which(rowSums(B) != 0))
  nbLV = ncol(B)
  
  # A = (Gamma Beta)
  A_FIM = B[(nbExo+1):nbLV, , drop = FALSE]
  Rlv = phi
  
  #-------------------------------------------------------------------------------
  # covariance matrix implied by the inner model
  # FIM method
  
  for(j in 1:nbEndo)
  {
    if(!is.null(R2))
    {
      temp = A_FIM[j, 1:(nbExo + j - 1),drop=F]
      tau = sqrt(R2[j] / (temp%*%
                            Rlv[1:(nbExo + j - 1), 1:(nbExo + j - 1), drop = F]%*%
                            t(temp)))
      A_FIM[j, 1:(nbExo + j - 1)] = as.numeric(tau) * temp
    }
    temp2 = A_FIM[j, 1:(nbExo + j - 1), drop = FALSE]%*%Rlv
    Rlv = rbind(Rlv, temp2)
    Rlv = cbind(Rlv, c(temp2, 1))
  }
  
  #-------------------------------------------------------------------------------
  # covariance matrix implied by the model
  
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
  S = H%*%Rlv%*%t(H)
  
  coeff = rbind(matrix(0, nrow = nbExo, ncol = nbLV), A_FIM)
  
  #-------------------------------------------------------------------------------
  # indicators generation
  
  allP = p+q
  
  if(DGP == 1){
    D = diag(sqrt(svd(S)$d))
    C = svd(S)$u%*%D
    A = matrix(rnorm(N*allP),N,allP)
    X = A%*%t(C)
    
  } else if(DGP == 2){
    C = chol(Rlv)
    A = matrix(rnorm(N*nbLV),N,nbLV)
    Z = A%*%C
    X = Z%*%t(H)
    
  } else if(DGP == 3)
  {
    D = diag(sqrt(svd(S)$d))
    C = svd(S)$u%*%D
    orth_A = svd(scale(matrix(rnorm(N*allP),N,allP)))$u
    X = orth_A%*%t(C)
    
  }

  #-------------------------------------------------------------------------------
  return(list(S = S, data = X, B = coeff, Rlv = Rlv))
}


#-------------------------------------------------------------------------------
# # example
# # param
# N = 100
# p1 = p2 = p3 = p4 = p5 = p6 = 4
# allP = p1+p2+p3+p4+p5+p6
# w = rep(0.4, allP) # weights all at 0.4 or 0.9
# 
# nbCoeff = 9
# B = matrix(c(0,0,0,0,0,0,
#              0.5,0,0,0,0,0,
#              0.15,0.5,0,0,0,0,
#              0.15,0.5,0.3,0,0,0,
#              0,0,0,0.5,0,0,
#              0,0,0,0.5,0.15,0), 6, 6, byrow = TRUE)
# nbLV = ncol(B)
# nbExo = length(which(rowSums(B) == 0))
# nbEndo = nbLV - nbExo
# 
# phi = diag(1)
# 
# blocks = list(X1=1:p1,X2=(p1+1):(p1+p2),X3=(p1+p2+1):(p1+p2+p3), X4=(p1+p2+p3+1):(p1+p2+p3+p4),
#               X5=(p1+p2+p3+p4+1):(p1+p2+p3+p4+p5),X6=(p1+p2+p3+p4+p5+1):(p1+p2+p3+p4+p5+p6))
# 
# indicatorx <- c(rep(1,p1))
# indicatory <- c(rep(1,p2),rep(2,p3),rep(3,p4),rep(4,p5),rep(5,p6))
# #-------------------------------------------------------------------------------
# 
# dataCOMP1 = generateCOMP(N, B, blocks, w, phi, R2 = NULL, DGP = 1)
# #-------------------------------------------------------------------------------
# library(cbsem)
# resCOMP1_A = plspath(dataCOMP1$data, B = matrix(c(0,0,0,0,0,0,
#                                                  1,0,0,0,0,0,
#                                                  1,1,0,0,0,0,
#                                                  1,1,1,0,0,0,
#                                                  0,0,0,1,0,0,
#                                                  0,0,0,1,1,0), 6, 6, byrow=TRUE),
#                      indicatorx = indicatorx, indicatory = indicatory,
#                      modex = "A", modey = "A")
# #-------------------------------------------------------------------------------
# 
# as.vector(resCOMP1_A$Bhat[which(resCOMP1_A$Bhat!=0)])
# resCOMP1_A$R2
# 
# as.vector(B[which(B!=0)])
# 
# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
# out = gscmcov(B = B, indicatorx = indicatorx, indicatory = indicatory,
#               lambdax = NULL, lambday = NULL, wx = w[1:(p1)], 
#               wy = w[(p1+p2+1):allP], Sxixi = phi, R2 = NULL)
# C = chol(out$S)
# dataSchlitt = matrix(rnorm(N*allP),N,allP)%*%C
# 
# resCOMP1_A = plspath(dataSchlitt, B = matrix(c(0,0,0,0,0,0,
#                                                   1,0,0,0,0,0,
#                                                   1,1,0,0,0,0,
#                                                   1,1,1,0,0,0,
#                                                   0,0,0,1,0,0,
#                                                   0,0,0,1,1,0), 6, 6, byrow=TRUE),
#                      indicatorx = indicatorx, indicatory = indicatory,
#                      modex = "A", modey = "A")
# #-------------------------------------------------------------------------------
# 
# as.vector(resCOMP1_A$Bhat[which(resCOMP1_A$Bhat!=0)])
# resCOMP1_A$R2
# 
# as.vector(B[which(B!=0)])
