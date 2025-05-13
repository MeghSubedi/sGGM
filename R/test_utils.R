
#' Cauchy combination of Pvalues from individual association tests.
#' @param Zs The matrix of Z scores (rows are SNPs and column are traits)
#' @returns A numeric vector of Pvalues rowwise(ie SNP)
#' @export
#'
#' @examples
#' \dontrun{
#' Zs<-sapply(1:5,function(x){rnorm(10)})
#' ACAT_Pval(Zs)}
#' @importFrom ACAT ACAT
#' @importFrom stats pchisq cor hclust cutree as.dist
#' @importFrom dplyr select
#' @importFrom stats lag
#' @importFrom MASS ginv
#' @importFrom dplyr filter

ACAT_Pval<-function(Zs){
  # Convert Zscores to pvalues
  pValues<-pchisq(Zs^2,df=1,lower.tail = FALSE)
  pValues<-as.matrix(pValues)
  K=ncol(pValues)
  if(K==1){
    return(as.vector(pValues))
  }
  else{
    # Add or subtract small values if P value is 0 or 1
    pValues[pValues==1]<-1-1e-16
    pValues[pValues==0]<-1e-16
    compute_ACAT<-apply(pValues,1,function(x){ACAT::ACAT(x)})
    return(compute_ACAT)
  }
}

#------------------------
# Wald
#-----------------------

#' Wald Test for testing association between SNPs and multiple traits.
#'
#' Computes the Wald test p-values for rows of Z-scores matrix, accounting for correlation among traits.
#' @param Zs The matrix of Z scores (rows are SNPs and column are traits)
#'
#' @returns A numeric vector vector of pvalues rowwise (per SNP)
#' @export
#'
#' @examples \dontrun{Zs<-sapply(1:5,function(x){rnorm(10)})
#' Wald_Pval(Zs)}

Wald_Pval<-function(Zs){
  Zs<-as.matrix(Zs)
  K=ncol(Zs)
  if(K==1){
    pval<-pchisq(Zs^2,df=1,lower.tail = FALSE)
    return(as.vector(pval))
  }
  else{
    R<-cor(Zs)
    Rinv<-MASS::ginv(R)
    W_stat<-apply(Zs,1,function(x){t(x)%*%Rinv%*%x})
    k=dim(Zs)[2]
    pval<-pchisq(W_stat,k,lower.tail = FALSE)
    return(pval)
  }
}


#' PCFisher test for testing an association between SNPs and multiple correlated traits
#' Computes the PCFisher test p-values for rows of Z-scores matrix, accounting for correlation among traits.
#' @param Zs a Matrix of Z scores (rows are SNPs and columns are traits)
#' @returns A numeric vector of pvalues for each row (i.e SNPs)
#' @export
#'
#' @examples
#' \dontrun{Zs<-sapply(1:5,function(x){rnorm(10)})
#' PCFisher_Pval(Zs)}


PCFisher_Pval<-function(Zs){
  Zs<-as.matrix(Zs)
  K=ncol(Zs)
  if(K==1){
    pval<-pchisq(Zs^2,1,lower.tail = FALSE)
    return(as.vector(pval))
  }
  else{
    R=cor(Zs)
    k=length(Zs)
    lambdas<-eigen(R)$values
    PCs<-eigen(R)$vectors
    # Function to compute pvalues rowwise
    compute_pval<-function(z_vector){
      lc_pc<-(t(PCs)%*%z_vector)
      PC_stat<-(lc_pc)^2/lambdas
      PC_p<-pchisq(PC_stat,1,lower.tail = FALSE)
      if (sum(PC_p == 0) != 0){
        PC_p <- PC_p[-which(PC_p == 0)]
      }
      PC_fisher<-sum((-2)*log(PC_p))
      PC_pval<-pchisq(PC_fisher,df=2*K,lower.tail = FALSE)
      return(PC_pval)
    }
    pval<-apply(Zs,1,function(x){compute_pval(x)})
    return(pval)
  }
}


#' HCLC Test for testing an association between SNPs and multiple correlated traits.
#' Computes the HCLC test p-values for rows of Z-scores matrix, accounting for correlation among traits.
#' @param Zs The matrix of Z scores (rows are SNPs and column are traits)
#' @returns A numeric vector of pvalues row-wise (SNP wise)
#' @export
#' @examples \dontrun{Zs<-sapply(1:5,function(x){rnorm(10)})
#' HCLC_Pval(Zs)}

HCLC_Pval<-function(Zs){
  K=ncol(as.matrix(Zs))
  if (K == 1) {
    pv0 <- pchisq(Zs^2, 1, lower.tail = FALSE)
    return(pv0)
  }

  else if(K==2){
    Sigma<-cor(Zs)
    W <- MASS::ginv(Sigma)
    compute_pval<-function(z_vector){
      CLC=t(z_vector)%*%W%*%z_vector
      pv0=pchisq(CLC,2,lower.tail = FALSE)
      return(pv0)
    }
    Pval<-apply(Zs,1,function(x){compute_pval(x)})
    return(Pval)
  }

  else{
    Sigma=cor(Zs)
    dist<-1-Sigma
    hc=hclust(as.dist(1-Sigma))
    H <- which.max(diff(hc$height)) + 1
    index <- cutree(hc, H)
    L <- max(index)
    B=sapply(1:L, function(t) as.numeric(index==t))
    W=t(B)%*%ginv(Sigma)
    U=t(W)%*%ginv(W%*%Sigma%*%t(W))%*%W
    # Function to compute Pvalues
    compute_pval<-function(z_vector){
      CLC=t(z_vector)%*%U%*%z_vector
      pv0=pchisq(CLC,L,lower.tail = FALSE)
      return(pv0)
    }
    Pval<-apply(Zs,1,function(x){compute_pval(x)})
    return(Pval)
  }

}


#' OBrien test for testing an association between SNPs and multiple correlated traits.
#' Computes the OBrien test p-values for rows of Z-scores matrix , accounting for correlation among traits.
#' @param Zs The matrix of Z scores (rows are SNPs and column are traits)
#' @returns A numeric vector of pvalues row-wise (SNP)
#' @export
#' @examples
#' \dontrun{Zs<-sapply(1:5,function(x){rnorm(10)})
#' OBrien_Pval(Zs)}

OBrien_Pval<-function(Zs){
  # Ensure Zs is a matrix
  Zs<-as.matrix(Zs)
  K=ncol(Zs)
  if(K==1){
    pv <- pchisq(Zs^2, 1, lower.tail = FALSE)
    return(as.vector(pv))
  }
  else{
    Sigma<-cor(Zs)
    compute_pval<-function(z_vector){
      k=length(z_vector)
      one_vec <- rep(1, k)
      Sigma_inv <- MASS::ginv(Sigma)
      num <- t(one_vec) %*% Sigma_inv %*% z_vector
      denom <- sqrt(t(one_vec) %*% Sigma_inv %*% one_vec)
      U1<-(num/denom)^2
      pv <- pchisq(U1, 1, lower.tail = FALSE)
      return(pv)
    }
    pval<-apply(Zs,1,function(x){compute_pval(x)})
    return(pval)
  }
}



