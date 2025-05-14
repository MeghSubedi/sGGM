
#' Function to run the Multiple Phenotype Association Test
#'
#' @param Zs The Matrix of Z scores (M by K) M is SNPs and K is traits
#' @param B Matrix of Community structure (K by L)
#' @param method The name of multiple phenotype test method. The default is ACAT.
#' @returns The Matrix of Pvalues, where each column is the Pvalues from each cluster
#' @import rlist
#'
#' @examples
#' #load the example data
#' data<-example_data
#' Cluster<-Get_Clusters(data)
#' pval_out<-Test_With_Network(Zs=data,B=Cluster$B_cor,method="ACAT")
#' head(pval_out)
#' @export

Test_With_Network<-function(Zs,B,method="ACAT"){
  # Ensure both Zs and B are matrices
  Zs<-as.matrix(Zs)
  B<-as.matrix(B)
  if(ncol(Zs)!=nrow(B)){
    stop("\n Mismatch between the Zs and cluster. Make sure they are correct. \n")
  }

  test_methods <- c("ACAT","Wald","PCFisher","HCLC","OBrien")

  # Normalize input (case-insensitive, partial match)
  method_norm <- toupper(method)
  method_index <- pmatch(method_norm, toupper(test_methods), nomatch = 0)

  if (method_index == 0) {
    stop(paste0("Error: 'method' must be one of: ", paste(test_methods, collapse = ", ")))
  }

  method_selected <- test_methods[method_index]

  # Make a matrix to store Pvalues

  Pval_list<-list()

  for(cols in 1:ncol(B)){
    pheno_group<-which(B[,cols]!=0)
    Z_Matrix<-Zs[,pheno_group]
    Z_Matrix<-as.matrix(Z_Matrix)
    # Apply multiple phenotype test
    pval <- switch(method_selected,
                   ACAT = ACAT_Pval(Z_Matrix),
                   Wald = Wald_Pval(Z_Matrix),
                   PCFisher = PCFisher_Pval(Z_Matrix),
                   HCLC = HCLC_Pval(Z_Matrix),
                   OBrien = OBrien_Pval(Z_Matrix))

    # Append Pvalues list
    Pval_list<-list.append(Pval_list,pval)
  }
  # Make a dara frame of Pvalues  where columns are Modules
  df<-do.call("cbind",Pval_list)
  colnames(df)<-paste0("Module ",1:ncol(B))
  return(df)
}






#' Run an association test of SNPs with multiple phenotype without considering network modules
#'
#' @param Zs  The matrix of Z scores; rows are SNPs and columns are traits
#' @param method method is the test method, There are five methods to check and default is ACAT
#' We plan to include other methods in the future
#'
#' @returns A numeric vector of pvalues of length equal to number of SNPs
#' @export
#'
#' @examples
#' #load the example data
#' data<-example_data
#' Cluster<-Get_Clusters(data)
#' pval_out<-Test_With_Network(Zs=data,B=Cluster$B_cor,method="ACAT")
#' head(pval_out)
#'
Test_Without_Network<-function(Zs,method="ACAT"){

  Zs<-as.matrix(Zs)

  test_methods <- c("ACAT","Wald","PCFisher","HCLC","OBrien")

  # Normalize input (case-insensitive, partial match)
  method_norm <- toupper(method)
  method_index <- pmatch(method_norm, toupper(test_methods), nomatch = 0)

  if (method_index == 0) {
    stop(paste0("Error: 'method' must be one of: ", paste(test_methods, collapse = ", ")))
  }

  method_selected <- test_methods[method_index]
  # Apply multiple phenotype test
  pval <- switch(method_selected,
                 ACAT = ACAT_Pval(Zs),
                 Wald = Wald_Pval(Zs),
                 PCFisher = PCFisher_Pval(Zs),
                 HCLC = HCLC_Pval(Zs),
                 OBrien = OBrien_Pval(Zs))
  return(pval)
}
