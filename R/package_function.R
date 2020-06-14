################

#' This function is to corect the bias of cross-trait PRS according to Theorem 1 and 4
#'
#' 
#' @param raw_estimator the raw genetic correlation estimator in PRS
#' @param n_train the sample size of training GWAS
#' @param p_indep the number of independent genetic variants
#' @param h2_trait1 heritability of the trait in training GWAS
#' @param h2_trait2 heritability of the trait in tetsing GWAS
#' @param n_overlap  the number of overlaped samples between training and testing GWAS, default=0
#' @param n_test the sample size of testing GWAS (only useful when n_overlap >0,default=0)
#' @keywords bc_prs
#' @return bias_corrected_estimator: raw bias-corrected genetic correlation estimator
#' @examples
#' 
#' ##toy example 
#' bc_prs(raw_estimator=0.1,n_train=50000,p_indep=500000,
#' h2_trait1=0.5,h2_trait2=0.5,n_overlap=0,n_test=0)
#' @export
bc_prs<-function(raw_estimator,n_train,p_indep,h2_trait1,h2_trait2,n_overlap=0,n_test){
  ####read in ####
  raw_estimator<-as.numeric(raw_estimator)
  n_train<-as.numeric(n_train)
  p_indep<-as.numeric(p_indep)
  h2_trait1<-as.numeric(h2_trait1)
  h2_trait2<-as.numeric(h2_trait2)
  n_overlap<-as.numeric(n_overlap)
  n_test<-as.numeric(n_test)
  if (any(is.na(c(raw_estimator,n_train,p_indep,h2_trait1,h2_trait2,n_overlap,n_test)))) stop("Some data are missing, please check.")
  if (min(c(n_train,p_indep,h2_trait1,h2_trait2,n_overlap,n_test))<0) stop("Some data are negative, please check.")
  if (n_overlap>min(n_train,n_test)) stop("n_overlap is larger than min(n_train,n_test), please check.")
  ###no overlaps
  if(n_overlap==0){
    corrected_estimator<-raw_estimator*sqrt((n_train+p_indep/h2_trait1)/(n_train*h2_trait2))
  }
  ###there are overlaps
  if(n_overlap>0){
    part1<-(1+n_overlap*p_indep/((n_train)*(n_test)*1))*sqrt(h2_trait2)
    part2<-1+p_indep/(n_train*h2_trait1)+2*n_overlap*p_indep/((n_train)*(n_test))+n_overlap*p_indep^2/((n_train)^2*(n_test)*h2_trait1)
    factor<-sqrt(part2)/part1
    ###
    corrected_estimator<-raw_estimator*factor
  }
  if(abs(corrected_estimator)>1.5){
    cat("Warning: the corrected estimator might be too big; you may consider a smaller R^2 threshold in LD-based prunning","\n")
  }
  result <- list(corrected_estimator=corrected_estimator)
  return(result)
}





#' This function is to corect the bias of cross-trait PRS according to Theorem 2,3 and Proposition S4
#'
#' 
#' @param raw_estimator the raw genetic correlation estimator in PRS
#' @param n_train1 the sample size of training GWAS1
#' @param n_train2 the sample size of training GWAS2
#' @param p_indep the number of independent genetic variants
#' @param h2_trait1 heritability of the trait in training GWAS1
#' @param h2_trait2 heritability of the trait in training GWAS2
#' @param n_overlap the number of overlaped samples between the two training GWAS, default=0
#' @keywords bc_prs2
#' @return bias_corrected_estimator: raw bias-corrected genetic correlation estimator
#' @examples
#' 
#' ##toy example 
#' bc_prs2(raw_estimator=0.1,n_train1=100000,
#' n_train2=100000,p_indep=300000,h2_trait1=0.5,h2_trait2=0.5,n_overlap=0)
#' @export
bc_prs2<-function(raw_estimator,n_train1,n_train2,p_indep,h2_trait1,h2_trait2,n_overlap=0){
  ####read in ####
  raw_estimator<-as.numeric(raw_estimator)
  n_train1<-as.numeric(n_train1)
  n_train2<-as.numeric(n_train2)
  p_indep<-as.numeric(p_indep)
  h2_trait1<-as.numeric(h2_trait1)
  h2_trait2<-as.numeric(h2_trait2)
  n_overlap<-as.numeric(n_overlap)
  if (any(is.na(c(raw_estimator,n_train1,n_train2,p_indep,h2_trait1,h2_trait2,n_overlap)))) stop("Some data are missing, please check.")
  if (min(c(n_train1,n_train2,p_indep,h2_trait1,h2_trait2,n_overlap))<0) stop("Some data are negative, please check.")
  if (n_overlap>min(n_train1,n_train2)) stop("n_overlap is larger than min(n_train1,n_train2), please check.")
  ###no overlaps
  if(n_overlap==0){
    corrected_estimator<-raw_estimator*sqrt((n_train1+p_indep/h2_trait1)*(n_train2+p_indep/h2_trait2)/(n_train1*n_train2))
  }
  ###there are overlaps
  if(n_overlap>0){
    part1<-sqrt(n_train1*n_train2)+n_overlap*p_indep/(sqrt(n_train1*n_train2)*1)
    part2<-sqrt((n_train1+p_indep/h2_trait1)*(n_train2+p_indep/h2_trait2))
    factor<-part2/part1
    ###
    corrected_estimator<-raw_estimator*factor
  }
  if(abs(corrected_estimator)>1.5){
    cat("Warning: the corrected estimator might be too big; you may consider a smaller R^2 threshold in LD-based prunning","\n")
  }
  result <- list(corrected_estimator=corrected_estimator)
  return(result)
}
