
2023.07.20 (목)
library(rTensor)   # 35p
tnsr <- rand_tensor(c(4,4,4,4))
tuckerD <- tucker(tnsr,ranks=c(2,2,2,2))
tuckerD$conv
tuckerD$Z
tuckerD$U

tnsr <- rand_tensor(c(4,4,4))
tuckerD <- tucker(tnsr,ranks=c(2,2,2))
tuckerD$Z
tuckerD$U


2023.07.16(일)
Zoom tucker matlab 코드를 R로 열어보자.

library(R.matlab)
path <- system.file("mat-files", package = "R.matlab")
path1 <- system.file("C://user/download/zoomtucker_code/demo", package = "R.matlab")

pathname <- file.path(path1, "run_demo_zoom_sample_preprocessing.mat")
data <- readMat("run_demo_zoom_sample_preprocessing.mat")
data <- readMat(pathname)

print(data)

# 안됨...matlab..

#다시 TRES + ManifoldOptim..에 
# envelope 에 regularization 있는 부분...코드실행

library(
  
  
  20230702 (일)
  # TRES() 패키지에 square() data x tensor data만 저장하고 싶으면
  
  library(TRES)
  data("square")
  str(square)
  head (square$x@data)  # tensor X data...
  
  
  #######################################
  # TPR model with envelope structure  (Square data)
  # Under the TPR model, we use cross-validation to select dimension based on the prediction mean squared error
  # X : 32*32*200 # 2-way tensor (matrix) predictor X
  # Y : 1*200  (# scalar response)
  # coef : 3-way tensor with dim 32*32*1  # 0.1 or 1 (=B)
  # Gamma : two basis matrices, Gamma1, Gamma2 with dim 32*2  #(SVD of B)
  
  data("square")
  fit_ols2 <- TPR.fit(square$x, square$y, method='standard')
  fit_1d2 <- TPR.fit(square$x, square$y, u=c(2,2), method='1D')
  fit_pls2 <- TPR.fit(square$x, square$y, u=c(2,2), method='PLS')
  fit_1d2
  coef(fit_1d2)
  dim( fitted(fit_1d2) )
  dim ( residuals(fit_1d2) )
  
  plot(fit_ols2)
  plot(fit_1d2)
  plot(fit_pls2)
  
  Stroke (뇌졸중), Brain Ischemia (뇌경색),Brain Hemorrhage  stroke, 
  
  230614 (수)
  
  # https://www.r-bloggers.com/2017/04/understanding-the-tucker-decomposition-and-# compressing-tensor-valued-data-with-r-code/   
  # Alexej Gossmann
  #Understanding the Tucker decomposition, and compressing tensor-valued data (with R code) 
  
  library(rTensor)
  tnsr <- as.tensor(array(1:12, dim = c(2, 2, 3)))
  tnsr 
  mat <- matrix(1:6, 3, 2)
  # 1-mode product performed via the function ttm in rTensor
  tnsr_times_mat <- ttm(tnsr = tnsr, mat = mat, m = 1)
  tnsr_times_mat
  tnsr_times_mat [ , ,1]   # vs. tnsr_times_mat@data [ , ,1]
  mat %*% as.matrix (tnsr@data [ , ,1])  # 왠일..@표시하면 data만 뽑히네....헐
  
  ####   Higher-order SVD (HOSVD) ####
  #how do you compute the Tucker decomposition?
  
  tnsr <- as.tensor(array(1:12, dim = c(2, 3, 2)))
  tnsr@data
  # mode-1 unfolding:
  k_unfold(tnsr, 1)@data
  
  # mode-2 unfolding:
  k_unfold(tnsr, 2)@data
  
  # mode-3 unfolding:
  k_unfold(tnsr, 3)@data
  
  tnsr <- rand_tensor(modes = c(30, 40, 50))
  hosv_decomp <- hosvd(tnsr)   # 6*40*50
  
  hosv_decomp$Z  # our matrix $G$  dim : 30*40*50
  hosv_decomp$U #list containing all the matrices $A^{(k)}$, 
  # dim [[1]] 는 30*30 [[2]]는 40*40 [[3]] 50*50
  
  ###############################
  use the function ttl, which performs multiple k-mode products on multiple modes successively given a tensor and a list of matrices, to check that up to numerical error the equation 
  X = G x_1 A1 x_2 A2 x_3... x_N  AN
  is satisfied:
    ################################
  HOSVD_prod <- ttl(hosv_decomp$Z, hosv_decomp$U, 1:3)
  error <- tnsr - HOSVD_prod
  table(abs(error@data) < 1e-12)
  
  # Neuroimaging 20230605 Youtube....Elizabeth Sweeney
  source("http://neuroconductor.org/neurocLite.R")
  # Image analysis and processing wiht R ....Dr. Bharatendra Rai youtube
  
  install.packages("BiocManager") 
  BiocManager::install("EBImage")
  
  biocLite()
  biocLite("EBImage")
  
  # Read Image file
  library(EBImage)
  
  
  # Neuroimaging 20230604 in R
  # 
  library(oro.nifti)
  
  ############################
  # Read in the baseline data
  ############################
  # Baseline data
  MPRAGE_base <- readNIfTI('SUBJ0001-02-MPRAGE.nii.gz', reorient =T) 
  # reorient : put the image into a standard space when you read them in.
  # 안하면 tons of different way로 표현될 것 . 항상 reorient=T로 !!!
  
  dim (MPRAGE_base ) 
  MPRAGE_base
  slotNames(MPRAGE_base)
  data <- slot(MPRAGE_base, ".Data")    
  # python에서는 t1_img = nib.load("/ nii ")이면 t1_imag.get_fdata() 이미지 숫자 자료 읽는법
  
  dim(MPRAGE_base [, , 170] )  # 170* 256
  dim(MPRAGE_base [, , 256] )  # 170* 256
  dim(MPRAGE_base [170, , ] )  # 256* 256
  
  #############################
  # Visualize the data 
  ##############################
  
  ## axials slice ###
  image(MPRAGE_base [, , 128])  # T1 
  
  ## coronal slice ###
  image(MPRAGE_base [, 128, ], col = rainbow(12) )
  
  ## sagittal slice ##
  image(MPRAGE_base [85, , ], col = topo.colors(12) )
  
  ###########################################
  # Using orthographic function (From oro.nifti package)
  ###########################################
  
  orthographic(MPRAGE_base ) 
  orthographic(MPRAGE_base,  xyz =c( 90,100,15) )   # z가 way down here.
  
  
  ############################################
  # Calculations on the nifti object 
  ############################################
  
  mean(MPRAGE_base)
  
  sd(MPRAGE_base)
  
  min(MPRAGE_base)
  
  max( MPRAGE_base) 
  
  ###############################################3
  ## Plot baseline and followup data
  #################################################
  
  library(fslr)
  
  MPRAGE_follow = readNIfTI ( 'SUBJ0001-02-MPRAGE.nii.gz', reorient =T)
  # plot the data together
  
  double_ortho(MPRAGE_base , MPRAGE_follow)
  
  ## plot the difference between the baseline and follow up data
  
  MPRAGE_diff <- MPRAGE_base - MPRAGE_follow 
  orthographic (MPRAGE_diff )
  
  ###########################################
  ## PRocess the data using fslr 
  ##########################################
  # set path to FSL
  
  options(fsl.path = =“/path/to/fsl/”)
  options(fsl.path= "C:/User/Gram/Download/fslinstaller")
  
  ############################################
  ## Inhomogeneity corretion woiht fslr 'fsk_biascorrect' 
  ############################################
  
  ## this takes awhile, so I went ahead and did it for you!
  ## MPRAGE_base_bias_corrected <- fsl_biascorrect(MPRAGE_base )  # 10분정도가 걸린다.
  
  ## writeNIfTI( MPRAGE_base_bias_corrected, filename = '   
  
  MPRAGE_base_bias_corrected <- readNIfTI( 'SUBJ0001-01-MPRAGE_bias_corr.nii.gz', reorient =T)
  
  bias_diff <-  MPRAGE_base_bias_corrected - MPRAGE_base
  
  orthographic( bias_diff )
  
  ##############################################
  ## Skull strip with fslr using 'fslbet' 
  ##############################################
  
  MPRAGE_base_bias_corrected_stripped <- fslbet( MPRAGE_base_bias_corrected, reorient =T) 
  
  #plot
  double_ortho(MPRAGE_base , MPRAGE_base_bias_corrected_stripped )
  
  # make a brain mask from the skull stripped image
  
  
  
  #################################################
  
  library(TRES)
  
  ################################
  # TRR model with envelope structure  (bat data) - 2 way tensor Y(64*64) and binary predictor X( 0,1)/// Y_i : 2 dimensional image, X_i: binary predict (two different groups)
  # 1D-BIC criterion to select the dimension under the TRR mode
  
  # X : 1*20
  # Y : conti 3-way tensor with dim 64*64*20
  # coef : conti 3-way tensor with dim 64*64*1
  # Gamma : two basis matrices, Gamma1, Gamma2 with dim 64*14 
  
  
  data("bat")
  fit_ols1 <- TRR.fit(bat$x, bat$y, method='standard')
  fit_1d1 <- TRR.fit(bat$x, bat$y, u=c(14,14), method='1D')
  fit_pls1 <- TRR.fit(bat$x, bat$y, u=c(14,14), method='PLS')
  fit_1d1
  coef(fit_1d1)
  fitted(fit_1d1)
  residuals(fit_1d1) 
  summary(fit_1d1)
  predict(fit_1d1, bat$x)
  
  plot(fit_ols1)
  plot(fit_1d1)
  plot(fit_pls1)
  
  dist_ols1 <- rTensor :: fnorm( coef( fit_ols1) - bat$coefficients ) 
  dist_1d1 <- rTensor :: fnorm (coef ( fit_1d1) - bat$coefficient ) 
  dist_pls1 <- rTensor :: fnorm( coef (fit_pls1) -bat$coefficient )
  
  c(dist_ols1 , dist_1d1, dist_pls1) 
  
  Pdist_1d1 <- rep(NA_real_, 2)
  Pdist_pls1 <- rep(NA_real_ , 2)
  # subspace ()  : compute subspace distance .... 
  for (i in 1:2) {
    Pdist_1d1 [i] <- subspace ( bat$ Gamma [[i]] , fit_1d1$Gamma [[i]] ) 
    Pdist_pls1[i] <- subspace ( bat $ Gamma [[i]] , fit_pls1 $Gamma [[i]] ) 
  }
  
  Pdist_1d1 <- sum( Pdist_1d1)
  Pdist_pls1 <- sum (Pdist_pls1)
  
  c(Pdist_1d1, Pdist_pls1) 
  
  #######################################
  # TPR model with envelope structure  (Square data)
  # Under the TPR model, we use cross-validation to select dimension based on the prediction mean squared error
  # X : 32*32*200 # 2-way tensor (matrix) predictor X
  # Y : 1*200  (# scalar response)
  # coef : 3-way tensor with dim 32*32*1  # 0.1 or 1 (=B)
  # Gamma : two basis matrices, Gamma1, Gamma2 with dim 32*2  #(SVD of B)
  
  data("square")
  fit_ols2 <- TPR.fit(square$x, square$y, method='standard')
  fit_1d2 <- TPR.fit(square$x, square$y, u=c(2,2), method='1D')
  fit_pls2 <- TPR.fit(square$x, square$y, u=c(2,2), method='PLS')
  fit_1d2
  coef(fit_1d2)
  dim( fitted(fit_1d2) )
  dim ( residuals(fit_1d2) )
  
  plot(fit_ols2)
  plot(fit_1d2)
  plot(fit_pls2)
  
  dist_ols2 <- rTensor :: fnorm( coef( fit_ols2) - square$coefficients ) 
  dist_1d2 <- rTensor :: fnorm (coef ( fit_1d2) - square$coefficient ) 
  dist_pls2 <- rTensor :: fnorm( coef (fit_pls2) -square$coefficient )
  
  c(dist_1d2, dist_pls2) 
  
  Pdist_1d2 <- rep(NA_real_, 2)
  Pdist_pls2 <- rep(NA_real_ , 2) 
  # subspace ()  : compute subspace distance .... 
  for (i in 1:2) {
    Pdist_1d2 [i] <- subspace ( square$ Gamma [[i]] , fit_1d2$Gamma [[i]] ) 
    Pdist_pls2[i] <- subspace ( square$ Gamma [[i]] , fit_pls2$Gamma [[i]] ) 
  }
  
  Pdist_1d2 <- sum( Pdist_1d2)
  Pdist_pls2 <- sum (Pdist_pls2)
  
  c(Pdist_1d2, Pdist_pls2) 
  
  ###################################
  # Envelope dimension selelction
  # TRRdim (x, y, C=Null, maxdim =10 , ...) # C: dim of the predictor. 
  # TPRdim (x, y, maxdim =10, nfolds = 5)  # CV
  # maxdim : half of the maximum dimension of the tensor. 
  
  # 3.4 p value in the TRR model
  # The coef B in the TRR model is sparse, which means that most of the elements are zero. 
  summary( fit_ols1) $p_val
  
  summary(fit_1d1)$p-val
  
  summary(fit_pls1) $p_val
  
  # We present the p value plots at significance level 0.05 for each estimator. 
  plot( fit_ols1, level =0.05)
  plot( fit_1d1, level =0.05)
  plot( fit_pls1, level =0.05)
  
  # 3.5 EEG data analysis.(electro encephalography - 뇌파검사)
  # They focus on detecting the difference of EEG images between the two groups, which helps expose the EEG correlates of genetic predisposition to alcoholism
  
  #For envelope-based methods, we first estimate the envelope dimension with the function TRRdim().
  
  data("EEG", package = "TRES")
  u_eeg <- TRRdim(EEG$x, EEG$y)
  u_eeg
  
  fit_eeg_ols <- TRR.fit ( EEG$x, EEG$y , u_eeg$u, method = "standard")
  fit_eeg_1D<- TRR.fit ( EEG$x, EEG$y , u_eeg$u, method = "1D")
  fit_eeg_pls <- TRR.fit ( EEG$x, EEG$y , u_eeg$u, method = "PLS")
  
  plot(fit_eeg_ols, xlab = "Time", ylab = "Channels", cex.main = 2, cex.axis = 2, cex.lab = 1.5)
  plot(fit_eeg_1D, xlab = "Time", ylab = "Channels", cex.main = 2, cex.axis = 2, cex.lab = 1.5)
  plot(fit_eeg_pls, xlab = "Time", ylab = "Channels", cex.main = 2, cex.axis = 2, cex.lab = 1.5)
  
  
  # 4. Envelope algorithms and dimension selection
  
  p <-20
  u <-5
  data <- MenvU_sim( p, u, jitter = 1e-5)
  Gamma <- data$ Gamma
  M <- data$M
  U  <- data$U
  G <- vector ("list", 8)
  
  # compare core function of each algorithm 
  G[[1]] <- simplsMU(M, U, u)
  G[[2]] <- ECD(M,U, u)
  G[[3]] <- manifold1D( M, U, u)
  G[[4]] <- OptM1D(M,U, u)
  G[[5]] <- manifoldFG(M ,U, u)   # initial value impacts
  G[[6]] <- OptMFG(M,U,u)
  
  # compare estimation error of the envelope basis using subspace() : 
  # true envelope basis Gamma
  
  d <- rep(NA_real_, 8)
  for (i in 1:6) {
    d[i] <- subspace( G[[i]] , Gamma)
  }
  res = d[1:6] 
  order(res)
  
  # 초기값에 따라 달라지는 FG 알고리즘을 보고자.....( manifoldFG(), OptMFG() ) 
  # 무작위로 생성된 p×u 행렬 A를 초기 값으로 전달함으로써, 두 가지 결과
  
  A <- matrix(runif( p * u), p, u)
  G[[7]] <- manifoldFG(M, U, u, Gamma_init = A)
  G[[8]] <- OptMFG(M, U, Gamma_init = A)
  for (i in 7:8) {
    d[i] <- subspace(G[[i]], Gamma)
  }
  d[5:8]   # 확실히 poor estimation of envelope that is far from the true envelope. 
  
  #####################
  #oneD_bic(M,U, n, C=1, maxdim =10 , ...) 
  
  # How the sample size influcences the behavior of the 1D-BIC criterion,
  # sample M, U from Wishart distribution.. (p*p)
  # fix p=50, envelope dimension u=5
  # vary the sample size n= 50,70,100,200,400,800
  
  p=50
  u <- 5
  n0 <- c(50,70,100,200,400,800)
  uhat3  <- rep(NA_integer_, length(n0) )
  for (i in seq_along(n0)) {
    n=n0[i]
    data = MenvU_sim( p, u, jitter = 1e-5, wishart =T, n=n)
    M <- data$M
    U <- data$U
    output <- oneD_bic (M,U, n, maxdim = p/2 ) 
    uhat3 [i] <- output$u
  }
  uhat3 
  # result indicate that as n increases, the function oneD_bic() consistently select the correct dimension. 
  
  