
HertzianFactor <- as.data.frame(cbind(c(seq(1, 0.950, length.out=51), seq(0.948, 0.900, length.out=25), seq(0.895, 0.850, length.out=10), seq(0.840, 0.800, length.out=5),
                                        seq(0.750, 0, length.out = 16)), 
                                      
                                      c(Inf, 18.53 , 14.25 , 12.26 , 11.02, 10.15 , 9.46 , 8.92 , 8.47 , 8.10 , 7.76 , 7.49, 7.25, 7.02 , 6.84 , 6.64 , 
                                        6.47 , 6.33 , 6.19 , 6.06 , 5.94 , 5.83 , 5.72 , 5.63 , 5.53 , 5.44 , 5.35 , 5.28 , 5.20 , 5.13 , 5.05, 4.98 , 4.92 , 4.86 , 
                                        4.81 , 4.76 , 4.70 , 4.65, 4.61 , 4.56 , 4.51 , 4.47 , 4.42 , 4.38 , 4.34 , 4.30 , 4.26 , 4.22, 4.19 , 4.15 , 4.12 , 4.05, 
                                        3.99 , 3.94, 3.88 , 3.83 , 3.78 , 3.73 , 3.68 , 3.63 , 3.59 , 3.55 , 3.51 , 3.47 , 3.43 , 3.40 , 3.36 , 3.33 , 3.30 , 3.27 , 
                                        3.23 , 3.20 , 3.17 , 3.15 , 3.12 , 3.09 , 3.03 , 2.97 , 2.92 , 2.86 , 2.82 , 2.77 , 2.72 , 2.68 , 2.64 , 2.60 , 2.53 , 2.46 , 
                                        2.40 , 2.35 , 2.30 , 2.07 , 1.91 , 1.77 , 1.66, 1.57 , 1.48 , 1.41 , 1.35, 1.29, 1.24 , 1.19 , 1.15 , 1.11 , 1.07, 1.03 , 1 ), 
                                      c(0, 0.185, 0.212 , 0.228 , 0.241 , 0.251, 0.260 , 0.268 , 0.275 , 0.281 , 0.287 , 0.292 , 0.297 , 0.301 , 0.305 , 0.310 , 0.314, 
                                        0.317 , 0.321 , 0.325 , 0.328 , 0.332 , 0.335 , 0.338 , 0.340 , 0.343 , 0.346 , 0.349 , 0.351 , 0.354 , 0.357 , 0.359 , 0.361 , 0.363 , 
                                        0.365 , 0.367 , 0.369 , 0.371 , 0.374 , 0.376 , 0.378 , 0.380, 0.382 , 0.384 , 0.386 , 0.388 , 0.390, 0.391 , 0.393 , 0.394 , 0.396 , 
                                        0.399 , 0.403, 0.406 , 0.409 , 0.412 , 0.415 , 0.418 , 0.420 , 0.423 , 0.426 , 0.428 , 0.431 , 0.433 , 0.436 , 0.438 , 0.441 , 0.443 , 
                                        0.445 , 0.448 , 0.450 , 0.452 , 0.454 , 0.456 , 0.459 , 0.461 , 0.466 , 0.471 , 0.476 , 0.481 , 0.485 , 0.490 , 0.494 , 0.498 , 0.502 , 
                                        0.507 , 0.515 , 0.523 , 0.530 , 0.537 , 0.544 , 0.577 , 0.607 , 0.637 , 0.664 , 0.690 , 0.718 , 0.745 , 0.771 , 0.796 , 0.824 , 0.850 , 
                                        0.879 , 0.908 , 0.938 , 0.969 , 1.0 ), 
                                      c(0, 0.207 ,0.249 ,0.279 ,0.302 ,0.320 ,0.336 ,0.350 ,0.362 ,0.373 ,0.384 ,0.393 ,0.402 ,0.411 ,0.420 ,0.427 ,0.433 ,0.440 ,
                                        0.447 ,0.453 ,0.459 ,0.465 ,0.470 ,0.476 ,0.481 ,0.486 ,0.491 ,0.495 ,0.500 ,0.505 ,0.509 ,0.513 ,0.518 ,0.522 ,0.525 ,0.530,0.533 ,
                                        0.536 ,0.540,0.543,0.546 ,0.550 ,0.553 ,0.556 ,0.559 ,0.562 ,0.565 ,0.568 ,0.571 ,0.574 ,0.577 ,0.583 ,0.588, 0.593 , 0.598 , 0.603 , 
                                        0.608 , 0.613 , 0.618 , 0.622 , 0.626 , 0.630 , 0.634 , 0.638 , 0.642 , 0.646 , 0.650 , 0.653 , 0.657 , 0.660, 0.664 , 0.667 , 0.671 , 
                                        0.674 , 0.677 , 0.678 , 0.688 , 0.695 , 0.702 , 0.709 , 0.715 , 0.721 , 0.727 , 0.733 , 0.739 , 0.745 , 0.755 , 0.765 , 0.774 , 0.783 , 
                                        0.792 , 0.829 , 0.855 , 0.884 , 0.904 , 0.922 , 0.938 , 0.951 , 0.962 , 0.971 , 0.979 , 0.986 , 0.991, 0.994 , 0.997 , 0.999 , 1.0  )))
colnames(HertzianFactor) <- c("cos_t", "mu", "v", "HertzValue")



StaticLoadFactor <- as.data.frame(cbind(c(0, 0.05 , 0.06 , 0.07 , 0.08 , 0.09 , 0.10, 0.12 , 0.14 , 0.16 , 0.18 , 0.20, 0.22 , 0.24 , 0.26 , 0.28 , 0.30, 0.32 , 0.34 , 0.36 , 0.38, 0.40), 
                                        c(14.7 , 15.7 , 15.9 , 16.1 , 16.3 , 16.5 , 16.4 , 15.9 , 15.4 , 14.9 , 14.4 , 14.0, 13.5 , 13 , 12.5 , 12.1 , 11.6 , 11.2 , 10.7 , 10.3, 9.8, 9.4)))
colnames(StaticLoadFactor) <- c("Dw_Cos(a)/Dpw", "f0")



LoadFactors <- as.data.frame(cbind(c(20, 25, 30, 35, 40, 45), 
                                   c(0.57, 0.68, 0.80, 0.95, 1.14, 1.33), 
                                   rep(1, 6), rep(0, 6),
                                   c(0.43, 0.41, 0.39, 0.37, 0.35, 0.33),
                                   c(1, 0.87, 0.76, 0.66, 0.57, 0.50))) 
colnames(LoadFactors) <- c("angles", "e", "lessEX", "lessEY", "moreEX", "moreEY")

#####---------------------Make numeric values-----------------------------------------------------------
Make.numeric <- function(InputSignals, Numeric.column){
  for (c in 1:ncol(InputSignals)){
    InputSignals[,c] <- as.character(InputSignals[,c])
    if (c %in% Numeric.column){
      InputSignals[,c] <- as.numeric(InputSignals[,c])
    }
  }
  return(InputSignals)
}
#####---------------------Function to calculate all bearing frequencies---------------------------------
Module_1 <- function(Speed, Pois, ElasMod, NumBalls, PCD, DInnerX, RInnerY, DOuterX, 
                     ROuterY, DBall, RClear, AxialLoad, RadLoad){
  
  i <- 1
  
  p11 <- +1/(DInnerX/2)
  p12 <- -1/RInnerY
  p21 <- -1/(DOuterX/2)
  p22 <- -1/ROuterY
  pb1 <- 1/(DBall/2)
  pb2 <- 1/(DBall/2)
  
  AuxNumber <- RInnerY+ROuterY-DBall
  
  Ep_1 <- sum(p11, p12, pb1, pb2) 
  Ep_2 <- sum(p21, p22, pb1, pb2)
  
  cosT_1 <- (p11-p12+pb1-pb2)/Ep_1
  cosT_2 <- (p21-p22+pb1-pb2)/Ep_2
  
  t_1 <- acos(cosT_1)
  t_2 <- acos(cosT_2)
  
  alpha_n <- acos(1-(RClear/(2*AuxNumber)))
  alpha_n_degrees <- alpha_n*180/pi
  
  f0_identifier <- (DBall*cos(alpha_n))/PCD
  
  f0 <- as.numeric(StaticLoadFactor[which(abs(StaticLoadFactor$`Dw_Cos(a)/Dpw`-f0_identifier)==min(abs(StaticLoadFactor$`Dw_Cos(a)/Dpw`-f0_identifier)))[1],2])
  C0 <- f0*i*NumBalls*(DBall^2)*cos(alpha_n)
  Q <- (5*C0)/(i*NumBalls*cos(alpha_n))
  
  TempLoadFactors <- as.numeric(LoadFactors[which(abs(LoadFactors$angles-alpha_n_degrees)==min(abs(LoadFactors$angles-alpha_n_degrees)))[1],])
  
  if ((AxialLoad/RadLoad)<=TempLoadFactors[2]){
    P0 <- TempLoadFactors[3] * RadLoad + TempLoadFactors[4] * AxialLoad
  } else if ((AxialLoad/RadLoad)>TempLoadFactors[2]){
    P0 <- TempLoadFactors[5] * RadLoad + TempLoadFactors[6] * AxialLoad
  }
  
  def_1 <- 1.5*as.numeric(HertzianFactor[which(abs(HertzianFactor$cos_t-cosT_1)==min(abs(HertzianFactor$cos_t-cosT_1)))[1],4])*((((1-Pois^2)^2)/ElasMod^2)*(Ep_1/3)*(Q^2))^(1/3)
  def_2 <- 1.5*as.numeric(HertzianFactor[which(abs(HertzianFactor$cos_t-cosT_2)==min(abs(HertzianFactor$cos_t-cosT_2)))[1],4])*((((1-Pois^2)^2)/ElasMod^2)*(Ep_2/3)*(Q^2))^(1/3)
  
  alpha_o_1 <- asin((AuxNumber*sin(alpha_n)+def_1)/((((AuxNumber^2)*cos(alpha_n)^2)+(((AuxNumber*sin(alpha_n))+def_1)^2))^(1/2)))
  alpha_o_2 <- asin((AuxNumber*sin(alpha_n)+def_2)/((((AuxNumber^2)*cos(alpha_n)^2)+(((AuxNumber*sin(alpha_n))+def_2)^2))^(1/2)))
  
  alpha_o_1_degrees <- alpha_o_1*180/pi
  alpha_o_2_degrees <- alpha_o_2*180/pi
  
  ShaftFreq <- Speed/60
  
  CageSpeed <- (Speed/2)*(1-((DBall*cos(alpha_o_1))/PCD))
  
  BPFO <- (NumBalls/2)*ShaftFreq*(1-((DBall/PCD)*cos(alpha_n)))
  BPFI <- (NumBalls/2)*ShaftFreq*(1+((DBall/PCD)*cos(alpha_n)))
  BCF=(ShaftFreq/2)*((PCD/DBall)-((DBall/PCD)*((cos(alpha_n))^2)))
  BK=(ShaftFreq/2)*(1-((DBall/PCD)*(cos(alpha_n))))
  
  Output <- as.data.frame(cbind(c("Auxiliary value", 
                                  "Inner ring curvature ratio (acos())", "Inner ring operating contact angle", "Deflection at inner ring", 
                                  "Outer ring curvature ratio (acos())", "Outer ring operating contact angle", "Deflection at outer ring",
                                  
                                  "Nominal contact angle", "Nominal contact angle", 
                                  
                                  "Equivalent load (P0)", "Static load factor (f0)", "Static load rating (C0)", "Maximum load on a ball (Qmax)", "Shaft speed", 
                                  "Ball passing frequency on outer ring (BPFO)", "Ball passing frequency on inner ring (BPFI)", "Ball circulation frequency (BSF)", "Fundamental train frequency (FTF)",
                                  "Cage speed"),
                                
                                c("mm", 
                                  "-", "°", "mm", 
                                  "-", "°", "mm", 
                                  
                                  "rad", "°", 
                                  
                                  "N", "-", "N", "N", "Hz",
                                  rep("Hz", 4), "rpm"),
                                round(c(AuxNumber, t_1, alpha_o_1_degrees, def_1, t_2, alpha_o_2_degrees, def_2, alpha_n, alpha_n_degrees,
                                        P0, f0, C0, Q, ShaftFreq, BPFO, BPFI, BCF, BK, CageSpeed), digits = 3)))
  
  colnames(Output) <- c("Parameters", "Units", "Values")
  rownames(Output) <- 1:nrow(Output)
  return(Output)
}
#####---------------------Analytical model to calculate all combinations of bearing frequencies---------
Module_2 <- function(ShRF, BSF, BPFO, BPFI, ShRF_Harmonics, MaxFreq, ArbAmp, ExempShRF_Harmonics,
                     BSF_Harmonics, BPFO_Harmonics, BPFI_Harmonics){
  
  
  Output <- as.data.frame(cbind(rep(ArbAmp, length(ShRF_Harmonics)), 
                                outer(c(ShRF_Harmonics), c(ShRF)),
                                
                                outer(ShRF*ShRF_Harmonics, c(BSF*(1:BSF_Harmonics), 
                                                             BPFO*(1:BPFO_Harmonics),
                                                             BPFI*(1:BPFI_Harmonics),
                                                             
                                                             BSF*(1:BSF_Harmonics)+BPFO,
                                                             BSF*(1:BSF_Harmonics)-BPFO,
                                                             BSF*(1:BSF_Harmonics)+BPFI,
                                                             BSF*(1:BSF_Harmonics)-BPFI,
                                                             
                                                             BPFO*(1:BPFO_Harmonics)+BSF,
                                                             BPFO*(1:BPFO_Harmonics)-BSF,
                                                             BPFO*(1:BPFO_Harmonics)+BPFI,
                                                             BPFO*(1:BPFO_Harmonics)-BPFI,
                                                             
                                                             BPFI*(1:BPFI_Harmonics)+BSF,
                                                             BPFI*(1:BPFI_Harmonics)-BSF,
                                                             BPFI*(1:BPFI_Harmonics)+BPFO,
                                                             BPFI*(1:BPFI_Harmonics)-BPFO,
                                                             
                                                             BSF*(1:BSF_Harmonics)+BPFO+BPFI,
                                                             BSF*(1:BSF_Harmonics)-BPFO+BPFI,
                                                             BSF*(1:BSF_Harmonics)+BPFO-BPFI,
                                                             BSF*(1:BSF_Harmonics)-BPFO-BPFI,
                                                             
                                                             BPFO*(1:BPFO_Harmonics)+BSF+BPFI,
                                                             BPFO*(1:BPFO_Harmonics)-BSF+BPFI,
                                                             BPFO*(1:BPFO_Harmonics)+BSF-BPFI,
                                                             BPFO*(1:BPFO_Harmonics)-BSF-BPFI,
                                                             
                                                             BPFI*(1:BPFI_Harmonics)+BSF+BPFO,
                                                             BPFI*(1:BPFI_Harmonics)-BSF+BPFO,
                                                             BPFI*(1:BPFI_Harmonics)+BSF-BPFO,
                                                             BPFI*(1:BPFI_Harmonics)-BSF-BPFO), 
                                      FUN = "+"),
                                
                                outer(ShRF*ShRF_Harmonics, c(BSF*(1:BSF_Harmonics), 
                                                             BPFO*(1:BPFO_Harmonics),
                                                             BPFI*(1:BPFI_Harmonics),
                                                             
                                                             BSF*(1:BSF_Harmonics)+BPFO,
                                                             BSF*(1:BSF_Harmonics)-BPFO,
                                                             BSF*(1:BSF_Harmonics)+BPFI,
                                                             BSF*(1:BSF_Harmonics)-BPFI,
                                                             
                                                             BPFO*(1:BPFO_Harmonics)+BSF,
                                                             BPFO*(1:BPFO_Harmonics)-BSF,
                                                             BPFO*(1:BPFO_Harmonics)+BPFI,
                                                             BPFO*(1:BPFO_Harmonics)-BPFI,
                                                             
                                                             BPFI*(1:BPFI_Harmonics)+BSF,
                                                             BPFI*(1:BPFI_Harmonics)-BSF,
                                                             BPFI*(1:BPFI_Harmonics)+BPFO,
                                                             BPFI*(1:BPFI_Harmonics)-BPFO,
                                                             
                                                             BSF*(1:BSF_Harmonics)+BPFO+BPFI,
                                                             BSF*(1:BSF_Harmonics)-BPFO+BPFI,
                                                             BSF*(1:BSF_Harmonics)+BPFO-BPFI,
                                                             BSF*(1:BSF_Harmonics)-BPFO-BPFI,
                                                             
                                                             BPFO*(1:BPFO_Harmonics)+BSF+BPFI,
                                                             BPFO*(1:BPFO_Harmonics)-BSF+BPFI,
                                                             BPFO*(1:BPFO_Harmonics)+BSF-BPFI,
                                                             BPFO*(1:BPFO_Harmonics)-BSF-BPFI,
                                                             
                                                             BPFI*(1:BPFI_Harmonics)+BSF+BPFO,
                                                             BPFI*(1:BPFI_Harmonics)-BSF+BPFO,
                                                             BPFI*(1:BPFI_Harmonics)+BSF-BPFO,
                                                             BPFI*(1:BPFI_Harmonics)-BSF-BPFO),
                                      FUN = "-")*(-1)))
  
  colnames(Output) <- c("YMag", "ShRF", 
                        paste0("+", c(paste0(1:BSF_Harmonics, "xBSF"), paste0(1:BPFO_Harmonics, "xBPFO"), paste0(1:BPFI_Harmonics, "xBPFI"), 
                                      paste0(1:BSF_Harmonics, "xBSF+BPFO"), paste0(1:BSF_Harmonics, "xBSF-BPFO"), paste0(1:BSF_Harmonics, "xBSF+BPFI"), paste0(1:BSF_Harmonics, "xBSF-BPFI"),
                                      paste0(1:BPFO_Harmonics, "xBPFO+BSF"), paste0(1:BPFO_Harmonics, "xBPFO-BSF"), paste0(1:BPFO_Harmonics, "xBPFO+BPFI"), paste0(1:BPFO_Harmonics, "xBPFO-BPFI"),
                                      paste0(1:BPFI_Harmonics, "xBPFI+BSF"), paste0(1:BPFI_Harmonics, "xBPFI-BSF"), paste0(1:BPFI_Harmonics, "xBPFI+BPFO"), paste0(1:BPFI_Harmonics, "xBPFI-BPFO"),
                                      paste0(1:BSF_Harmonics, "xBSF+BPFO+BPFI"), paste0(1:BSF_Harmonics, "xBSF-BPFO+BPFI"), paste0(1:BSF_Harmonics, "xBSF+BPFO-BPFI"), paste0(1:BSF_Harmonics, "xBSF-BPFO-BPFI"),
                                      paste0(1:BPFO_Harmonics, "xBPFO+BSF+BPFI"), paste0(1:BPFO_Harmonics, "xBPFO-BSF+BPFI"), paste0(1:BPFO_Harmonics, "xBPFO+BSF-BPFI"), paste0(1:BPFO_Harmonics, "xBPFO-BSF-BPFI"),
                                      paste0(1:BPFI_Harmonics, "xBPFI+BSF+BPFO"), paste0(1:BPFI_Harmonics, "xBPFI-BSF+BPFO"), paste0(1:BPFI_Harmonics, "xBPFI+BSF-BPFO"), paste0(1:BPFI_Harmonics, "xBPFI-BSF-BPFO"))),
                        
                        paste0("-", c(paste0(1:BSF_Harmonics, "xBSF"), paste0(1:BPFO_Harmonics, "xBPFO"), paste0(1:BPFI_Harmonics, "xBPFI"), 
                                      paste0(1:BSF_Harmonics, "xBSF+BPFO"), paste0(1:BSF_Harmonics, "xBSF-BPFO"), paste0(1:BSF_Harmonics, "xBSF+BPFI"), paste0(1:BSF_Harmonics, "xBSF-BPFI"),
                                      paste0(1:BPFO_Harmonics, "xBPFO+BSF"), paste0(1:BPFO_Harmonics, "xBPFO-BSF"), paste0(1:BPFO_Harmonics, "xBPFO+BPFI"), paste0(1:BPFO_Harmonics, "xBPFO-BPFI"),
                                      paste0(1:BPFI_Harmonics, "xBPFI+BSF"), paste0(1:BPFI_Harmonics, "xBPFI-BSF"), paste0(1:BPFI_Harmonics, "xBPFI+BPFO"), paste0(1:BPFI_Harmonics, "xBPFI-BPFO"),
                                      paste0(1:BSF_Harmonics, "xBSF+BPFO+BPFI"), paste0(1:BSF_Harmonics, "xBSF-BPFO+BPFI"), paste0(1:BSF_Harmonics, "xBSF+BPFO-BPFI"), paste0(1:BSF_Harmonics, "xBSF-BPFO-BPFI"),
                                      paste0(1:BPFO_Harmonics, "xBPFO+BSF+BPFI"), paste0(1:BPFO_Harmonics, "xBPFO-BSF+BPFI"), paste0(1:BPFO_Harmonics, "xBPFO+BSF-BPFI"), paste0(1:BPFO_Harmonics, "xBPFO-BSF-BPFI"),
                                      paste0(1:BPFI_Harmonics, "xBPFI+BSF+BPFO"), paste0(1:BPFI_Harmonics, "xBPFI-BSF+BPFO"), paste0(1:BPFI_Harmonics, "xBPFI+BSF-BPFO"), paste0(1:BPFI_Harmonics, "xBPFI-BSF-BPFO"))))
  
  Output <- as.data.frame(cbind(paste0(ShRF_Harmonics, "ShRF"), Output))
  colnames(Output)[1] <- "ShRF_Harmonics"
  Output <- reshape2::melt(Output, id=c("YMag", "ShRF_Harmonics"))
  
  Output$variable <- paste0(Output$ShRF_Harmonics, Output$variable)
  Output <- subset(Output, select=-c(ShRF_Harmonics))
  Output$variable <- gsub("ShRFShRF", "ShRF", Output$variable)
  
  
  id.minus <- which(grepl(paste0("ShRF-", collapse = "|"), Output$variable))
  Output$variable[id.minus] <- paste0("-", gsub("ShRF-", "ShRF\\+", Output$variable[id.minus]))
  Output <- Output[which(!grepl(paste0(paste0(ExempShRF_Harmonics, "ShRF-", collapse = "|"), "|", paste0(ExempShRF_Harmonics, "ShRF\\+", collapse = "|")), Output$variable)),]
  Output$variable <- gsub("^0ShRF", "", Output$variable)
  
  return(Output)
  
}
#####---------------------Function to convert signals to frequency form---------------------------------
Module_3 <- function(InputSignals, SamplingRate, Method, RangeOfBPF, OrderOfBPF){
  nfft <- length(InputSignals)
  wnd <- hamming.window(nfft)
  df  <- nfft/SamplingRate
  Frequency <- (1:(nfft/2))*(as.numeric(SamplingRate)/nfft)
  
  if (Method == "Fast-Fourier transform"){
    Mag.label <- "Amplitude (g)"
    Amplitude <- abs(fft(InputSignals))[(1:(nfft/2))]
    
  } else if (Method == "Power-spectral density"){
    Mag.label <- expression(paste("Amplitude (g"^2, "/Hz)"))
    Amplitude <- pwelch(InputSignals,wnd,0,nfft,SamplingRate, plot = FALSE)$spec * df
    
  } else if (Method == "Envelope"){
    Mag.label <- expression(paste("Amplitude (g"^2, "/Hz)"))
    Amplitude <- pwelch(abs(hilbert(InputSignals, f=SamplingRate)),wnd,0,nfft,
                        SamplingRate, plot = FALSE)$spec * df
    
  } else if (Method == "High-frequency resonance technique"){
    Mag.label <- expression(paste("Amplitude (g"^2, "/Hz)"))
    Amplitude <- pwelch(abs(hilbert(bwfilter(InputSignals, f=SamplingRate, 
                                             n=OrderOfBPF, from=RangeOfBPF[1], to=RangeOfBPF[2]), 
                                    f=SamplingRate)),wnd,0,nfft,SamplingRate, plot = FALSE)$spec * df
    
  } else {
    stop("Invalid 'Method' argument.")
  }
  
  Output <- as.data.frame(cbind(Frequency, Method, Amplitude))
  Output <- Make.numeric(Output, c(1,3))
  colnames(Output) <- c("Frequency", "Method", "Amplitude")
  
  return(Output)
}
#####---------------------Find threshold and return the peaks-------------------------------------------
Module_4 <- function(InputSignals, WinSize, FreqColumn, MethodColumn, AmpColumn, 
                     QuantileFactor, ThresholdFactor){
  
  InputSignals <- Defect.data
  WinSize <- 500
  FreqColumn <- 1
  MethodColumn <- 2
  AmpColumn <- 3
  QuantileFactor <- 0.99
  ThresholdFactor <- 8
  
  TempWin <- floor((max(InputSignals[,FreqColumn])-min(InputSignals[,FreqColumn]))/WinSize)
  OutputTemp <- NULL
  
  if (TempWin == 0){
    y.threshold <- max(InputSignals[AmpColumn])
    OutputTemp <- rbind(OutputTemp, InputSignals[InputSignals[AmpColumn] > (y.threshold/ThresholdFactor),])
  } else{
    for (t in 1:TempWin){
      Temp.dataset <- InputSignals[InputSignals[,FreqColumn] >= (((t-1)*WinSize)+1) & InputSignals[,FreqColumn] < (t*WinSize),]
      y.threshold <- max(c(max(Temp.dataset[AmpColumn])/ThresholdFactor, quantile(Temp.dataset[,AmpColumn], QuantileFactor)))
      OutputTemp <- rbind(OutputTemp, Temp.dataset[Temp.dataset[AmpColumn] > y.threshold,])
    }
  }
  Output <- NULL
  Breaks <- c(0, which(diff(OutputTemp[,FreqColumn]) != 1), length(OutputTemp[,FreqColumn]))
  Temp.list <- sapply(seq(length(Breaks) - 1), function(i) OutputTemp[(Breaks[i] + 1):Breaks[i+1], FreqColumn])
  for (t in 1:(length(Breaks) - 1)){
    Temp <- OutputTemp[OutputTemp[,FreqColumn] %in% unlist(Temp.list[[t]]),]
    Output <- rbind(Output, Temp[Temp[,AmpColumn]==max(Temp[,AmpColumn]),])
  }
  return(Output)
}
#####---------------------Peak frequency identifying and labelling--------------------------------------
Module_5 <- function(DefectDataset, FrequencyLibrary,  MultiplierFactor, OffsetFactor){
  Output <- FrequencyLibrary
  Output[,4] <- ""
  colnames(Output)[4] <- "Magnitude"
  for (z in 1:nrow(DefectDataset)){
    Temp.row <- DefectDataset[z,]
    Temp.data <- Output[(Output$value <= ((Temp.row$Frequency * (1+(MultiplierFactor/1000)))+OffsetFactor)) & 
                          (Output$value >= ((Temp.row$Frequency * (1-(MultiplierFactor/1000)))-OffsetFactor)),]
    id.magnitude <- which((Output$value <= ((Temp.row$Frequency * (1+(MultiplierFactor/1000)))+OffsetFactor)) & 
                            (Output$value >= ((Temp.row$Frequency * (1-(MultiplierFactor/1000)))-OffsetFactor)))
    if (nrow(Temp.data) > 0){
      Output[id.magnitude, 4] <- paste0(Output[id.magnitude, 4], 
                                        Temp.row[,1], ":", Temp.row[,3], ", ")
    }
  }
  return(Output)
}
