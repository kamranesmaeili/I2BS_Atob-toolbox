
##### Input section #####

library(plotly)
library(reshape2)
library(ggplot2)
library(grid)
library(signal)
library(seewave)
library(oce)
library(e1071)
library(stringr)

damping <- 0.05
Noise.items <- c(1)
Operation <- c(TRUE, TRUE, TRUE, TRUE)
Side.bands <- TRUE
Plots.Signal <- c(Operation, TRUE)
Plots.method <- c(TRUE, TRUE, TRUE, TRUE)

n <- 1
Speed.shaft <- 13940

Plots <- TRUE
Prime.Har.factor <- 0.2
Resonance <- c(NA, 25000, 30000)
theta <- c(NA, 0, 0)
T <- c(1/233, 1/2585, 1/2064)
Sampling.rate <- 100000
Pres.shaft.values <- c(3, 4, 2, 4, 6, 1, as.numeric(Sampling.rate)/2, 1)
Over.sampling <- 1
BPF.order <- 1
BPF.values <- c(22500, 32500)
Split.x <- 10000
Primary_frequencies <- c(233, 2083.5, 2632, 0)

Defect.list <- c(paste0(1:7, "ShRF"), 
                 # paste0("+", 1:3, "xBSP"),
                 paste0("+", 1:20, "xBPFO"),
                 paste0("+", 1:20, "xBPFI"),
                 as.vector(outer(paste0(-4:4, "ShRF"), paste0("+", 1:20, "xBPFI"), paste, sep="")))

Plot.directory <- paste0(getwd(), "/Images/2.3 Toolbox/", Prime.Har.factor, " Prime/")



##### First section, creating simulated datasets #####

options(scipen=0)
windowsFonts(Times=windowsFont("Times New Roman"))
source("Atob modules.R")
Frequency <- seq(1, Sampling.rate, 1)[1:(Sampling.rate/2)]
k <- 1
name.items <- paste(as.integer(as.logical(Operation)),collapse="-")
Noise.SD <- Noise.items[k]

t=seq(0,1, length.out = Sampling.rate*Over.sampling)

signal.A <- 0
signal.1 <- 0
signal.2 <- 0
signal.3 <- 0
signal.e <- 0


if (Operation[1] == TRUE){
  for (n in 1:(floor(Sampling.rate/2/(1/T[1]))-2)){
    signal.1 <- signal.1+((0.5/n)*sin(2*pi*(1/T[1])*t*n))
  }
  
} else{
  signal.1 <- rep(0, length(t))
}
signal.1 <- signal.1[(seq(1, length(t), length.out = length(t)/Over.sampling))]
ShRF.fft <- abs(fft(signal.1))[1:(Sampling.rate/2)]
plot(Frequency, ShRF.fft, type = "l", ylab = "Amplitude", xlab = "Frequency (Hz)")

if (Operation[2] == TRUE){
  
  X_piece <- (exp((-damping/sqrt(1-damping^2))*(2*pi*Resonance[2]*t))*
                sin(2*pi*Resonance[2]*(t-theta[2]-(T[2]))))[1:floor(T[2]*(length(t)/max(t)))]
  while(length(signal.2) < length(t)){signal.2 <- c(signal.2, X_piece)}
  signal.2 <- signal.2[1:length(t)]
  
} else{
  signal.2 <- rep(0, length(t))
}

signal.2 <- signal.2[(seq(1, length(t), length.out = length(t)/Over.sampling))]
BPFI.fft <- abs(fft(signal.2))[1:(Sampling.rate/2)]
plot(Frequency, BPFI.fft, type = "l", ylab = "Amplitude", xlab = "Frequency (Hz)")

if (Operation[3] == TRUE){
  X_piece <- (exp((-damping/sqrt(1-damping^2))*(2*pi*Resonance[3]*t))*
                sin(2*pi*Resonance[3]*(t-theta[3]-(T[3]))))[1:round(T[3]*(length(t)/max(t)))]
  
  while(length(signal.3) < length(t)){signal.3 <- c(signal.3, X_piece)}
  signal.3 <- signal.3[1:length(t)]
  
} else{
  signal.3 <- rep(0, length(t))
}
signal.3 <- signal.3[(seq(1, length(t), length.out = length(t)/Over.sampling))]
BPFO.fft <- abs(fft(signal.3))[1:(Sampling.rate/2)]
plot(Frequency, BPFO.fft, type = "l", ylab = "Amplitude", xlab = "Frequency (Hz)")

if (Operation[4] == TRUE){
  signal.e <- rnorm(Sampling.rate, mean = 0, sd = Noise.SD)
} else{
  signal.e <- rep(0, Sampling.rate)
}
e.fft <- abs(fft(signal.e))[1:(Sampling.rate/2)]
plot(Frequency, e.fft, type = "l", ylab = "Amplitude", xlab = "Frequency (Hz)")

Noise.level <- max(e.fft)
if (Side.bands == FALSE){
  signal.A <- signal.1+signal.2+signal.3+signal.e
} else {
  if (Operation[3] == TRUE){
    print("FAULT: OR & IR")
    signal.A <- (signal.1*signal.2)+(signal.2*signal.3)+(signal.1*signal.3)+ ((signal.1+signal.2+signal.3)*Prime.Har.factor)+signal.e
  } else if (Operation[3] == FALSE){
    print("FAULT: IR")
    signal.A <- (signal.1*signal.2)+((signal.1+signal.2+signal.3)*Prime.Har.factor)+signal.e
  }
}

Signal.A.fft <- abs(fft(signal.A))[1:(Sampling.rate/2)]
t <- t[seq(1, length(t), length.out = Sampling.rate)]

plot(t, signal.A, type = "l")
plot(Frequency, Signal.A.fft, type = "l", ylab = "Amplitude", xlab = "Frequency (Hz)")


SNR <- 10*log10(sum((signal.2+signal.3)^2)/sum((signal.e)^2))
SRR <- 10*log10(sum((signal.2+signal.3)^2)/sum((signal.e+signal.1)^2))

print(paste0("Signal-to-noise ratio: ", round(SNR, 3), "     Signal-to-rest ratio: ", round(SRR, 3)))


###### Module 1 ######
Defect.values <- Module_1(Speed.shaft, 0.5, 203000, 20, 75, 32.7375*2, 4.955, 
                          42.2625*2, 4.905, 4.7625*2, 0.08, 250, 250)

###### Module 2 ######
Frequency.data <- Module_2(Primary_frequencies[1], Primary_frequencies[4], 
                           Primary_frequencies[2], Primary_frequencies[3], 
                           0:38, Split.x, 1, 11:40, 
                           1, 22, 22)

###### Some customised cleaning or sorting of the library #####
if (!is.null(Defect.list)){
  Frequency.data <- Frequency.data[Frequency.data$variable %in% Defect.list,]
} else{
  Frequency.data <- Frequency.data[Frequency.data$variable %in% c(paste0(1:7, "ShRF"), as.vector(outer(paste0("+", 1:25), c("xBPFO", "xBPFI", "xBSF"), paste, sep=""))),]
}

Frequency.data <- Frequency.data[!Frequency.data$variable %in% c(paste0(8:40, "ShRF")),]
Frequency.data <- Frequency.data[Frequency.data$value > 0,]
Frequency.data <- Frequency.data[!grepl("BSF", Frequency.data$variable),]

###### Module 3 ######
Freq.signals <- Module_3(signal.A, Sampling.rate, Method="Fast-Fourier transform", 
                   RangeOfBPF=c(20000, 25000), OrderOfBPF=1)
Defect.data <- Freq.signals[Freq.signals$`Frequency` > 5,]

###### Module 4 ######
Peaks.data <- Module_4(Defect.data, 500, 1, 2, 3, 0.99,  8)

###### Module 5 ######
Frequency.data <- Module_5(Peaks.data, Frequency.data, 0, 2)
Frequency.data <- Frequency.data[!duplicated(Frequency.data$value) & 
                                   !grepl(paste0("^", 9:round(Split.x/Primary_frequencies[1]), "ShRF$", collapse = "|"), 
                                          Frequency.data$variable),]

Frequency.data$variable <- gsub(paste0("^", 0:round(Split.x/Primary_frequencies[1]), "ShRF$", collapse = "|"), 
                                "ShRF", Frequency.data$variable)

Frequency.data$variable <- gsub(paste0("^", -round(Split.x/Primary_frequencies[1]):round(Split.x/Primary_frequencies[1]), 
                                       "ShRF", collapse = "|"), "XShRF", Frequency.data$variable)

