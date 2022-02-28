###### SCRIPT ARTICLE


# Packages
pacman::p_load(fBasics, quantmod, tseries, MSGARCH, rugarch, reshape2)




# Getting Data
tickers <- "^BVSP"
suppressWarnings(getSymbols(tickers,to = "2022-02-01"))

bvsp <- `BVSP`
bvsp <- na.omit(diff(log(bvsp$BVSP.Adjusted)))
bvsp <- bvsp*100



###############     FIGURE 1

p <- ggplot(bvsp, aes(x = index(bvsp), y = bvsp)) +
  geom_line(color = "deepskyblue4") +
  ggtitle("") + xlab(" ") + ylab(" ") +
  theme(axis.title.x = element_text(angle = 45),plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")+
  theme_linedraw()

p

## saving image as .jpg in high resolution
ggsave("return_linedraw.jpg", width = 7, height = 5, units = c("in"), dpi = 900)




###############     FORCASTING MEASURES


# Fitting ARMA(1,1) model 
m1=arima(bvsp,order=c(1,0,1),include.mean=T)
resid <- m1$residuals





#### Creating variables
N <- dim(bvsp)[1]
WE <- 1500
VaR.SR <- matrix(nrow = (N-WE), ncol=6*2*2)# 6 dists, two models and 2 alphas (0.01 and 0.025)
VaR.MS <- matrix(nrow = (N-WE), ncol=12*2)
ES.SR <- matrix(nrow = (N-WE), ncol=12*2)
ES.MS <- matrix(nrow = (N-WE), ncol=12*2)



## computing forecastss
for (t in (WE+1):N) {
  
  t1 <- t - WE
  t2 <- t - 1
  window <- resid[t1:t2] # I will use the de-meaned time series
  # That is, I will use the residuals time series.
  print(paste0("Im on : ",t1))
  
  tempo <- Sys.time()
  
  
  
  if(t1%%10 == 0){ # day to estimate the parameters
    
    
    
    ################################################################################################
    ###### SINGLE REGIME
    ################################################################################################
    
    
    # Garch(1,1) single regime norm
    sr.garch.n <- CreateSpec(variance.spec = list(model=c("sGARCH")),
                             distribution.spec = list(distribution = c("norm")))
    fit.ml.n <- try(FitML(spec = sr.garch.n, data = window), silent=TRUE)
    
    risk <- try(Risk(fit.ml.n, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.SR[t1,c(1, 2)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(1, 2)] <- try(risk$ES, silent=TRUE)
    
    # Garch(1,1) single regime s-norm
    sr.garch.sn <- CreateSpec(variance.spec = list(model=c("sGARCH")),
                              distribution.spec = list(distribution = c("snorm")))
    fit.ml.sn <- try(FitML(spec = sr.garch.sn, data = window), silent=TRUE)
    
    risk <- try(Risk(fit.ml.sn, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.SR[t1,c(3, 4)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(3, 4)] <- try(risk$ES, silent=TRUE)    
    
    # Garch(1,1) single regime std
    sr.garch.std <- CreateSpec(variance.spec = list(model=c("sGARCH")),
                               distribution.spec = list(distribution = c("std")))
    fit.ml.std <- try(FitML(spec = sr.garch.std, data = window), silent=TRUE)
    
    risk <- try(Risk(fit.ml.std, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.SR[t1,c(5, 6)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(5, 6)] <- try(risk$ES, silent=TRUE)        
    
    
    # Garch(1,1) single regime sstd
    sr.garch.sstd <- CreateSpec(variance.spec = list(model=c("sGARCH")),
                                distribution.spec = list(distribution = c("sstd")))
    fit.ml.sstd <- try(FitML(spec = sr.garch.sstd, data = window), silent=TRUE)
    
    risk <- try(Risk(fit.ml.sstd, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.SR[t1,c(7, 8)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(7, 8)] <- try(risk$ES, silent=TRUE)   
    
    
    
    # Garch(1,1) single regime ged
    sr.garch.ged <- CreateSpec(variance.spec = list(model=c("sGARCH")),
                               distribution.spec = list(distribution = c("ged")))
    fit.ml.ged <- try(FitML(spec = sr.garch.ged, data = window), silent=TRUE)
    
    risk <- try(Risk(fit.ml.ged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.SR[t1,c(9,10)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(9,10)] <- try(risk$ES, silent=TRUE)   
    
    
    
    # Garch(1,1) single regime sged
    sr.garch.sged <- CreateSpec(variance.spec = list(model=c("sGARCH")),
                                distribution.spec = list(distribution = c("sged")))
    fit.ml.sged <- try(FitML(spec = sr.garch.sged, data = window), silent=TRUE)
    
    risk <- try(Risk(fit.ml.sged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.SR[t1,c(11,12)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(11,12)] <- try(risk$ES, silent=TRUE)  
    
    
    ###################################################################################################
    ###  GJR-Models
    
    # GJR(1,1) single regime norm
    sr.gjr.n <- CreateSpec(variance.spec = list(model=c("gjrGARCH")),
                           distribution.spec = list(distribution = c("norm")))
    fit.ml.gjr.n <- try(FitML(spec = sr.gjr.n, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.gjr.n, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.SR[t1,c(13,14)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(13,14)] <- try(risk$ES, silent=TRUE)
    
    # GJR(1,1) single regime snorm
    sr.gjr.sn <- CreateSpec(variance.spec = list(model=c("gjrGARCH")),
                            distribution.spec = list(distribution = c("snorm")))
    fit.ml.gjr.sn <- try(FitML(spec = sr.gjr.sn, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.gjr.sn, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.SR[t1,c(15,16)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(15,16)] <- try(risk$ES, silent=TRUE)
    
    # GJR(1,1) single regime std
    sr.gjr.std <- CreateSpec(variance.spec = list(model=c("gjrGARCH")),
                             distribution.spec = list(distribution = c("std")))
    fit.ml.gjr.std <- try(FitML(spec = sr.gjr.std, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.gjr.std, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.SR[t1,c(17,18)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(17,18)] <- try(risk$ES, silent=TRUE)
    
    
    # GJR(1,1) single regime sstd
    sr.gjr.sstd <- CreateSpec(variance.spec = list(model=c("gjrGARCH")),
                              distribution.spec = list(distribution = c("sstd")))
    fit.ml.gjr.sstd <- try(FitML(spec = sr.gjr.sstd, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.gjr.sstd, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.SR[t1,c(19,20)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(19,20)] <- try(risk$ES, silent=TRUE)
    
    
    # GJR(1,1) single regime ged
    sr.gjr.ged <- CreateSpec(variance.spec = list(model=c("gjrGARCH")),
                             distribution.spec = list(distribution = c("ged")))
    fit.ml.gjr.ged <- try(FitML(spec = sr.gjr.ged, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.gjr.ged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.SR[t1,c(21,22)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(21,22)] <- try(risk$ES, silent=TRUE)
    
    
    # GJR(1,1) single regime sged
    sr.gjr.sged <- CreateSpec(variance.spec = list(model=c("gjrGARCH")),
                              distribution.spec = list(distribution = c("sged")))
    fit.ml.gjr.sged <- try(FitML(spec = sr.gjr.sged, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.gjr.sged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.SR[t1,c(23,24)] <- try(risk$VaR, silent=TRUE)
    ES.SR[t1,c(23,24)] <- try(risk$ES, silent=TRUE)
    
    
    
    ########################################################################################
    ###### MARKOV SWITCHING
    ########################################################################################
    
    # MS-Garch(1,1) norm
    ms.garch.n <- CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                             distribution.spec = list(distribution = c("norm","norm")))
    fit.ml.ms.n <- try(FitML(spec = ms.garch.n, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.n, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.MS[t1,c(1, 2)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(1, 2)] <- try(risk$VaR, silent=TRUE)
    
    
    # MS-Garch(1,1) snorm
    ms.garch.sn <- CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                              distribution.spec = list(distribution = c("snorm","snorm")))
    fit.ml.ms.sn <- try(FitML(spec = ms.garch.sn, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.sn, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.MS[t1,c(3, 4)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(3, 4)] <- try(risk$VaR, silent=TRUE)
    
    # MS-Garch(1,1) std
    ms.garch.std <- CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                               distribution.spec = list(distribution = c("std","std")))
    fit.ml.ms.std <- try(FitML(spec = ms.garch.std, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.std, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.MS[t1,c(5, 6)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(5, 6)] <- try(risk$VaR, silent=TRUE)
    
    
    # MS-Garch(1,1) sstd
    ms.garch.sstd <- CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                                distribution.spec = list(distribution = c("sstd","sstd")))
    fit.ml.ms.sstd <- try(FitML(spec = ms.garch.sstd, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.sstd, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.MS[t1,c(7, 8)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(7, 8)] <- try(risk$VaR, silent=TRUE)
    
    
    # MS-Garch(1,1) ged
    ms.garch.ged <- CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                               distribution.spec = list(distribution = c("ged","ged")))
    fit.ml.ms.ged <- try(FitML(spec = ms.garch.ged, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.ged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.MS[t1,c(9, 10)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(9, 10)] <- try(risk$VaR, silent=TRUE)
    
    # MS-Garch(1,1) sged
    ms.garch.sged <- CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                                distribution.spec = list(distribution = c("sged","sged")))
    fit.ml.ms.sged <- try(FitML(spec = ms.garch.sged, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.sged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)# os alphas default são 1% e 5%
    VaR.MS[t1,c(11, 12)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(11, 12)] <- try(risk$VaR, silent=TRUE)
    
    
    ########################################################################################
    ### GJR-Models
    
    # MS-GJR(1,1) norm
    ms.gjr.n <- CreateSpec(variance.spec = list(model=c("gjrGARCH","gjrGARCH")),
                           distribution.spec = list(distribution = c("norm","norm")))
    fit.ml.ms.gjr.n <- try(FitML(spec = ms.gjr.n, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.gjr.n, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.MS[t1,c(13, 14)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(13, 14)] <- try(risk$VaR, silent=TRUE)
    
    # MS-GJR(1,1) snorm
    ms.gjr.sn <- CreateSpec(variance.spec = list(model=c("gjrGARCH","gjrGARCH")),
                            distribution.spec = list(distribution = c("snorm","snorm")))
    fit.ml.ms.gjr.sn <- try(FitML(spec = ms.gjr.sn, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.gjr.sn, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.MS[t1,c(15, 16)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(15, 16)] <- try(risk$VaR, silent=TRUE)
    
    
    # MS-GJR(1,1) std
    ms.gjr.std <- CreateSpec(variance.spec = list(model=c("gjrGARCH","gjrGARCH")),
                             distribution.spec = list(distribution = c("std","std")))
    fit.ml.ms.gjr.std <- try(FitML(spec = ms.gjr.std, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.gjr.std, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.MS[t1,c(17, 18)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(17, 18)] <- try(risk$VaR, silent=TRUE)
    
    # MS-GJR(1,1) sstd
    ms.gjr.sstd <- CreateSpec(variance.spec = list(model=c("gjrGARCH","gjrGARCH")),
                              distribution.spec = list(distribution = c("sstd","sstd")))
    fit.ml.ms.gjr.sstd <- try(FitML(spec = ms.gjr.sstd, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.gjr.sstd, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.MS[t1,c(19, 20)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(19, 20)] <- try(risk$VaR, silent=TRUE)
    
    
    # MS-GJR(1,1) ged
    ms.gjr.ged <- CreateSpec(variance.spec = list(model=c("gjrGARCH","gjrGARCH")),
                             distribution.spec = list(distribution = c("ged","ged")))
    fit.ml.ms.gjr.ged <- try(FitML(spec = ms.gjr.ged, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.gjr.ged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.MS[t1,c(21, 22)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(21, 22)] <- try(risk$VaR, silent=TRUE)
    
    # MS-GJR(1,1) sged
    ms.gjr.sged <- CreateSpec(variance.spec = list(model=c("gjrGARCH","gjrGARCH")),
                              distribution.spec = list(distribution = c("sged","sged")))
    fit.ml.ms.gjr.sged <- try(FitML(spec = ms.gjr.sged, data = window), silent=TRUE)
    risk <- try(Risk(fit.ml.ms.gjr.sged, do.its= FALSE, alpha = c(0.01, 0.025)), silent=TRUE)
    VaR.MS[t1,c(23, 24)] <- try(risk$VaR, silent=TRUE)
    ES.MS[t1,c(23, 24)] <- try(risk$VaR, silent=TRUE)
    
    
  
} else{ # I will use the obtained parameters
  
  print("Computing VaR and ES")
  
  resto <- t1%%10
  
  
  ########################################################################################
  ###### SINGLE REGIME
  ########################################################################################
  
  risk <- try(Risk(fit.ml.n, do.its= FALSE, alpha = c(0.01, 0.025), newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(1,2)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(1,2)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.sn, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(3,4)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(3,4)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.std, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(5,6)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(5,6)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.sstd, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(7,8)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(7,8)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(9,10)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(9,10)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.sged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(11,12)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(11,12)] <- try(risk$ES, silent=TRUE)
  
  ###############################################################################################
  ###  GJR- MODELS
  
  
  risk <- try(Risk(fit.ml.gjr.n, do.its= FALSE, alpha = c(0.01, 0.025), newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(13,14)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(13,14)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.gjr.sn, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(15,16)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(15,16)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.gjr.std, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(17,18)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(17,18)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.gjr.sstd, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(19,20)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(19,20)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.gjr.ged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(21,22)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(21,22)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.gjr.sged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(23,24)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(23,24)] <- try(risk$ES, silent=TRUE)
  
  
  
  ########################################################################################
  ###### MARKOV SWITCHING
  ########################################################################################
  
  
  risk <- try(Risk(fit.ml.ms.n, do.its= FALSE, alpha = c(0.01, 0.025), newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(1,2)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(1,2)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.sn, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(3,4)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(3,4)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.std, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(5,6)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(5,6)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.sstd, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(7,8)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(7,8)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.ged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(9,10)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(9,10)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.sged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(11,12)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(11,12)] <- try(risk$ES, silent=TRUE)
  
  ###############################################################################################
  #### GJR-Models
  
  
  risk <- try(Risk(fit.ml.ms.gjr.n, do.its= FALSE, alpha = c(0.01, 0.025), newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(13,14)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(13,14)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.gjr.sn, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(15,16)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(15,16)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.gjr.std, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(17,18)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(17,18)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.gjr.sstd, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(19,20)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(19,20)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.gjr.ged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(21,22)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(21,22)] <- try(risk$ES, silent=TRUE)
  
  risk <- try(Risk(fit.ml.ms.gjr.sged, do.its= FALSE, alpha = c(0.01, 0.025),newdata = resid[(t2-resto+1),(t2+1)]), silent=TRUE)
  VaR.SR[t1,c(23,24)] <- try(risk$VaR, silent=TRUE)
  ES.SR[t1,c(23,24)] <- try(risk$ES, silent=TRUE)
  
  
  
}
  print(Sys.time() - tempo)
}





###############     TABLE DESCRIPTIVE STATISTICS

## Whole Sample

ds1 <- basicStats(resid)%>% 
  slice(c(3,4,7,8,14,15,16))

round(ds1,3)

## Out-of-Sample
N <- dim(bvsp)[1]
WE <- 1500

ds2 <- basicStats(resid[WE+1:N]) %>% 
  slice(c(3,4,7,8,14,15,16))

round(ds2,3)






###############     TABLE SCORE FUNCTIONS


# Conditional mean ( y_t - e_t [resid])
cond_mean <- bvsp-resid # length 3715
cond_mean_os <- cond_mean[1501:3715,1] # os: out of sample

## We will build the risk measures here (that is, we will add the conditional mean)

# Loading my data
load("risk_measures.RData")

## addiing conditional mean in all 'measures'
VaR.SR.f <- VaR.SR + cond_mean_os
ES.SR.f <- ES.SR + cond_mean_os

### transforming the matrix in numeric
class(VaR.MS) <- "numeric"
class(ES.MS) <- "numeric"

#### editting the NA's
#which(is.na(VaR.MS[,17:18])!= FALSE) # NA's position
VaR.MS <- na.locf(VaR.MS)
ES.MS <- na.locf(ES.MS)

## Keeping adding the conditional mean
ES.MS.f <- ES.MS + cond_mean_os
VaR.MS.f <- VaR.MS + cond_mean_os


## Computing the mean, sd and score function of the measures
### loading the score functions
source("VaR_score_function.R")
source("ES_score_function.R")

#creating variables
scoresVaR <- rep(0,24)
scoresES <- rep(0,24)
med.VaR <- rep(0,24)
med.ES <- rep(0,24)
dp.VaR <- rep(0,24)
dp.ES <- rep(0,24)


for (i in 1:24) { # computing score function
  
  if (i<=12) {
    print(paste0(i," Single Regime Models"))
    
    #score functions
    scoresVaR[i] <- ScoreVaR(0.01,bvsp[1501:3715],-VaR.SR.f[,(i*2-1)])
    scoresES[i] <- ScoreES(0.025,bvsp[1501:3715],-VaR.SR.f[,(i*2)]
                           ,-ES.SR.f[,(i*2)])
    
    #VaR mean
    med.VaR[i] <- mean(VaR.SR.f[,(i*2-1)])
    
    #ES mean
    med.ES[i] <- mean(ES.SR.f[,(i*2)])
    
    #VaR standard deviation
    dp.VaR[i] <- sd(VaR.SR.f[,(i*2-1)])
    
    #ES standard deviation
    dp.ES[i] <- sd(ES.SR.f[,(i*2)])
    
  } 
  else{
    j <- i-12
    print(paste0(i," Markov Switching Models, column: ",j))
    
    #score functions
    scoresVaR[i] <- ScoreVaR(0.01,bvsp[1501:3715],-VaR.MS.f[,(j*2-1)])
    scoresES[i] <- ScoreES(0.025,bvsp[1501:3715],-VaR.MS.f[,(j*2)]
                           ,-ES.MS.f[,(j*2)])
    
    #VaR mean
    med.VaR[i] <- mean(VaR.MS.f[,(j*2-1)])
    
    #ES mean
    med.ES[i] <- mean(ES.MS.f[,(j*2)])
    
    #VaR standard deviation
    dp.VaR[i] <- sd(VaR.MS.f[,(j*2-1)])
    
    #ES standard deviation
    dp.ES[i] <- sd(ES.MS.f[,(j*2)])
  }
  
}

## Creating Latex table

#Extracted from latex

dists <- c( "GARCH_{norm} 
            GARCH_{snorm} 
            GARCH_{std} 
            GARCH_{sstd} 
            GARCH_{ged} 
            GARCH_{sged} 
            GJR_{norm} 
            GJR_{snorm} 
            GJR_{std} 
            GJR_{sstd} 
            GJR_{ged} 
            GJR_{sged} 
            MS-GARCH_{norm}
            MS-GARCH_{snorm} 
            MS-GARCH_{std} 
            MS-GARCH_{sstd}
            MS-GARCH_{ged} 
            MS-GARCH_{sged} 
            MS-GJR_{norm} 
            MS-GJR_{snorm} 
            MS-GJR_{std}
            MS-GJR_{sstd} 
            MS-GJR_{ged} 
            MS-GJR_{sged}")



a <- str_split(dists,"\n")
dists <- str_trim(a[[1]])


## creating the latex code - Table Score Functions 
for( i in 1:24){
  cat(dists[i]," & " ,round(med.VaR[i],4)," & ",
      round(dp.VaR[i],4)," & ", round(scoresVaR[i],4)," & ",
      round(med.ES[i],4)," & ", round(dp.ES[i],4)," & ",
      round(scoresES[i],4), " \\\\", "\n")
}





###############     TABLE TESTS

options(scipen = 999)## to remove the e+ number


#creating variables.
nomes_VaR <- c("UC Stat", "CC Stat","UC P-value", "CC P-value")
nomes_ES <- c("P-Value","Expc Exceed","Actual Exceed")
VaR_testes_critical <- matrix(ncol=6,nrow = 24)# uc, cc, critical and test
ES_testes_critical <- matrix(ncol=3,nrow = 24)# pvalue, expected and actual exceed
VaR_testes_decision <- matrix(ncol=2,nrow = 24)# 2 character columns
ES_testes_decision <- matrix(ncol=1,nrow = 24)# 1 character column
colnames(VaR_testes_critical) <- nomes_VaR
colnames(ES_testes_critical) <- nomes_ES



for (i in 1:24) {
  
  if (i<=12) {
    print(paste0(i," Single Regime Models"))
    
    #sVaR test
    a <- VaRTest(0.01,bvsp[1501:3715],VaR.SR.f[,(i*2-1)],conf.level = 0.99)
    VaR_testes_critical[i,1:6] <- c(a$uc.LRstat,a$cc.LRstat,
                                    a$uc.LRp,a$cc.LRp,a$uc.critical,
                                    a$cc.critical)
    VaR_testes_decision[i,1:2] <- c(a$uc.Decision,a$cc.Decision)
    
    
    
    b <- ESTest(alpha = 0.025, bvsp[1501:3715],
                VaR.SR.f[,(i*2)], ES.SR.f[,(i*2)],conf.level = 0.99)
    ES_testes_critical[i,1:3] <- c(b$p.value,b$expected.exceed,
                                   b$actual.exceed)
    ES_testes_decision[i,1] <- c(b$Decision)
    
    
  } 
  else{
    j <- i-12
    print(paste0(i," Markov Switching Models, column: ",j))
    
    a <- VaRTest(0.01,bvsp[1501:3715],VaR.MS.f[,(j*2-1)],conf.level = 0.99)
    VaR_testes_critical[i,1:6] <-c(a$uc.LRstat,a$cc.LRstat,
                                   a$uc.LRp,a$cc.LRp,a$uc.critical,
                                   a$cc.critical)
    VaR_testes_decision[i,1:2] <- c(a$uc.Decision,a$cc.Decision)
    
    
    
    b <- ESTest(alpha = 0.025, bvsp[1501:3715],
                VaR.MS.f[,(j*2)], ES.MS.f[,(j*2)],conf.level = 0.99)
    ES_testes_critical[i,1:3] <-  c(b$p.value,b$expected.exceed,
                                    b$actual.exceed)
    ES_testes_decision[i,1] <- c(b$Decision)
  }
  
}

VaR_testes_decision[,1] <- gsub("Reject", "Rej",VaR_testes_decision[,1])# shorthand


## creating the latex code - Table Tests
for( i in 1:24){
  cat(dists[i]," & " ,
      round(VaR_testes_critical[i,1],4)," & ", # UC Stat
      round(VaR_testes_critical[i,3],4)," & ", # UC p-value
      round(VaR_testes_critical[i,2],4), " & ", # CC Stat
      round(VaR_testes_critical[i,4],4), " & ", # CC p-value
      #VaR_testes_decision[i,1], " & ",
      #round(ES_testes_critical[i,2],4)," & ",
      #round(ES_testes_critical[i,3],4)," & ", # actual exceed
      round(ES_testes_critical[i,1],4), # p-value
      " \\\\", "\n")
}




###############     VaR PLOT


# building dfs
VaR.SR.df <- data.frame(index(bvsp[1501:3715]),bvsp[1501:3715],VaR.SR.f)
VaR.MS.df <- data.frame(index(bvsp[1501:3715]),bvsp[1501:3715],VaR.MS.f)
colnames(VaR.SR.df) <- c("Data","Bovespa",dists)
colnames(VaR.MS.df) <- c("Data","Bovespa",dists)


# Selecting the best models
df.sr.garch <-  VaR.SR.df %>% 
  select(Data, Bovespa, `GARCH_{norm}`, `GARCH_{snorm}`,`GARCH_{std}`)
df.ms.garch <-  VaR.MS.df %>% 
  select(`MS-GARCH_{sstd}`,`MS-GJR_{sstd}`) 

df <- cbind(df.sr.garch,df.ms.garch)%>%
  gather(key = "Models", value = "value", -Data)



#vectos of colors

colors <- c("Bovespa"= "black","GARCH_{norm}" = "red", "GARCH_{snorm}" = "green",
            "GARCH_{std}" = "purple", "MS_GARCH_{sstd}"= "gray", 
            "MS-GJR_{sstd}"= "orange")


VaR_plot <- ggplot(df, aes(x = Data, y = value)) + 
  geom_point(aes(color =Models),size =0.5)+ 
  scale_fill_manual(values = c("black","red","green","purple","gray",
                               "orange"))+
  labs(x= " ",y = "Return",title= " ",
       color = "Models")+
  scale_color_manual(values = colors)

#plot
VaR_plot + theme_linedraw()




# saving plot
ggsave("VaR_plot_linedraw.jpg", width = 7, height = 5, units = c("in"), dpi = 900)





###############     THE END     ###############





