#Script Dl2019_Cmdr_mcds_mrds.R...Megan C. Ferguson...1 August 2025
#
#  NOTES
#    1. This script was based on MCF scripts Dl2019_Cmdr_ddf.R and DL2019_quatre_mrds.R.
#
#    2. Required input files:
#       a. cmdr_dat.rds
#       b. no_ir.rds
#
#    3. Filters for ddf estimation:
      #       a. size >= 1 & is.na(size) == FALSE
      #       b. distance >= 0 & is.na(distance) == FALSE
      #       c. is.na(iBeauf) == FALSE & iBeauf <= 4
      #       d. PrimObs > 0 & is.na(PrimObs) == FALSE
      #       e. Behavior != "dead"
      #       f. Entry == "s on transect" | Entry == "s on search"
      #       g. Clino >= 0 & Clino <= 90 
      #       h. VisX.km > 0 & is.na(VisX.km) == FALSE
      #       i. Species == "beluga"
      #       j. is.na(best.Alt) == FALSE & best.Alt > 0 & best.Alt <= 2.0
      #       k. Also only used data from cmdr during: 2019-08-05 to 2019-08-27, excluding FOV-only
      #          flights (14, 17) or effort (portion of flt 409).
#    4. Data are left-truncated to account for the lower sighting probability very close to 
#       the aircraft. Ltrnc = 100 m (0.1 km). See HistDlPrimObs_MML_cmdr_Dl_sightings_by_25m.png.
#
#    5. MCDS Detection probability models were built below. Sea ice pct was not considered as a
#      covariate because only 24 out of 436 sightings had catIcePct > 10%. SkyCon was not considered
#      as a covariate because only 20 out of 436 sightings had "clear". Covariates considered 
#      included:
#                                                      size
#                                                      catsize
#                                                      catsizeGT2
#                                                      log10z
#                                                      catZ
#                                                      iBeauf
#                                                      f4Beauf
#                                                      best.Alt
#    6. Best mcds model based on Cmdr data, VisX.km > 0.0, max.w.Dl=1.140511: Dl.cmdr.hn.alt
#                   
#    7. The dataset used to build the mrds model consists of all imagery and aerial
#       observer data from all 2018 & 2019 Turbo Commander flights with the belly
#       port camera (bpc), plus non-bpc Turbo Commander flights during the period 
#       from 2019-08-05 to 2019-08-27. The dataset excludes imagery sightings that
#       could not be conclusively categorized as either a match or a mismatch (i.e.,
#       it excludes the "IR" sightings). For the mrds model, observer 1 = aerial 
#       observers, observer 2 = BPC imagery. 
#
#    8. The lowest AIC mrds model was no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf. 
#       However, during the parametric bootstrap procedure used to estimate 
#       variance in estimated BS beluga abundance (Supplement 6), 5 485 out of 
#       10 000 iterations resulted in estimates of transect detection probability 
#       equal to zero for the largest group size category. This was due to the 
#       large estimated standard error for the MR parameter for groups of $>$ 2 
#       belugas ($\beta_{S>2}=17.2$, $\hat{SE}(\beta_{S>2})=1217.1$). In the data 
#       subset used to build the MR models, there were relatively few sightings 
#       in the largest group size category: 662 sightings of individuals; 
#       169 sightings of two belugas; and 118 sightings of $>$ 2 belugas. We 
#       considered it outside the scope of reality for all of the largest beluga 
#       groups to be undetected. Therefore, we ultimately selected the MR model 
#       with distance and the 2-category group size variable (catsize) that distinguished 
#       sightings of single belugas from sightings with $>$ 1 beluga. This model 
#       was 2.89 AIC units away from the model with the lowest AIC value, 
#       and it resulted in an overall estimated transect detection probability, 
#       $\hat{p}_1(0)$, of 0.745. The name of the selected model is
#       no.ir.aba.bpc.hn.mrds.best.Alt.catsize.

  #setwd("C:\\Users\\megan.ferguson\\OneDrive - Biodiversity Research Institute\\Documents\\Belugas\\2019_BeaufortSeaBelugas\\Analysis\\FergusonEtal_BeaufortSea_Beluga_Nhat_2019")

  library(mrds)
  library(mgcv)
  library(lubridate)
  library(Distance) #needed for fcn gof_ds()

  #Input required data
    cmdr.dat <- readRDS("data//cmdr_dat.rds")

    #Examine histogram of perpendicular distances of beluga sightings from
    #primary observers to determine whether and how much to left truncate. Plot 
    #a variety of bin widths. 
        
      #Identify primary observer sightings of belugas
        idx <- which(cmdr.dat$PrimObs > 0 & cmdr.dat$Species == "beluga")
        length(idx) #518 sightings
      
      #Determine maximum sighting distance for creating histogram
        max(cmdr.dat$Pdist[which(cmdr.dat$PrimObs > 0 & cmdr.dat$Species == "beluga")], na.rm=TRUE) #3.873978 km is max sighting distance
        
      #Create and output histograms to figures 
        png(file="figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_15m.png",height=1000,width=2000,pointsize=24)
          hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.015), 
               main="Beluga Primary Observer Sightings from Turbo Commander",
               xlab="Pdist\n (15-m bins)")
        dev.off()
                
        png(file="figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_25m.png",height=1000,width=2000,pointsize=24)
          hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.025), 
               main="Beluga Primary Observer Sightings from Turbo Commander",
               xlab="Pdist\n (25-m bins)")
        dev.off()
        
        png(file="figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_30m.png",height=1000,width=2000,pointsize=24)
          hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.03), 
               main="Beluga Primary Observer Sightings from cmdr",
               xlab="Pdist\n (30-m bins)")
        dev.off()
        
        png(file="figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_50m.png",height=1000,width=2000,pointsize=24)
        hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.05), 
               main="Beluga Primary Observer Sightings from Turbo Commander",
               xlab="Pdist\n (50-m bins)")
        dev.off()

    #Create truncated distance variable.
       
      #Delete 100 m from perpendicular distances of all sightings   
        Pdist100 <- cmdr.dat$Pdist - 0.1 
      
      #Set distance of sightings < 100 m to NA. Will omit these sightings later.
        Pdist100[which(Pdist100 < 0)] <- NA 
          
      cmdr.dat$distance <- Pdist100
      #Ck
        summary(cmdr.dat$Pdist)
        summary(cmdr.dat$distance) # >= 0 or NA
    
        max(cmdr.dat$Pdist, na.rm=TRUE)
        max(cmdr.dat$distance, na.rm=TRUE)   
        max(cmdr.dat$Pdist, na.rm=TRUE) - max(cmdr.dat$distance, na.rm=TRUE) #0.1
            
    #Create detection function models using fixed strip width 
      #       a. size >= 1 & is.na(size) == FALSE
      #       b. distance >= 0 & is.na(distance) == FALSE
      #       c. is.na(iBeauf) == FALSE & iBeauf <= 4
      #       d. PrimObs > 0 & is.na(PrimObs) == FALSE
      #       e. Behavior != "dead"
      #       f. Entry == "s on transect" | Entry == "s on search"
      #       g. Clino >= 0 & Clino <= 90 
      #       h. VisX.km > 0 & is.na(VisX.km) == FALSE
      #       i. Species == "beluga"
      #       j. is.na(best.Alt) == FALSE & best.Alt > 0 & best.Alt <= 2.0
        
        Dl.idx <- which(cmdr.dat$size > 0 & 
                     is.na(cmdr.dat$size) == FALSE &
                     cmdr.dat$distance >= 0 & 
                     is.na(cmdr.dat$distance) == FALSE &
                     is.na(cmdr.dat$iBeauf) == FALSE &
                     cmdr.dat$iBeauf <= 4 &
                     cmdr.dat$PrimObs > 0 &
                     is.na(cmdr.dat$PrimObs) == FALSE &                     
                     (cmdr.dat$Behavior != "dead"  |
                     is.na(cmdr.dat$Behavior)==TRUE) &
                     (cmdr.dat$Entry == "s on transect" |
                     cmdr.dat$Entry == "s on search") &
                     cmdr.dat$Species == "beluga" &
                     cmdr.dat$Clino >= 0 &
                     cmdr.dat$Clino <= 90 &
                     cmdr.dat$VisX.km > 0 & 
                     is.na(cmdr.dat$VisX.km) == FALSE &
                     is.na(cmdr.dat$best.Alt) == FALSE & 
                     cmdr.dat$best.Alt > 0 & 
                     cmdr.dat$best.Alt <= 2.0)
        cmdr.Dl.dx.dat <- cmdr.dat[Dl.idx,] 
        cmdr.Dl.dx.dat$SkyCon <- as.vector(cmdr.Dl.dx.dat$SkyCon)
        #Ck
          length(Dl.idx)
          summary(cmdr.Dl.dx.dat$size) #24 max
          summary(cmdr.Dl.dx.dat$distance) #3.77398 max
          summary(cmdr.Dl.dx.dat$iBeauf)
          summary(cmdr.Dl.dx.dat$PrimObs)
          summary(as.factor(as.vector(cmdr.Dl.dx.dat$Behavior)))
          summary(as.factor(as.vector(cmdr.Dl.dx.dat$Entry)))
          summary(as.factor(as.vector(cmdr.Dl.dx.dat$Species))) #459 beluga
          summary(cmdr.Dl.dx.dat$Clino)
          summary(cmdr.Dl.dx.dat$VisX.km) # >0 and no NA
          summary(as.factor(cmdr.Dl.dx.dat$VisX.km))
          summary(cmdr.Dl.dx.dat$Dayt)
          summary(cmdr.Dl.dx.dat$best.Alt)
          summary(as.factor(cmdr.Dl.dx.dat$SkyCon)) #No NAs
          summary(cmdr.Dl.dx.dat$ArcLong) #>= -156.0
          summary(cmdr.Dl.dx.dat$Long100) #Should be no NAs
          
          #histogram of distance values in cmdr.Dl.dx.dat
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.025), by=0.025)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_25m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="magenta")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.05), by=0.05)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_50m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="pink")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.1), by=0.1)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_100m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.15), by=0.15)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_150m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="rosybrown")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.2), by=0.2)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_200m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink2")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.25), by=0.25)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_250m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink3")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.275), by=0.275)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_275m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink4")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.28), by=0.28)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_distance_280m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="rosybrown2")
            dev.off()

          #histogram of clino values in cmdr.Dl.dx.dat. Slight evidence of rounding.
            hist(cmdr.Dl.dx.dat$Clino, breaks=81)
            png(file="figures//cmdr//Hist_cmdr_Dldxdat_clino.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$Clino, breaks=81, main="MML cmdr")
            dev.off()
            
      #Omit farthest ~5% of sightings, which is approx. 23 sightings. In this case, that's 
      #anything with distance > 1.3 km
        n.omit <- .05*nrow(cmdr.Dl.dx.dat)
        n.omit
        
        dx.max <- sort(cmdr.Dl.dx.dat$distance, decreasing=TRUE)[round(n.omit)]
        dx.max #1.149897

        idx95 <- which(cmdr.Dl.dx.dat$distance < dx.max) 

        Dl.dx.x95 <- cmdr.Dl.dx.dat[idx95,]
          #CK
            length(idx95) #should be nrow(cmdr.Dl.dx.dat) - 23
            nrow(cmdr.Dl.dx.dat) - 23
        
            summary(as.factor(as.vector(Dl.dx.x95$Species))) #436
            summary(Dl.dx.x95$distance) #1.14051
            summary(Dl.dx.x95$Yr)
            summary(Dl.dx.x95$observer)
            summary(Dl.dx.x95$iBeauf)
            summary(as.factor(Dl.dx.x95$f4Beauf))
            summary(Dl.dx.x95$VisX.km) # >0
            summary(Dl.dx.x95$best.Alt)
            summary(as.factor(Dl.dx.x95$SkyCon))
            summary(as.factor(Dl.dx.x95$catIcePct)) #0=412; 1=24
            summary(as.factor(Dl.dx.x95$catsize))   #1=323; 2=113
            summary(as.factor(Dl.dx.x95$catsizeGT2)) #1=323; 2=79; gt2=34
            summary(as.factor(Dl.dx.x95$catsize5))   #1=323; 2=79; 3to5=23; gt5=11
            summary(as.factor(Dl.dx.x95$catsizeGT5)) #lteq5=425; gt5=11
            nrow(Dl.dx.x95) #436 
            class(Dl.dx.x95$Sample.Label) #character
            class(Dl.dx.x95$Region.Label) #character
            names(Dl.dx.x95)
            dim(Dl.dx.x95) #436 x 101
  
      #Multiple covariate distance sampling analysis of primary observer data.  
                    
        max.w.Dl <- max(Dl.dx.x95$distance) #should be 1.140511
        max.w.Dl
  
        #Create ddf
            
          #Determine key function for null model
            
              Dl.cmdr.hr  <-
                  ddf(dsmodel = ~mcds(key = "hr", formula = ~1),
                      data = Dl.dx.x95,
                      method = "ds") 
              
              Dl.cmdr.hn  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
                      data = Dl.dx.x95,
                      method = "ds") 
              
              Dl.cmdr.gamma  <-
                  ddf(dsmodel = ~mcds(key = "gamma", formula = ~1),
                      data = Dl.dx.x95,
                      method = "ds") 
              #CK
                plot(Dl.cmdr.hr)
                plot(Dl.cmdr.hn)
                plot(Dl.cmdr.gamma)
                  png(file="figures//cmdr//Dl_cmdr_null_modelfit.png", height=3000, width=1500,
                      pointsize=64)
                    par(mfrow=c(3,1))
                      plot(Dl.cmdr.hr)
                      plot(Dl.cmdr.hn)
                      plot(Dl.cmdr.gamma)
                  dev.off()
  
                summary(Dl.cmdr.hr)$aic
                summary(Dl.cmdr.hn)$aic
                summary(Dl.cmdr.gamma)$aic
                
              #Go with hn based on aic and plot

          #Proceed with covariate selection:      
          #            size
          #            catsize
          #            catsizeGT2
          #            log10z
          #            catZ
          #            iBeauf
          #            f4Beauf
          #            best.Alt

            #Round 1 hn: 

              Dl.cmdr.hn.size  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~size),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.catsize  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~catsize),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.catsizeGT2  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~catsizeGT2),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.log10z  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~log10z),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.catZ  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~catZ),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.iBeauf  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~iBeauf),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.f4Beauf  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~f4Beauf),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.alt  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                      data = Dl.dx.x95,
                      method = "ds") 

              #List Model Names and AIC Values
                Dl.dx.cmdr.names <- c("Dl.cmdr.hn.size",
                                     "Dl.cmdr.hn.catsize",
                                     "Dl.cmdr.hn.catsizeGT2",
                                     "Dl.cmdr.hn.log10z",
                                     "Dl.cmdr.hn.catZ",
                                     "Dl.cmdr.hn.iBeauf",
                                     "Dl.cmdr.hn.f4Beauf",
                                     "Dl.cmdr.hn.alt")

                Dl.dx.cmdr.aic <- c(summary(Dl.cmdr.hn.size)$aic,
                                                    summary(Dl.cmdr.hn.catsize)$aic,
                                                    summary(Dl.cmdr.hn.catsizeGT2)$aic,
                                                    summary(Dl.cmdr.hn.log10z)$aic,
                                                    summary(Dl.cmdr.hn.catZ)$aic,
                                                    summary(Dl.cmdr.hn.iBeauf)$aic,
                                                    summary(Dl.cmdr.hn.f4Beauf)$aic,
                                                    summary(Dl.cmdr.hn.alt)$aic)
                Dl.dx.cmdr.aic.sort <- sort(Dl.dx.cmdr.aic)
                Dl.dx.cmdr.min.aic <- Dl.dx.cmdr.aic.sort[1]
                Dl.dx.cmdr.model.idx.sort <- order(Dl.dx.cmdr.aic)
                Dl.dx.cmdr.best.model <- Dl.dx.cmdr.names[Dl.dx.cmdr.model.idx.sort[1]]
                #
                Dl.dx.cmdr.best.model    #Dl.cmdr.hn.alt
                Dl.dx.cmdr.min.aic       #-12.07678
                length(Dl.dx.cmdr.names)
                length(Dl.dx.cmdr.aic)
                
                Dl.dx.cmdr.aic.sort
                Dl.dx.cmdr.names[Dl.dx.cmdr.model.idx.sort]
                  #size, log10z, f4Beauf are the best covariates of each type

                x <- summary(Dl.cmdr.hn.alt)
                x$average.p * x$width #0.6560647 ESW
                
                plot(Dl.cmdr.hn.alt)
                summary(Dl.dx.x95$best.Alt)
                hist(Dl.dx.x95$best.Alt)
                
            #Round 2 hn: 

              Dl.cmdr.hn.alt.size  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + size),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.alt.log10z  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + log10z),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.cmdr.hn.alt.f4Beauf  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + f4Beauf),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              #List Model Names and AIC Values
                Dl.dx.cmdr.names <- c(Dl.dx.cmdr.names,
                                     "Dl.cmdr.hn.alt.size",
                                     "Dl.cmdr.hn.alt.log10z",
                                     "Dl.cmdr.hn.alt.f4Beauf")

                Dl.dx.cmdr.aic <- c(Dl.dx.cmdr.aic,
                                   summary(Dl.cmdr.hn.alt.size)$aic,
                                                    summary(Dl.cmdr.hn.alt.log10z)$aic,
                                                    summary(Dl.cmdr.hn.alt.f4Beauf)$aic)
                Dl.dx.cmdr.aic.sort <- sort(Dl.dx.cmdr.aic)
                Dl.dx.cmdr.min.aic <- Dl.dx.cmdr.aic.sort[1]
                Dl.dx.cmdr.model.idx.sort <- order(Dl.dx.cmdr.aic)
                Dl.dx.cmdr.best.model <- Dl.dx.cmdr.names[Dl.dx.cmdr.model.idx.sort[1]]
                #
                Dl.dx.cmdr.best.model    #Dl.cmdr.hn.alt is the best overall MCDS model for the cmdr
                Dl.dx.cmdr.min.aic       #-12.07678
                length(Dl.dx.cmdr.names)
                length(Dl.dx.cmdr.aic)
                
                Dl.dx.cmdr.aic.sort
                Dl.dx.cmdr.names[Dl.dx.cmdr.model.idx.sort]

                x <- summary(Dl.cmdr.hn.alt)
                x$average.p * x$width #0.6560647 ESW
                
                plot(Dl.cmdr.hn.alt)
                
                #output
                  saveRDS(Dl.cmdr.hn.alt, file="data//Dlcmdrhnalt.rds")
                
################
################
#  MRDS MODEL  #
################
################
                
  #Input required data
    no.ir <- readRDS("data//no_ir.rds")

    #Build mrds detection fcn models                

      no.ir.aba.bpc.hn.mrds  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
                                          mrmodel = ~glm(link = "logit", formula = ~ distance),
                                          data = no.ir,
                                          method = "trial",
                                          meta.data=list(width=max.w.Dl))
      
      no.ir.aba.bpc.hr.mrds  <- ddf(dsmodel = ~mcds(key = "hr", formula = ~1),
                                          mrmodel = ~glm(link = "logit", formula = ~ distance),
                                          data = no.ir,
                                          method = "trial",
                                          meta.data=list(width=max.w.Dl))
      #CK
        summary(no.ir.aba.bpc.hn.mrds)$AIC #243.7454
        summary(no.ir.aba.bpc.hr.mrds)$AIC #250.8129
        #Use hn model based on AIC and plots.
        
        gof_ds(list("ddf"=no.ir.aba.bpc.hn.mrds))
          #Distance sampling Cramer-von Mises test (unweighted)
          #Test statistic = 0.0641781 p-value = 0.787709
        gof_ds(list("ddf"=no.ir.aba.bpc.hr.mrds))
          #Distance sampling Cramer-von Mises test (unweighted)
          #Test statistic = 0.214347 p-value = 0.24132
        #No evidence of lack of fit
      
        x <- summary(no.ir.aba.bpc.hn.mrds)
        x
          #Estimate transect detection probability (p0) two ways to confirm
                
            #First, use methods from internal mrds function
                
              newdat <- no.ir
              newdat <- newdat[newdat$observer == 1 & newdat$detected == 1, ]
              nrow(newdat) #905 use only observer 1 data
              newdat$distance <- rep(0, length(newdat$distance)) #predict to y=0
              pred.at0 <- predict(no.ir.aba.bpc.hn.mrds$mr$mr, newdat, type = "response")
              summary(pred.at0) #0.7449  
              length(pred.at0) #905
                      
          ####
          #### MCF Ck pred.at0
          ####
                  
            b0 <- as.vector(no.ir.aba.bpc.hn.mrds$mr$mr$coefficients[1])
            b1 <- as.vector(no.ir.aba.bpc.hn.mrds$mr$mr$coefficients[2])
    
            distance.y <- 0
                              
            #Compute linear predictor
              lpred <- b0 + (b1*distance.y)
                        
            #Define logistic function for p_1|2 
            #See eq. 6.48, p. 154 of Laake and Borchers 2004
              logistic.p.fun <- function(lpred){
                p.val <- exp(lpred)/(1+exp(lpred))
                return(p.val)
              }
                            
            #Compute p(0). 
            #This is p(0)_1|2 = p(0)_1 under the assumption of point independence. See
            #p. 118 Laake and Borchers (2004) for definition of point independence.
              p0 <- logistic.p.fun(lpred)
              #CK
                summary(no.ir.aba.bpc.hn.mrds$mr)$average.p0.1
                summary(no.ir.aba.bpc.hn.mrds)$mr.summary$average.p0.1 #should equal the value above bc
                                                                       #it's the null model.
                p0 #0.7449
                summary(pred.at0)
                summary(p0 - pred.at0) #should be zero
              
          ####
          #### End MCF Ck pred.at0
          ####

      #Try including covariates f4Beauf, iBeauf, and size in mrmodel, and best.Alt 
      #in dsmodel. Setting dsmodel constant at best.Alt because that 
      #was the best mcds model selected above, and the dsmodel and mrmodel are fit
      #independently by mrds.
                      
        no.ir.aba.bpc.hn.mrds.best.Alt  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))
          
        no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + f4Beauf),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))
          
        no.ir.aba.bpc.hn.mrds.best.Alt.iBeauf  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + iBeauf),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))
          
        #  The mrmodel with size would not run. See note at top of script for full explanation.
          no.ir.aba.bpc.hn.mrds.best.Alt.size  <- try(ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + size),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl)))
          
        no.ir.aba.bpc.hn.mrds.best.Alt.catsize  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + catsize),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))
          
        no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + catsizeGT2),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))
          
        #The catsize10 model gives the following error:
        #Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
        #  contrasts can be applied only to factors with 2 or more levels
        #and I think it's because the BPC had no sightings gt10
          summary(as.factor(no.ir$catsizeGT10))
          idx <- which(no.ir$observer == 1 & no.ir$detected == 1)
          summary(as.factor(no.ir$catsizeGT10[idx]))
          idx <- which(no.ir$observer == 2 & no.ir$detected == 1)
          summary(as.factor(no.ir$catsizeGT10[idx]))
          
          no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT10  <- try(ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + catsizeGT10),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl)))
          
        no.ir.aba.bpc.hn.mrds.best.Alt.loggs  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + loggs),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))

      #Examine the model components.  
      #Here's the summary for the best mcds model found above:
            summary(Dl.cmdr.hn.alt)
#              
#            Summary for ds object
#            Number of observations :  436 
#            Distance range         :  0  -  1.140511 
#            AIC                    :  -12.07678 
#            Optimisation           :  mrds (nlminb) 
#            
#            Detection function:
#             Half-normal key function 
#            
#            Detection function parameters 
#            Scale coefficient(s): 
#                          estimate        se
#            (Intercept) -1.6454113 0.5520893
#            best.Alt     0.7976918 0.4318209
#            
#                                   Estimate         SE         CV
#            Average p             0.5752376  0.0232511 0.04041999
#            N in covered region 757.9476182 38.7559285 0.05113273
               
          summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf) #best.Alt coef is positive; f4Beauf coef is negative
          summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf)$mr.summary$average.p0.1 #0.734644 for the full model
          summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr)$average.p0.1 #0.7382614 for the mrmodel
                  
            #Work through sample computations
          
              #1. Derive summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf)$mr.summary$average.p0.1
              #   using Nhat_ds and Nhat_full
            
                x.full <- summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf)
                x.ds <- summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$ds)
                #The following three values are equal: 
                  x.ds$Nhat/x.full$Nhat 
                  sum(1/no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$ds$fitted)/x.full$Nhat
                  summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf)$mr.summary$average.p0.1
              
              #2. Derive summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr)$average.p0.1
              #   using Nhat_ds and Nhat_full
              #     NOTE: For method="io.fi" or method="trial.fi" if integrate=FALSE, 
              #     predict returns the value of the conditional detection probability and if 
              #     integrate=TRUE, it returns the average conditional detection probability by 
              #     integrating over x (distance) with respect to a uniform distribution. 
              #     The fitted values for full$mr (object classes "trial.fi" and "ddf") are calculated 
              #     using integrate=TRUE.
                  
                x.mr <- summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr)  
                newdat <- no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$data  
                newdat$distance <- 0
                newdat <- newdat[newdat$observer == 1 & newdat$detected == 1, ] #extract only data detected by observer 1
                pred.at0 <- predict(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr, newdat, type = "response")$fitted
                #The following two values are equal:
                  sum(pred.at0/no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr$fitted[,1])/x.mr$Nhat
                  summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr)$average.p0.1
                  
          #Try building no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf and see how it compares to
          #no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2, which is the current best model based on AIC.
                  
            no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf  <- ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                                        mrmodel = ~glm(link = "logit", formula = ~ distance + catsizeGT2 + f4Beauf),
                                        data = no.ir,
                                        method = "trial",
                                        meta.data=list(width=max.w.Dl))
            
          #mrds model summary  
          
            mrds.model.name <- c("no.ir.aba.bpc.hn.mrds.best.Alt",
                                   "no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf",
                                   "no.ir.aba.bpc.hn.mrds.best.Alt.iBeauf",
                                   "no.ir.aba.bpc.hn.mrds.best.Alt.catsize",
                                   "no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2",
                                   "no.ir.aba.bpc.hn.mrds.best.Alt.loggs",
                                   "no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf")
              
            mrds.model.AIC <- c(summary(no.ir.aba.bpc.hn.mrds.best.Alt)$AIC,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf)$AIC,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.iBeauf)$AIC,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsize)$AIC,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2)$AIC,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.loggs)$AIC,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf)$AIC)
            
          #Compute AIC weights 
            #Anderson et al. 2000. J. Wildlife Management 64(4): 912-923.
            #Burnham and Anderson. 2004. DOI: 10.1177/0049124104268644.
            #From Burnham and Anderson (2004, pg. 271): "The simple transformation 
            #exp(???delta.i/2), for i = 1, 2, . . . , R,
            #provides the likelihood of the model (Akaike 1981) given the data: L(gi|data). 
            #The relative likelihood of model i versus model j is L(gi|data)/L(gj|data); 
            #this is termed the evidence ratio, and it does not depend
            #on any of the other models under consideration.Without loss of generality,
            #we may assume that model gi is more likely than gj . Then, if
            #this evidence ratio is large (e.g., > 150 is quite large), model gj is a
            #poor model relative to model gi , based on the data."
                  
              aic.wt.fun <- function(aic.vec){
                      min.aic <- min(aic.vec)
                      delta.aic <- aic.vec - min.aic
                      x <- exp(-0.5*delta.aic)
                      sum.x <- sum(x)
                      aic.w <- x/sum.x
                      return(aic.w)
              }
            
              mrds.model.aic.wt <- aic.wt.fun(mrds.model.AIC)
              
            #Extract log-likelihood (lnl) from each model
               mrds.model.lnl <- c(no.ir.aba.bpc.hn.mrds.best.Alt$lnl,
                                            no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$lnl,
                                            no.ir.aba.bpc.hn.mrds.best.Alt.iBeauf$lnl,
                                            no.ir.aba.bpc.hn.mrds.best.Alt.catsize$lnl,
                                            no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2$lnl,
                                            no.ir.aba.bpc.hn.mrds.best.Alt.loggs$lnl,
                                            no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf$lnl)
                    
            #Compute delta AIC for each model
              mrds.model.min.aic <- min(mrds.model.AIC)
              mrds.model.delta.aic <- mrds.model.AIC - mrds.model.min.aic
    
            mrds.model.average.p0.1.full <- c(summary(no.ir.aba.bpc.hn.mrds.best.Alt)$mr.summary$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf)$mr.summary$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.iBeauf)$mr.summary$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsize)$mr.summary$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2)$mr.summary$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.loggs)$mr.summary$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf)$mr.summary$average.p0.1)
              
            mrds.model.average.p0.1.mr <- c(summary(no.ir.aba.bpc.hn.mrds.best.Alt$mr)$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.f4Beauf$mr)$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.iBeauf$mr)$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsize$mr)$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2$mr)$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.loggs$mr)$average.p0.1,
                                    summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf$mr)$average.p0.1)
                
            mrds.model.smry <- cbind.data.frame("model.name"=mrds.model.name,
                                                  "AIC"=mrds.model.AIC,
                                                  "delta.AIC"=mrds.model.delta.aic,
                                                  "AIC.wt"=mrds.model.aic.wt,
                                                  "lnl"=mrds.model.lnl,
                                                  "avg.p0.1.full"=mrds.model.average.p0.1.full,
                                                  "avg.p0.1.mr"=mrds.model.average.p0.1.mr)
            
            mrds.model.sort <- order(mrds.model.smry$AIC)
            mrds.model.smry <- mrds.model.smry[mrds.model.sort,]
            #Ck
              mrds.model.smry 
              sum(mrds.model.smry$AIC.wt) #should equal 1
              
#       The lowest AIC mrds model was no.ir.aba.bpc.hn.mrds.best.Alt.catsizeGT2.f4Beauf. 
#       However, during the parametric bootstrap procedure used to estimate 
#       variance in estimated BS beluga abundance (Supplement 6), 5 485 out of 
#       10 000 iterations resulted in estimates of transect detection probability 
#       equal to zero for the largest group size category. This was due to the 
#       large estimated standard error for the MR parameter for groups of $>$ 2 
#       belugas ($\beta_{S>2}=17.2$, $\hat{SE}(\beta_{S>2})=1217.1$). In the data 
#       subset used to build the MR models, there were relatively few sightings 
#       in the largest group size category: 662 sightings of individuals; 
#       169 sightings of two belugas; and 118 sightings of $>$ 2 belugas. We 
#       considered it outside the scope of reality for all of the largest beluga 
#       groups to be undetected. Therefore, we ultimately selected the MR model 
#       with distance and the 2-category group size variable (catsize) that distinguished 
#       sightings of single belugas from sightings with $>$ 1 beluga. This model 
#       was 2.89 AIC units away from the model with the lowest AIC value, 
#       and it resulted in an overall estimated transect detection probability, 
#       $\hat{p}_1(0)$, of 0.745. The name of the selected model is
#       no.ir.aba.bpc.hn.mrds.best.Alt.catsize.

  #Output mrds detection function to use in spatial modeling and abundance estimation
                
    saveRDS(no.ir.aba.bpc.hn.mrds.best.Alt.catsize,
            file="data//noirababpchnmrdsbestAltcatsize.rds")        
            
    
    
    
    
    
                  