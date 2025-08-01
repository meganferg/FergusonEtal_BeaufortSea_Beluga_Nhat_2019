#Script Dl2019_Otter_mcds.R...Megan C. Ferguson...1 August 2025
#
#  NOTES
#    1. This script was based on MCF script Dl2019_Otter_ddf_Deux.R
#    2. Required input files: ott_dat.rds
#    3. This script uses an equidistant conic map projection for Beaufort Sea study area:
#         aba.prj <- CRS('+proj=eqdc +lat_1=69.42d
#                              +lat_2=71.48d
#                              +lat_0=69.75d
#                              +lon_0=-129.2d
#                              +x_0=0
#                              +y_0=0')
#    4. VisX.km was defined as the farthest endpoint in the VisX range
#    5. Filters for ddf estimation (variables are designated below notes):
#       Otter: 
      #       a. size >= 1 & is.na(size) == FALSE
      #       b. distance >= 0 & is.na(distance) == FALSE
      #       c. is.na(iBeauf) == FALSE & iBeauf <= 4
      #       d. PrimObs > 0 & is.na(PrimObs) == FALSE
      #       e. Behavior != "dead"
      #       f. Entry == "s on transect" | Entry == "s on search"
      #       g. Clino >= 0 & Clino <= 90 
      #       h. VisX.km > 0 & is.na(VisX.km) == FALSE
      #       i. Species == "beluga"
      #       j. is.na(best.Alt) == FALSE & best.Alt > 0 & best.Alt <= 2.0 #That's 2000 ft.
#    6. Data are left-truncated to account for the lower sighting probability very close to the aircraft. 
#       otter Ltrnc = 75 m (0.075 km). See HistDlPrimObs_ott_Dl_sightings5km_by_25m.png.
#       I looked at histograms with various bin widths; liked the 25-min bin the best.
#    7. MCDS Detection probability models were built below. Sea ice pct was not considered as a
#      covariate because only 2 out of 64 sightings had catIcePct > 10%. SkyCon was not considered
#      as a covariate because only 2 out of 64 sightings had "clear". Also considered using Long100
#      and Region.Label, but (see below) decided to omit them because they are somewhat 
#      correlated with group size. Covariates considered 
#      included:
#                                                      Observer
#                                                      size
#                                                      catsize
#                                                      log10z
#                                                      catZ
#                                                      iBeauf
#                                                      f4Beauf
#                                                      best.Alt
#    8. Best model based on Otter data, VisX.km > 0.0, max.w.Dl=1.095209:
#                   Dl.ott.hn.alt

  #setwd("C:\\Users\\megan.ferguson\\OneDrive - Biodiversity Research Institute\\Documents\\Belugas\\2019_BeaufortSeaBelugas\\Analysis\\FergusonEtal_BeaufortSea_Beluga_Nhat_2019")

  library(mrds)

  #Designate file paths
    fig.dir <- "figures"
  
  #Input required functions, model objects, and data
    ott.dat <- readRDS("data//ott_dat.rds")
    
    #Edit variable names for later use in DSM script
      idx <- which(names(ott.dat) == "Sample.Label")
      names(ott.dat)[idx] <- "Strat.Sample.Label"
      
      idx <- which(names(ott.dat) == "seg.ID2.B4")
      names(ott.dat)[idx] <- "Sample.Label"

    #Examine histogram of perpendicular distances of Dl sightings in ott.dat to figure out ltrnc.

        length(which(ott.dat$PrimObs > 0 & ott.dat$Species == "beluga")) #102
        max(ott.dat$Pdist[which(ott.dat$PrimObs > 0 & ott.dat$Species == "beluga")], na.rm=TRUE) #2.33
        
        idx <- which(ott.dat$Pdist<=3.0 & ott.dat$PrimObs > 0 & 
                             ott.dat$Species == "beluga")
        length(idx) #99 sightings
        #
          png(file=file.path(fig.dir,"ott","HistDlPrimObs_ott_Dl_sightings_by_15m.png"),height=1000,width=2000,pointsize=24)
            hist(ott.dat$Pdist[idx], breaks=seq(from=0.0, to=2.5, by=0.015), 
                 main="Beluga Primary Observer Sightings from Otter",
                 xlab="Pdist\n (15-m bins)")
          dev.off()
                
          png(file=file.path(fig.dir,"ott","HistDlPrimObs_ott_Dl_sightings_by_25m.png"),height=1000,width=2000,pointsize=24)
            hist(ott.dat$Pdist[idx], breaks=seq(from=0.0, to=2.5, by=0.025), 
                 main="Beluga Primary Observer Sightings from Otter",
                 xlab="Pdist\n (25-m bins)")
          dev.off()

          png(file=file.path(fig.dir,"ott","HistDlPrimObs_ott_Dl_sightings_by_30m.png"),height=1000,width=2000,pointsize=24)
            hist(ott.dat$Pdist[idx], breaks=seq(from=0.0, to=2.5, by=0.03), 
                 main="Beluga Primary Observer Sightings from Otter",
                 xlab="Pdist\n (30-m bins)")
          dev.off()

          png(file=file.path(fig.dir,"ott","HistDlPrimObs_ott_Dl_sightings_by_50m.png"),height=1000,width=2000,pointsize=24)
            hist(ott.dat$Pdist[idx], breaks=seq(from=0.0, to=2.5, by=0.05), 
                 main="Beluga Primary Observer Sightings from Otter",
                 xlab="Pdist\n (50-m bins)")
          dev.off()

    #Create truncated distance variable.
      Pdist075 <- ott.dat$Pdist - 0.075
      Pdist075[which(Pdist075 < 0)] <- NA
          
      ott.dat$distance <- Pdist075
      #Ck
        summary(ott.dat$Pdist)
        summary(ott.dat$distance)                            # >= 0 or NA
    
        max(ott.dat$Pdist, na.rm=TRUE)
        max(ott.dat$distance, na.rm=TRUE)   
        max(ott.dat$Pdist, na.rm=TRUE) - max(ott.dat$distance, na.rm=TRUE) #0.075
            
    #Create detection function models
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
        
      #Dl            
        Dl.idx <- which(ott.dat$size > 0 & 
                     is.na(ott.dat$size) == FALSE &
                     ott.dat$distance >= 0 & 
                     is.na(ott.dat$distance) == FALSE &
                     is.na(ott.dat$iBeauf) == FALSE &
                     ott.dat$iBeauf <= 4 &
                     ott.dat$PrimObs > 0 &
                     is.na(ott.dat$PrimObs) == FALSE &                     
                     (ott.dat$Behavior != "dead"  |
                     is.na(ott.dat$Behavior)==TRUE) &
                     (ott.dat$Entry == "s on transect" |
                     ott.dat$Entry == "s on search") &
                     ott.dat$Species == "beluga" &
                     ott.dat$Clino >= 0 &
                     ott.dat$Clino <= 90 &
                     ott.dat$VisX.km > 0 & is.na(ott.dat$VisX.km) == FALSE &
                     is.na(ott.dat$best.Alt) == FALSE & ott.dat$best.Alt > 0 & ott.dat$best.Alt <= 2.0)
        ott.Dl.dx.dat <- ott.dat[Dl.idx,] 
        ott.Dl.dx.dat$SkyCon <- as.vector(ott.Dl.dx.dat$SkyCon)
        #Ck
          summary(ott.Dl.dx.dat$size) #30 max
          summary(ott.Dl.dx.dat$distance) #2.25758 max
          summary(ott.Dl.dx.dat$iBeauf)
          summary(ott.Dl.dx.dat$PrimObs)
          summary(as.factor(as.vector(ott.Dl.dx.dat$Behavior)))
          summary(as.factor(as.vector(ott.Dl.dx.dat$Entry)))
          summary(as.factor(as.vector(ott.Dl.dx.dat$Species))) #79 beluga
          summary(ott.Dl.dx.dat$Clino)
          summary(ott.Dl.dx.dat$VisX.km) # >0 and no NA
          summary(as.factor(ott.Dl.dx.dat$VisX.km))
          summary(ott.Dl.dx.dat$Dayt)
          summary(ott.Dl.dx.dat$best.Alt)
          summary(as.factor(ott.Dl.dx.dat$SkyCon)) #No NAs
          summary(ott.Dl.dx.dat$ArcLong) #>= -156.0
          summary(ott.Dl.dx.dat$Long100) #Should be no NAs
          
          #histogram of distance values in ott.Dl.dx.dat
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.025), by=0.025)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_25m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="magenta")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.05), by=0.05)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_50m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="pink")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.1), by=0.1)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_100m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="hotpink")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.15), by=0.15)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_150m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="rosybrown")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.2), by=0.2)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_200m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="hotpink2")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.25), by=0.25)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_250m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="hotpink3")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.275), by=0.275)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_275m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="hotpink4")
            dev.off()
            
            brk <- seq(from=0, to=(max(ott.Dl.dx.dat$distance)+0.28), by=0.28)
            hist(ott.Dl.dx.dat$distance, breaks=brk)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_distance_280m_bins.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$distance, breaks=brk, main="MML Otter", col="rosybrown2")
            dev.off()

          #histogram of clino values in ott.Dl.dx.dat. Slight evidence of rounding.
            hist(ott.Dl.dx.dat$Clino, breaks=81)
            png(file=file.path(fig.dir,"ott","Hist_ott_Dldxdat_clino.png"), height=500, width=1500,
                pointsize=20)
              hist(ott.Dl.dx.dat$Clino, breaks=81, main="MML Otter")
            dev.off()
            
      #Omit farthest ~5% of sightings, which is approx. 3 sightings. In this case, that's 
      #anything with distance > 1.3 km
        n.omit <- .05*nrow(ott.Dl.dx.dat)
        n.omit
        idx95 <- which(ott.Dl.dx.dat$distance <= 1.3) #I determined the 1.3 km cutpoint manually

        Dl.dx.x95 <- ott.Dl.dx.dat[idx95,]
          #CK
            length(idx95) #should be nrow(ott.Dl.dx.dat) - 3
            nrow(ott.Dl.dx.dat) - 3
        
            summary(as.factor(as.vector(Dl.dx.x95$Species))) #76
            summary(Dl.dx.x95$SurvNam)
            summary(Dl.dx.x95$distance) #1.09521 
            summary(Dl.dx.x95$Yr)
            summary(Dl.dx.x95$observer)
            summary(Dl.dx.x95$iBeauf)
            summary(as.factor(Dl.dx.x95$f4Beauf))
            summary(Dl.dx.x95$VisX.km) # >0
            summary(Dl.dx.x95$best.Alt)
            summary(as.factor(Dl.dx.x95$SkyCon))
            summary(as.factor(Dl.dx.x95$catIcePct)) #0=74; 1=2
            summary(as.factor(Dl.dx.x95$catsize))   #1=55; 2=21
            summary(as.factor(Dl.dx.x95$catsizeGT2)) #1=55; 2=15; gt2=6
            summary(as.factor(Dl.dx.x95$catsize5))   #1=55; 2=15; 3to5=5; gt5=1
            summary(as.factor(Dl.dx.x95$catsizeGT5)) #lteq5=75; gt5=1
            summary(as.factor(Dl.dx.x95$Observer))
            nrow(Dl.dx.x95) #76 
            class(Dl.dx.x95$Sample.Label) #character
            class(Dl.dx.x95$Region.Label) #character
            summary(as.factor(Dl.dx.x95$Sample.Label))
            names(Dl.dx.x95)
            dim(Dl.dx.x95) #76 x 101

      #MCDS analysis.  
                    
        max.w.Dl <- max(Dl.dx.x95$distance) #should be 1.095209
        max.w.Dl
  
          #Determine key function for null model
            
              Dl.ott.hr  <-
                  ddf(dsmodel = ~mcds(key = "hr", formula = ~1),
                      data = Dl.dx.x95,
                      method = "ds") 
              
              Dl.ott.hn  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
                      data = Dl.dx.x95,
                      method = "ds") 
              
              Dl.ott.gamma  <-
                  ddf(dsmodel = ~mcds(key = "gamma", formula = ~1),
                      data = Dl.dx.x95,
                      method = "ds") 
              #CK
                plot(Dl.ott.hr)
                plot(Dl.ott.hn)
                plot(Dl.ott.gamma)
                  png(file=file.path(fig.dir,"ott","Dl_ott_null_modelfit.png"), height=3000, width=1500,
                      pointsize=64)
                    par(mfrow=c(3,1))
                      plot(Dl.ott.hr)
                      plot(Dl.ott.hn)
                      plot(Dl.ott.gamma)
                  dev.off()
  
                summary(Dl.ott.hr)$aic
                summary(Dl.ott.hn)$aic
                summary(Dl.ott.gamma)$aic
                
              #hn has lowest AIC

          #Examine a few potential covariates   
                
            hist(Dl.dx.x95$Long100)
            summary(as.factor(Dl.dx.x95$Region.Label)) #east=19; west=57
            summary(as.factor(Dl.dx.x95$SkyCon)) #Only 2 "clear". Do not use.
            summary(as.factor(Dl.dx.x95$catIcePct)) #Only 2 >10%
                  
          #Proceed with covariate selection:      
          #            Observer
          #            size
          #            catsize
          #            log10z
          #            catZ
          #            iBeauf
          #            f4Beauf
          #            Long100
          #            Region.Label
          #            best.Alt

            #Round 1 hn: 

              Dl.ott.hn.Obs  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~Observer),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.size  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~size),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.catsize  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~catsize),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.log10z  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~log10z),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.catZ  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~catZ),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.iBeauf  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~iBeauf),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.f4Beauf  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~f4Beauf),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.Long100  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~Long100),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.Reg  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~Region.Label),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.alt  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt),
                      data = Dl.dx.x95,
                      method = "ds") 

              #List Model Names and AIC Values
                Dl.dx.ott.names <- c("Dl.ott.hn.Obs",
                                     "Dl.ott.hn.size",
                                     "Dl.ott.hn.catsize",
                                     "Dl.ott.hn.log10z",
                                     "Dl.ott.hn.catZ",
                                     "Dl.ott.hn.iBeauf",
                                     "Dl.ott.hn.f4Beauf",
                                     "Dl.ott.hn.Long100",
                                     "Dl.ott.hn.Reg",
                                     "Dl.ott.hn.alt")

                Dl.dx.ott.aic <- c(summary(Dl.ott.hn.Obs)$aic,
                                                    summary(Dl.ott.hn.size)$aic,
                                                    summary(Dl.ott.hn.catsize)$aic,
                                                    summary(Dl.ott.hn.log10z)$aic,
                                                    summary(Dl.ott.hn.catZ)$aic,
                                                    summary(Dl.ott.hn.iBeauf)$aic,
                                                    summary(Dl.ott.hn.f4Beauf)$aic,
                                                    summary(Dl.ott.hn.Long100)$aic,
                                                    summary(Dl.ott.hn.Reg)$aic,
                                                    summary(Dl.ott.hn.alt)$aic)
                Dl.dx.ott.aic.sort <- sort(Dl.dx.ott.aic)
                Dl.dx.ott.min.aic <- Dl.dx.ott.aic.sort[1]
                Dl.dx.ott.model.idx.sort <- order(Dl.dx.ott.aic)
                Dl.dx.ott.best.model <- Dl.dx.ott.names[Dl.dx.ott.model.idx.sort[1]]
                #
                Dl.dx.ott.best.model    #Dl.ott.hn.Long100
                Dl.dx.ott.min.aic       #-13.39712
                length(Dl.dx.ott.names)
                length(Dl.dx.ott.aic)
                
                Dl.dx.ott.aic.sort
                Dl.dx.ott.names[Dl.dx.ott.model.idx.sort]

                x <- summary(Dl.ott.hn.Long100)
                x$average.p * x$width #0.5935299 ESW
                
                plot(Dl.dx.x95$size, Dl.dx.x95$Long100)
                  #Largest group size was in the east
                  summary(Dl.ott.hn.Long100)
                  #Coefficient for Long100 is -7.5, meaning detection probability
                  #decreases towards the west. This is probably really a group size
                  #issue.
                  plot(Dl.ott.hn.Long100)
                  plot(Dl.ott.hn.size)
                  #Omit Long100 and Reg from analysis.
                
              #List Model Names and AIC Values, omitting Long100 and Region.Label
                Dl.dx.ott.names <- c("Dl.ott.hn.Obs",
                                     "Dl.ott.hn.size",
                                     "Dl.ott.hn.catsize",
                                     "Dl.ott.hn.log10z",
                                     "Dl.ott.hn.catZ",
                                     "Dl.ott.hn.iBeauf",
                                     "Dl.ott.hn.f4Beauf",
                                     "Dl.ott.hn.alt")

                Dl.dx.ott.aic <- c(summary(Dl.ott.hn.Obs)$aic,
                                                    summary(Dl.ott.hn.size)$aic,
                                                    summary(Dl.ott.hn.catsize)$aic,
                                                    summary(Dl.ott.hn.log10z)$aic,
                                                    summary(Dl.ott.hn.catZ)$aic,
                                                    summary(Dl.ott.hn.iBeauf)$aic,
                                                    summary(Dl.ott.hn.f4Beauf)$aic,
                                                    summary(Dl.ott.hn.alt)$aic)
                Dl.dx.ott.aic.sort <- sort(Dl.dx.ott.aic)
                Dl.dx.ott.min.aic <- Dl.dx.ott.aic.sort[1]
                Dl.dx.ott.model.idx.sort <- order(Dl.dx.ott.aic)
                Dl.dx.ott.best.model <- Dl.dx.ott.names[Dl.dx.ott.model.idx.sort[1]]
                #
                Dl.dx.ott.best.model    #Dl.ott.hn.alt
                Dl.dx.ott.min.aic       #-9.810406
                length(Dl.dx.ott.names)
                length(Dl.dx.ott.aic)
                
                Dl.dx.ott.aic.sort
                Dl.dx.ott.names[Dl.dx.ott.model.idx.sort]
                  #catsize, f4Beauf, log10z were the best performers out of 
                  #the options for size, beauf, and depth

                x <- summary(Dl.ott.hn.alt)
                x$average.p * x$width #0.5983228 ESW
                
                plot(Dl.ott.hn.alt)
                summary(Dl.dx.x95$best.Alt)
                hist(Dl.dx.x95$best.Alt)
                
            #Round 2 hn: 

              Dl.ott.hn.alt.Obs  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + Observer),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.alt.catsize  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + catsize),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.alt.log10z  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + log10z),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              Dl.ott.hn.alt.f4Beauf  <-
                  ddf(dsmodel = ~mcds(key = "hn", formula = ~best.Alt + f4Beauf),
                      data = Dl.dx.x95,
                      method = "ds") 
                
              #List Model Names and AIC Values
                Dl.dx.ott.names <- c(Dl.dx.ott.names,
                                     "Dl.ott.hn.alt.Obs",
                                     "Dl.ott.hn.alt.catsize",
                                     "Dl.ott.hn.alt.log10z",
                                     "Dl.ott.hn.alt.f4Beauf")

                Dl.dx.ott.aic <- c(Dl.dx.ott.aic,
                                   summary(Dl.ott.hn.alt.Obs)$aic,
                                                    summary(Dl.ott.hn.alt.catsize)$aic,
                                                    summary(Dl.ott.hn.alt.log10z)$aic,
                                                    summary(Dl.ott.hn.alt.f4Beauf)$aic)
                Dl.dx.ott.aic.sort <- sort(Dl.dx.ott.aic)
                Dl.dx.ott.min.aic <- Dl.dx.ott.aic.sort[1]
                Dl.dx.ott.model.idx.sort <- order(Dl.dx.ott.aic)
                Dl.dx.ott.best.model <- Dl.dx.ott.names[Dl.dx.ott.model.idx.sort[1]]
                #
                Dl.dx.ott.best.model    #Dl.ott.hn.alt.log10z
                Dl.dx.ott.min.aic       #-10.35483
                length(Dl.dx.ott.names)
                length(Dl.dx.ott.aic)
                
                Dl.dx.ott.aic.sort #<2 AIC units between Dl.ott.hn.alt.log10z and Dl.ott.hn.alt
                Dl.dx.ott.names[Dl.dx.ott.model.idx.sort]

                x <- summary(Dl.ott.hn.alt)
                x$average.p * x$width #0.5983228 ESW
                
                plot(Dl.ott.hn.alt)
                
  #Save to .rds file
    saveRDS(Dl.dx.x95, file="data//Dl_ott_dat.rds")
    saveRDS(Dl.ott.hn.alt, file=("data//Dlotthnalt.rds"))
                
                
                
                
                

                