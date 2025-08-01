#Script Dl2019_Cmdr_ddf.r...Megan C. Ferguson...16 March 2022
#
#  NOTES
#    1. This script was based on Dl2019_Otter_ddf.r
#    2. Required input files:
#       a. workspace from DL2019_Segs_Deux.r
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
      #       e. SurvNam == "MML_cmdr" 
      #       f. Behavior != "dead"
      #       g. Entry == "s on transect" | Entry == "s on search"
      #       h. Clino >= 0 & Clino <= 90 
      #       i. VisX.km > 0 & is.na(VisX.km) == FALSE
      #       j. Species == "beluga"
      #       k. is.na(best.Alt) == FALSE & best.Alt > 0 & best.Alt <= 2.0 #That's 2000 ft. Zero
      #                                                                     sightings were > 2.0.
      #       l. Also only used data from ABA survey period: 2019-08-05 to 2019-08-27, excluding FOV-only
      #          flights (14, 17) or effort (portion of flt 409).
#    6. Data are left-truncated to account for the lower sighting probability very close to 
#       the aircraft. Ltrnc = 100 m (0.1 km). See HistDlPrimObs_MML_cmdr_Dl_sightings_by_25m.png.
#       I looked at histograms with various bin widths; liked the 25-min bin the best.
#    7. Variable strip width models (i.e., with an int.range matrix) were not used because the 
#       width of the "beluga zone" is very narrow.
#    8. MCDS Detection probability models were built below. Sea ice pct was not considered as a
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
#        a. Best model based on Cmdr data, VisX.km > 0.0, max.w.Dl=1.140511: Dl.cmdr.hn.alt
#                   

  #rm(list=ls())
  
  library(sp)
  library(maptools)
  library(rgeos)
  library(rgdal)
  library(raster)   
  library(mrds)
  library(mgcv)

  #Input required functions, model objects, and data
    load("Functions//xtract_flightlines_SL_byBeaufVis_ReturnList_20181223.Rdata")
    load("Functions//LTsamplFCNS.Rdata")
    load("Functions//ClinoArcDist.Rdata")

    #Examine histogram of perpendicular distances of Dl sightings in ott.dat to figure out ltrnc.
    #  NOTE: NSDL analysis used 50-m truncation; the hists used to determine that 
    #        truncation distance had only 50-m bins.

      #hist of PrimObs sightings out to 4.0 km. Use the 25-m bin width and choose to ltrnc 100 m.
        
        idx <- which(cmdr.dat$PrimObs > 0 & cmdr.dat$Species == "beluga")
        length(idx) #518 sightings
        max(cmdr.dat$Pdist[which(cmdr.dat$PrimObs > 0 & cmdr.dat$Species == "beluga")], na.rm=TRUE) #3.873978
        
        hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.025))
        dev.off()
        #
          png(file="Figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_15m.png",height=1000,width=2000,pointsize=24)
            hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.015), 
                 main="Beluga Primary Observer Sightings from cmdr",
                 xlab="Pdist\n (15-m bins)")
          dev.off()
                
          png(file="Figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_25m.png",height=1000,width=2000,pointsize=24)
            hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.025), 
                 main="Beluga Primary Observer Sightings from cmdr",
                 xlab="Pdist\n (25-m bins)")
          dev.off()

          png(file="Figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_30m.png",height=1000,width=2000,pointsize=24)
            hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.03), 
                 main="Beluga Primary Observer Sightings from cmdr",
                 xlab="Pdist\n (30-m bins)")
          dev.off()

          png(file="Figures//cmdr//HistDlPrimObs_MML_cmdr_Dl_sightings_by_50m.png",height=1000,width=2000,pointsize=24)
            hist(cmdr.dat$Pdist[idx], breaks=seq(from=0.0, to=4.0, by=0.05), 
                 main="Beluga Primary Observer Sightings from cmdr",
                 xlab="Pdist\n (50-m bins)")
          dev.off()

    #Create truncated distance variable.
      Pdist100 <- cmdr.dat$Pdist - 0.1
      Pdist100[which(Pdist100 < 0)] <- NA
          
      cmdr.dat$distance <- Pdist100
      #Ck
        summary(cmdr.dat$Pdist)
        summary(cmdr.dat$distance)                            # >= 0 or NA
    
        max(cmdr.dat$Pdist, na.rm=TRUE)
        max(cmdr.dat$distance, na.rm=TRUE)   
        max(cmdr.dat$Pdist, na.rm=TRUE) - max(cmdr.dat$distance, na.rm=TRUE) #0.1
            
        nrow(cmdr.dat)
        length(unique(cmdr.dat$object)) #same as above
        summary(cmdr.dat$observer) #1
        summary(cmdr.dat$detected) #1
        names(cmdr.dat) #gs should be "size" and Pdist100 should be "distance"

    #Create detection function models using fixed strip width 
      #       a. size >= 1 & is.na(size) == FALSE
      #       b. distance >= 0 & is.na(distance) == FALSE
      #       c. is.na(iBeauf) == FALSE & iBeauf <= 4
      #       d. PrimObs > 0 & is.na(PrimObs) == FALSE
      #       e. SurvNam == "MML_Cmdr" 
      #       f. Behavior != "dead"
      #       g. Entry == "s on transect" | Entry == "s on search"
      #       h. Clino >= 0 & Clino <= 90 
      #       i. VisX.km > 0 & is.na(VisX.km) == FALSE
      #       j. Species == "beluga"
      #       k. is.na(best.Alt) == FALSE & best.Alt > 0 & best.Alt <= 2000
        
      #Check how many Dl sightings have best.Alt > 2.0  
        test.idx <- which(cmdr.dat$size > 0 & 
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
                     cmdr.dat$VisX.km > 0 & is.na(cmdr.dat$VisX.km) == FALSE &
                     is.na(cmdr.dat$best.Alt) == FALSE & cmdr.dat$best.Alt > 2.0)

      #Dl            
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
                     cmdr.dat$VisX.km > 0 & is.na(cmdr.dat$VisX.km) == FALSE &
                     is.na(cmdr.dat$best.Alt) == FALSE & cmdr.dat$best.Alt > 0 & cmdr.dat$best.Alt <= 2.0)
        cmdr.Dl.dx.dat <- cmdr.dat[Dl.idx,] 
        cmdr.Dl.dx.dat$SkyCon <- as.vector(cmdr.Dl.dx.dat$SkyCon)
        #Ck
          length(test.idx) #0
          length(Dl.idx)
          summary(cmdr.Dl.dx.dat$size) #24 max
          summary(cmdr.Dl.dx.dat$distance) #3.77398 max
          summary(cmdr.Dl.dx.dat$iBeauf)
          summary(cmdr.Dl.dx.dat$Depth)
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
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_25m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="magenta")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.05), by=0.05)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_50m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="pink")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.1), by=0.1)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_100m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.15), by=0.15)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_150m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="rosybrown")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.2), by=0.2)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_200m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink2")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.25), by=0.25)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_250m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink3")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.275), by=0.275)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_275m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="hotpink4")
            dev.off()
            
            brk <- seq(from=0, to=(max(cmdr.Dl.dx.dat$distance)+0.28), by=0.28)
            hist(cmdr.Dl.dx.dat$distance, breaks=brk)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_distance_280m_bins.png", height=500, width=1500,
                pointsize=20)
              hist(cmdr.Dl.dx.dat$distance, breaks=brk, main="MML cmdr", col="rosybrown2")
            dev.off()

          #histogram of clino values in cmdr.Dl.dx.dat. Slight evidence of rounding.
            hist(cmdr.Dl.dx.dat$Clino, breaks=81)
            png(file="Figures//cmdr//Hist_cmdr_Dldxdat_clino.png", height=500, width=1500,
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
            summary(as.factor(Dl.dx.x95$SurvNam))
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
            summary(as.factor(Dl.dx.x95$Observer))
            nrow(Dl.dx.x95) #436 
            class(Dl.dx.x95$Sample.Label) #character
            class(Dl.dx.x95$Region.Label) #character
            summary(as.factor(Dl.dx.x95$Sample.Label))
            names(Dl.dx.x95)
            dim(Dl.dx.x95) #436 x 101
  
      #max.w analysis.  
                    
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
                  png(file="Figures//cmdr//Dl_cmdr_null_modelfit.png", height=3000, width=1500,
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
                  #size, log10z, f4Beauf

                x <- summary(Dl.cmdr.hn.alt)
                x$average.p * x$width #0.6560647
                
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
                Dl.dx.cmdr.best.model    #Dl.cmdr.hn.alt
                Dl.dx.cmdr.min.aic       #-12.07678
                length(Dl.dx.cmdr.names)
                length(Dl.dx.cmdr.aic)
                
                Dl.dx.cmdr.aic.sort
                Dl.dx.cmdr.names[Dl.dx.cmdr.model.idx.sort]

                x <- summary(Dl.cmdr.hn.alt)
                x$average.p * x$width #0.6560647
                
                plot(Dl.cmdr.hn.alt)
                
  #Save to .Rdata file
    
    Dl.cmdr.dat <- Dl.dx.x95
    save(Dl.cmdr.dat, Dl.cmdr.hn.alt, file="Output//ddf//cmdr_ddf.Rdata")
                
                
                
                
                
                

                