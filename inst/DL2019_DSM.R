#Script DL2019_DSM.R...Megan C. Ferguson...21 August 2025
#
#  0. This script was based on MCF scripts DL2019_dsm_trois.R and 
#     DL2019_dsm_trois_stats.R.
#
#  1. Required input files:
#     a. aba_prj.rds
#     b. shor_SP_prj.rds
#     c. tx_10k4.rds
#     d. noirababpchnmrdsbestAltcatsize.rds
#     e. no_ir.rds
#     f. inla_mesh2sp_sf_fun.R
#     g. MCF_MuVarAvg_rr_plot_sdmTMB_tweedie.R
#     h. Dl_ott_dat.rds
#     i. seg_dat.rds
#     j. dl2019_mrds_dat.rds
#
#  2. Uses an equidistant conic map projection for for eastern Beaufort and Amundsen 
#     68.9 (min coastal Long in actual ABA tx flown) to 72 (max offshore Long in 
#     actual ABA tx flown) N, 118.9 (T.ID BCB401) to 139.5 (T.ID BCB442) W. 
#
#    aba.prj <- CRS('+proj=eqdc +lat_1=69.42d
#                              +lat_2=71.48d
#                              +lat_0=69.75d
#                              +lon_0=-129.2d
#                              +x_0=0
#                              +y_0=0')
#
#   Use projection with units=km to build spde.
#
#    aba.prj.km <- CRS('+proj=eqdc +lat_1=69.42d
#                              +lat_2=71.48d
#                              +lat_0=69.75d
#                              +lon_0=-129.2d
#                              +x_0=0
#                              +y_0=0
#                              +units=km')
#
#  3. Used Horvitz-Thompson estimator as response variable in DSMs because the 
#     mr model to estimate transect detection probability included observation
#     level covariate catsize
#
#  4. Availability bias correction factors were computed in DL2019_pAvail.R.
#     For the DL2019 analysis, I used p.avail = 0.56 (CF=1.8). All aircraft
#     surveyed at the same target speed of 115 kts = 212.98 km/hr = 59.16111 m/s.
#     Their aircraft-specific estimates of p.avail are shown below; note that they
#     are the same out to 1 decimal place (2 decimal places after rounding), which 
#     is probably sufficient given the amount of uncertainty and expected variability 
#     in the p.avail estimate.
#
#       w.cmdr = 1140 m
#       tiv.cmdr = 19.3 sec
#       p.avail.cmdr = 1.752267
#
#       w.ott = 1100 m
#       tiv.ott = 18.6 sec
#       p.avail.ott = 1.798703

  #setwd("C:\\Users\\megan.ferguson\\OneDrive - Biodiversity Research Institute\\Documents\\Belugas\\2019_BeaufortSeaBelugas\\Analysis\\FergusonEtal_BeaufortSea_Beluga_Nhat_2019")

  library(sp)
  library(sf)
  library(tidyverse)
  library(ggplot2)
  library(dsm)
  library(mrds)
  library(mgcv)
  library(sdmTMB)
  library(sdmTMBextra)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
  library(classInt)

  #Define file paths
    dat.dir <- "data"
    fig.dir <- "figures"

  #Input objects 
    aba.prj <- readRDS("data//aba_prj.rds")
    shor.SP.prj <- readRDS("data//shor_SP_prj.rds")
    tx.10k4 <- readRDS("data//tx_10k4.rds")
    no.ir <- readRDS("data//no_ir.rds")
    seg.dat <- readRDS("data//seg_dat.rds")
    dl2019.mrds.dat <- readRDS("data//dl2019_mrds_dat.rds")

  #Define aba.prj in km
    aba.prj.km <- CRS('+proj=eqdc +lat_1=69.42d
                              +lat_2=71.48d
                              +lat_0=69.75d
                              +lon_0=-129.2d
                              +x_0=0
                              +y_0=0
                              +units=km')
    
  #Convert shoreline to sf
    shor.sf.prj <- st_as_sf(shor.SP.prj) %>% st_transform(aba.prj.km)
    
  #Input map of North America from Rnaturalearth
    NoAm <- ne_countries(continent = "North America", scale = "large")
    NoAm.prj <- st_transform(NoAm, aba.prj.km)
    
    #Preliminary crop to shor.sf.prj, for plotting 
      NoAm.crop <- sf::st_crop(NoAm.prj, st_bbox(shor.sf.prj))
      #CK
        ggplot() + geom_sf(data=NoAm.prj) +
          geom_sf(data=NoAm.crop, fill="violet")
   
###############
###############
# Set p.avail #
###############
###############
  #See script header notes for further details. 
    CF.avail <- 1.8    
    p.avail <- 1/CF.avail
    #CK
      p.avail #Should be 0.555....

  #Input objects from mrds scripts and define average primary p(0), which is the 
  #trackline detection probability that will be used in the offset for the Otter
  #in order to build the dsm.
    no.ir.aba.bpc.hn.mrds.best.Alt.catsize <- readRDS(file.path(dat.dir,"noirababpchnmrdsbestAltcatsize.rds"))
    p0.1 <- summary(no.ir.aba.bpc.hn.mrds.best.Alt.catsize$mr)$average.p0.1
    p0.1 #0.7451675
      #Work through some calculations to get insight into model components
    
        full.model <- no.ir.aba.bpc.hn.mrds.best.Alt.catsize
        newdat <- no.ir
        newdat <- newdat[newdat$observer == 1 & newdat$detected == 1, ]
        nrow(newdat) #905 use only observer 1 data
        summary(newdat$distance) #nonzero
        newdat$distance <- rep(0, length(newdat$distance)) #predict to y=0
        
        #Nhat_full
          summary(full.model)$Nhat #2078.131
          sum(1/full.model$fitted) #2078.131
          summary(full.model$ds)$Nhat/summary(full.model)$mr.summary$average.p0.1 #2078.131
          full.model$Nhat #2078.131
          
        #Nhat_ds
          summary(full.model$ds)$Nhat #1527.235
          sum(1/full.model$ds$fitted) #1527.235
          full.model$ds$Nhat #1527.235

        #Nhat_mr
          summary(full.model$mr)$Nhat #949.8361
          sum(1/full.model$mr$fitted) #949.8361
          full.model$mr$Nhat #949.8361
          
        #average.p_full
          summary(full.model)$mr.summary$n1/summary(full.model)$Nhat #0.4354874
          summary(full.model)$mr.summary$n1/full.model$Nhat #0.4354874
          summary(full.model)$average.p #0.4354874
          summary(full.model)$mr.summary$average.p0.1*summary(full.model$ds)$average.p #0.4354874

        #average.p_ds
          summary(full.model)$mr.summary$n1/summary(full.model$ds)$Nhat #0.5925743
          summary(full.model)$mr.summary$n1/full.model$ds$Nhat #0.5925743
          summary(full.model$ds)$average.p #0.5925743
          
        #average.p_mr  
          summary(full.model)$mr.summary$n1/summary(full.model$mr)$Nhat #0.952796
          summary(full.model)$mr.summary$n1/full.model$mr$Nhat #0.952796
          summary(full.model$mr)$average.p #0.952796
        
        #average.p0.1_full   
          summary(full.model)$mr.summary$average.p0.1 #0.7349078
          full.model$ds$Nhat/full.model$Nhat #0.7349078

        #average.p0.1_mr  
          summary(full.model$mr)$average.p0.1 #0.7451675
          
          pred.mr.at0 <- predict(full.model$mr, newdat, type = "response")
          length(pred.mr.at0$fitted) #905
          unique(pred.mr.at0$fitted) #2 values
          summary(pred.mr.at0$p1 - pred.mr.at0$fitted) #identical
          
          pdot.mr <- full.model$mr$fitted
          length(pdot.mr) #905
          unique(pdot.mr) #2 values
          test <- predict(full.model$mr, integrate=TRUE)
          test.new <- predict(full.model$mr, data=newdat, integrate=TRUE)
          summary(test$fitted - pdot.mr) #identical
          summary(test.new$fitted - pdot.mr) #identical
          
          (1/full.model$mr$Nhat)*sum(pred.mr.at0$fitted/pdot.mr) #0.7451675
          
            ####
            #### MCF Ck pred.mr.at0
            ####
              
              b0 <- as.vector(full.model$mr$par[1])
              b1 <- as.vector(full.model$mr$par[2])
              b2 <- as.vector(full.model$mr$par[3])

              distance.y <- 0
                          
              #Compute linear predictor
              
                #For catsize = 1
                  lpred.sz1 <- b0 + (b1*distance.y)
                  
                #For catsize = 2
                  lpred.sz2 <- b0 + (b1*distance.y) + b2
                  
              #Define logistic function for p_1|2 
              #See eq. 6.48, p. 154 of Laake and Borchers 2004
                logistic.p.fun <- function(lpred){
                  p.val <- exp(lpred)/(1+exp(lpred))
                  return(p.val)
                }
                            
              #Compute p(0). 
              #This is p(0)_1|2 = p(0)_1 under the assumption of point independence. See
              #p. 118 Laake and Borchers (2004) for definition of point independence.
                p0.sz1 <- logistic.p.fun(lpred.sz1)
                p0.sz2 <- logistic.p.fun(lpred.sz2)
                #CK
                  unique(pred.mr.at0$fitted)
                  p0.sz1 #0.6778702
                  p0.sz2 #0.9008536

            ####
            #### End MCF Ck pred.mr.at0
            ####
          
        #full.model_fitted
          full.fitted <- pred.mr.at0$fitted*predict(full.model$ds, esw=FALSE)$fitted
          summary(full.model$fitted - full.fitted) #identical
          
        #class  
          class(full.model) #"trial" "ddf"
          class(full.model$mr$mr) #"glm" "lm" 
          
  #Input objects from Otter detection function code Dl2019_Otter_mcds.R
    Dl.ott.dat <- readRDS(file.path(dat.dir, "Dl_ott_dat.rds"))
    Dl.ott.hn.alt <- readRDS(file.path(dat.dir, "Dlotthnalt.rds"))
    #Ck
      summary(Dl.ott.hn.alt)

      x <- summary(Dl.ott.hn.alt)
      x$average.p * x$width #0.5983228

  #Fcn to create sf object from inla mesh
    source("inst/inla_mesh2sp_sf_fun.R")       
    
  #Define necessary functions
    
    #Fcn to create plot of tweedie variance vs. mean
      source("inst/MCF_MuVarAvg_rr_plot_sdmTMB_tweedie.R")

    #Define function to plot sdmTMB model predictions, observations, counts per sample,
    #and effort. Return area-integrated abundance (Nhat).

      plot.sdmTMB.pred.obs.map.fun <- function(m, preddat, grd, obs, counts, 
                                               effort.sf, coast, pth, title, NSIM){
        
        #This function outputs m, obs, counts, effort.sf, and grd to .Rdata object.  
        #It outputs a map of predictions, obs, and effort to a .png file.
        #It returns epsilon-corrected estimate of area-integrated Nhat.
        #
        # m: sdmTMB model object
        # preddat: Data.frame with covariates for each cell in grd, 
        #          in the order of grd
        # grd: sf object with prediction grid
        # obs: sf object with locations of observations to plot
        # counts: sf object with number of individuals on each sample unit, saved
        #         to seg.ind column
        # effort.sf: sf object with effort
        # pth: Filepath for outputting figure
        # title: Title for map and base name for output file.
        # NSIM: number of simulations from the posterior covariance mtx
        #
        # The first estimate returned is the area-integrated abundance estimate
        # without epsilon detransformation bias correction factor. The second
        # estimate returned is the area-integrated abundance estimate with
        # epsilon detransformation bias correction factor.
        
        #Predict number of individuals across grid 
              
          grd$pred.ind <- predict(m, newdata=preddat, type="response",
                                  offset=log(preddat$a))$est 
          Nhat <- sum(grd$pred.ind, na.rm=TRUE)
          #CK
            grd #should be an sf object with pred.ind as one column
            Nhat #area-integrated abundance w/o epsilon detransformation bias correction
        
        #Plot. 
  
          ggplot(data=grd) + geom_sf(aes(fill=pred.ind)) +
            geom_sf(data=effort.sf, color="darkgray", linewidth=0.2) +
            geom_sf(data=counts[which(counts$seg.ind > 0),], pch=16, color="white", aes(size=seg.ind), alpha=0.5) +
            geom_sf(data=obs, pch=16, color="orange", size=0.5) +
            ggtitle(title) +
            theme(text = element_text(size = 10),
                  axis.text.x = element_text(angle = 45, hjust=1),
                  legend.key.size = unit(0.5, 'cm'), 
                  legend.text = element_text(size=8),
                  plot.title = element_text(size = 10))
                
          ggsave(filename=file.path(pth, paste0(title,".png")), dpi="retina", height=10,
                 width=15, units="cm")
        
        #Pretty plot    
          cst.crop <- sf::st_crop(coast, st_bbox(grd))
          
          #Divide pred.ind using jenks natural breaks
            class.int <- classIntervals(grd$pred.ind, 
                                            n = 10, 
                                            style = "jenks")
            #CK
              summary(grd$pred.ind)
              class.int$brks
        
          ggplot(data=grd) + geom_sf(data = cst.crop, fill = 'darkgray') +
            geom_sf(aes(fill=pred.ind)) +
            scale_fill_continuous(breaks = class.int$brks, 
                                  guide=guide_colourbar(barheight = unit(8, "cm"),
                                                        title="Predicted \nNumber of Belugas"),
                                  labels=function(x) format(round(x, 0), nsmall = 0)) + #zero decimal places
            geom_sf(data=effort.sf, color="darkgray", linewidth=0.2) +
            geom_sf(data=obs, pch=16, color="orange", size=0.25, alpha=0.5) +
            ggtitle(title) +
            theme(text = element_text(size = 10),
                  axis.text.x = element_text(angle = 45, hjust=1),
                  legend.key.size = unit(0.5, 'cm'), 
                  legend.text = element_text(size=8),
                  plot.title = element_text(size = 10))

          ggsave(filename=file.path(pth, paste0(title,"_jenks.png")), dpi="retina", 
                 height=10, width=15, units="cm")
          
        #Apply epsilon detransformation bias correction to Nhat estimate. Note
        #that get_index returns a data frame with a columns for time, estimate, 
        #lower and upper confidence intervals, log estimate, and standard error 
        #of the ******log estimate******.
          
          pred2 <- predict(m, newdata=preddat, return_tmb_object=TRUE) 
          Nhat2.vec <- try(get_index(pred2, bias_correct=TRUE, area=preddat$a,
                                 level = 0.95))
            Nhat2.vec
            
          if(class(Nhat2.vec)[1] == "try-error"){
            Nhat2 <- NA
          } else {
            Nhat2 <- Nhat2.vec$est
          }
          #CK
            Nhat2
            
        #Simulate from the fitted model. Draw new random effects; fix fixed
        #effects at MLE; no observation error. Returns a matrix of nrow(data)
        #by nsim representing the estimates of the response variable.
        #Use to derive uncertainty on predictions.
          pred3 <- try(simulate(m, 
                               newdata=preddat, 
                               offset=log(preddat$a), 
                               nsim=NSIM,
                               type = "mle-mvn", # fixed effects at MLE values and random effect MVN draws
                               mle_mvn_samples = "multiple", # take an MVN draw for each sample
                               observation_error = FALSE, # do not include observation error
                               model = NA
                               ))  
            
        #Output stuff to .Rdata
          save(m, grd, obs, counts, effort.sf, Nhat2.vec, pred3,
               file=file.path(pth, paste0(title,".Rdata")))
          
        #Return the area-integrated predicted number of individuals
          return(c("Nhat"=Nhat,"Nhat.bc"=Nhat2))
      
      }
      
  #Import analytical study area boundary and create 1.25-km buffer s.t. it includes all Dl
  #sightings on transect. (See Dl2019_Otter_mcds.R and Dl2019_Cmdr_mcds_mrds.R for max.w.Dl.)
  #Omit land.
    dl2019.sf <- st_read(file.path(dat.dir, "Dl2019_surveyed_polygon.shp"))
    
    #Transform to aba.prj.km
      dl2019.sf.prj <- st_transform(dl2019.sf, aba.prj.km) 
    
    #1.25-km buffer, bc max.w.Dl is ~1 km for both cmdr & ott    
      buff.sf <- st_buffer(dl2019.sf.prj, dist = 1.250) 

    #Omit land from buffered study area

      st_crs(NoAm.prj) == st_crs(buff.sf)

      dl2019.buff.sf <- st_difference(buff.sf, st_union(st_combine(NoAm.prj)))
      #Ck
        ggplot() + geom_sf(data=dl2019.buff.sf, fill="cyan") +
          geom_sf(data=NoAm.crop, fill="purple") +
          geom_sf(data=buff.sf, color="red", fill=NA, linewidth=1.0) 

  #Compute Horvitz-Thompson estimator for each observation and each segment. 
  #For each ddf, use predict() with new data corresponding to the valid sighting 
  #locations.
      
    dl2019.mrds.dat$avg.p <- NA  
      
    ddf1 <- no.ir.aba.bpc.hn.mrds.best.Alt.catsize
    ddf2 <- Dl.ott.hn.alt

    i1 <- which(dl2019.mrds.dat$ddfobj == 1)
    dat1 <- dl2019.mrds.dat[i1,]
    dl2019.mrds.dat$avg.p[i1] <- predict(ddf1, newdata = dat1, esw=FALSE)$fitted
    
    i2 <- which(dl2019.mrds.dat$ddfobj == 2)
    dat2 <- dl2019.mrds.dat[i2,]
    dl2019.mrds.dat$avg.p[i2] <- predict(ddf2, newdata = dat2, esw=FALSE)$fitted
    #CK
      summary(dl2019.mrds.dat$avg.p[i2])
      summary(predict(ddf2, newdata = dat2, esw=FALSE, integrate=TRUE)$fitted) #Identical
    
    #Multiply avg.p2 by catsz-specific value of p0.1 from ddf1
    
      b0 <- as.vector(ddf1$mr$par[1])
      b1 <- as.vector(ddf1$mr$par[2])
      b2 <- as.vector(ddf1$mr$par[3])

      distance.y <- 0
                    
      #Compute linear predictor
              
        #For catsize = 1
          lpred.1 <- b0 + (b1*distance.y)
                  
        #For catsize = 2
          lpred.2 <- b0 + (b1*distance.y) + b2
          
      #Define logistic function for p_1|2 
      #See eq. 6.48, p. 154 of Laake and Borchers 2004
        logistic.p.fun <- function(lpred){
          p.val <- exp(lpred)/(1+exp(lpred))
                  return(p.val)
        }
                            
      #Compute p(0) and multiply it by $avg.p[i2] computed above. 
      #This is p(0)_1|2 = p(0)_1 under the assumption of point independence. See
      #p. 118 Laake and Borchers (2004) for definition of point independence.
        p0.gs1 <- logistic.p.fun(lpred.1)
        p0.gs2 <- logistic.p.fun(lpred.2)

        idx <- which(dl2019.mrds.dat$ddf == 2 & dl2019.mrds.dat$catsize == "1")
        dl2019.mrds.dat$avg.p[idx] <- p0.gs1*dl2019.mrds.dat$avg.p[idx]

        idx <- which(dl2019.mrds.dat$ddf == 2 & dl2019.mrds.dat$catsize == "2")
        dl2019.mrds.dat$avg.p[idx] <- p0.gs2*dl2019.mrds.dat$avg.p[idx]
        #CK
          p0.gs1
          p0.sz1 #same
          
          p0.gs2
          p0.sz2 #same
          
          p0.1
          
          summary(dl2019.mrds.dat$avg.p[which(dl2019.mrds.dat$ddfobj == 1)])
          #CK
            summary(predict(ddf1, newdata = dat1, esw=FALSE, integrate=TRUE)$fitted) #Identical
            
          summary(dl2019.mrds.dat$avg.p[which(dl2019.mrds.dat$ddfobj == 2)])
          
        #Calc HT estimator, Nht  
          dl2019.mrds.dat$Nht <- dl2019.mrds.dat$size/(dl2019.mrds.dat$avg.p*p.avail)
          
  #Compute Nht per segment. 
    seg.Nht <- sapply(1:nrow(seg.dat), function(i){
                idx <- which(dl2019.mrds.dat$Sample.Label == seg.dat$Sample.Label[i])
                Nht.i <- sum(dl2019.mrds.dat$Nht[idx])
                return(Nht.i)
    })
    seg.dat$seg.Nht <- seg.Nht
    #Ck
      length(seg.Nht) #1170
      nrow(seg.dat) #1170
      summary(seg.dat$seg.Nht)
      sum(seg.dat$seg.Nht) #2462.496
      sum(dl2019.mrds.dat$Nht) #2462.496

   #Define offsets as area in km2, using:
   # i. aircraft-specific w
   # ii. seg.km

    w1 <- ddf1$meta.data$width
    w2 <- ddf2$meta.data$width
  
    i1 <- which(seg.dat$ddfobj == 1)
    i2 <- which(seg.dat$ddfobj == 2)
    
    seg.dat$a <- NA
    
    seg.dat$a[i1] <- seg.dat$seg.km[i1]*w1
    seg.dat$a[i2] <- seg.dat$seg.km[i2]*w2
    #CK
      w1
      w2
      summary(seg.dat$seg.km)
      summary(seg.dat$a)

   #Plot stuff. 
      
       dl2019.mrds.sf.LL <- st_as_sf(x=dl2019.mrds.dat,
                                     coords=c("ArcLong", "ArcLat"),
                                     crs="EPSG:4326")
       dl2019.mrds.sf <- st_transform(dl2019.mrds.sf.LL, aba.prj.km)
       
       seg.dat.sf <- st_as_sf(x=seg.dat, coords=c("x","y"), crs=aba.prj) %>% 
         st_transform(crs=aba.prj.km)
       
       tx.10k4.sf <- st_as_sf(x=tx.10k4, crs=aba.prj) %>% 
         st_transform(crs=aba.prj.km)
       
        #Extract tx that are located within dl2019.buff.sf      
          tx.10k4.in <- st_intersection(dl2019.buff.sf, tx.10k4.sf)           

       ggplot() + geom_sf(data=dl2019.buff.sf, fill="cyan") +
          geom_sf(data=NoAm.crop, fill="purple") + 
          geom_sf(data=tx.10k4.sf, linewidth=1, color="white") +
          geom_sf(data=tx.10k4.in, linewidth=1) +
          geom_sf(data=seg.dat.sf, col="red", size=0.25) +
          geom_sf(data=dl2019.mrds.sf, col="yellow", size=0.25) 

  #Create prediction grid
  #
  #  x: x-coord of cell midpoint
  #  y: y-coord of cell midpoint

     #Create hex grid with 10-km cells
      set.seed(1865068920) #To get the same grid every time
      hex.grid <- st_make_grid(dl2019.buff.sf, 
                               cellsize = 10, #10 km
                               what = "polygons",
                               square = FALSE)
       #CK
         ggplot() + geom_sf(data=dl2019.buff.sf, fill="cyan") +
                    geom_sf(data=NoAm.crop, fill="purple") + 
                    geom_sf(data=hex.grid, fill=NA, color="red")
      
    #Intersect cells with dl2019.buff.sf
       hex4pred <- st_intersection(dl2019.buff.sf, hex.grid)

       hex4pred.pts <- st_centroid(hex4pred) #extract centroids of cells
       hex4pred$x <- st_coordinates(hex4pred.pts)[,1]
       hex4pred$y <- st_coordinates(hex4pred.pts)[,2]
       hex4pred$a <- as.numeric(st_area(hex4pred)) #units = km^2

      #Check for and remove cells with ~zero area
        idx <- which(hex4pred$a < 0.001)
        if(length(idx) > 0){
          print(nrow(hex4pred))
          
          hex4pred <- hex4pred[-idx,]
          
          print(nrow(hex4pred))
        }
        #CK
          ggplot(data=dl2019.buff.sf) + geom_sf(fill="pink") +
            geom_sf(data=hex.grid, fill=NA, color="red") +
            geom_sf(data=hex4pred, fill="lightblue") +
            geom_sf(data=hex4pred.pts, cex=0.25)
          
          summary(hex4pred)
          class(hex4pred) #sf
          
          dim(hex4pred) #1901 x 5
          
      #Convert to df for prediction
        hex4pred.df <- st_drop_geometry(hex4pred)[,-1] #Omit initial "NA_" column
        summary(hex4pred.df)

  #Create mesh for spde knots
          
    #First, need to extract seg.dat that are located within dl2019.buff.sf      
      seg.dat.in <- st_intersection(dl2019.buff.sf, seg.dat.sf)    
      names(seg.dat.in)
      seg.dat.in <- seg.dat.in[,-1]
      #CK
        dim(seg.dat.sf) #1170 x 11
        dim(seg.dat.in) #823 x 11
        
        names(seg.dat.sf)
        names(seg.dat.in)
          
         ggplot(data=dl2019.buff.sf) + geom_sf(color="pink", linewidth=2) +
           geom_sf(data=hex4pred, fill=NA) +
           geom_sf(data=seg.dat.in, cex=1)
         
         ggplot() + geom_sf(data=tx.10k4.sf) + 
           geom_sf(data=dl2019.buff.sf, color="pink", linewidth=2, fill=NA) +
           geom_sf(data=seg.dat.in, cex=1, color="red") 

    #Set up initial knot locations and inner and outer boundaries for mesh   
      loc <- st_coordinates(seg.dat.in)
      loc.df <- as.data.frame(loc)
      names(loc.df) <- c("x","y")
      boundary = INLA::inla.nonconvex.hull(loc)      
      boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35) 
          
    #Mesh 70: max.edge 70
            
      max.edg70 <- 70
            
      mesh70 <- INLA::inla.mesh.2d(
        loc=loc,
        boundary=list(boundary,boundary2),
        max.edge=c(1,2)*max.edg70,
        cutoff = max.edg70/5
      )
      
      mesh70.sdmTMB <- make_mesh(loc.df, c("x","y"), mesh=mesh70) #for sdmTMB 
      mesh70.sf <- inla.mesh2sp.sf.fun(mesh70, crs=aba.prj.km) #for ggplot
      
      #Crop North America to mesh70 for plotting and to later create spde.bar
        st_crs(mesh70.sf[[2]]) == st_crs(NoAm.prj)
      
        NoAm.coast <- sf::st_crop(NoAm.prj, st_bbox(mesh70.sf[[2]]))
        #CK
            summary(mesh70)
            plot(mesh70)
            
            ggplot() + geom_sf(data=dl2019.buff.sf, color="magenta", linewidth=2) +
              geom_sf(data=NoAm.coast, fill="purple") + 
              geom_sf(data=mesh70.sf[[2]], color="cyan", fill=NA) +
              geom_sf(data=mesh70.sf[[4]], color="blue", fill=NA) +
              geom_sf(data=hex4pred, fill=NA) +
              geom_sf(data=tx.10k4.in, linewidth=1) +
              geom_sf(data=seg.dat.in, size=0.25, color="gray") +
              geom_sf(data=seg.dat.in, col="red", size=0.25) +
              geom_sf(data=dl2019.mrds.sf, col="yellow", size=0.25)
  
  #Create dataframe for saving model highlights
    # df <- data.frame(Doubles=double(),
    #               Ints=integer(),
    #               Factors=factor(),
    #               Logicals=logical(),
    #               Characters=character(),
    #               stringsAsFactors=FALSE)    
      M.df <- data.frame(typ = character(),
                                 BS = character(),
                                 Max.Edg = double(),
                                 n.re = integer(),
                                 edf = double(),
                                 cAIC = double(),
                                 n.flags = integer(),
                                 n.no0.segs = integer(),
                                 n.obs = integer(),
                                 n.ind = integer(),
                                 p.avail = double(),
                                 Nhat.plugin = double(),
                                 Nhat.bc = double(),
                               stringsAsFactors=FALSE)
        
    #Initialize with fake values
      M.df[1,] <- cbind.data.frame(typ="ZZ",
                                           BS="ZZ",
                                           Max.Edg=-99.9,
                                           n.re=0,
                                           edf=-99.9,
                                           cAIC=-99.9,
                                           n.flags=-1,
                                           n.no0.segs=0,
                                           n.obs=0,
                                           n.ind=0,
                                           p.avail,
                                           Nhat.plugin=-99.9,
                                           Nhat.bc=-99.9)

  #Build dsms
    
    #sdmTMB requires data be class data.frame    
      seg.dat.in.df <- st_drop_geometry(seg.dat.in)  
      seg.dat.in.df$x <- st_coordinates(seg.dat.in)[,1]
      seg.dat.in.df$y <- st_coordinates(seg.dat.in)[,2]
      #CK
        summary(seg.dat.in.df)
        
###################    
###################    
###################   
                
M.idx <- 1
typ <- "tweedie"
BS <- "spde.bar"
Max.Edg <- max.edg70

###################
      
    #Purely spatial models using spde with barriers and evaluating different 
    #observation likelihoods. Anisotropy not allowed with barrier mesh.

    #Add barrier mesh component. Setting range_fraction to 0.1 means that the
    #spatial range will be assumed to be 0.1 over land compared to over water.
  
      bspde.70 <- sdmTMBextra::add_barrier_mesh(
                              mesh70.sdmTMB, NoAm.coast, 
                              range_fraction = 0.1,
                              proj_scaling = 1, plot = TRUE
                            )
      #CK
        mesh.df.water <- bspde.70$mesh_sf[bspde.70$normal_triangles, ]
        mesh.df.land <- bspde.70$mesh_sf[bspde.70$barrier_triangles, ]
        
        ggplot(NoAm.coast) +
          geom_sf() +
          geom_sf(data = mesh.df.water, size = 1, color = "blue") +
          geom_sf(data = mesh.df.land, size = 1, color = "green") +
          geom_sf(data = dl2019.buff.sf, color = "purple", fill=NA)

    #Fit barrier spde model
      M <- try(sdmTMB(
                  seg.Nht ~ 1,
                  data = seg.dat.in.df,
                  mesh = bspde.70,
                  family = tweedie(),
                  offset = log(seg.dat.in.df$a),
                  spatial = "on",
                  ))
      
    #Evaluate model fit 

      TITLE <- paste(typ,BS,Max.Edg,"catsize",sep="_")
      Fnam <- file.path(fig.dir, "sdmTMB", TITLE)   
      
      if(class(M)[1] == "try-error"){
                  
        sink(file=file.path(fig.dir, "sdmTMB", paste0(TITLE,"_TryError.txt")))
          print(M)
        sink()  
        
      } else {  
      
        #Plot predictions, observations, and effort; output model and predictions
          Nhat <- plot.sdmTMB.pred.obs.map.fun(m=M, 
                                        preddat=hex4pred.df, 
                                        grd=hex4pred, 
                                        obs=dl2019.mrds.sf, #observations as points
                                        counts=seg.dat.in, #seg.ind
                                        effort.sf=tx.10k4.in, 
                                        coast=NoAm.coast,
                                        pth=file.path(fig.dir, "sdmTMB"),
                                        title=TITLE,
                                        NSIM=1000)
          
        #Check simulation-based randomized quantile residuals using DHARMa  
  
          M.sim <- simulate(M, nsim=1000, type = "mle-mvn")
          
          png(paste(Fnam, "_DHARMa_sdmTMB.png",sep=""),bg="white",
                      height=640, width=960, units="px", pointsize=16)
          
            dharma_residuals(M.sim, M, test_uniformity = TRUE)
            
          dev.off()  
          
        #Plot tweedie variance vs. mean
          sdmTMB.tweedie.mu.var.plot(m=M, dat=seg.dat.in.df, 
                                     OFFSET=log(seg.dat.in.df$a), 
                                     fnam=Fnam)
          
        #Output stuff
            
          M.cAIC <- try(cAIC(M, what="cAIC"))
          if(class(M.cAIC)[1] == "try-error"){
            M.cAIC <- NA
          }
          
          M.edf <- try(cAIC(M, what="EDF"))
          if(class(M.edf)[1] == "try-error"){
            M.edf <- NA
          }
          
          M.sanity <- sanity(M)
          
          sink(file=file.path(fig.dir, "sdmTMB", paste0(TITLE,".txt")))
            print(paste0("M.idx=",M.idx))
            print(M.sanity)
            summary(M)
          sink()  
    
          M.df[M.idx,] <- cbind.data.frame(typ,
                                               BS,
                                               Max.Edg,
                                               n.re=(length(M$tmb_params$omega_s)), #Assuming a purely spatial spde model
                                               edf=sum(M.edf),
                                               cAIC=M.cAIC,
                                               n.flags=length(which(as.vector(unlist(M.sanity)) == "FALSE")),
                                               n.no0.segs=length(which(seg.dat.in$seg.ind > 0)),
                                               n.obs=nrow(dl2019.mrds.sf),
                                               n.ind=sum(seg.dat.in$seg.ind),
                                               p.avail,
                                               Nhat.plugin=as.numeric(Nhat[1]),
                                               Nhat.bc=as.numeric(Nhat[2]))
  
        #Evaluate other observation likelihoods
          
          M.dg <- try(update(M, family = delta_gamma()))
          M.dl <- try(update(M, family = delta_lognormal()))
          M.dpg <- try(update(M, family = delta_gamma(type="poisson-link")))
          M.dpl <- try(update(M, family = delta_lognormal(type="poisson-link")))
          M.nb2 <- try(update(M, family = nbinom2(link = "log")))
          M.pois <- try(update(M, family = poisson(link = "log")))
          
          M.list <- list(M.dg, M.dl, M.dpg, M.dpl, M.nb2, M.pois)
          M.typ <- c("dg", "dl", "dpg", "dpl", "nb2", "pois")
          
          for(i in 1:length(M.list)){
            
            M.idx.i <- M.idx + 1 
            
            TITLE <- paste(M.typ[i],BS,Max.Edg,"catsize",sep="_")
            Fnam <- file.path(fig.dir, "sdmTMB", TITLE)
            
            if(class(M.list[[i]])[1] == "try-error"){
                        
              sink(file=file.path(fig.dir, "sdmTMB", paste0(TITLE,"_TryError.txt")))
                print(M)
              sink()  
              
            } else {  
            
                        
              #Plot predictions, observations, and effort
                Nhat <- plot.sdmTMB.pred.obs.map.fun(m=M.list[[i]], 
                                        preddat=hex4pred.df, 
                                        grd=hex4pred, 
                                        obs=dl2019.mrds.sf, #observations as points
                                        counts=seg.dat.in, #seg.ind
                                        effort.sf=tx.10k4.in, 
                                        coast=NoAm.coast,
                                        pth=file.path(fig.dir, "sdmTMB"),
                                        title=TITLE,
                                        NSIM=1000)
                
              #Check simulation-based randomized quantile residuals using DHARMa  
        
                M.sim <- simulate(M.list[[i]], nsim=1000, type = "mle-mvn")
                
                png(paste(Fnam, "_DHARMa_sdmTMB.png",sep=""),bg="white",
                            height=640, width=960, units="px", pointsize=16)
                
                  dharma_residuals(M.sim, M.list[[i]], test_uniformity = TRUE)
                  
                dev.off()  

              #Output stuff
                  
                M.cAIC <- try(cAIC(M.list[[i]], what="cAIC"))
                if(class(M.cAIC)[1] == "try-error"){
                  M.cAIC <- NA
                }
                
                M.edf <- try(cAIC(M.list[[i]], what="EDF"))
                if(class(M.edf)[1] == "try-error"){
                  M.edf <- NA
                }
                
                M.sanity <- sanity(M.list[[i]])
                
                sink(file=file.path(fig.dir, "sdmTMB", paste0(TITLE,".txt")))
                  print(paste0("M.idx=",M.idx.i))
                  print(M.sanity)
                  summary(M.list[[i]])
                sink()  
          
                M.df[M.idx.i,] <- cbind.data.frame(M.typ[i],
                                                     BS,
                                                     Max.Edg,
                                                     n.re=(length(M.list[[i]]$tmb_params$omega_s)), #Assuming a purely spatial spde model
                                                     edf=sum(M.edf),
                                                     cAIC=M.cAIC,
                                                     n.flags=length(which(as.vector(unlist(M.sanity)) == "FALSE")),
                                                     n.no0.segs=length(which(seg.dat.in$seg.ind > 0)),
                                                     n.obs=nrow(dl2019.mrds.sf),
                                                     n.ind=sum(seg.dat.in$seg.ind),
                                                     p.avail,
                                                     Nhat.plugin=as.numeric(Nhat[1]),
                                                     Nhat.bc=as.numeric(Nhat[2]))
                
                M.idx <- M.idx.i
            }    
    
          }
      }

write.csv(M.df, file.path(fig.dir, "sdmTMB", "M_df.csv"))      
                       