#Script DL2019_boot.R...Megan C. Ferguson...22 August 2025

# 1. This script is based on DL2019_boot_trois.R (v. 7.18.25)
#
# 2. a. cmdr detection function: no.ir.aba.bpc.hn.mrds.best.Alt.catsize
#    b. ott detection function: Dl.ott.hn.alt
#    c. sdmTMB model: tweedie_spde.bar_70_catsize (tweedie_spde.bar_70_catsize.Rdata)
#
# 3. Computing SE and CV on sdmTMB index values: This response from Sean 
#    Anderson was copied from
#    https://github.com/pbs-assess/sdmTMB/discussions/175#discussioncomment-4865196
#      "sdmTMB uses TMB's generalized delta method (see ?TMB::sdreport) on the 
#       area-weighted sum of abundance/biomass to calculate the SE. The upper and 
#       lower 95% CIs are then +/- 1.96 the SE from log abundance or biomass 
#       followed by exp(). Where 1.96 is from qnorm(1 - (1 - 0.95)/2).
#       Those are available in the output of get_index() in the se column. 
#       Those are on log abundance. Assuming a lognormal variable, which we are, 
#       you can get a CV as: sqrt(exp(se^2) - 1). For small values, the SE of the 
#       log variable and the CV will be very similar."
#
#  4. Availability bias correction factors were computed in DL2019_pAvail.R.
#
#  5. Statistical computations for total variance were based on 
#    NSDL_dsm_tmb_mgcv_smry_create_No_osc.R. 
#
      #Use law of total variance to compute total variance. In the case of DL2019,
      #THIS VAR DOES INCLUDE UNCERTAINTY FROM p0.1.
      #
      #The Law of Total Variance says:
      #
      # Var(Nhat) = E(Var(Nhat|p)) + Var(E(Nhat|p))
      #
      #where E(Var(Nhat|p)) = tmb.sd^2 and Var(E(Nhat|p)) is the 
      #variance in the bootstrap estimates of Nhat
      #
      #  VAR <- tmb.sd^2 + sd^2
      #  CV <- sqrt(VAR)/tmb.Nhat

  library(sdmTMB)
  library(sdmTMBextra)
  library(mgcv)
  library(Matrix)
  library(mrds)
  library(tidyverse)
  library(sf)
  library(ggbreak)

  #Set number of bootstrap iterations
    nboot <- 100
    #nboot <- 10000

  #Set random seed
    set.seed(20250327)

  ###############
  ###############
  # Set p.avail #
  ###############
  ###############
    #See script DL2019_pAvail.R for details
      p.avail <- 1/1.8 #should be 0.555....

  #Define file paths
    proj.dir <- "C:\\Users\\megan.ferguson\\OneDrive - Biodiversity Research Institute\\Documents\\Belugas\\2019_BeaufortSeaBelugas\\Analysis\\FergusonEtal_BeaufortSea_Beluga_Nhat_2019"
    out.dir <- file.path(proj.dir, "figures")
    
  #Define logistic function for p_1|2 
  #See eq. 6.48, p. 154 of Laake and Borchers 2004
    logistic.p.fun <- function(lpred){
      p.val <- exp(lpred)/(1+exp(lpred))
              return(p.val)
    }
    
  #Input stuff 
  # i. no.ir.aba.bpc.hn.mrds.best.Alt.catsize: mrds model for Cmdr
  # ii. Dl.ott.hn.alt: mcds model for Ott
  # iii. hex4pred.df: Prediction grid, as a data.frame. 
  # iv. seg.dat.in: Segments that are located within dl2019.buff.sf. seg.dat.in$a 
  #     is the strip searched, computed as seg.dat$seg.km[i1]*w1
  #     and seg.dat$seg.km[i2]*w2 for ddf1 and 2, respectively. These $a values
  #     remain constant.
  # v. dl2019.mrds.dat: Valid beluga observations
  # vi. sdmTMB model stuff in tweedie_spde.bar_70.Rdata
  # vii. bspde.70
    
    no.ir.aba.bpc.hn.mrds.best.Alt.catsize <- readRDS(file.path(proj.dir,"data","noirababpchnmrdsbestAltcatsize.rds"))
    Dl.ott.hn.alt <- readRDS(file.path(proj.dir,"data","Dlotthnalt.rds"))
    hex4pred.df <- readRDS(file.path(proj.dir,"data","hex4preddf.rds"))
    seg.dat.in <- readRDS(file.path(proj.dir,"data","segdatin.rds"))
    dl2019.mrds.dat <- readRDS(file.path(proj.dir,"data","dl2019_mrds_dat_Nhat.rds"))
    load(file.path(out.dir, "sdmTMB", "tweedie_spde.bar_70_catsize.Rdata"))
    bspde.70 <- readRDS(file.path(proj.dir, "data", "bspde70.rds"))
    
  #Identify which observations in dl2019.mrds.dat belong to ddf1 vs. ddf2
    i1 <- which(dl2019.mrds.dat$ddfobj == 1)
    i2 <- which(dl2019.mrds.dat$ddfobj == 2)
    
    n.obs1 <- length(i1)
    n.obs2 <- length(i2)
    
    obs1 <- dl2019.mrds.dat[i1,]
    obs2 <- dl2019.mrds.dat[i2,]
    #CK
      n.obs1
      n.obs2

  #####################    
  #Detection Functions#
  #####################    
    
    ddf1 <- no.ir.aba.bpc.hn.mrds.best.Alt.catsize
    ddf2 <- Dl.ott.hn.alt
      
    #Extract the MLEs and covariance matrix from fitted detection function model
      
      #Cmdr
        ddf1.beta <- ddf1$par
        ddf1.Sigma <- solve(ddf1$hessian)
        #CK
          summary(ddf1)
          
          ddf1.beta
            
          ddf1$hessian
          class(ddf1$hessian) #matrix
          dim(ddf1$hessian)   #5 x 5. Upper left 3x3 submatrix is for mr.
                              #       Lower right 2x2 submatrix is for ds.
                              #       The mr terms are independent of the ds
                              #       terms because the pairwise covariances are zero.            
          sqrt(diag(ddf1.Sigma)) #These values should equal the SEs for the
                                 #parameter estimates in summary(ddf).
          coefficients(ddf1)

      #Otter
        ddf2.beta <- ddf2$par
        ddf2.Sigma <- solve(ddf2$hessian)
        #CK
          summary(ddf2)
          
          ddf2.beta #These values should equal the parameter estimates in the 
                    #summary call above. First value is for shape coeff.
          
          sqrt(diag(ddf2.Sigma)) #These values should equal the SEs for the 
                                #parameter estimates in the summary(ddf) line above.
                                #First value is for shape coeff.
          coefficients(ddf2)

    #Implement ddf bootstrap
  
        DEBUGG <- FALSE
        #DEBUGG <- TRUE
        
        #Create matrices to hold bootstrapped values of average.p. Rows
        #correspond to sightings and columns correspond to bootstrap iterations.
          bs.p1 <- matrix(data=rep(0, nboot*n.obs1),
                               nrow=n.obs1)
          
          bs.p2 <- matrix(data=rep(0, nboot*n.obs2),
                               nrow=n.obs2)
          
          bs.Nht <- matrix(data=rep(NA, nboot*nrow(seg.dat.in)),
                           nrow=nrow(seg.dat.in))
          #CK
            n.obs1
            dim(bs.p1)
            
            n.obs2
            dim(bs.p2)
          
        for(ddf.boot in 1:nboot){
          
          #ddf.boot <- 1
          
          bs.obs1 <- obs1
          bs.obs2 <- obs2
  
          cat(paste0("ddf.boot = ",ddf.boot,"\n"))
          
          #Sample ddf betas, assuming a multivariate normal distribution
            ddf1.beta.i <- mgcv::rmvn(1, ddf1.beta, ddf1.Sigma)
            ddf2.beta.i <- mgcv::rmvn(1, ddf2.beta, ddf2.Sigma)
            #CK
              if(DEBUGG){
                print(ddf1.beta)
                print(ddf1.beta.i)
                print(ddf2.beta)
                print(ddf2.beta.i)
              }
            
          #Generate new bootstrap detection probabilities.
            
            #Compute bootstrap p(0) from ddf1.beta.i. This is p(0)_1|2 = p(0)_1  
            #under the assumption of point independence. See
            #p. 118 Laake and Borchers (2004) for definition of point independence.
              
              b0.i <- ddf1.beta.i[1]
              b1.i <- ddf1.beta.i[2]
              b2.i <- ddf1.beta.i[3]

              distance.y <- 0
                            
              #Compute linear predictor
                      
                #For catsize = 1
                  lpred.1.i <- b0.i + (b1.i*distance.y)
                          
                #For catsize = 2
                  lpred.2.i <- b0.i + (b1.i*distance.y) + b2.i
                  
              #Compute p(0). 
                p0.gs1.i <- logistic.p.fun(lpred.1.i)
                p0.gs2.i <- logistic.p.fun(lpred.2.i)
                #CK
                  p0.gs1.i
                  p0.gs2.i

          #Generate bootstrap predictions for the ddf1 and ddf2 ds models.
          #As a workaround, I apply the p0.gs... values computed above 
          #to the bootstrapped ds model estimates for both ddf1 and ddf2.
            
            bs.ddf1 <- ddf1
            bs.ddf1$par <- ddf1.beta.i
            bs.ddf1$mr$par <- ddf1.beta.i[1:3]
            bs.ddf1$ds$par <- ddf1.beta.i[4:5]
            bs.p1[,ddf.boot] <- predict(bs.ddf1$ds, newdata=bs.obs1, esw=FALSE)$fitted
            #CK
              test <- predict(ddf1$ds, newdata=obs1, esw=FALSE)$fitted
              summary(test - bs.p1[,ddf.boot]) #These should be different
              
              summary(test) 
              summary(bs.p1[,ddf.boot]) #Should be similar but not identical to test
            
            bs.ddf2 <- ddf2
            bs.ddf2$par <- ddf2.beta.i
            bs.p2[,ddf.boot] <- predict(bs.ddf2, newdata=bs.obs2, esw=FALSE)$fitted 
            #CK
              test <- predict(bs.ddf2, newdata=obs2, compute=TRUE, esw=FALSE)$fitted 
              summary(test - bs.p2[,ddf.boot]) #These should be identical
              
              test <- predict(ddf2, newdata=obs2, compute=TRUE, esw=FALSE)$fitted
              summary(test - bs.p2[,ddf.boot]) #These should be different
              
              summary(test) 
              summary(bs.p2[,ddf.boot]) #Should be similar but not identical to test

          #Multiply bs.p1[,ddf.boot] and bs.p2[,ddf.boot] by the p0.gs... values
          #computed above. 
            
            if(DEBUGG){
              print(summary(bs.p1[,ddf.boot]))
                
              print(summary(bs.p2[,ddf.boot]))
            }  
            
            #Cmdr
            
              idx <- which(bs.obs1$catsize == "1")
              bs.p1[idx,ddf.boot] <- p0.gs1.i*bs.p1[idx,ddf.boot]
        
              idx <- which(bs.obs1$catsize == "2")
              bs.p1[idx,ddf.boot] <- p0.gs2.i*bs.p1[idx,ddf.boot]
        
            #Ott

              idx <- which(bs.obs2$catsize == "1")
              bs.p2[idx,ddf.boot] <- p0.gs1.i*bs.p2[idx,ddf.boot]
        
              idx <- which(bs.obs2$catsize == "2")
              bs.p2[idx,ddf.boot] <- p0.gs2.i*bs.p2[idx,ddf.boot]
        
              #CK
                if(DEBUGG){
                  print(p0.gs1.i)
                  print(p0.gs2.i)
                  print(summary(bs.p1[,ddf.boot]))
                  print(summary(bs.p2[,ddf.boot]))
                }  

          #Calc HT estimator, Nht. 
              
            bs.obs1$Nht <- NA
            bs.obs2$Nht <- NA
            seg.Nht <- NA

            #If bootstrapped p(0) = 0, then keep Nht at NA and later remove this iteration 
            #from the analysis because it is outside the realm of reality.
            
            print(p0.gs1.i)
            print(p0.gs2.i)
            print(p0.gs1.i * p0.gs2.i)
            
            if((p0.gs1.i * p0.gs2.i > 0.0) &
               (!is.nan(p0.gs1.i * p0.gs2.i))){ #All p0 > 0.0 and none are NaN 

              bs.obs1$Nht <- bs.obs1$size/(bs.p1[,ddf.boot]*p.avail)
              bs.obs2$Nht <- bs.obs2$size/(bs.p2[,ddf.boot]*p.avail)
              #CK
                if(DEBUGG){
                  print(summary(bs.obs1$Nht))
                  print(summary(bs.obs2$Nht))
                }  
                  
              #Compute Nht per segment. 
                  
                all.obs <- rbind.data.frame(bs.obs1, bs.obs2)
                  
                seg.Nht <- sapply(1:nrow(seg.dat.in), function(i){
                            idx <- which(all.obs$Sample.Label == seg.dat.in$Sample.Label[i])
                            Nht.i <- sum(all.obs$Nht[idx])
                            return(Nht.i)
                })
                
                bs.Nht[,ddf.boot] <- seg.Nht
                
            }    
            #Ck
              if(DEBUGG){
                print(length(seg.Nht)) #823 if all p0 > 0.0
                print(summary(bs.Nht[,ddf.boot]))
                print(sum(bs.Nht[,ddf.boot])) #bootstrapped obs, pooled into segs
                print(sum(bs.obs1$Nht) + sum(bs.obs2$Nht)) #bootstrapped obs
                print(sum(dl2019.mrds.dat$Nht)) #original obs 2462.496
              }
              
        } #End ddf bootstrap            
  ############
  ############
  ############
            
  ##############      
  #sdmTMB model#
  ##############      
        
    M.name <- "tweedie_spde.bar_70_catsize"
            
    #sdmTMB requires data be class data.frame    
      seg.dat.in.df <- st_drop_geometry(seg.dat.in)  
      seg.dat.in.df$x <- st_coordinates(seg.dat.in)[,1]
      seg.dat.in.df$y <- st_coordinates(seg.dat.in)[,2]
      #CK
        summary(seg.dat.in.df)
        
    #Prep dataframe to save bootstrap Nhat and se values
      M.df <- data.frame(boot.idx = integer(),
                         M.name = character(),
                         n.re = integer(),
                         edf = double(),
                         n.flags = integer(),
                         Nhat = double(),
                         Nhat.bc.lo = double(),
                         Nhat.bc = double(),
                         Nhat.bc.hi = double(),
                         Nhat.bc.se = double(),
                         Nhat.bc.cv = double(),
                         stringsAsFactors=FALSE)
        
      #Initialize with fake values
        M.df[1,] <- cbind.data.frame(boot.idx=-1,
                                     M.name="ZZ",
                                     n.re=0,
                                     edf=-99.9,
                                     n.flags=-1,
                                     Nhat=-99.9,
                                     Nhat.bc.lo=-99.9,
                                     Nhat.bc=-99.9,
                                     Nhat.bc.hi=-99.9,
                                     Nhat.bc.se=-99.0,
                                     Nhat.bc.cv=-99.0)

    #Bootstrap new sdmTMB models, abundance estimates, and associated
    #uncertainty. Use bs.Nht instead of seg.Nht in seg.dat.in.df.

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
        
      set.seed(20250414)
        
      #Set the number of sdmTMB bootstrap iterations to the number of valid
      #ddf iterations
        
        sum.bs.Nht <- apply(bs.Nht, 2, sum)
        good.bs.Nht <- bs.Nht[,which(is.na(sum.bs.Nht) == FALSE)]  
        n.sdm.boot <- ncol(good.bs.Nht)
        #CK
          dim(bs.Nht) #ncol should equal nboot
          length(sum.bs.Nht) #should equal nboot
          nboot
          
          summary(sum.bs.Nht) #may have NA
          
          summary(apply(good.bs.Nht,2,sum)) #should have no NA
          n.sdm.boot #will be < nboot if sum.bs.Nhat had NAs

      for(sdm.boot in 1:n.sdm.boot){

        #sdm.boot <- 2 #For debugging

        print(paste0(" sdm.boot = ",sdm.boot))
        
        #Generate bootstrap data.frame using valid bootstrapped Nht in good.bs.Nht
          seg.dat.bs.df <- seg.dat.in.df
          seg.dat.bs.df$seg.Nht <- good.bs.Nht[,sdm.boot]
          
        #Try building bs sdmTMB model
          bs.M <- try(sdmTMB(
            seg.Nht ~ 1, #bootstrapped value
            data = seg.dat.bs.df, #use bootstrapped value of seg.Nht
            mesh = bspde.70,
            family = tweedie(),
            offset = log(seg.dat.in.df$a),
            spatial = "on",
          ))
          
        #Save stuff to M.df
          
          if(class(bs.M)[1] == "try-error"){
            M.df[sdm.boot,] <- cbind.data.frame(boot.idx=sdm.boot,
                                                M.name=M.name,
                                                n.re=0,
                                                edf=-99.9,
                                                n.flags=-1,
                                                Nhat=-99.9,
                                                Nhat.bc.lo=-99.9,
                                                Nhat.bc=-99.9,
                                                Nhat.bc.hi=-99.9,
                                                Nhat.bc.se=-99.0,
                                                Nhat.bc.cv=-99.0)
          } else {
            
            M.edf <- try(cAIC(bs.M, what="EDF"))
            if(class(M.edf)[1] == "try-error"){
              M.edf <- NA
            }
            
            M.sanity <- sanity(bs.M)
            
            #Plug-in estimator of predicted number of individuals across grid 
              pred.ind <- predict(bs.M, newdata=hex4pred.df, type="response",
                                      offset=log(hex4pred.df$a))$est 
              Nhat <- sum(pred.ind, na.rm=TRUE)
              #CK
                Nhat #area-integrated abundance w/o epsilon detransformation bias correction
            
            #Apply epsilon detransformation bias correction to Nhat estimate. Note
            #that get_index returns a data frame with a columns for time, estimate, 
            #lower and upper confidence intervals, log estimate, and standard error 
            #of the ******log estimate******.
            #
            #Computing SE and CV on sdmTMB index values: This response from Sean 
            #Anderson was copied from
            #    https://github.com/pbs-assess/sdmTMB/discussions/175#discussioncomment-4865196
            #      "sdmTMB uses TMB's generalized delta method (see ?TMB::sdreport) on the 
            #       area-weighted sum of abundance/biomass to calculate the SE. The upper and 
            #       lower 95% CIs are then +/- 1.96 the SE from log abundance or biomass 
            #       followed by exp(). Where 1.96 is from qnorm(1 - (1 - 0.95)/2).
            #       Those are available in the output of get_index() in the se column. 
            #       Those are on log abundance. Assuming a lognormal variable, which we are, 
            #       you can get a CV as: sqrt(exp(se^2) - 1). For small values, the SE of the 
            #       log variable and the CV will be very similar."        

              pred2 <- predict(bs.M, newdata=hex4pred.df, return_tmb_object=TRUE) 
              Nhat2.vec <- try(get_index(pred2, bias_correct=TRUE, area=hex4pred.df$a,
                                         level = 0.95))
              
              if(class(Nhat2.vec)[1] == "try-error"){
                Nhat2.vec$est <- NA
                Nhat2.vec$lwr <- NA
                Nhat2.vec$upr <- NA
                se <- NA
                cv <- NA
              } else {
                cv <- sqrt(exp(Nhat2.vec$se^2) - 1)
                se <- cv*Nhat2.vec$est
              }

              if(DEBUGG){
                print(Nhat2.vec) 
                print(se)
                print(cv)
              }
              
              M.df[sdm.boot,] <- cbind.data.frame(boot.idx=sdm.boot,
                                                  M.name=M.name,
                                                  n.re=(length(bs.M$tmb_params$omega_s)), #Assuming a purely spatial spde model
                                                  edf=sum(M.edf),
                                                  n.flags=length(which(as.vector(unlist(M.sanity)) == "FALSE")),
                                                  Nhat=Nhat,
                                                  Nhat.bc.lo=Nhat2.vec$lwr,
                                                  Nhat.bc=Nhat2.vec$est,
                                                  Nhat.bc.hi=Nhat2.vec$upr,
                                                  Nhat.bc.se=se,
                                                  Nhat.bc.cv=cv)
              rm(M.edf, M.sanity, Nhat, Nhat2.vec, se, cv)
          }
        
          rm(bs.M)
      } #end sdm.boot loop
          
    #Identify the number of valid ddf iterations
      
        sum.bs.Nht <- apply(bs.Nht, 2, sum)
        good.bs.Nht <- bs.Nht[,which(is.na(sum.bs.Nht) == FALSE)]  
        n.sdm.boot <- ncol(good.bs.Nht)
        #CK
          dim(bs.Nht) #823 x nboot; nrow = num hex cells; ncol = num ddf iterations.
          length(sum.bs.Nht) #should equal num ddf iterations
          summary(sum.bs.Nht) #may have NA
          
          summary(apply(good.bs.Nht,2,sum)) #should have no NA
          n.sdm.boot #will be < nboot if sum.bs.Nhat had NAs
      
    #Examine M.df
      dim(M.df)   #nrow = num bs.Nht iterations without NA; may have Inf
      summary(M.df)
      
      M.df$bad.sdmTMB <- 0
      M.df$bad.cAIC <- 0
      M.df$bad.bc <- 0
      M.df$bad.sanity <- 0
      
      #Identify sdmTMB BS iterations that generated "try-error" from function
      #sdmTMB. 
        M.df$bad.sdmTMB[which(M.df$n.re == 0)] <- 1

      #Identify sdmTMB BS iterations that generated "try-error" from function 
      #cAIC
        M.df$bad.cAIC[which(is.na(M.df$edf))] <- 1

      #Identify sdmTMB BS iterations that generated "try-error" from function
      #get_index
        M.df$bad.bc[which(is.na(M.df$Nhat.bc.lo))] <- 1

      #Identify sdmTMB BS iterations that generated flags during sanity check
        M.df$bad.sanity[which(M.df$n.flags > 0)] <- 1
        
      #Create a dataframe of "good" iterations
        M.df$test.idx <- M.df$bad.sdmTMB + M.df$bad.cAIC + M.df$bad.bc + M.df$bad.sanity
        M.df.good <- M.df[which(M.df$test.idx == 0),]
        #CK
          summary(M.df.good) 
          nrow(M.df.good)    
          
          summary(M.df$test.idx)
          test.tbl <- table(M.df$test.idx)
          test.tbl
          test.df <- data.frame(n.failed.tests=dimnames(test.tbl)[[1]],
                                n.iterations=as.vector(test.tbl))
          ggplot(test.df, aes(x=n.failed.tests, y=n.iterations)) + 
            geom_bar(stat = "identity", fill="darkturquoise", color="darkturquoise")
          ggsave(filename=file.path(out.dir, "Boot", "DL2019_dsm_boot_failed_test_barplot.png"),
                 dpi = "retina")
          
          bad.df <- M.df[,c(12:16)]
          apply(bad.df, 2, sum)
          
          test.tbl
          
          bad1 <- bad.df[which(bad.df$test.idx ==1),]
          apply(bad1, 2, sum)

          bad2 <- bad.df[which(bad.df$test.idx ==2),]
          apply(bad2, 2, sum) 
          
          bad3 <- bad.df[which(bad.df$test.idx ==3),]
          apply(bad3, 2, sum) 
          
      #Create histogram of resulting Nhat estimates
        sort(M.df.good$Nhat)[nrow(M.df.good):1]  
        sort(M.df.good$Nhat.bc)[nrow(M.df.good):1]
          
        d1 <- ggplot(M.df.good) +
          geom_density(fill="navy", color="navy", aes(Nhat.bc))
        d1
        ggsave(plot=d1, filename=file.path(out.dir, "Boot", "DL2019_dsm_boot_M_df_good_dens.png"),
               dpi = "retina")
        
        h1 <- ggplot(M.df.good) +
          geom_histogram(fill="navy", color="navy", aes(Nhat.bc), bins=100)
        h1
        ggsave(plot=h1, filename=file.path(out.dir, "Boot", "DL2019_dsm_boot_M_df_good_hist.png"),
               dpi = "retina")
        
      #Use law of total variance to compute total variance. The Law of Total Variance says:
      #
      # Var(Nhat) = E(Var(Nhat|p)) + Var(E(Nhat|p))
      #
      #where E(Var(Nhat|p)) = tmb.sd^2 and Var(E(Nhat|p)) is the 
      #variance in the bootstrap estimates of Nhat. 
      #
      
        #First, compute tmb.se based on output from sdmTMB's get_index function that
        #is saved to Nhat2.vec: This response from Sean Anderson was copied from
        #    https://github.com/pbs-assess/sdmTMB/discussions/175#discussioncomment-4865196
        #      "sdmTMB uses TMB's generalized delta method (see ?TMB::sdreport) on the 
        #       area-weighted sum of abundance/biomass to calculate the SE. The upper and 
        #       lower 95% CIs are then +/- 1.96 the SE from log abundance or biomass 
        #       followed by exp(). Where 1.96 is from qnorm(1 - (1 - 0.95)/2).
        #       Those are available in the output of get_index() in the se column. 
        #       Those are on log abundance. Assuming a lognormal variable, which we are, 
        #       you can get a CV as: sqrt(exp(se^2) - 1). For small values, the SE of the 
        #       log variable and the CV will be very similar."   
          
          tmb.cv <- sqrt(exp(Nhat2.vec$se^2) - 1)
          tmb.se <- tmb.cv*Nhat2.vec$est
          #CK
            Nhat2.vec$se
            tmb.cv
            tmb.se
      
          tot.var.good <- tmb.se^2 + var(M.df.good$Nhat.bc)
          tot.se.good <- sqrt(tot.var.good)
          tot.cv.good <- sqrt(tot.var.good)/Nhat2.vec$est
          #CK
            tot.var.good

            tot.se.good 

            tot.cv.good

  #Output stuff
    write.csv(M.df, file=file.path(out.dir, "boot", "DL2019_dsm_boot_M_df.csv"))
    save(bs.Nht, M.df, file=file.path(out.dir, "boot", "DL2019_dsm_boot_M_df.Rdata"))
    saveRDS(M.df.good, file=file.path(out.dir, "boot", "DL2019_dsm_boot_M_df_good.rds"))
    save.image(file=file.path(out.dir, "boot", "DL2019_dsm_boot_stats.Rdata"))
        
        sink(file=file.path(out.dir, "boot", "DL2019_dsm_boot_M_df_good.txt"))
          print(summary(M.df.good))
          print(paste0("num hex cells:", dim(bs.Nht)[1]))
          print(paste0("num ddf iterations:", dim(bs.Nht)[2]))
          print(paste0("num valid ddf iterations; may have Nhat=Inf:", nrow(M.df))) 
          print(paste0("num sdmTMB BS iterations that generated try-error from fcn sdmTMB:", sum(M.df$bad.sdmTMB))) 
          print(paste0("num sdmTMB BS iterations that generated try-error from function cAIC:", sum(M.df$bad.cAIC)))
          print(paste0("num sdmTMB BS iterations that generated try-error from fcn get_index:", sum(M.df$bad.bc)))   
          print(paste0("num sdmTMB BS iterations that generated flags during sanity check:", sum(M.df$bad.sanity))) 
          print(paste0("num iterations without p0=0, run-time errors, or flags:", nrow(M.df.good))) 
          print(paste0("tmb.se based on output from sdmTMB's get_index function (bias-corrected):", tmb.se))
          print(paste0("tmb.cv based on output from sdmTMB's get_index function (bias-corrected):", tmb.cv))
          print(paste0("standard deviation of M.df.good$Nhat.bc (with outliers):", sd(M.df.good$Nhat.bc)))
          print(paste0("variance of M.df.good$Nhat.bc (with outliers):", var(M.df.good$Nhat.bc)))
          print(paste0("CV of M.df.good$Nhat.bc (with outliers):", sd(M.df.good$Nhat.bc)/mean(M.df.good$Nhat.bc)))
          print(paste0("TOTAL variance of M.df.good$Nhat.bc:", tot.var.good))
          print(paste0("TOTAL SE of M.df.good$Nhat.bc:", tot.se.good))
          print(paste0("TOTAL CV of M.df.good$Nhat.bc:", tot.cv.good))
          print(paste0("Number of bootstrap iterations with failed tests:", paste(names(test.tbl), test.tbl, sep="=")))
          print(paste0("Breakdown of bootstrap iterations with only 1 failed test:", paste(names(bad1), apply(bad1, 2, sum), sep="=")))
        sink()
          
      
      
 