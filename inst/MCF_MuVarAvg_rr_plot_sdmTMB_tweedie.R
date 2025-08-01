#Script MCF_MuVarAvg_rr_plot_sdmTMB_tweedie.R...Megan C. Ferguson...10 Februrary 2025

# 1. This script is based on mcf_mod_eval_plots.R and NSDL17dsm_tmb_spde_tw.R.
#
# 2. This script takes an sdmTMB model with family = tweedie() and generates 
#    a plot of: 
#    i) mean squared response residuals vs. predicted values
#   ii) theoretical variance vs. mean for a Tweedie distribution with 
#           var = disp*mu^pwr
#       where
#           var = variance
#           mu = mean
#           disp = dispersion parameter, estimated by sdmTMB
#           pwr = power parameter, estimated by sdmTMB
#
# 3. Variables:
#       m: sdmTMB tweedie model object
#       dat: original data used to build m
#       OFFSET: offset term used in model formula for m
#       fnam: file path and the base of the filename for outputting figure

      sdmTMB.tweedie.mu.var.plot <- function(m, dat, OFFSET, fnam){

          #squared response residuals vs. mean
            
            #Compute predictions based on original data
              p <- predict(m, newdata=dat, type="response",
                                  offset=OFFSET)$est
            
            #Divide predictions into bins
              phist <- hist(p, breaks=20, plot=FALSE)
              pint <- findInterval(p, phist$breaks)
            
            #Sequence of means for plotting
              mu <- seq(from=range(p)[1], to=range(p)[2], length.out=100)
            
            #Compute squared residuals and average within bins
              #Response residuals
                modl.rr <- residuals(m, type="response")
                modl.rr2 <- (modl.rr)^2
                mean.modl.rr2 <- sapply(1:max(pint), function(i){
                  rr2.i <- modl.rr2[pint == i]
                  mean.rr2.i <- mean(rr2.i)
                  return(mean.rr2.i)
                })
                
            #Compute estimated variances based on mu or p and estimated dispersion parameter
                disp.idx <- which(names(m$sd_report$par.fixed) == "ln_phi")
                disp <- as.numeric(exp(m$sd_report$par.fixed[disp.idx]))
                disp
                
                finv.pwr.idx <- which(names(m$sd_report$par.fixed) == "thetaf")
                finv.pwr <- as.numeric(m$sd_report$par.fixed[finv.pwr.idx])
                pwr <- 1.0 + (exp(finv.pwr) / (1.0 + exp(finv.pwr)))
                pwr
                
                v <- disp*mu^pwr
                v.all <- disp*p^pwr
                  
            #Plot binned results
              png(paste0(fnam,"MuVarAvg_rr.png"),bg="white")
                plot(phist$mids, mean.modl.rr2, pch=19, xlab="mu",
                     ylab="Mean Squared Response Residuals")
                points(mu, v, pch="*", col="blue")
              dev.off()
              
      }        
