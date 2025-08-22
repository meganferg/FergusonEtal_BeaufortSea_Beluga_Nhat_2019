#Script DL2019_AvailBias.r...Megan C. Ferguson...12 July 2025

# NOTES
#
# 1. From Citta et al. 2013. Dive behavior of Eastern Chukchi 
#    beluga whales (Delphinapterus leucas), 1998-2008:
#    *[p. 402] Dive Duration
#      In general, belugas dove 5.1 - 9.8 times per hour. Fewer 
#      dives occurred during periods when belugas were diving to 
#      deeper depths (Table 5). The average duration of dives was 
#      very consistent across dive types (i.e., shallow, 
#      intermediate, and deep), although the longest dives (> 21 
#      min) were observed during periods with the deepest dives.
#
#      Average dive duration ranged from 3.4 to 5.1 min (Table 5). 
#      After accounting for region, we found no significant relationship
#      between sex or age and average or maximum dive durations. 
#      While there were some significant interactions between 
#      region and age, the effects on dive duration were slight. 
#      For example, although average duration had a significant 
#      interaction with region*age, it varied by less than a minute
#      within regions (Fig. 8a). In general, average dive duration 
#      was fairly consistent across regions, sexes, and ages (Fig. 
#      8a), while the maximum dive duration increased as depth 
#      increased (Fig. 8b).
#
#      Our findings generally agree with other studies of beluga 
#      dive behavior. Kingsley et al. (2001) found that most dives 
#      were shorter than 10 min in duration and that dive duration 
#      did not vary greatly with dive depth in Hudson Bay, where 
#      depth ranged from 40 to 400 m. Near Somerset Island, Canada,
#      Martin and Smith (1999) observed a longer mean dive time 
#      (13.1 min) than we observed, but a similar maximum dive 
#      duration (22.9 min). For beluga whales near Devon Island, 
#      Canada, Heide-J?rgenson et al. (1998) found that most dives
#      ranged from 9 to 18 min, which is also similar to the
#      durations we observed (Table 5).
#
# 2. Availability bias estimator is from Laake et al. 1997. Journal of 
#    Wildlife Management 61(1):63-75.
#
# 3. See also p_avail_fun.r for more examples using p.avail.fun().
#
# 4. Estimate of time in view (tiv) for belugas for ASAMM surveys
#    on the Cmdr and Otter were computed as
#          tiv = div/speed
#    where div = distance in view. Because the estimated field of 
#    view for both the Cmdr and Otter (see FB_FOV_Analysis_PartTrois.r 
#    and Inuvik_FOV_Analysis.r, respectively) is greater than the 
#    right-truncation distance for this analysis (Cmdr max.w.Dl=1.140511, 
#    from Dl2019_Cmdr_ddf.r; Otter max.w.Dl=1.095209, from Dl2019_Otter_ddf.r),
#    use the max.w.Dl to inform div. --> div = 1.1 km 
#
# 5. DFO estimates of avg s and avg d are from Table 3 in:
#
#    Marcoux, Marianne, Luke Storrie, Shannon MacPhee, and Lisa Loseto. 2025. 
#    “Availability Bias Adjustment for Calculating Aerial Survey Abundance Estimates 
#    for Belugas (Delphinapterus Leucas) in the Eastern Beaufort Sea.” 2025/002. 
#    Canadian Science Advisory Secretariat Research Document. Winnipeg, Manitoba, 
#    Canada: Fisheries and Oceans Canada.
#
#    Marcoux et al. recommended using the 5-m values for the offshore area.


  #Define p.avail.fun. (See p_avail_fun.r)
  # s: avg surface interval
  # d: avg dive interval
  # tiv: time in view
    p.avail.fun <- function(s, d, tiv){
      p.a <- s/(s+d) + d*( 1 - exp(-(tiv/d)) )/(s + d)
      return(p.a)
    }  
    
  #Compute tiv from div and speed
    div <- 1100 #meters
    
    speed115.kmh <- 212.98 #115 kts = 212.98 km/hr
    speed115.ms <- speed115.kmh*(1000/(60*60)) #115 kts = 59.1611 m/s
    
    tiv <- div/speed115.ms #seconds
    
  #Estimate ballpark surface interval based on Citta et al. (2013)
    avg.d <- c(3.4, 5.1) #average dive duration in mins
    num.d <- c(5.1, 9.8) #average number of dives per hour
    
    #Function to compute average surface interval (minutes) from average dive interval (d, in minutes)
    #and average number of dives per hour (num)
      avg.s.fun <- function(d, num){
        s <- (60/num) - d
        return(s)
      }
      
      #Set up d and num to get all combinations
        d.vec <- c(avg.d, avg.d[2:1])
        n.vec <- rep(num.d, 2)
        
      #Compute avg.s for all combinations of d and num
        avg.s <- avg.s.fun(d.vec, n.vec)
        
      #Convert avg dive and surface intervals to secs
        avg.d.secs <- d.vec*60
        avg.s.secs <- avg.s*60
        #CK
          d.vec
          n.vec
          avg.s #ranges from 1.0 to 8.4 mins
          
          avg.d.secs
          avg.s.secs
      
    #Compute probability of availability (pA) for different surface and dive intervals
      
      pA <- p.avail.fun(avg.s.secs, avg.d.secs, tiv)
      pA.df <- cbind.data.frame(avg.s.secs, 
                                "avg.s.mins"=avg.s.secs/60,
                                avg.d.secs, 
                                "avg.d.mins"=avg.d.secs/60,
                                "n.dive.per.hr"=n.vec,
                                pA)
      #CK
        pA
        pA.df
        p.avail.fun(avg.s.secs[1], avg.d.secs[1], tiv) #should equal pA[1]
        p.avail.fun(avg.s.secs[2], avg.d.secs[2], tiv) #should equal pA[2]
        p.avail.fun(avg.s.secs[3], avg.d.secs[3], tiv) #should equal pA[3]
        p.avail.fun(avg.s.secs[4], avg.d.secs[4], tiv) #should equal pA[4]  
        
    #Output pA and pA.df
#      save(pA, pA.df, file="Output//DL2019_AvailBias.Rdata")
#      write.csv(pA.df, file="Output//DL2019_AvailBias.csv", row.names=FALSE)
      
####
####
####

  #Compute pA using DFO August data from Marcoux et al. (2025) Table 3.
  #Marcoux et al. recommended using the 5-m values for the offshore area.
    
    dfo.jul.z <- c(1,2,5) #meters
    dfo.jul.s.min <- c(4.01, 4.59, 5.38) #avg. suface time (min)
    dfo.jul.d.min <- c(5.47, 5.63, 6.14) #avg. dive time (min)
    
    dfo.jul.s.sec <- dfo.jul.s.min*60
    dfo.jul.d.sec <- dfo.jul.d.min*60
    
    test.tiv <- 13.87 #DFO's estimate for their surveys (pg. 4)
    
    tiv.cmdr.Dl <- 1140/speed115.ms #seconds. width searched from DL2019_DetectionFcns_S4.pdf
    tiv.ott.Dl <- 1100/speed115.ms  #seconds. width searched from DL2019_DetectionFcns_S4.pdf
    #CK
      tiv.cmdr.Dl
      tiv.ott.Dl
    
    pA.dfo.jul.cmdr <- p.avail.fun(dfo.jul.s.sec, dfo.jul.d.sec, tiv.cmdr.Dl)
    pA.dfo.jul.cmdr.df <- cbind.data.frame(dfo.jul.s.sec, 
                                dfo.jul.s.min,
                                dfo.jul.d.sec, 
                                dfo.jul.d.min,
                                "tiv.cmdr.sec"=rep(tiv.cmdr.Dl,3),
                                pA.dfo.jul.cmdr)
    
    pA.dfo.jul.ott <- p.avail.fun(dfo.jul.s.sec, dfo.jul.d.sec, tiv.ott.Dl)
    pA.dfo.jul.ott.df <- cbind.data.frame(dfo.jul.s.sec, 
                                dfo.jul.s.min,
                                dfo.jul.d.sec, 
                                dfo.jul.d.min,
                                "tiv.ott.sec"=rep(tiv.ott.Dl,3),
                                pA.dfo.jul.ott)
    
    #CK
      pA.dfo.jul.cmdr
      1/pA.dfo.jul.cmdr
      pA.dfo.jul.cmdr.df
      
      pA.dfo.jul.ott
      1/pA.dfo.jul.ott
      pA.dfo.jul.ott.df
      
      p.avail.fun(dfo.jul.s.sec, dfo.jul.d.sec, test.tiv)
      
  #Compute pA using DFO August data from Marcoux et al. (2025) Table 3.
  #Marcoux et al. recommended using the 5-m values for the offshore area.
    
    dfo.aug.z <- c(1,2,5) #meters
    dfo.aug.s.min <- c(4.52, 5.40, 6.37) #avg. suface time (min)
    dfo.aug.d.min <- c(7.64, 7.88, 8.40) #avg. dive time (min)
    
    dfo.aug.s.sec <- dfo.aug.s.min*60
    dfo.aug.d.sec <- dfo.aug.d.min*60
    
    tiv.cmdr.Dl <- 1140/speed115.ms #seconds. width searched from DL2019_DetectionFcns_S4.pdf
    tiv.ott.Dl <- 1100/speed115.ms  #seconds. width searched from DL2019_DetectionFcns_S4.pdf
    #CK
      tiv.cmdr.Dl
      tiv.ott.Dl
    
    pA.dfo.aug.cmdr <- p.avail.fun(dfo.aug.s.sec, dfo.aug.d.sec, tiv.cmdr.Dl)
    pA.dfo.aug.cmdr.df <- cbind.data.frame(dfo.aug.s.sec, 
                                dfo.aug.s.min,
                                dfo.aug.d.sec, 
                                dfo.aug.d.min,
                                "tiv.cmdr.sec"=rep(tiv.cmdr.Dl,3),
                                pA.dfo.aug.cmdr)
    
    pA.dfo.aug.ott <- p.avail.fun(dfo.aug.s.sec, dfo.aug.d.sec, tiv.ott.Dl)
    pA.dfo.aug.ott.df <- cbind.data.frame(dfo.aug.s.sec, 
                                dfo.aug.s.min,
                                dfo.aug.d.sec, 
                                dfo.aug.d.min,
                                "tiv.ott.sec"=rep(tiv.ott.Dl,3),
                                pA.dfo.aug.ott)
    
    #CK
      pA.dfo.aug.cmdr
      1/pA.dfo.aug.cmdr
      pA.dfo.aug.cmdr.df
      
      pA.dfo.aug.ott
      1/pA.dfo.aug.ott
      pA.dfo.aug.ott.df
      
      p.avail.fun(dfo.aug.s.sec, dfo.aug.d.sec, test.tiv)
    
    
    
    
    
    
    
    
    
    

