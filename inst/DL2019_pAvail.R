#Script DL2019_pAvail.R...Megan C. Ferguson...21 August 2025

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
# 2. Availability bias estimator is from:
#    McLaren, I., 1961. Methods of determining the numbers and availability of 
#    ring seals in the eastern Canadian Arctic. Arctic 14, 162–175.
#
#    Frost, K., Lowry, L., 1995. Radio tag based correction factors for use in 
#    beluga whale population estimates. Working paper for Alaska Beluga Whale 
#    Committee Scientific Workshop, Anchorage, AK, 5-7 April 1995.
#
# 3. Estimate of time in view (tiv) for belugas for ASAMM surveys
#    on the Cmdr and Otter were computed as
#          tiv = div/speed
#    where div = distance in view. Because the estimated field of 
#    view for both the Cmdr and Otter is greater than the 
#    right-truncation distance for this analysis (Cmdr max.w.Dl=1.140511; Otter 
#    max.w.Dl=1.095209), use the max.w.Dl to inform div. --> div = 1.1 km 

  #Define p.avail.fun from McLaren (1961)
  # s: avg surface interval
  # d: avg dive interval
  # tiv: time in view
    Mc.p.avail <- function(s, d, tiv){
      p.a <- (s + tiv)/(s + d)
      return(p.a)
    }  
    
  #Define Frost and Lowry's (1995) availability probability for individual beluga
  # s.short: total duration of surfacings associated with short dives
  # s.long: total duration of surfacings associated with long dives
  # d.short: total duration of short dives
  # d.long: total duration of long dives
  # pA: availability probability from McLaren's estimator
    FL.p.avail <- function(s.short, s.long, d.short, d.long, pA){
      p.a <- ( (s.short + d.short) + (s.long + d.long)*pA )/(s.short + s.long + d.short + d.long)
      return(p.a)
    }
    
  #Compute tiv from div and speed
    cmdr.div <- 1140 #meters
    ott.div <- 1100 #meters
    
    speed115.kmh <- 212.98 #115 kts = 212.98 km/hr
    speed115.ms <- speed115.kmh*(1000/(60*60)) #115 kts = 59.1611 m/s
    
    cmdr.tiv <- cmdr.div/speed115.ms #seconds
    ott.tiv <- ott.div/speed115.ms #seconds
    
  #Define s.avg, d.avg, s.short, s.long, d.short, and d.long for belugas BB-585,  
  #CI-657, and CI-457 in Table 1 from Frost and Lowry(1995). All vectors present 
  #data for individual belugas in the following order: BB-585, CI-657, CI-457. 
  #These are data associated with the breakpoint of 10.3 sec.
    s.avg <- c(4, 2.66, 1.65) #sec
    d.avg <- c(37.01, 42.51, 40.08) #sec
    s.short <- c(181, 1553, 156) #sec 
    d.short <- c(195, 3783, 678) #sec
    s.long <- c(560, 1879, 262) #sec
    d.long <- c(5182, 30057, 6372) #sec
    
  #Calc McLaren's p.avail for each beluga using cmdr.tiv, ott.tiv, s.avg,
  #and d.avg
    Mc.pa.cmdr <- Mc.p.avail(s=s.avg,
                             d=d.avg,
                             tiv=cmdr.tiv)
    
    Mc.pa.ott <- Mc.p.avail(s=s.avg,
                             d=d.avg,
                             tiv=ott.tiv)
    #CK
      Mc.pa.cmdr
      Mc.pa.ott
      
  #Calc Frost and Lowry's p.avail for each beluga
    FL.pa.cmdr <- FL.p.avail(s.short, s.long, d.short, d.long, Mc.pa.cmdr)
    FL.pa.ott <- FL.p.avail(s.short, s.long, d.short, d.long, Mc.pa.ott)
    #CK
      FL.pa.cmdr
      FL.pa.ott
      
  #Calc Frost and Lowry's average p.avail
    FL.pa.cmdr.avg <- mean(FL.pa.cmdr)
    FL.pa.ott.avg <- mean(FL.pa.ott)
    #CK
      FL.pa.cmdr.avg #0.570
      FL.pa.ott.avg  #0.556
      
      1/FL.pa.cmdr.avg #1.754
      1/FL.pa.ott.avg  #1.799
    
    
