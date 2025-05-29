# CUSTOM FUNCTIONS

# Effect size (lnVR and lnCVR) for 2 main effects and interaction effect ----
effect_setV <- function(CC_n, CC_mean, CC_SD,
                       EC_n, EC_mean, EC_SD,
                       CS_n, CS_mean, CS_SD,
                       ES_n, ES_mean, ES_SD,
                       percent){
  
  if(percent == "No"){
  # lnRR----
  # main effect Environmental enrichment----
  lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                         log(0.5*(CS_mean+ CC_mean))
  
  lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
    (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
  
  # main effect Stress----
  lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                         log(0.5*(EC_mean+ CC_mean))
  
  lnRRV_S <- lnRRV_E
  
  # interaction----
  
  lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                            (log(EC_mean) - log(CC_mean))
  
  
  lnRRV_ES <- 
    (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
     ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
      ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
       ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
  
  # SMD
  SD_pool <- sqrt(((ES_n-1)*ES_SD^2 + 
                                (EC_n-1)*EC_SD^2 + 
                                (CS_n-1)*CS_SD^2 +
                                (CC_n-1)*CC_SD^2) / 
                               (ES_n + EC_n + CS_n + CC_n - 4))
  
  
  
  # lnVR
  # main effect Environmental enrichment----
  lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
  # main effect Stress----
  lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnVRV_S <- lnVRV_E
  
  # interaction----
  
  lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  # lnCVR
  # main effect Environmental enrichment----
  ES_CV <- ES_SD/ES_mean
  EC_CV <- EC_SD/EC_mean
  CS_CV <- CS_SD/CS_mean
  CC_CV <- CS_SD/CS_mean
  
  lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnCVRV_E <- lnRRV_E + lnVRV_E 
  
  # main effect Stress----
  lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnCVRV_S <- lnRRV_S + lnVRV_S
  
  # interaction----
  lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
  effect <- tibble(
    # lnVR
    lnVR_E = lnVR_E,
    lnVRV_E = lnVRV_E, 
    lnVR_S = lnVR_S, 
    lnVRV_S = lnVRV_S,
    lnVR_ES =lnVR_ES, 
    lnVRV_ES = lnVRV_ES,
    # lnCVR
    lnCVR_E = lnCVR_E,
    lnCVRV_E = lnCVRV_E, 
    lnCVR_S = lnCVR_S, 
    lnCVRV_S = lnCVRV_S,
    lnCVR_ES =lnCVR_ES, 
    lnCVRV_ES = lnCVRV_ES
  )
  effect
  }
  
  else {
    
    asin_trans <- function(percent) { asin(sqrt(percent/100)) }
    
    
    # transforming SD 
    ES_SD <- sqrt((ES_SD/100)^2/(4*(ES_mean/100)*(1-(ES_mean/100))))
    EC_SD <- sqrt((EC_SD/100)^2/(4*(EC_mean/100)*(1-(EC_mean/100))))
    CS_SD <- sqrt((CS_SD/100)^2/(4*(CS_mean/100)*(1-(CS_mean/100))))
    CC_SD <- sqrt((CC_SD/100)^2/(4*(CC_mean/100)*(1-(CC_mean/100))))
    
    # transformaing mean
    ES_mean <- asin_trans(ES_mean)
    EC_mean <- asin_trans(EC_mean)
    CS_mean <- asin_trans(CS_mean)
    CC_mean <- asin_trans(CC_mean)
    
    # lnRR
    # main effect Enrichment
    lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                           log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) +  
                             (1/(CS_mean + CC_mean))^2*(CS_SD^2 /CS_n + CC_SD^2 / CC_n) 
    
    # main effect Stress
    lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                           log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_S <- lnRRV_E
    
    # interaction----
    
    lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                              (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ES <- (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
                    ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
                    ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
                    ((CC_SD)^2 / ((CC_mean)^2*CC_n)))    
     
     
    # lnVR
    # main effect Environmental enrichment----              
    lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
    # main effect Stress----
    lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnVRV_S <- lnVRV_E
  
    # interaction----
  
    lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    # lnCVR
    # main effect Environmental enrichment----
    ES_CV = ES_SD/ES_mean
    EC_CV = EC_SD/EC_mean
    CS_CV = CS_SD/CS_mean
    CC_CV = CS_SD/CS_mean
  
    lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnCVRV_E <- lnRRV_E + lnVRV_E 
  
    # main effect Stress----
    lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnCVRV_S <- lnRRV_S + lnVRV_S
  
    # interaction----
  
    lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
    effect <- tibble(
    # lnVR
    lnVR_E = lnVR_E,
    lnVRV_E = lnVRV_E, 
    lnVR_S = lnVR_S, 
    lnVRV_S = lnVRV_S,
    lnVR_ES =lnVR_ES, 
    lnVRV_ES = lnVRV_ES,
    # lnCVR
    lnCVR_E = lnCVR_E,
    lnCVRV_E = lnCVRV_E, 
    lnCVR_S = lnCVR_S, 
    lnCVRV_S = lnCVRV_S,
    lnCVR_ES = lnCVR_ES, 
    lnCVRV_ES = lnCVRV_ES
  )
    effect
  }
  
}



# Effect size (lnRR and SMD) for 2 main effects and interaction effect ----
effect_set <- function(CC_n, CC_mean, CC_SD,
                       EC_n, EC_mean, EC_SD,
                       CS_n, CS_mean, CS_SD,
                       ES_n, ES_mean, ES_SD,
                       percent){
  
  if(percent == "No"){
  
  # lnRR----
  # main effect Environmental enrichment----
  lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                         log(0.5*(CS_mean+ CC_mean))
  
  lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
    (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
  
  # main effect Stress----
  lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                         log(0.5*(EC_mean+ CC_mean))
  
  lnRRV_S <- lnRRV_E
  
  # interaction----
  
  lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                            (log(EC_mean) - log(CC_mean))
  
  
  lnRRV_ES <- 
    (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
     ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
      ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
       ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
  
  # SMD
  SD_pool <- sqrt(((ES_n-1)*ES_SD^2 + 
                                (EC_n-1)*EC_SD^2 + 
                                (CS_n-1)*CS_SD^2 +
                                (CC_n-1)*CC_SD^2) / 
                               (ES_n + EC_n + CS_n + CC_n - 4))
  
  
  # main effect Environment enrichment
  SMD_E <- ((ES_mean + EC_mean) - (CS_mean + CC_mean))/ (2*SD_pool)
  
  
  SMDV_E <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + 
                    (SMD_E^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
  
  
  
  # main effect Stress
  SMD_S <- ((ES_mean + CS_mean) - (EC_mean + CC_mean)) / (2*SD_pool)
  
  SMDV_S <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + 
                    (SMD_S^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
  
  # interaction
  SMD_ES <- ((ES_mean - EC_mean) - (CS_mean - CC_mean)) / SD_pool
  
  SMDV_ES <- (1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_ES^2 / (2*(ES_n + EC_n + CS_n + CC_n)))
  
  
  # lnVR
  # main effect Environmental enrichment----
  lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
  # main effect Stress----
  lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnVRV_S <- lnVRV_E
  
  # interaction----
  
  lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  # lnCVR
  # main effect Environmental enrichment----
  ES_CV <- ES_SD/ES_mean
  EC_CV <- EC_SD/EC_mean
  CS_CV <- CS_SD/CS_mean
  CC_CV <- CS_SD/CS_mean
  
  lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnCVRV_E <- lnRRV_E + lnVRV_E 
  
  # main effect Stress----
  lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnCVRV_S <- lnRRV_S + lnVRV_S
  
  # interaction----
  lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
  effect <- tibble(
    # lnRR
    lnRR_E = lnRR_E,
    lnRRV_E = lnRRV_E, 
    lnRR_S = lnRR_S, 
    lnRRV_S = lnRRV_S,
    lnRR_ES =lnRR_ES, 
    lnRRV_ES = lnRRV_ES,
    #SMD
    SMD_E = SMD_E,
    SMDV_E = SMDV_E, 
    SMD_S = SMD_S, 
    SMDV_S = SMDV_S, 
    SMD_ES = SMD_ES, 
    SMDV_ES = SMDV_ES
    # lnVR
    #lnVR_E = lnVR_E,
    #lnVRV_E = lnVRV_E, 
    #lnVR_S = lnVR_S, 
    #lnVRV_S = lnVRV_S,
    #lnVR_ES =lnVR_ES, 
    #lnVRV_ES = lnVRV_ES,
    # lnCVR
    #lnCVR_E = lnCVR_E,
    #lnCVRV_E = lnCVRV_E, 
    #lnCVR_S = lnCVR_S, 
    #lnCVRV_S = lnCVRV_S,
    #lnCVR_ES =lnCVR_ES, 
    #lnCVRV_ES = lnCVRV_ES
  )
  effect
  }
  
  else {
    
    asin_trans <- function(percent) { asin(sqrt(percent/100)) }
    
    
    # transforming SD 
    ES_SD <- sqrt((ES_SD/100)^2/(4*(ES_mean/100)*(1-(ES_mean/100))))
    EC_SD <- sqrt((EC_SD/100)^2/(4*(EC_mean/100)*(1-(EC_mean/100))))
    CS_SD <- sqrt((CS_SD/100)^2/(4*(CS_mean/100)*(1-(CS_mean/100))))
    CC_SD <- sqrt((CC_SD/100)^2/(4*(CC_mean/100)*(1-(CC_mean/100))))
    
    # transformaing mean
    ES_mean <- asin_trans(ES_mean)
    EC_mean <- asin_trans(EC_mean)
    CS_mean <- asin_trans(CS_mean)
    CC_mean <- asin_trans(CC_mean)
     
    # lnRR
    # main effect Enrichment
    lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                           log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) +  
                             (1/(CS_mean + CC_mean))^2*(CS_SD^2 /CS_n + CC_SD^2 / CC_n) 
    
    # main effect Stress
    lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                           log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_S <- lnRRV_E
    
    # interaction----
    
    lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                              (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ES <- (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
                    ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
                    ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
                    ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
    
    # SMD
    SD_pool <- sqrt(((ES_n-1)*ES_SD^2 + 
                    (EC_n-1)*EC_SD^2 + 
                    (CS_n-1)*CS_SD^2 +
                    (CC_n-1)*CC_SD^2) / 
                    (ES_n + EC_n + CS_n + CC_n - 4))
    
    
    # main effect Environment enrichment
    SMD_E <- ((ES_mean + EC_mean) - (CS_mean + CC_mean))/ (2*SD_pool)
    
    
    SMDV_E <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_E^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
    
    # main effect Stress
    SMD_S <- ((ES_mean + CS_mean) - (EC_mean + CC_mean)) / (2*SD_pool)
    
    SMDV_S <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_S^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
    
    # interaction
    SMD_ES <- ((ES_mean - EC_mean) - (CS_mean - CC_mean)) / SD_pool
    
    SMDV_ES <- (1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_ES^2 / (2*(ES_n + EC_n + CS_n + CC_n)))
    
    # lnVR
    # main effect Environmental enrichment----              
    lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
    # main effect Stress----
    lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnVRV_S <- lnVRV_E
  
    # interaction----
  
    lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    # lnCVR
    # main effect Environmental enrichment----
    ES_CV = ES_SD/ES_mean
    EC_CV = EC_SD/EC_mean
    CS_CV = CS_SD/CS_mean
    CC_CV = CS_SD/CS_mean
  
    lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnCVRV_E <- lnRRV_E + lnVRV_E 
  
    # main effect Stress----
    lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnCVRV_S <- lnRRV_S + lnVRV_S
  
    # interaction----
  
    lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
    effect <- tibble(
    # lnRR
    lnRR_E = lnRR_E,
    lnRRV_E = lnRRV_E, 
    lnRR_S = lnRR_S, 
    lnRRV_S = lnRRV_S,
    lnRR_ES =lnRR_ES, 
    lnRRV_ES = lnRRV_ES,
    #SMD
    SMD_E = SMD_E,
    SMDV_E = SMDV_E, 
    SMD_S = SMD_S, 
    SMDV_S = SMDV_S, 
    SMD_ES = SMD_ES, 
    SMDV_ES = SMDV_ES
    # lnVR
    #lnVR_E = lnVR_E,
    #lnVRV_E = lnVRV_E, 
    #lnVR_S = lnVR_S, 
    #lnVRV_S = lnVRV_S,
    #lnVR_ES =lnVR_ES, 
    #lnVRV_ES = lnVRV_ES,
    # lnCVR
    #lnCVR_E = lnCVR_E,
    #lnCVRV_E = lnCVRV_E, 
    #lnCVR_S = lnCVR_S, 
    #lnCVRV_S = lnCVRV_S,
    #lnCVR_ES = lnCVR_ES, 
    #lnCVRV_ES = lnCVRV_ES
  )
    effect
  }
  
}


# Removing asin_trans for sensitivity analysis----

effect_setb <- function(CC_n, CC_mean, CC_SD,
                        EC_n, EC_mean, EC_SD,
                        CS_n, CS_mean, CS_SD,
                        ES_n, ES_mean, ES_SD)
  {
    # lnRR----
    # main effect Environmental enrichment----
    lnRR_Eb <- log(0.5*(ES_mean + EC_mean)) - 
      log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_Eb <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
      (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
    
    # main effect Stress----
    lnRR_Sb <- log(0.5*(ES_mean + CS_mean)) - 
      log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_Sb <- lnRRV_Eb
    
    # interaction----
    
    lnRR_ESb <-   (log(ES_mean) - log(CS_mean)) - 
      (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ESb <- 
      (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
         ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
         ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
         ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
         
         
  # lnVR
  # main effect Environmental enrichment----
                         
  lnVR_Eb <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnVRV_Eb <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
  # main effect Stress----
  lnVR_Sb <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnVRV_Sb <- lnVRV_Eb
  
  # interaction----
  
  lnVR_ESb <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnVRV_ESb <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  # lnCVR
  # main effect Environmental enrichment----
  ES_CV = ES_SD/ES_mean
  EC_CV = EC_SD/EC_mean
  CS_CV = CS_SD/CS_mean
  CC_CV = CS_SD/CS_mean
  
  lnCVR_Eb <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnCVRV_Eb <- lnRRV_Eb + lnVRV_Eb 
  
  # main effect Stress----
  lnCVR_Sb <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnCVRV_Sb <- lnRRV_Sb + lnVRV_Sb
  
  # interaction----
  
  lnCVR_ESb <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnCVRV_ESb <- lnRRV_ESb + lnVRV_ESb
  
        
    effectb <- tibble(
      # lnRR
      lnRR_Eb = lnRR_Eb,
      lnRRV_Eb = lnRRV_Eb, 
      lnRR_Sb = lnRR_Sb, 
      lnRRV_Sb = lnRRV_Sb,
      lnRR_ESb =lnRR_ESb, 
      lnRRV_ESb = lnRRV_ESb,
      # lnVR
      lnVR_Eb = lnVR_Eb,
      lnVRV_Eb = lnVRV_Eb, 
      lnVR_Sb = lnVR_Sb, 
      lnVRV_Sb = lnVRV_Sb,
      lnVR_ESb =lnVR_ESb, 
      lnVRV_ESb = lnVRV_ESb,
      # lnCVR
      lnCVR_Eb = lnCVR_Eb,
      lnCVRV_Eb = lnCVRV_Eb, 
      lnCVR_Sb = lnCVR_Sb, 
      lnCVRV_Sb = lnCVRV_Sb,
      lnCVR_ESb =lnCVR_ESb, 
      lnCVRV_ESb = lnCVRV_ESb
    )
    effectb
}


# Pairwise comparisons lnRR (not for SMD) -----

effect_set2 <- function(CC_n, CC_mean, CC_SD,
                        EC_n, EC_mean, EC_SD,
                        CS_n, CS_mean, CS_SD,
                        ES_n, ES_mean, ES_SD,
                        percent){
  
  if(percent == "No"){
  
  # EE vs control
  lnRR_E2 <- log(EC_mean) - log(CC_mean)
  
  
  lnRRV_E2 <-  (EC_SD^2 / (EC_mean^2*EC_n)) + 
                            (CC_SD^2 / (CC_mean^2*CC_n))
  
  
  # Stress vs control
  lnRR_S2 <- log(CS_mean) - log(CC_mean)
  
  lnRRV_S2 <- (CS_SD^2 / (CS_mean^2*CS_n)) + 
                           (CC_SD^2 / (CC_mean^2*CC_n))
  
  # EE + stress vs control
  lnRR_ES2 <- log(ES_mean) - log(CC_mean)
  
  lnRRV_ES2 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                            (CC_SD^2 / (CC_mean^2*CC_n))
  
  # EE + stress vs stress (the effect of E in the presence of S)
  lnRR_E3 <- log(ES_mean) - log(CS_mean)
  
  lnRRV_E3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                           (CS_SD^2 / (CS_mean^2*CS_n))
  
  # EE + stress vs EE (the effect of S in the presence of E)
  lnRR_S3 <- log(ES_mean) - log(EC_mean)
  
  lnRRV_S3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                           (EC_SD^2 / (EC_mean^2*EC_n))
  
  effect2 <- tibble(
    lnRR_E2 = lnRR_E2,
    lnRRV_E2 = lnRRV_E2, 
    lnRR_S2 = lnRR_S2, 
    lnRRV_S2 = lnRRV_S2, 
    lnRR_ES2 =lnRR_ES2, 
    lnRRV_ES2 = lnRRV_ES2,
    lnRR_E3 =lnRR_E3, 
    lnRRV_E3 = lnRRV_E3,
    lnRR_S3 = lnRR_S3,
    lnRRV_S3 = lnRRV_S3
  )
  effect2
  }
  else {
    asin_trans <- function(percent) { asin(sqrt(percent/100)) }
    
    # transforming SD 
    ES_SD <- sqrt((ES_SD/100)^2/(4*(ES_mean/100)*(1-(ES_mean/100))))
    EC_SD <- sqrt((EC_SD/100)^2/(4*(EC_mean/100)*(1-(EC_mean/100))))
    CS_SD <- sqrt((CS_SD/100)^2/(4*(CS_mean/100)*(1-(CS_mean/100))))
    CC_SD <- sqrt((CC_SD/100)^2/(4*(CC_mean/100)*(1-(CC_mean/100))))
    
    # transformaing mean
    ES_mean <- asin_trans(ES_mean)
    EC_mean <- asin_trans(EC_mean)
    CS_mean <- asin_trans(CS_mean)
    CC_mean <- asin_trans(CC_mean)
    
    # EE vs control
    lnRR_E2 <- log(EC_mean) - log(CC_mean)
    
    
    lnRRV_E2 <- (EC_SD^2 / (EC_mean^2*EC_n)) + 
                              (CC_SD^2 / (CC_mean^2*CC_n))
    
    # Stress vs control
    lnRR_S2 <- log(CS_mean) - log(CC_mean)
    
    lnRRV_S2 <- (CS_SD^2 / (CS_mean^2*CS_n)) + 
                             (CC_SD^2 / (CC_mean^2*CC_n))
    
    # EE + stress vs control
    lnRR_ES2 <- log(ES_mean) - log(CC_mean)
    
    lnRRV_ES2 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                              (CC_SD^2 / (CC_mean^2*CC_n))
    
    # EE + stress vs stress (the effect of E in the presence of S)
    lnRR_E3 <-log(ES_mean) - log(CS_mean)
    
    lnRRV_E3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                             (CS_SD^2 / (CS_mean^2*CS_n))
    
    # EE + stress vs EE (the effect of S in the presence of E)
    lnRR_S3 <- log(ES_mean) - log(EC_mean)
    
    lnRRV_S3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                             (EC_SD^2 / (EC_mean^2*EC_n))
    
    
    effect2 <- tibble(
      lnRR_E2 = lnRR_E2,
      lnRRV_E2 = lnRRV_E2, 
      lnRR_S2 = lnRR_S2, 
      lnRRV_S2 = lnRRV_S2, 
      lnRR_ES2 =lnRR_ES2, 
      lnRRV_ES2 = lnRRV_ES2,
      lnRR_E3 =lnRR_E3, 
      lnRRV_E3 = lnRRV_E3,
      lnRR_S3 = lnRR_S3,
      lnRRV_S3 = lnRRV_S3
    )
    effect2
  }
  
}

# making take from metafor models
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}

# pairwise comparisons
### Step 1 ####
# contrast function
# takes: 1) data frame, 2) moderator of interest, and 3) VCV matrix

contrast_fun <- function(data, response, moderator, VCV){
  
  mod <- deparse(substitute(moderator))
  resp <- deparse(substitute(response))
  
  # get names of different levels of a categorical variable 
  names <- levels(as.factor(data[[mod]]))
  moderator <- as.factor(data[[mod]])
  
  # function for running a meta-regression
  run_rma <- function(name) {
    rma.mv(yi = data[[resp]], 
           V = VCV, 
           mods = ~ relevel(moderator, ref = name),
           random = list(~1|Study_ID, ~1|ES_ID,~1|Strain),
           test = "t",
           data = data,
           sparse = TRUE)
  }
  
  
  # running all the models to get all the contrasts
  contra <- map(names[-length(names)], run_rma)
}

### Step 2 ####
# function to get relevant estimates from normal and contast models
# takes: 1) normal model output 2) contrast model outputs, and 3) moderator of interest, 
get_estimate <- function(model, contra, moderator) {
  
  # making into the name
  mod <- deparse(substitute(moderator))
  
  # getting estimations 
  get_pred1 <- function (model, mod =  mod) {
    name <- name <- firstup(as.character(stringr::str_replace(row.names(model$beta), mod, "")))
    len <- length(name)
    newdata <- matrix(NA, ncol = len, nrow = len)
    for (i in 1:len) {
      pos <- which(model$X[, i] == 1)[[1]]
      newdata[, i] <- model$X[pos, ]
    }
    pred <- metafor::predict.rma(model, newmods = newdata)
    
    estimate <- pred$pred
    lowerCL <- pred$ci.lb
    upperCL <- pred$ci.ub 
    lowerPR <- pred$cr.lb
    upperPR <- pred$cr.ub 
    
    table <- tibble(name = factor(name, levels = name, labels = name), estimate = estimate,
                    lowerCL = lowerCL, upperCL = upperCL,
                    pval = model$pval,
                    lowerPR = lowerPR, upperPR = upperPR)
    table
  }
  # something has changed here
  # get estiamtions from normal model
  get_pred2 <- function (model) {
    name <- as.factor(str_replace(row.names(model$beta), 
                                  paste0("relevel", "\\(", "moderator, ref = name", "\\)"),
                                  ""))
    len <- length(name)
    newdata <- diag(len)
    pred1 <- data.frame(predict.rma(model, newmods = newdata[1,-1]))
    pred2 <- data.frame(predict.rma(model, intercept = FALSE, newmods = newdata[-1,-1]))
    pred <- rbind(pred1, pred2)
    
    estimate <- pred$pred
    lowerCL <- pred$ci.lb
    upperCL <- pred$ci.ub 
    lowerPR <- pred$cr.lb
    upperPR <- pred$cr.ub 
    
    table <- tibble(name = factor(name, levels = name, labels = name), 
                    estimate = estimate,
                    lowerCL = lowerCL, 
                    upperCL = upperCL,
                    pval = model$pval,
                    lowerPR = lowerPR, upperPR = upperPR)
    table
  }
  
  # getting contrasts name
  cont_gen <- function(name) {
    combination <- combn(name, 2)
    name_dat <- t(combination)
    names <- paste(name_dat[, 1], name_dat[, 2], sep = "-")
    return(names)
  }
  
  # put into a proper tale
  mr_results <- function(res1, res2) {
    restuls <-tibble(
      Levels = c(as.character(res1$name), cont_gen(res1$name)),
      Estimate = c(res1$estimate, res2$estimate),
      Lower_CI = c(res1$lowerCL, res2$lowerCL),
      Upper_CI = c(res1$upperCL, res2$upperCL),
      P_value = c(res1$pval, res2$pval),
      Lower_PI = c(res1$lowerPR, res2$lowerPR),
      Upper_PI = c(res1$upperPR, res2$upperPR),
    )
    restuls
  }
  
  res1 <- get_pred1(model, mod = mod)
  
  # applying the estimation function for each model output
  estiamtes <- map(contra, ~ get_pred2(.x))
  
  # a list of the numbers to take out unnecessary contrasts
  contra_list <- Map(seq, from=1, to =1:length(contra))
  
  # you need to flatten twice: first to make it a list and make it a vector
  # I guess we can make it a loop (or vectorise)
  res2 <- map2_dfr(estiamtes, contra_list, ~.x[-(.y), ])   
  
  # putting two results together
  
  results <-  mr_results(res1, res2)
}


# Functions for plotting
############################################################################
orchard_plot <- function(object, mod = "1", group, data, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         trunk.size = 3, branch.size = 1.2, twig.size = 0.5, errorbar = 0.03, precision.size = 3.5, circle.col = "#999999",
                         transfm = c("none", "tanh", "exp"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  mod_results(object, mod, group, data, N,
                              by = by, at = at, weights = weights, upper = upper)
    } else {
      results <-  mod_results(object, mod = "1", group, data, N,
                              by = by, at = at, weights = weights, upper = upper)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  data_trim <- results$data
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision"
  
  if(any(N != "none")){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
    #latex2exp::TeX()
  }
  
  if(transfm == "tanh"){
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- Zr_to_r(mod_table[,cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }else{
    label <- xlab
  }
  
  
  if(transfm == "exp"){
    lnRR_to_perc <- function(df){
      EXP <- function(d) {1-exp(d)}
      return(sapply(df, EXP))
    }
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- lnRR_to_perc(mod_table[,cols])
    data_trim$yi <- lnRR_to_perc(data_trim$yi)
    label <- xlab
  }else{
    label <- xlab
  }
  
  # Add in total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))
  
  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])
  
  # the number of groups in a moderator & data points
  group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  # colour blind friendly colours with grey
  cbpl <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  # setting fruit colour
  if(colour == TRUE){
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  }else{
    color <- data_trim$mod
    color2 <- mod_table$name
  }
  
  # whether we fill fruit or not
  if(fill == TRUE){
    fill <- color
  }else{
    fill <- NULL
  }
  
  # whether marginal
  if(names(mod_table)[2] == "condition"){
    
    # the number of levels in the condition
    condition_no <- length(unique(mod_table[, "condition"]))
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  
                                                              shape = as.factor(condition), 
                                                              #shape = 5,
                                                              fill = color2), 
                               size = twig.size, position = ggplot2::position_dodge2(width = 0.3), fatten = trunk.size) +
      # this will only work for up to 5 different conditions
      # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
      ggplot2::scale_shape_manual(values =  20 + (1:condition_no)) + ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::labs(shape = condition.lab) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle))
    
  } else {
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21, col="#999999") + 
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.3) +
      # creating CI # geom_linerange
      ggplot2::geom_errorbar(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size, width = errorbar) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name,  ymin = lowerPR, ymax = upperPR, fill = color2), size = twig.size, fatten = trunk.size, shape = 23) +
      # change the point estimate to white color
      ggplot2::geom_point(data = mod_table,ggplot2::aes(y = estimate, x = name),shape = 23, fill = "white") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle)) + #+
      #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
      ggplot2::theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line.x = element_line(colour = "black"),
                     axis.ticks.y = element_blank())
    
  }
  
  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting colors in
  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values = cbpl) +
      ggplot2::scale_colour_manual(values = cbpl)
  }
  
  # putting k and g in
  if(k == TRUE && g == FALSE && k.pos == "right"){
    plot <- plot +
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = precision.size)
  } else if(k == TRUE && g == FALSE && k.pos == "left") {
    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = precision.size)
  } else if (k == TRUE && g == TRUE && k.pos == "right"){
    # get group numbers for moderator
    plot <- plot + #ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
      #label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
      #parse = TRUE, hjust = "right", size = 3.5)  
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(N)[observation]==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = precision.size) +
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.2),
                        label= paste("italic(N)['rodent strain']==", mod_table$g[1:group_no]), parse = TRUE, hjust = "right", size = precision.size)
    
  } else if (k == TRUE && g == TRUE && k.pos == "left"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "left", size = 3.5)
  }
  
  
  return(plot)
}



i2_ml <- function(model, method = c("ratio", "matrix"), data, boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment i2_ml cannot take models with heterogeneous variance.")}
  
  ## evaluate choices
  method <- match.arg(method)
  
  if (method == "matrix") {
    # Wolfgang Viechtbauer's method
    I2s <- matrix_i2(model)
  } else {
    # Nakagawa & Santos (2012)
    I2s <- ratio_i2(model)
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot)
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Paramatric bootsrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    
    I2_each <- sapply(sim, function(ysim) {
      
      # The model
      tmp <- metafor::rma.mv( ysim, vi,
                              mods = mods_formula,
                              random = random_formula,
                              data = data)
      pb$tick()
      Sys.sleep(1 / boot)
      
      if(method == "matrix"){
        I2 <- matrix_i2(tmp)
      } else {
        I2 <- ratio_i2(tmp)
      }
      
      return(I2) })
    
    # Summarise the bootstrapped distribution.
    I2s_each_95 <- data.frame(t(apply(I2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    I2s <-  round(I2s_each_95, digits = 3)
    colnames(I2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(I2s)
}

#' @title matrix_i2
#' @description Calculated I2 (I-squared) for mulilevel meta-analytic models, based on a matrix method proposed by Wolfgang Viechtbauer.
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence. The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim <- i2_ml(lim_MR, data=lim, method = "matrix")
#' }
#' @export
matrix_i2 <- function(model){
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  W <- solve(model$V)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total <- 100* (sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
  I2_each <- 100* (model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"
  I2s <- c(I2_total, I2_each)
  return(I2s)
}


#' @title ratio_i2
#' @description I2 (I-squared) for mulilevel meta-analytic models based on Nakagawa & Santos (2012). Under multilevel models, we can have a multiple I2 (see also Senior et al. 2016).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt,
#' sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 1000)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' }
#' @export
ratio_i2 <- function(model){
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / model$vi) * (model$k - 1) /
    (sum(1 / model$vi)^2 - sum((1 / model$vi)^2))
  
  # s^2_t = total variance
  I2_total <- 100 * (sum(model$sigma2) / (sum(model$sigma2) + sigma2_v))
  I2_each <- 100 * (model$sigma2 / (sum(model$sigma2) + sigma2_v))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"
  
  I2s <- c(I2_total, I2_each)
  return(I2s)
}

r2_ml <- function(model, data, boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment r2_ml cannot take models with heterogeneous variance.")}
  
  R2 <- R2_calc(model)
  
  if(!is.null(boot)){
    
    if(any(class(model) %in% c("robust.rma")) == TRUE){stop("Sorry, bootstrapping currently doesn't work for robust.rma objects. Please use rma.mv instead.")}
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot)
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    # Paramatric bootsrap
    R2 <- sapply(sim, function(ysim) {
      # The model
      tmp <- metafor::rma.mv( ysim, vi,
                              mods = mods_formula,
                              random = random_formula,
                              data = data)
      R2s <- R2_calc(tmp)
      pb$tick()
      Sys.sleep(1 / boot)
      return(R2s)
    })
    
    # Summarise the bootstrapped distribution.
    R2 <- data.frame(t(apply(R2, 1, stats::quantile, probs=c(0.5, .025, .975))))
    R2 <-  round(R2, digits = 3)
    colnames(R2) = c("Est.", "2.5%", "97.5%")
  }
  
  return(R2)
  
}

#' @title R2_calc
#' @description Calculated R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- R2_calc(lim_MR)
#' }
#' @export

R2_calc <- function(model){
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  # fixed effect variance
  fix <- stats::var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  
  # conditional. Need to remove 'residual' variance; assume this is the sigma level with the largest k. Really the only way we can get that to work.
  R2c <- (fix + sum(model$sigma2) - model$sigma2[which(model$s.nlevels.f == max(model$s.nlevels.f))]) /
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}

bubble_plot <- function(object, mod, group = NULL, data,
                        xlab = "Moderator", ylab = "Effect size", N = "none",
                        alpha = 0.5, cb = TRUE, k = TRUE, g = FALSE,
                        est.lwd = 1, ci.lwd = 0.5, pi.lwd = 0.5, precision.size = 2.5,
                        est.col = "black", ci.col = "black", pi.col = "black", bubble.col = "grey90",
                        legend.pos = c("top.left", "top.right",
                                       "bottom.right", "bottom.left",
                                       "top.out", "bottom.out",
                                       "none"),
                        k.pos = c("top.right", "top.left",
                                  "bottom.right", "bottom.left",
                                  "none"),
                        condition.nrow = 2,
                        #condition.lab = "Condition",
                        weights = "prop", by = NULL, at = NULL)
{
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  #facet <- match.arg(NULL, choices = facet)
  
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?bubble_plot")
  }
  
  if(any(grepl(mod, colnames(data))) == FALSE){
    error("The moderator specified is not found in your data. Did you transform the variable in the model and forget to add it to your dataframe?")
  }
  
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?bubble_plot")
  }
  
  if(is.numeric(by)){
    k = FALSE
    g = FALSE
  }
  
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  mod_results(object, mod, group, data,
                              by = by, at = at, weights = weights)
    } else {
      results <-  mod_results(object, mod = "1", group, data,
                              by = by, at = at, weights = weights)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  data_trim <- results$data
  #data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision"
  
  if(any(N != "none")){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
  }
  
  label <- xlab
  # if(transfm == "tanh"){
  #   cols <- sapply(mod_table, is.numeric)
  #   mod_table[,cols] <- Zr_to_r(mod_table[,cols])
  #   data_trim$yi <- Zr_to_r(data_trim$yi)
  #   label <- xlab
  # }else{
  #   label <- xlab
  # }
  
  
  if(is.null(data_trim$condition) == TRUE){
    # the number of effect sizes
    effect_num <- nrow(data_trim)
    
    # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
    group_num <- length(unique(data_trim$stdy))
    
    dat_text <- data.frame(K = effect_num, G = group_num)
    
  }else{
    effect_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(x[,"yi"])))
    
    # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
    #group_num <- c(2,4)
    group_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(base::unique(x[,"stdy"]))))
    
    dat_text <- data.frame(K = effect_num, G = group_num, condition = as.vector(base::unique(data_trim$condition)))
  }
  # the number of groups in a moderator & data points
  #group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  # # colour blind friendly colours with grey
  # cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  

    plot <- ggplot2::ggplot() +
      # putting bubbles
      #geom_point(data = data, aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() # +
    #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
    # theme(legend.direction="horizontal") +
    # #theme(legend.background = element_rect(fill = "white", colour = "black")) +
    # theme(legend.background = element_blank()) +
    # theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))

  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting k and g in
  # c("top.right", "top.left", "bottom.right", "bottom.left","none")
  if(k == TRUE && g == FALSE && k.pos == "top.right"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = 2,
                         vjust   = 2.5
      )
    
  } else if(k == TRUE && g == FALSE && k.pos == "top.left") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2.5
      )
  } else if(k == TRUE && g == FALSE && k.pos == "bottom.right") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = 2,
                         vjust   = -1.5
      )
  } else if (k == TRUE && g == FALSE && k.pos == "bottom.left"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -1.5
      )
    # below get g ----
    
  } else if (k == TRUE && g == TRUE && k.pos == "top.right"){
    # get group numbers for moderator
    plot <- plot +
      #ggplot2::geom_text(data = dat_text,
      #                   mapping = ggplot2::aes(x = Inf, y = Inf),
      #                   label =  paste("italic(k)==",
      #                                  rev(dat_text$K),
      #                                  "~","(", rev(dat_text$G), ")"),
      #                   parse = TRUE,
      #                   hjust   = 1.5,
      #                   vjust   = 2) +
      
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = Inf),
                         label =  paste("italic(N[observation])==",
                                        rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = 1.5,
                         vjust   = 2,
                         size = precision.size) +
      
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = Inf),
                         label =  paste("italic(N['rodent strain'])==",
                                        rev(dat_text$G)),
                         parse = TRUE,
                         hjust   = 1.35,
                         vjust   = 3,
                         size = precision.size)
    
    
    
    
    
  } else if (k == TRUE && g == TRUE && k.pos == "top.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==",
                                        rev(dat_text$K),
                                        "~","(", rev(dat_text$G), ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==",
                                        rev(dat_text$K),
                                        "~","(", rev(dat_text$G), ")"),
                         parse = TRUE,
                         hjust   = 1.5,
                         vjust   = -0.5)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==",
                                        rev(dat_text$K),
                                        "~","(", rev(dat_text$G), ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -0.5)
  }
  
  # # putting colors in
  # if(cb == TRUE){
  #   plot <- plot +
  #     ggplot2::scale_fill_manual(values=cbpl) +
  #     ggplot2::scale_colour_manual(values=cbpl)
  # }
  
  return(plot)
}
mod_results <- function(model, mod = "1", group, data, N = NULL,  weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...){
  
  if(any(grepl("-1|0", as.character(model$formula.mods)))){
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }
  
  if(missing(model)){
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")}
  
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  
  if(is.null(stats::formula(model))){
    #model <- stats::update(model, "~1")
    model$formula.mods <- ~ 1
    dat_tmp <- data$`1` <- "Intrcpt"
    model$data <- dat_tmp
  } else{
    model$data <- data
  }
  
  if(model$test == "t"){
    df_mod = as.numeric(model$ddf[[1]])
  } else{
    df_mod = 1.0e6 # almost identical to z value
  }
  
  if(is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1) ## NOTE: Added data argument emmeans >vers 1.7.4. Object is unstable so feeding in the relevant arguments from model object directly. Note, we should think about df!
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    
    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    
    if(is.null(by)){
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
      
    } else{
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              condition = mm_pi[,2], estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }
    
    # Extract data
    data2 <- get_data_raw(model, mod, group, N, data, at = at, subset)
    
    mod_table$name <- factor(mod_table$name,
                             levels = mod_table$name,
                             labels = mod_table$name)
    
  } else{
    at2 <- list(mod = seq(min(data[,mod], na.rm = TRUE), max(data[,mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula =  stats::formula(model), data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1, at = c(at2, at))  # getting 100 points. Fixing this to make it more general
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)
    
    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if(is.null(by)){
      mod_table <- data.frame(moderator = mm_pi[,1],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    } else{
      mod_table <- data.frame(moderator = mm_pi[,1],
                              condition = mm_pi[,2],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }
    
    # extract data
    data2 <- get_data_raw_cont(model, mod, group, N, data, by = by)
    
  }
  
  
  output <- list(mod_table = mod_table,
                 data = data2)
  
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}




############# Key Sub-functions #############

#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (credibility intervals) from \code{esmeans} objects (\pkg{metafor}).
#' @param model \code{rma.mv} object.
#' @param mm result from \code{emmeans::emmeans} object.
#' @param mod Moderator of interest.
#' @param ... other arguments passed to function.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export


pred_interval_esmeans <- function(model, mm, mod, ...){
  
  tmp <- summary(mm)
  tmp <- tmp[ , ]
  test.stat <- stats::qt(0.975, tmp$df[[1]])
  
  if(length(model$tau2) <= 1){ # including gamma2
    sigmas <- sum(model$sigma2)
    PI <- test.stat * base::sqrt(tmp$SE^2 + sigmas)
  } else {
    sigmas <- sum(model$sigma2)
    taus   <- model$tau2
    gammas <- model$gamma2
    w_tau <- model$g.levels.k
    w_gamma <- model$g.levels.k
    
    if(mod == "1"){
      tau <- weighted_var(taus, weights = w_tau)
      gamma <- weighted_var(gamma, weights = w_gamma)
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau + gamma)
      
    } else {
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus + gammas)
    }
  }
  
  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI
  
  # renaming "overall" to ""
  if(tmp[1,1] == "overall"){tmp[,1] <- "intrcpt"}
  
  return(tmp)
}

#' @title get_data_raw
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor}.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or whatever other grouping variable one wishes to present sample sizes.
#' @param data The data frame used to fit the \code{rma.mv} model object.
#' @param N The name of the column in the data specifying the sample size, N. Defaults to \code{NULL}, so precision is plotted instead of sample size.
#' @param at List of moderators. If \code{at} is equal to \code{mod} then levels specified within \code{at} will be used to subset levels when \code{subset = TRUE}. Otherwise, it will marginalise over the moderators at the specified levels.
#' @param subset Whether or not to subset levels within the \code{mod} argument. Defaults to \code{FALSE}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' @examples \dontrun{
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'  test <- get_data_raw(model, mod = "trait.type", group = "group_ID", data = warm_dat, at = list(trait.type = c("physiology", "morphology")))
#'  test2 <- get_data_raw(model, mod = "1", group = "group_ID", data = warm_dat)
#'
#'  data(english)
#'  # We need to calculate the effect sizes, in this case d
#'  english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"), data = english)
#'  model <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#'  test3 <-  get_data_raw(model, mod = "1", group = "StudyNo", data = english)}

get_data_raw <- function(model, mod, group, N = NULL, data, at = NULL, subset = TRUE){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  # Extract data
  # Check first if missing data exists
  if(length(attr(model$X, "dimnames")[[1]]) > 0){
    # full model delete missing values so need to adjust
    position <- attr(model$X, "dimnames")[[1]]
    data <- data[position, ] }
  if(!is.null(at) & subset){
    # Find the at slot in list that pertains to the moderator and extract levels
    at_mod <- at[[mod]]
    position2 <- which(data[,mod] %in% at_mod)
    # Subset the data to only the levels in the moderator
    data <- data[position2,]
    yi <- model$yi[position2]
    vi <- model$vi[position2]
    type <- attr(model$yi, "measure")
  } else {
    # Extract effect sizes
    yi <- model$yi
    vi <- model$vi
    type <- attr(model$yi, "measure")
  }
  if(mod == "1"){
    moderator <- "Intrcpt"
  }else{
    # Get moderator
    moderator <- as.character(data[[mod]]) # Could default to base instead of tidy
    moderator <- firstup(moderator)
  }
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  #names(data_reorg)[4] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  
  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }
  
  return(data_reorg)
}

#' @title get_data_raw_cont
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor} when a continuous variable is fit within a model object.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param N  The name of the column in the data specifying the sample size, N. Defaults to \code{NULL} so that precision is plotted instead of sample size.
#' @param data The data frame used to fit the \code{rma.mv} model object.
#' @param by Character name(s) of the 'condition' variables to use for grouping into separate tables.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export


get_data_raw_cont <- function(model, mod, group, N = NULL, data, by){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  # Extract data
  # Check first if missing data exists
  if(length(attr(model$X, "dimnames")[[1]]) > 0){
    # full model delete missing values so need to adjust
    position <- attr(model$X, "dimnames")[[1]]
    data <- data[position, ] }
  # Extract effect sizes
  yi <- model$yi
  vi <- model$vi
  type <- attr(model$yi, "measure")
  # Get moderator
  moderator <- data[[mod]] # Could default to base instead of tidy
  #names(moderator) <  "moderator"
  if(is.null(by)){
    condition <- data[ , by]
  }else{
    condition <- data[[by]]
  }
  #names(condition) <  "condition"
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, condition, stdy, type)
  # if(!is.na(names(data_reorg)[names(data_reorg) == by]) == TRUE) {  ## FAILING HERE
  #   names(data_reorg)[names(data_reorg) == by] <- "condition"
  # }
  #names(data_reorg)[5] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  
  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }
  
  return(data_reorg)
}

############# Helper-functions #############

#' @title firstup
#' @description Uppercase moderator names
#' @param x a character string
#' @param upper logical indicating if the first letter of the character string should be capitalized. Defaults to TRUE.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a character string with all combinations of the moderator level names with upper case first letters
#' @export

firstup <- function(x, upper = TRUE) {
  if(upper){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  } else{ x }
}


#' @title print.orchard
#' @description Print method for class 'orchard'
#' @param x an R object of class orchard
#' @param ... Other arguments passed to print
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#'

print.orchard <- function(x, ...){
  return(print(x$mod_table))
}

#' @title weighted_var
#' @description Calculate weighted variance
#' @param x A vector of tau2s to be averaged
#' @param weights Weights, or sample sizes, used to average the variance
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a vector with a single weighted variance
#' @export
#'

weighted_var <- function(x, weights){
  weight_var <- sum(x * weights) / sum(weights)
  return(weight_var)
}


#' @title num_studies
#' @description Computes how many studies are in each level of categorical moderators of a \code{rma.mv} model object.
#' @param mod Character string describing the moderator of interest.
#' @param data Raw data from object of class "orchard"
#' @param group A character string specifying the column name of the study ID grouping variable.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a table with the number of studies in each level of all parameters within a \code{rma.mv} or \code{rma} object.
#' @export
#' @examples
#' \dontrun{data(fish)
#'warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list( ~1 | es_ID,~1 | group_ID), mods = ~experimental_design-1, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' num_studies(model$data, experimental_design, group_ID)
#' }

num_studies <- function(data, mod, group){
  
  # Summarize the number of studies within each level of moderator
  table <- data        %>%
    dplyr::group_by({{mod}}) %>%
    dplyr::summarise(stdy = length(unique({{group}})))
  
  table <- table[!is.na(table$moderator),]
  # Rename, and return
  colnames(table) <- c("Parameter", "Num_Studies")
  return(data.frame(table))
  
}


#' @title submerge
#' @description Merge two model results tables (orchard objects).
#' @param object1  object of class \code{orchard}.
#' @param object2  object of class \code{orchard}.
#' @param ... Other arguments passed to submerge.
#' @param mix If \code{TRUE}, it will add the number to the moderator name.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame.
#' @export
#'
submerge <- function(object1, object2, ..., mix = FALSE){
  orchard_list <- list(object1, object2, ...)
  
  len <- length(orchard_list)
  # merging tables
  tables <- lapply(orchard_list, function(x) x$mod_table)
  tables <- do.call("rbind", tables)
  
  # merging data
  ## checking moderator names are the same or not
  datas <- lapply(orchard_list, function(x) x$data)
  datas <- do.call("rbind", datas)
  
  # renaming
  if(mix == TRUE){
    names <- lapply(orchard_list, function(x) x$data$moderator)
    names <- as.vector(unlist(mapply(function(x, y) paste0(x, y), x = names, y = 1:len)))
    datas$moderator <- factor(names)
    tables$name <- levels(factor(names))
  }
  
  model_results <- list(mod_table = tables, data = datas)
  
  class(model_results) <- "orchard"
  
  return(model_results)
  
}
############################################################################
calc_evi <- function(m1, m2, sd1, sd2, n1, n2, rscale = 1) {
  # pooled standard deviation
  spooled <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Cohen's d
  d <- (m1 - m2) / spooled
  
  # Hedge's g
  J <- 1 - (3 / (4 * (n1 + n2 - 2) - 1))  # Correction factor
  g <- J * d
  
  # p-value
  sepooled <- sqrt((spooled^2 / n1) + (spooled^2 / n2))
  t_value <- (m1 - m2) / sepooled
  p_value <- 2 * pt(abs(t_value), df = n1 + n2 - 2, lower.tail = FALSE)
  
  # Bayes factor
  bayes_factor <- as.data.frame(meta.ttestBF(t = t_value, n1 = n1, n2 = n2, rscale = rscale))$bf
  
  # post-hoc power
  power <- pwr.t2n.test(n1 = n1, n2 = n2, d = d)$power
  
  # return a data frame with the results
  return(data.frame(
    d = d,
    g = g,
    p = p_value,
    bf = bayes_factor,
    pwr = power
  ))
}


# folded distribution
folddist <-  function(mu, sd){
  postfnorm <- stats::dnorm(mu, 0, sd)*2*sd^2 + mu*(2*stats::pnorm(mu, 0, sd) -1)
    est <- postfnorm
    var.fnorm <- mu^2 + sd^2 - (sd*sqrt(2/pi)*exp((-1*mu^2)/(2*sd^2)) + mu*(1-2*stats::pnorm(-1*mu/sd, 0, 1)))^2
    est <- data.frame(Mean=est, Variance = var.fnorm)
  return(est)
}

# colour blind friendly colours with grey
cbpl <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


#' @title cvh1_ml
#' @description CVH1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH1. TODO - we need to cite original CVH1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH1_eng_1 <- cvh1_ml(english_MA, boot = 10)
#' CVH1_eng_2 <- cvh1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH1_fish_1 <- cvh1_ml(model, boot = 10)
#' CVH1_fish_2 <- cvh1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH1_lim_1 <- cvh1_ml(lim_MR, boot = 10)
#' CVH1_lim_2 <- cvh1_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh1_ml <- function(model,
                    boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh1_ml cannot take models with heterogeneous variance.")}
  
  CVH1s <- ml_cvh1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CVH1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH1 <- ml_cvh1(tmp)
      })
    } else{
      CVH1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH1 <- ml_cvh1(tmp)
        return(CVH1) })
    }
    # Summarise the bootstrapped distribution.
    CVH1s_each_95 <- data.frame(t(apply(CVH1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH1s <-  round(CVH1s_each_95, digits = 3)
    colnames(CVH1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVH1s)
}


#' @title ml_cvh1
#' @description Calculated CVH1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh1 <- function(model){
  
  # total cvh1
  CVH1_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
  # cvh1 at different levels
  CVH1_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
  names(CVH1_each) <- paste0("CVH1_", model$s.names)
  names(CVH1_total) <- "CVH1_Total"
  
  CVH1s <- c(CVH1_total, CVH1_each)
  
  return(CVH1s)
}

#' @title cvh2_ml
#' @description CVH2 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH2. TODO - we need to cite original CVH2 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH2_eng_1 <- cvh2_ml(english_MA, boot = 10)
#' CVH2_eng_2 <- cvh2_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH2_fish_1 <- cvh2_ml(model, boot = 10)
#' CVH2_fish_2 <- cvh2_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH2_lim_1 <- cvh2_ml(lim_MR, boot = 10)
#' CVH2_lim_2 <- cvh2_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh2_ml <- function(model,
                    boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh2_ml cannot take models with heterogeneous variance.")}
  
  CVH2s <- ml_cvh2(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CVH2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH2 <- ml_cvh2(tmp)
      })
    } else{
      CVH2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH2 <- ml_cvh2(tmp)
        return(CVH2) })
    }
    # Summarise the bootstrapped distribution.
    CVH2s_each_95 <- data.frame(t(apply(CVH2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH2s <-  round(CVH2s_each_95, digits = 3)
    colnames(CVH2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVH2s)
}


#' @title ml_cvh2
#' @description Calculated CVH2 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh2 <- function(model){
  
  # total cvh2
  CVH2_total <- sum(model$sigma2) / (model$beta[[1]])^2
  # cvh2 at different levels
  CVH2_each <-  model$sigma2 / (model$beta[[1]])^2
  names(CVH2_each) <- paste0("CVH2_", model$s.names)
  names(CVH2_total) <- "CVH2_Total"
  
  CVH2s <- c(CVH2_total, CVH2_each)
  
  return(CVH2s)
}

#' @title m1_ml
#' @description M1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M1 - TODO - we need to cite original M1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M1. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M1_eng_1 <- m1_ml(english_MA, boot = 10)
#' M1_eng_2 <- m1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M1_fish_1 <- m1_ml(model, boot = 10)
#' M1_fish_2 <- m1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M1_lim_1 <- m1_ml(lim_MR, boot = 10)
#' M1_lim_2 <- m1_ml(lim_MR)
#' }
#' @references TODO
#' @export

m1_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m1_ml cannot take models with heterogeneous variance.")}
  
  M1s <- ml_m1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M1 <- ml_m1(tmp)
      })
    } else{
      M1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M1 <- ml_m1(tmp)
        return(M1) })
    }
    # Summarise the bootstrapped distribution.
    M1s_each_95 <- data.frame(t(apply(M1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M1s <-  round(M1s_each_95, digits = 3)
    colnames(M1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M1s)
}


#' @title ml_m1
#' @description Calculated M1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m1 <- function(model){
  
  # total m1
  M1_total <- sqrt(sum(model$sigma2)) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  # m1 at different levels
  M1_each <-  sqrt(model$sigma2) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  names(M1_each) <- paste0("M1_", model$s.names)
  names(M1_total) <- "M1_Total"
  
  M1s <- c(M1_total, M1_each)
  
  return(M1s)
}

#' @title m2_ml
#' @description M2 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M2 - TODO - we need to cite original M2 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M2_eng_1 <- m2_ml(english_MA, boot = 10)
#' M2_eng_2 <- m2_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M2_fish_1 <- m2_ml(model, boot = 10)
#' M2_fish_2 <- m2_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M2_lim_1 <- m2_ml(lim_MR, boot = 10)
#' M2_lim_2 <- m2_ml(lim_MR)
#' }
#' @references TODO
#' @export

m2_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m2_ml cannot take models with heterogeneous variance.")}
  
  M2s <- ml_m2(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M2 <- ml_m2(tmp)
      })
    } else{
      M2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M2 <- ml_m2(tmp)
        return(M2) })
    }
    # Summarise the bootstrapped distribution.
    M2s_each_95 <- data.frame(t(apply(M2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M2s <-  round(M2s_each_95, digits = 3)
    colnames(M2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M2s)
}


#' @title ml_m2
#' @description Calculated CV for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m2 <- function(model){
  
  # total m2
  M2_total <- sum(model$sigma2) / ( (model$beta[[1]])^2 + sum(model$sigma2) )
  # m2 at different levels
  M2_each <-  model$sigma2 / ( (model$beta[[1]])^2 + sum(model$sigma2) )
  names(M2_each) <- paste0("M2_", model$s.names)
  names(M2_total) <- "M2_Total"
  
  M2s <- c(M2_total, M2_each)
  
  return(M2s)
}

# prepare data
prepare_mod_data <- function(mod, moderator = "moderator", stdy = "stdy") {
  # 
  mod_table <- mod$mod_table
  data_trim <- mod$data
  
  # factor
  data_trim[[moderator]] <- factor(
    data_trim[[moderator]],
    levels = mod_table$name,
    labels = levels(mod_table$name)
  )
  
  # scale for point sizes
  data_trim$scale <- sqrt(data_trim[,"vi"])
  
  # total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[[moderator]], function(x) length(x[,"yi"])))
  
  # total number of studies (grouping variable) per moderator level
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])
  
  # addition
  legend <- "Uncertainty"
  group_no <- length(mod_table[["name"]])
  color <- data_trim[[moderator]]
  color2 <- mod_table$name
  fill <- color
  
  # a list of processed elements
  return(list(
    mod_table = mod_table,
    data_trim = data_trim,
    legend = legend,
    group_no = group_no,
    color = color,
    color2 = color2,
    fill = fill
  ))
}

# plot layer 
add_my_layers <- function(mod_table, cbpl) {
  list(
    geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.3),
    
    geom_errorbar(
      data = mod_table,
      aes(x = name, ymin = lowerCL, ymax = upperCL),
      size = 0.7, width = 0.1
    ),
    
    geom_pointrange(
      data = mod_table,
      aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR, fill = color2),
      size = 1.5, fatten = 1, shape = 23
    ),
    
    geom_point(
      data = mod_table,
      aes(y = estimate, x = name),
      shape = 23, fill = "white", size = 4
    ),
    
    coord_flip(),
    theme_bw(),
    guides(fill = "none", colour = "none"),
    
    theme(
      legend.title = element_text(size = 9),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      axis.text.y = element_text(size = 12, colour = "black", hjust = 0.5, angle = 90),
      axis.text.x = element_text(size = 10, colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
    ),
    
    scale_fill_manual(values = cbpl),
    scale_colour_manual(values = cbpl)
  )
}

# model selection
eval(metafor:::.MuMIn)
eval(metafor:::.glmulti)

# publication bias test
pub_bias_plot <- function(plot, fe_model, v_model = NULL, col = c("red", "blue"), plotadj = -0.05, textadj = 0.05, branch.size = 1.2, trunk.size = 3){
  
  # Add check to make sure it's an intercept ONLY model being added. Message to user if not.
  if(length(fe_model$b) > 1){
    stop("The model you are trying to add to the plot is not an intercept only model. Please ensure you have fit an intercept only meta-analysis. See vignette for details: https://daniel1noble.github.io/orchaRd/")
  }
  
  # Get the predictions from the final model and create a label for the plot
  pub_bias_data <- get_ints_dat(fe_model, type = "br")
  
  if(is.null(v_model)){
    
    # Add to Existing Orchard Plot
    plot + geom_pub_stats_yang(pub_bias_data, plotadj = plotadj, textadj = textadj, branch.size = branch.size, trunk.size = trunk.size) 
    
  } else{
    # Extract the corrected meta-analytic mean and CI
    pub_bias_data2 <- get_ints_dat(v_model, type = "bc")
    
    plot + geom_pub_stats_yang(pub_bias_data, plotadj = plotadj, textadj = textadj, branch.size = branch.size, trunk.size = trunk.size) + geom_pub_stats_naka(pub_bias_data2, plotadj = plotadj, textadj = textadj, branch.size = branch.size, trunk.size = trunk.size) 
  }		
}


#####################
## Helper functioons
#####################


#' @title geom_pub_stats_yang
#' @description This function adds a corrected meta-analytic mean, sensu Yang et al. 2023, confidence interval and text annotation to an intercept only orchard plot.
#' @param data The data frame containing the corrected meta-analytic mean and confidence intervals.
#' @param col The colour of the mean and confidence intervals.
#' @param plotadj The adjustment to the x-axis position of the mean and confidence intervals.
#' @param textadj The adjustment to the y-axis position of the mean and confidence intervals for the text displaying the type of correction.
#' @param branch.size Size of the confidence intervals.
#' @param trunk.size Size of the mean, or central point.
#' @return A list of ggplot2 objects to be added to the orchard plot.
#' @author Daniel Noble - daniel.noble@anu.edu.au

geom_pub_stats_yang <-  function(data, col = "red", plotadj = -0.05, textadj = 0.05, branch.size = 1.2, trunk.size = 3){
  list(ggplot2::geom_point(data = data[[1]], ggplot2::aes(x = name, y = pred), color = col, shape = "diamond", position = position_nudge(plotadj), size = trunk.size), 
       ggplot2::geom_linerange(data = data[[1]], ggplot2::aes(x = name, ymin = ci.lb, ymax = ci.ub), color = col, position = position_nudge(plotadj), size = branch.size),
       ggplot2::annotate("text", x = 1+plotadj-textadj, y = data[[1]]$pred+textadj, label = data[[2]], color = col, size = 4, hjust = data[[1]]$ci.ub -0.2)	
  )
}

#' @title geom_pub_stats_naka
#' @description This function adds a corrected meta-analytic mean, sensu Nakagawa et al. 2022, confidence interval and text annotation to an intercept only orchard plot.
#' @param data The data frame containing the corrected meta-analytic mean and confidence intervals.
#' @param col The colour of the mean and confidence intervals.
#' @param plotadj The adjustment to the x-axis position of the mean and confidence intervals.
#' @param textadj The adjustment to the y-axis position of the mean and confidence intervals for the text displaying the type of correction.
#' @param branch.size Size of the confidence intervals.
#' @param trunk.size Size of the mean, or central point.
#' @return A list of ggplot2 objects to be added to the orchard plot.
#' @author Daniel Noble - daniel.noble@anu.edu.au

geom_pub_stats_naka <- function(data, col = "blue", plotadj = -0.05, textadj = 0.05, branch.size = 1.2, trunk.size = 3) {
  list(ggplot2::geom_point(data = data[[1]], ggplot2::aes(x = name, y = pred), color = col, shape = "diamond", position = position_nudge(abs(plotadj)), size = trunk.size), 
       ggplot2::geom_linerange(data = data[[1]], ggplot2::aes(x = name, ymin = ci.lb, ymax = ci.ub), color = col, position = position_nudge(abs(plotadj)), size = branch.size), 
       ggplot2::annotate("text", x = 1+abs(plotadj)+textadj, y = data[[1]]$pred-textadj, label = data[[2]], color = col, size = 4, hjust = data[[1]]$ci.ub +0.2))
}

#' @title get_ints_dat
#' @description This function extracts the corrected meta-analytic mean and confidence intervals from a model object.
#' @param model The rma model object containing the corrected meta-analytic mean and confidence intervals.
#' @param type The type of correction to extract the corrected meta-analytic mean and confidence intervals from. "br" (i.e., Bias Robust) for Yang et al. 2023, "bc" (i.e., Bias-Corrected) for Nakagawa et al. 2023.
#' @return A list containing the corrected meta-analytic mean and confidence intervals, and a label for the plot.
#' @author Daniel Noble - daniel.noble@anu.edu.au

get_ints_dat <- function(model, type = c("bc", "br")){
  # Extract the corrected meta-analytic mean and CI
  type = match.arg(type)
  
  dat <- data.frame(name  =  "Intrcpt", 
                    pred = model$b["intrcpt",], 
                    ci.lb = model$ci.lb[1], 
                    ci.ub = model$ci.ub[1])
  if(type == "bc"){
    lab <- paste0("Bias Corrected Estimate: ", round(dat$pred, 2), 
                  ", 95% CI [", round(dat$ci.lb, 2), ",", round(dat$ci.ub, 2), "]")
  }
  
  if(type == "br"){
    lab <- paste0("Bias Robust Estimate: ", round(dat$pred, 2), 
                  ", 95% CI [", round(dat$ci.lb, 2), ",", round(dat$ci.ub, 2), "]")
  }
  
  return(list(dat, lab))
}

# leave one out
leave_one_out <- function(model, group, vcalc_args = NULL, robust_args = NULL, phylo_args = NULL) {

  # Check if we have at least 2 groups
  if (length(unique(model$data[[group]])) < 2) {
    stop("Need at least 2 groups for leave-one-out analysis", call. = FALSE)
  }
  
  if (!is.null(vcalc_args)) {
    .validate_vcalc_args(model$data, vcalc_args)
  } 
  
  if (!is.null(robust_args)) {
    .validate_robust_args(model$data, robust_args)
  }
  
  if (!is.null(phylo_args)) {
    .validate_phylo_args(model, phylo_args) 
  }
  
  # Run leave-one-out analysis
  models_outputs <- .run_leave1out(model, group, vcalc_args, robust_args, phylo_args)
  estimates      <- .get_estimates(models_outputs, group)
  effect_sizes   <- .get_effectsizes(models_outputs, group)
  
  # Get the original model results
  orig_mod_results <- mod_results(model, group = group)
  
  # Immitates the output of mod_results().
  #   - mod_table: In this case, the estimates from each model ran
  #   - data:  Effect sizes and sampling variances from each model
  # Plus one more element:
  #   - orig_mod_results: mod_results() output from the original model
  output <- list(mod_table = estimates,
                 data = effect_sizes,
                 orig_mod_results = orig_mod_results)
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}


#' Fit Multiple Meta-Analytic Models For Leave-One-Out Analysis
#'
#' Iteratively refits a meta-analytic model, leaving out one level of a specified 
#' grouping variable in each iteration. This internal function handles the actual model
#' refitting process for the leave-one-out analysis.
#'
#' @param model A fitted metafor model object containing a \code{data} element and \code{call} object.
#' @param group A character string specifying the column in \code{model$data} that 
#'        defines the groups to be omitted one at a time.
#' @param vcalc_args Optional list of arguments for the variance-covariance calculation using 
#'        metafor's vcalc function. See \code{leave_one_out()} for details.
#' @param robust_args Optional list of arguments for robust variance estimation using
#'        metafor's robust function. See \code{leave_one_out()} for details.
#' @param phylo_args Optional list of arguments for phylogenetic matrix calculation.
#'        See \code{leave_one_out()} for details.
#'
#' @details The function creates a subset of the original data for each unique value in the
#'          grouping variable, removing that group and refitting the model using the same
#'          specification as the original model. If variance-covariance matrices or 
#'          phylogenetic corrections were used in the original model, these are 
#'          recalculated for each subset. If a model fails to fit, NULL is returned for that group.
#'
#' @return A named list of models, each fitted after omitting one group. Names correspond 
#'         to the omitted group IDs. Any failed model fits will be NULL.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal

.run_leave1out <- function(model, group, vcalc_args = NULL, robust_args = NULL, phylo_args = NULL) {
  group_ids <- unique(model$data[[group]])
  
  models_outputs <- lapply(group_ids, function(id_left_out) {
    # Create a new call to fit the model. Modify the data to leave out the group
    # and change de VCV and phylo matrix if needed. Then evaluate the new call.
    
    tmp_model_call <- model$call
    tmp_model_call$data <- subset(model$data, model$data[[group]] != id_left_out)
    
    # If vcalc_args are provided, create a temporary VCV matrix
    if (!is.null(vcalc_args)) {
      tmp_model_call$V <- .create_tmp_vcv(tmp_model_call$data, vcalc_args)
    }
    
    # If the model uses phylogenetic matrix, recalculate it using the original tree
    # The model object contains the correlation matrices in 'R'. This is a list
    # where the names are the random effects and the elements are the correlation matrices.
    # So, first create the new matrix, then use it as the matrix linked to phylo_args$species_colname.
    if (!is.null(phylo_args)) {
      tmp_phylo_matrix <- .create_tmp_phylo_matrix(tmp_model_call$data, phylo_args) 
      tmp_model_call$R[[phylo_args$species_colname]] <- tmp_phylo_matrix
    }
    
    # Evaluate the new call. If something happens, return NULL.
    # In some cases the fixed or random effects are not represented
    # when one group is left out and the model fails to fit.
    tmp_res <- tryCatch({
      eval(tmp_model_call)
    }, error = function(e) {
      warning(sprintf("Error fitting model when leaving out '%s': %s", 
                      id_left_out, e$message))
      return(NULL)
    })
    
    if(!is.null(robust_args)) {
      # cluster_var has to be a vector, not a string. robust_args$cluster is a string.
      cluster_var <- tmp_model_call$data[[robust_args$cluster]]
      
      if(!is.null(robust_args$clubSandwich)) {
        clubSandwich_arg <- robust_args$clubSandwich
      } else {
        clubSandwich_arg <- FALSE
      }
      
      tmp_res <- metafor::robust(tmp_res, cluster = cluster_var, clubSandwich = clubSandwich_arg)
    }
    
    # Return the model output so it is saved in 'models_outputs' list
    tmp_res
  })
  
  names(models_outputs) <- group_ids
  
  return(models_outputs)
}

#' Create Temporary Variance-Covariance Matrix
#'
#' Creates a variance-covariance matrix for a subset of data using metafor's vcalc function.
#'
#' @param data A data frame containing the variables specified in vcalc_args.
#' @param vcalc_args A list of arguments for metafor::vcalc function, including:
#'   \itemize{
#'     \item vi: Name of the variance column
#'     \item cluster: Name of the clustering variable column
#'     \item obs: Name of the observation ID column
#'     \item rho: Correlation coefficient between effect sizes
#'   }
#'
#' @return A variance-covariance matrix for use in meta-analytic models.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.create_tmp_vcv <- function(data, vcalc_args) {
  tryCatch({
    metafor::vcalc(vi      = data[[vcalc_args$vi]],
                   cluster = data[[vcalc_args$cluster]],
                   obs     = data[[vcalc_args$obs]],
                   data    = data,
                   rho     = vcalc_args$rho)
  }, error = function(e) {
    stop(sprintf("Error creating VCV: %s", e$message))
  })
}


#' Validate Variance-Covariance Calculation Arguments
#'
#' Ensures that the arguments provided for variance-covariance calculation are
#' valid and refer to existing variables in the model data. Performs checks on
#' the structure and content of the vcalc_args list.
#'
#' @param model_data A data frame containing the variables used in the model.
#' @param vcalc_args A list of arguments for the metafor::vcalc function.
#'
#' @return The validated vcalc_args list if all checks pass.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.validate_vcalc_args <- function(model_data, vcalc_args) {
  if (!is.list(vcalc_args)) {
    stop("vcalc must be a list with the arguments for the 'vcalc' function: e.g., vcalc_args = list(vi = 'lnrr_vi', cluster = 'paper_ID', obs = 'es_ID', rho = 0.5)",
         call. = FALSE)
  }
  
  # Check if required arguments for vcalc are present
  if (!all(c("vi", "cluster", "obs", "rho") %in% names(vcalc_args))) {
    stop("vcalc_args must contain at least the following elements: 'vi', 'cluster', 'obs', 'rho'", call. = FALSE)
  }
  
  # Check if the vcalc arguments are present in the model data
  if (is.null(model_data[[vcalc_args$vi]]) || is.null(model_data[[vcalc_args$cluster]]) || is.null(model_data[[vcalc_args$obs]])) {
    stop("One or more of the vcalc arguments are not found in the model data", call. = FALSE)
  }
  
  return(vcalc_args)
}

#' Validate Robust Variance Estimation Arguments
#'
#' Validates that robust_args contains the required parameters and that they
#' reference valid columns in the data. Ensures that the cluster variable exists
#' in the provided model data.
#'
#' @param model_data A data frame containing the variables used in the model.
#' @param robust_args A list of arguments for the metafor::robust function.
#'
#' @return The validated robust_args list if all checks pass.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.validate_robust_args <- function(model_data, robust_args) {
  if (!is.list(robust_args)) {
    stop("robust_args must be a list with the arguments for the 'robust' function: e.g., robust_args = list(cluster = 'paper_ID')",
         call. = FALSE)
  }
  
  if (!("cluster" %in% names(robust_args))) {
    stop("robust_args must contain at least the following elements: 'cluster'", call. = FALSE)
  }
  
  # Check if the cluster variable is present in the model data
  if (is.null(model_data[[robust_args$cluster]])) {
    stop("The cluster variable specified in robust_args is not found in the model data", call. = FALSE)
  }
  
  return(robust_args)
}


#' Get Leave-One-Out Model Estimates
#'
#' Extracts and combines the model estimates from each leave-one-out iteration into a 
#' single data frame. Applies mod_results() to each model and combines the resulting
#' mod_table values.
#'
#' @param outputs A named list of model objects from leave-one-out analysis where each
#'        name corresponds to the omitted group ID.
#' @param group A character string specifying the grouping variable used in the analysis.
#'
#' @return A data frame of model estimates with an added column 'name' indicating the omitted group.
#'         Contains the same structure as the mod_table from mod_results() for each model iteration.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.get_estimates <- function(outputs, group) {
  # Call `mod_results` for each model ran in the leave-one-out,
  # transform its output to a dataframe, and then rbind()  
  # to create a long data frame with the estimates of all the models.
  estimates <- do.call(rbind, lapply(names(outputs), function(name) {
    res <- mod_results(outputs[[name]], group = group)
    df <- res$mod_table
    df$name <- name
    df
  }))
  
  row.names(estimates) <- NULL
  return(estimates)
}


#' Get Leave-One-Out Effect Sizes
#'
#' Extracts and aggregates effect size data from each leave-one-out iteration into a
#' single data frame. Applies mod_results() to each model and combines the resulting
#' data values.
#'
#' @param outputs A named list of model objects from leave-one-out analysis where each
#'        name corresponds to the omitted group ID.
#' @param group A character string specifying the grouping variable used in the analysis.
#'
#' @return A data frame of effect sizes with a column 'moderator' indicating the omitted group.
#'         Contains the same structure as the data element from mod_results() for each model iteration.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.get_effectsizes <- function(outputs, group) {
  effect_sizes <- do.call(rbind, lapply(names(outputs), function(name) {
    res <- mod_results(outputs[[name]], group = group)
    df <- res$data
    df$moderator <- name  
    df
  }))
  
  row.names(effect_sizes) <- NULL
  return(effect_sizes)
}

#' Prune Phylogenetic Tree For Leave-One-Out Analysis
#' 
#' Removes species from a phylogenetic tree that are not present in the current data subset
#' during leave-one-out analysis.
#'
#' @param tree A phylogenetic tree object of class "phylo".
#' @param species_names A vector of species names that should remain in the pruned tree.
#'        These are the species present in the current data subset.
#' 
#' @return A pruned phylogenetic tree containing only the tips for species in species_names.
#'
#' @author Daniel Noble  - daniel.noble@anu.edu.au
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.prune_tree <- function(tree, species_names) {
  tree_species <- tree$tip.label
  data_species <- unique(species_names)
  
  species_to_prune <- setdiff(tree_species, data_species)
  
  if (length(species_to_prune) > 0) {
    tree <- ape::drop.tip(tree, species_to_prune)
  }
  
  return(tree)
}

#' Create Temporary Phylogenetic Matrix For Leave-One-Out Analysis
#' 
#' Creates a correlation matrix from a phylogenetic tree for the species remaining
#' in the data subset during leave-one-out analysis.
#'
#' @param data A data frame containing the variables specified in phylo_args,
#'        representing the current data subset after removing one group.
#' @param phylo_args A list of arguments for the phylogenetic matrix calculation, including:
#'   \itemize{
#'     \item tree: A phylogenetic tree object of class "phylo"
#'     \item species_colname: Name of the column in the data that contains species names
#'   }
#'
#' @return A phylogenetic correlation matrix for the species in the data subset.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.create_tmp_phylo_matrix <- function(data, phylo_args) {
  orig_tree <- phylo_args$tree
  species_colname <- phylo_args$species_colname
  
  # Remove species that are left out
  pruned_tree <- .prune_tree(orig_tree, data[[species_colname]])
  pruned_tree <- ape::compute.brlen(pruned_tree)
  
  # Compute the phylo matrix 
  tmp_phylo_matrix <- ape::vcv(pruned_tree, corr = TRUE)
  
  return(tmp_phylo_matrix)
}

#' Validate Phylogenetic Arguments
#'
#' Checks the validity of the arguments provided for phylogenetic matrix calculations.
#' Ensures the tree is a proper phylogenetic object and that the species column name
#' is correctly linked to a random effect in the model.
#'
#' @param model A metafor model object with random effects.
#' @param phylo_args A list containing the arguments for the phylogenetic matrix calculation:
#'   \itemize{
#'     \item tree: A phylogenetic tree object of class "phylo"
#'     \item species_colname: Name of the column in the data corresponding to species names,
#'           which should be a random factor linked to a matrix in the model
#'   }
#'
#' @return The validated phylo_args list if all checks pass.
#'
#' @keywords internal
.validate_phylo_args <- function(model, phylo_args) {
  if (!is.list(phylo_args)) {
    stop("phylo_args must be a list with the arguments for the phylogenetic matrix calculation: e.g., phylo_args = list(tree = tree, species_colname = 'species')",
         call. = FALSE)
  }
  
  if (!all(c("tree", "species_colname") %in% names(phylo_args))) {
    stop("phylo_args must contain at least the following elements: 'tree', 'species_colname'", call. = FALSE)
  }
  
  if (class(phylo_args$tree) != "phylo") {
    stop("The 'tree' argument in phylo_args must be a phylogenetic tree object", call. = FALSE)
  }
  
  # Check if the species_colname is the name of a random factor linked to a matrix in the model
  linked_random_factors <- names(model$Rfix[model$Rfix == TRUE])
  if (!(phylo_args$species_colname %in% linked_random_factors)) {
    stop("The 'species_colname' argument in phylo_args must be the name of the random factor linked to the phylo matrix in the model",
         call. = FALSE)
  }
  
  return(phylo_args)
}

#--------define a function to run moderator analysis--------#
run_mod <- function(mod_var, data, V_mat, yi_var) {
  metafor::rma.mv(
    yi = get(yi_var, data),
    V = V_mat,
    random = list(~1 | Strain, ~1 | Study_ID, ~1 | ES_ID),
    mods = as.formula(paste0("~ ", mod_var, " - 1")),
    test = "t",
    method = "REML",
    data = data,
    subset = !is.na(data[[mod_var]]),
    sparse = TRUE
  )
}

#--------define a function to run null model--------#
run_null_mod <- function(data, V, yi) {
  metafor::rma.mv(
    yi = data[[yi]],
    V = V,
    random = list(~1 | Strain, ~1 | Study_ID, ~1 | ES_ID),
    test = "t",
    method = "REML",
    data = data,
    sparse = TRUE
  )
}

#--------define a function to summarize model--------#
sum_model <- function(model, model_label = NULL, cluster_var = "Study_ID", data = NULL) {
  # print summary
  cat("------ Model summary ------\n")
  print(summary(model))
  
  # robust error
  cat("\n------ Robust variance estimation ------\n")
  robust_out <- metafor::robust(model, cluster = model$data[[cluster_var]], adjust = TRUE)
  print(robust_out)
  
  # Heterogeneity
  cat("\n------ I^2 Statistics ------\n")
  print(i2_ml(model))
  
  # visualization
  if (!is.null(data)) {
    plot <- orchard_plot(
      model,
      mod = "1",
      xlab = "Effect size",
      group = "Strain",
      k = TRUE,
      g = TRUE,
      trunk.size = 1.5,
      data = data
    )
    
    if (!is.null(model_label)) {
      plot <- plot + scale_x_discrete(labels = model_label)
    }
    
    print(plot)
  } else {
    cat("\n[No data provided for plotting.]\n")
  }
}


#--------define a function for model selection-------#
run_model_selection_pipeline <- function(yi_name,
                                         V_matrix,
                                         fixed_effects,
                                         dat,
                                         model_rds_path,
                                         dredge_rds_path,
                                         best_models_rds_path,
                                         label = "Model",
                                         sparse = TRUE,
                                         rerun_model = FALSE,
                                         rerun_dredge = FALSE) {
  cat("\n=========== Running:", label, "===========\n")
  
  # step 1: Extract yi values from the data
  yi_vals <- dat[[yi_name]]
  
  # step 2: Fit or load full model
  model_rds_full <- here(model_rds_path)
  if (rerun_model || !file.exists(model_rds_full)) {
    cat("\n--- Fitting full model ---\n")
    full_model <- rma.mv(
      yi = yi_vals,
      V = V_matrix,
      random = list(~1 | Strain, ~1 | Study_ID, ~1 | ES_ID),
      mods = as.formula(fixed_effects),
      method = "ML",
      test = "t",
      data = dat,
      sparse = sparse
    )
    saveRDS(full_model, model_rds_full)
    cat(paste0("\n--- Saved full model to: ", model_rds_full, " ---\n"))
  } else {
    cat("\n--- Loading full model from cache ---\n")
    full_model <- readRDS(model_rds_full)
  }
  
  
  # step 3: Dredge or load
  dredge_rds_full <- here(dredge_rds_path)
  if (rerun_dredge || !file.exists(dredge_rds_full)) {
    cat("\n--- Running dredge (this may take time) ---\n")
    candidate_models <- dredge(full_model, beta = "none", evaluate = TRUE, rank = "AICc", trace = 2)
    saveRDS(candidate_models, dredge_rds_full)
    cat(paste0("\n--- Saved dredge results to: ", dredge_rds_full, " ---\n"))
  } else {
    cat("\n--- Loading cached candidate models ---\n")
    candidate_models <- readRDS(dredge_rds_full)
  }
  
  # step 4: Best models (ΔAICc ≤ 2)
  best_models <- subset(candidate_models, delta <= 2, recalc.weights = FALSE)
  best_rds <- here(best_models_rds_path)
  saveRDS(best_models, best_rds)
  
  # step 5: View in browser
  print(DT::datatable(best_models))
  
  # step 6: Variable importance
  cat("\n--- Variable Importance (SW) ---\n")
  var_imp <- sw(candidate_models)
  print(var_imp)
  
  # step 7: Model averaging summary
  cat("\n--- Model Averaged Coefficients ---\n")
  avg_summary <- summary(model.avg(candidate_models))
  coef_table <- avg_summary$coefmat.full
  coef_rounded <- dfround(coef_table, 3)
  print(coef_rounded)
  
  # return useful objects
  return(list(
    full_model = full_model,
    candidate_models = candidate_models,
    best_models = best_models,
    sw = var_imp,
    avg_coefmat = coef_rounded
  ))
}

#--------define a function for alternative model selection--------#
# custom fitting function for glmulti
custom_rma_mv <- function(formula, data, ...) {
  metafor::rma.mv(
    formula,
    V = VCV_E,
    random = list(~1 | Strain, ~1 | Study_ID, ~1 | ES_ID),
    data = data,
    test = "t",
    method = "ML",
    sparse = TRUE,
    ...
  )
}



analyze_mlmr_models <- function(candidate_models) {

  dfround <- function(x, digits) {
    x[] <- lapply(x, function(y) if (is.numeric(y)) round(y, digits) else y)
    return(x)
  }
  
  best_models <- weightable(candidate_models)
  best_models_subset <- best_models[best_models$aicc <= min(best_models$aicc) + 2, ]

  best_models_table <- DT::datatable(best_models_subset)
  
  ww <- exp(-(candidate_models@crits - candidate_models@crits[1]) / 2)
  ww <- ww / sum(ww)
  
  clartou <- function(x) {
    pieces <- sort(strsplit(x, ":")[[1]])
    if (length(pieces) > 1)
      paste(pieces[1], ":", pieces[2], sep = "")
    else x
  }
  
  tet <- lapply(candidate_models@formulas, function(x)
    sapply(attr(delete.response(terms(x)), "term.labels"), clartou))
  
  allt <- unique(unlist(tet))
  imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
  vip_predictors <- data.frame(predictor = names(imp), weight = as.numeric(imp))
  
  coef_table <- coef(candidate_models, varweighting = "Johnson")
  coef_df <- as.data.frame(coef_table)
  coef_df <- data.frame(
    Estimate = coef_df$Estimate,
    SE = sqrt(coef_df$`Uncond. variance`),
    Importance = coef_df$Importance,
    row.names = row.names(coef_df)
  ) %>%
    mutate(
      z = Estimate / SE,
      p = 2 * pnorm(abs(z), lower.tail = FALSE),
      ci.lb = Estimate - qnorm(1 - 0.05 / 2, lower.tail = TRUE) * SE,
      ci.ub = Estimate + qnorm(1 - 0.05 / 2, lower.tail = TRUE) * SE
    )
  

  final_table <- coef_df[order(coef_df$Importance, decreasing = TRUE),
                         c("Estimate", "SE", "z", "p", "ci.lb", "ci.ub", "Importance")] %>%
    dfround(3)
  
  return(list(
    best_models = best_models_table,
    variable_importance = vip_predictors,
    multimodel_inference = final_table
  ))
}


analyze_mlmr_models <- function(candidate_models) {

  dfround <- function(x, digits) {
    x[] <- lapply(x, function(y) if (is.numeric(y)) round(y, digits) else y)
    return(x)
  }
  
  best_models <- weightable(candidate_models)
  best_models_subset <- best_models[best_models$aicc <= min(best_models$aicc) + 2, ]
  
  best_models_table <- DT::datatable(best_models_subset)
  
  ww <- exp(-(candidate_models@crits - candidate_models@crits[1]) / 2)
  ww <- ww / sum(ww)
  
  clartou <- function(x) {
    pieces <- sort(strsplit(x, ":")[[1]])
    if (length(pieces) > 1)
      paste(pieces[1], ":", pieces[2], sep = "")
    else x
  }
  
  tet <- lapply(candidate_models@formulas, function(x)
    sapply(attr(delete.response(terms(x)), "term.labels"), clartou))
  
  allt <- unique(unlist(tet))
  imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
  vip_predictors <- data.frame(predictor = names(imp), weight = as.numeric(imp))
  
  coef_table <- coef(candidate_models, varweighting = "Johnson")
  coef_df <- as.data.frame(coef_table)
  coef_df <- data.frame(
    Estimate = coef_df$Estimate,
    SE = sqrt(coef_df$`Uncond. variance`),
    Importance = coef_df$Importance,
    row.names = row.names(coef_df)
  ) %>%
    mutate(
      z = Estimate / SE,
      p = 2 * pnorm(abs(z), lower.tail = FALSE),
      ci.lb = Estimate - qnorm(1 - 0.05 / 2, lower.tail = TRUE) * SE,
      ci.ub = Estimate + qnorm(1 - 0.05 / 2, lower.tail = TRUE) * SE
    )
  
  final_table <- coef_df[order(coef_df$Importance, decreasing = TRUE),
                         c("Estimate", "SE", "z", "p", "ci.lb", "ci.ub", "Importance")] %>%
    dfround(3)
  
  return(list(
    best_models = best_models_table,
    variable_importance = vip_predictors,
    multimodel_inference = final_table
  ))
}


#--------define a function for sensitivity analysis--------#
post_strat_analysis <- function(moderator, data, VCV_E, response = "lnRR_Ea") {
  res <- rma.mv(
    yi = as.formula(paste(response, "~", moderator)),
    V = VCV_E,
    random = list(~1|Study_ID, ~1|ES_ID),
    test = "t",
    method = "REML",
    data = data,
    subset = !is.na(data[[moderator]]),
    sparse = TRUE
  )
  
  sav <- emmprep(res)
  
  preds <- emmeans(sav, as.formula(paste("~", moderator)), type = "response")
  preds_df <- as.data.frame(preds)
  
  weights <- rep(1 / nrow(preds_df), nrow(preds_df))
  
  post_strat_mean <- sum(weights * preds_df$emmean)
  
  vcov_mat <- vcov(preds)
  post_strat_var <- as.numeric(t(weights) %*% vcov_mat %*% weights)
  post_strat_se <- sqrt(post_strat_var)
  
  t_stat <- post_strat_mean / post_strat_se
  df <- res$ddf[1]
  p <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  
  result <- data.frame(
    moderator = moderator,
    post_strat_mean = post_strat_mean,
    post_strat_se = post_strat_se,
    df = df,
    p = p
  )
  
  return(result)
}

#--------define a function for sensitivity analysis2--------#
leave_one_out_analysis <- function(data, moderator = "Strain", response = "lnRR_Ea", vi = "lnRRV_E", cluster = "Study_ID") {
  data$leave_out <- as.factor(data[[moderator]])
  levels_leave_out <- levels(data$leave_out)
  
  leave1out_mod <- list()
  VCV_leave1out <- list()
  
  for (i in 1:length(levels_leave_out)) {
    dat2 <- data %>% filter(leave_out != levels_leave_out[i])
    
    VCV_leave1out[[i]] <- impute_covariance_matrix(
      vi = dat2[[vi]],
      cluster = dat2[[cluster]],
      r = 0.5
    )
    
    leave1out_mod[[i]] <- rma.mv(
      yi = as.formula(paste(response, "~ 1")),
      V = VCV_leave1out[[i]],
      random = list(~1|Strain, ~1|Study_ID, ~1|ES_ID),
      method = "REML",
      test = "t",
      data = dat2,
      sparse = TRUE
    )
  }
  
  est.func <- function(mod) {
    data.frame(
      est = mod$b,
      lower = mod$ci.lb,
      upper = mod$ci.ub
    )
  }
  
  leave1out_results <- lapply(leave1out_mod, est.func) %>%
    bind_rows() %>%
    mutate(left_out = levels_leave_out)
  
  leave1out_results$left_out <- factor(leave1out_results$left_out, levels = leave1out_results$left_out)
  
  return(leave1out_results)
}

#--------plot function--------#

# Custom mod_results function (provided by user)
mod_results <- function(model, mod = "1", group, data, N = NULL, weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...) {
  if (any(grepl("-1|0", as.character(model$formula.mods)))) {
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }
  
  if (missing(model)) {
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }
  
  if (all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {
    stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")
  }
  
  if (missing(group)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  if (missing(data)) {
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  
  if (is.null(stats::formula(model))) {
    model$formula.mods <- ~ 1
    dat_tmp <- data
    dat_tmp$`1` <- "Intrcpt"
    model$data <- dat_tmp
  } else {
    model$data <- data
  }
  
  if (model$test == "t") {
    df_mod <- as.numeric(model$ddf[[1]])
  } else {
    df_mod <- 1.0e6
  }
  
  if (is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k - 1)
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if (is.null(by)) {
      mod_table <- data.frame(
        name = firstup(as.character(mm_pi[, 1]), upper = upper),
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    } else {
      mod_table <- data.frame(
        name = firstup(as.character(mm_pi[, 1]), upper = upper),
        condition = mm_pi[, 2],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    }
    
    data2 <- get_data_raw(model, mod, group, N, data, at = at, subset)
    
    mod_table$name <- factor(mod_table$name, levels = mod_table$name, labels = mod_table$name)
  } else {
    at2 <- list(mod = seq(min(data[, mod], na.rm = TRUE), max(data[, mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula = stats::formula(model), data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k - 1, at = c(at2, at))
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)
    
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if (is.null(by)) {
      mod_table <- data.frame(
        moderator = mm_pi[, 1],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    } else {
      mod_table <- data.frame(
        moderator = mm_pi[, 1],
        condition = mm_pi[, 2],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    }
    
    data2 <- get_data_raw_cont(model, mod, group, N, data, by = by)
  }
  
  output <- list(mod_table = mod_table, data = data2)
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}

# Function to prepare data for plotting
prepare_mod_results <- function(mod_E0, mod_S0, mod_ES0, data, vi_E = "lnRRV_E", vi_S = "lnRRV_S", vi_ES = "lnRRV_ES", group = "Strain") {
  library(dplyr)
  library(metafor)
  
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }
  
  if (!(vi_E %in% names(data) && vi_S %in% names(data) && vi_ES %in% names(data))) {
    stop("One or more variance columns (vi_E, vi_S, vi_ES) not found in data.")
  }
  
  mod_E0_results <- mod_results(mod_E0, mod = "1", group = group, data = data[!is.na(data[[vi_E]]), ])
  class(mod_E0_results) <- "list"
  mod_E0_results$data$moderator <- rep("E", length(mod_E0_results$data$moderator))
  
  mod_S0_results <- mod_results(mod_S0, mod = "1", group = group, data = data[!is.na(data[[vi_S]]), ])
  class(mod_S0_results) <- "list"
  mod_S0_results$data$moderator <- rep("S", length(mod_S0_results$data$moderator))
  
  mod_ES0_results <- mod_results(mod_ES0, mod = "1", group = group, data = data[!is.na(data[[vi_ES]]), ])
  class(mod_ES0_results) <- "list"
  mod_ES0_results$data$moderator <- rep("ES", length(mod_ES0_results$data$moderator))
  
  submerge <- function(...) {
    models <- list(...)
    list(
      data = do.call(rbind, lapply(models, function(x) x$data)),
      mod_table = do.call(rbind, lapply(models, function(x) x$mod_table))
    )
  }
  
  mod_all_results <- submerge(mod_ES0_results, mod_S0_results, mod_E0_results)
  mod_all_results$mod_table$name <- c("ES", "S", "E")
  
  # Calculate K (number of observations) and g (number of unique studies)
  num_studies <- function(data, moderator, stdy) {
    if (!stdy %in% names(data)) {
      warning("Column '", stdy, "' not found in data. Returning NA for number of studies.")
      return(data.frame(moderator = unique(data[[moderator]]), n = NA))
    }
    data %>%
      group_by(.data[[moderator]]) %>%
      summarise(n = length(unique(.data[[stdy]])))
  }
  
  mod_all_results$mod_table <- mod_all_results$mod_table %>%
    mutate(
      sd = c(
        sqrt(mod_ES0$se[[1]]^2 + sum(mod_ES0$sigma2)),
        sqrt(mod_S0$se[[1]]^2 + sum(mod_S0$sigma2)),
        sqrt(mod_E0$se[[1]]^2 + sum(mod_E0$sigma2))
      ),
      sd2 = c(
        sqrt(mod_ES0$se[[1]]^2 + sum(mod_ES0$sigma2[1])),
        sqrt(mod_S0$se[[1]]^2 + sum(mod_S0$sigma2[1])),
        sqrt(mod_E0$se[[1]]^2 + sum(mod_E0$sigma2[1]))
      ),
      sd3 = c(
        sqrt(mod_ES0$se[[1]]^2 + sum(mod_ES0$sigma2[2])),
        sqrt(mod_S0$se[[1]]^2 + sum(mod_S0$sigma2[2])),
        sqrt(mod_E0$se[[1]]^2 + sum(mod_E0$sigma2[2]))
      ),
      df = c(mod_ES0$ddf, mod_S0$ddf, mod_E0$ddf),
      K = as.vector(by(mod_all_results$data, mod_all_results$data$moderator, function(x) length(x$yi))),
      g = as.vector(num_studies(mod_all_results$data, "moderator", "stdy")[, "n"])
    )
  
  return(mod_all_results)
}

# Function to plot model results
plot_base_mod_results <- function(results, cbpl = c("#E69F00", "#56B4E9", "#009E73")) {
  mod_table <- results$mod_table
  mod_table$name <- factor(mod_table$name, levels = c("ES", "S", "E"), labels = c("ES", "S", "E"))
  data_trim <- results$data
  data_trim$moderator <- factor(data_trim$moderator, levels = c("ES", "S", "E"), labels = c("ES", "S", "E"))
  data_trim$scale <- sqrt(data_trim$vi)
  
  p <- ggplot() +
    stat_slab(
      data = mod_table,
      aes(x = name, ydist = dist_student_t(df = df - 2, mu = estimate, sigma = sd)),
      size = 2, stroke = 10, fill = "#FFFACD", alpha = 1
    ) +
    ggbeeswarm::geom_quasirandom(
      data = data_trim,
      aes(y = yi, x = moderator, size = scale, color = moderator),
      alpha = 1, shape = 21, col = "#999999", width = 0.5
    ) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.3) +
    geom_errorbar(
      data = mod_table,
      aes(x = name, ymin = lowerCL, ymax = upperCL),
      size = 0.7, width = 0.1
    ) +
    geom_pointrange(
      data = mod_table,
      aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR, fill = name),
      size = 1.5, fatten = 1, shape = 23
    ) +
    geom_point(
      data = mod_table,
      aes(y = estimate, x = name),
      shape = 23, fill = "white", size = 4
    ) +
    coord_flip() +
    theme_bw() +
    guides(fill = "none", colour = "none") +
    theme(
      legend.title = element_text(size = 9),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      axis.text.y = element_text(size = 12, colour = "black", hjust = 0.5, angle = 90),
      axis.text.x = element_text(size = 10, colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
    ) +
    scale_fill_manual(values = cbpl) +
    scale_colour_manual(values = cbpl) +
    scale_color_npg() +
    scale_fill_npg()
  
  return(p)
}

# Custom mod_results function (unchanged, included for completeness)
mod_results <- function(model, mod = "1", group, data, N = NULL, weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...) {
  if (any(grepl("-1|0", as.character(model$formula.mods)))) {
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }
  
  if (missing(model)) {
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }
  
  if (all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {
    stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")
  }
  
  if (missing(group)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  if (missing(data)) {
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  
  if (is.null(stats::formula(model))) {
    model$formula.mods <- ~ 1
    dat_tmp <- data
    dat_tmp$`1` <- "Intrcpt"
    model$data <- dat_tmp
  } else {
    model$data <- data
  }
  
  if (model$test == "t") {
    df_mod <- as.numeric(model$ddf[[1]])
  } else {
    df_mod <- 1.0e6
  }
  
  if (is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k - 1)
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if (is.null(by)) {
      mod_table <- data.frame(
        name = firstup(as.character(mm_pi[, 1]), upper = upper),
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    } else {
      mod_table <- data.frame(
        name = firstup(as.character(mm_pi[, 1]), upper = upper),
        condition = mm_pi[, 2],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    }
    
    data2 <- get_data_raw(model, mod, group, N, data, at = at, subset)
    
    mod_table$name <- factor(mod_table$name, levels = mod_table$name, labels = mod_table$name)
  } else {
    at2 <- list(mod = seq(min(data[, mod], na.rm = TRUE), max(data[, mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula = stats::formula(model), data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k - 1, at = c(at2, at))
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)
    
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if (is.null(by)) {
      mod_table <- data.frame(
        moderator = mm_pi[, 1],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    } else {
      mod_table <- data.frame(
        moderator = mm_pi[, 1],
        condition = mm_pi[, 2],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    }
    
    data2 <- get_data_raw_cont(model, mod, group, N, data, by = by)
  }
  
  output <- list(mod_table = mod_table, data = data2)
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}

# Updated function to prepare data for lnCVR plotting
# Custom mod_results function (unchanged, included for completeness)
mod_results <- function(model, mod = "1", group, data, N = NULL, weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...) {
  if (any(grepl("-1|0", as.character(model$formula.mods)))) {
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }
  
  if (missing(model)) {
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }
  
  if (all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {
    stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")
  }
  
  if (missing(group)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  if (missing(data)) {
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  
  if (is.null(stats::formula(model))) {
    model$formula.mods <- ~ 1
    dat_tmp <- data
    dat_tmp$`1` <- "Intrcpt"
    model$data <- dat_tmp
  } else {
    model$data <- data
  }
  
  if (model$test == "t") {
    df_mod <- as.numeric(model$ddf[[1]])
  } else {
    df_mod <- 1.0e6
  }
  
  if (is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k - 1)
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if (is.null(by)) {
      mod_table <- data.frame(
        name = firstup(as.character(mm_pi[, 1]), upper = upper),
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    } else {
      mod_table <- data.frame(
        name = firstup(as.character(mm_pi[, 1]), upper = upper),
        condition = mm_pi[, 2],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    }
    
    data2 <- get_data_raw(model, mod, group, N, data, at = at, subset)
    
    mod_table$name <- factor(mod_table$name, levels = mod_table$name, labels = mod_table$name)
  } else {
    at2 <- list(mod = seq(min(data[, mod], na.rm = TRUE), max(data[, mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula = stats::formula(model), data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k - 1, at = c(at2, at))
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)
    
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if (is.null(by)) {
      mod_table <- data.frame(
        moderator = mm_pi[, 1],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    } else {
      mod_table <- data.frame(
        moderator = mm_pi[, 1],
        condition = mm_pi[, 2],
        estimate = mm_pi[, "emmean"],
        lowerCL = mm_pi[, "lower.CL"],
        upperCL = mm_pi[, "upper.CL"],
        lowerPR = mm_pi[, "lower.PI"],
        upperPR = mm_pi[, "upper.PI"]
      )
    }
    
    data2 <- get_data_raw_cont(model, mod, group, N, data, by = by)
  }
  
  output <- list(mod_table = mod_table, data = data2)
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}

# Updated function to prepare data for lnCVR plotting
prepare_mod_results_lncvr <- function(mod_lnCVR_E0, mod_lnCVR_S0, mod_lnCVR_ES0, data, vi_E = "lnCVRV_E", vi_S = "lnCVRV_S", vi_ES = "lnCVRV_ES", group = "Strain") {
  library(dplyr)
  library(metafor)
  
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }
  
  if (!(vi_E %in% names(data) && vi_S %in% names(data) && vi_ES %in% names(data))) {
    stop("One or more variance columns (vi_E, vi_S, vi_ES) not found in data.")
  }
  
  # Process each model
  mod_E0_results <- mod_results(mod_lnCVR_E0, mod = "1", group = group, data = data[!is.na(data[[vi_E]]), ])
  class(mod_E0_results) <- "list"
  mod_E0_results$data$moderator <- rep("E", length(mod_E0_results$data$moderator))
  
  mod_S0_results <- mod_results(mod_lnCVR_S0, mod = "1", group = group, data = data[!is.na(data[[vi_S]]), ])
  class(mod_S0_results) <- "list"
  mod_S0_results$data$moderator <- rep("S", length(mod_S0_results$data$moderator))
  
  mod_ES0_results <- mod_results(mod_lnCVR_ES0, mod = "1", group = group, data = data[!is.na(data[[vi_ES]]), ])
  class(mod_ES0_results) <- "list"
  mod_ES0_results$data$moderator <- rep("ES", length(mod_ES0_results$data$moderator))
  
  # Combine results with consistent columns
  submerge <- function(...) {
    models <- list(...)
    # Ensure consistent columns in data
    data_list <- lapply(models, function(x) {
      df <- x$data
      # Add missing columns with NA if necessary
      expected_cols <- c("yi", "vi", "moderator", "stdy", group)
      for (col in expected_cols) {
        if (!col %in% names(df)) {
          df[[col]] <- NA
        }
      }
      df[, expected_cols, drop = FALSE]
    })
    # Combine data
    combined_data <- do.call(rbind, data_list)
    # Ensure consistent columns in mod_table
    mod_table_list <- lapply(models, function(x) {
      mt <- x$mod_table
      expected_mt_cols <- c("name", "estimate", "lowerCL", "upperCL", "lowerPR", "upperPR")
      for (col in expected_mt_cols) {
        if (!col %in% names(mt)) {
          mt[[col]] <- NA
        }
      }
      mt[, expected_mt_cols, drop = FALSE]
    })
    combined_mod_table <- do.call(rbind, mod_table_list)
    list(data = combined_data, mod_table = combined_mod_table)
  }
  
  mod_all_results <- submerge(mod_ES0_results, mod_S0_results, mod_E0_results)
  mod_all_results$mod_table$name <- c("ES", "S", "E")
  
  # Calculate K and g
  num_studies <- function(data, moderator, stdy) {
    if (!stdy %in% names(data)) {
      warning("Column '", stdy, "' not found in data. Returning NA for number of studies.")
      return(data.frame(moderator = unique(data[[moderator]]), n = NA))
    }
    data %>%
      group_by(.data[[moderator]]) %>%
      summarise(n = length(unique(.data[[stdy]])))
  }
  
  mod_all_results$mod_table <- mod_all_results$mod_table %>%
    mutate(
      sd = c(
        sqrt(mod_lnCVR_ES0$se[[1]]^2 + sum(mod_lnCVR_ES0$sigma2)),
        sqrt(mod_lnCVR_S0$se[[1]]^2 + sum(mod_lnCVR_S0$sigma2)),
        sqrt(mod_lnCVR_E0$se[[1]]^2 + sum(mod_lnCVR_E0$sigma2))
      ),
      sd2 = c(
        sqrt(mod_lnCVR_ES0$se[[1]]^2 + sum(mod_lnCVR_ES0$sigma2[1])),
        sqrt(mod_lnCVR_S0$se[[1]]^2 + sum(mod_lnCVR_S0$sigma2[1])),
        sqrt(mod_lnCVR_E0$se[[1]]^2 + sum(mod_lnCVR_E0$sigma2[1]))
      ),
      sd3 = c(
        sqrt(mod_lnCVR_ES0$se[[1]]^2 + sum(mod_lnCVR_ES0$sigma2[3])),
        sqrt(mod_lnCVR_S0$se[[1]]^2 + sum(mod_lnCVR_S0$sigma2[3])),
        sqrt(mod_lnCVR_E0$se[[1]]^2 + sum(mod_lnCVR_E0$sigma2[3]))
      ),
      df = c(mod_lnCVR_ES0$ddf, mod_lnCVR_S0$ddf, mod_lnCVR_E0$ddf),
      K = as.vector(by(mod_all_results$data, mod_all_results$data$moderator, function(x) length(x$yi))),
      g = as.vector(num_studies(mod_all_results$data, "moderator", "stdy")[, "n"])
    )
  
  return(mod_all_results)
}

# Updated function to create base ggplot object for lnCVR
plot_base_mod_results_lncvr <- function(results, cbpl = c("#E69F00", "#56B4E9", "#009E73")) {
  library(ggplot2)
  library(ggdist)
  library(ggbeeswarm)
  library(ggsci)
  library(dplyr)
  library(stringr)
  
  mod_table <- results$mod_table
  mod_table$name <- factor(mod_table$name, levels = c("ES", "S", "E"), labels = c("ES", "S", "E"))
  data_trim <- results$data
  data_trim$moderator <- factor(data_trim$moderator, levels = c("ES", "S", "E"), labels = c("ES", "S", "E"))
  data_trim$scale <- sqrt(data_trim$vi)
  data_trim$color <- data_trim$moderator
  
  p <- ggplot() +
    stat_slab(
      data = mod_table,
      aes(x = name, ydist = dist_student_t(df = df - 2, mu = estimate, sigma = sd)),
      size = 2, stroke = 10, fill = "#FFFACD", alpha = 1
    ) +
    ggbeeswarm::geom_quasirandom(
      data = data_trim,
      aes(y = yi, x = moderator, size = scale, colour = color),
      col = "#999999", alpha = 1, shape = 21, width = 0.35
    ) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.3) +
    geom_errorbar(
      data = mod_table,
      aes(x = name, ymin = lowerCL, ymax = upperCL),
      size = 0.8, width = 0.1
    ) +
    geom_pointrange(
      data = mod_table,
      aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR, fill = name),
      size = 1.5, fatten = 1, shape = 23
    ) +
    geom_point(
      data = mod_table,
      aes(y = estimate, x = name),
      shape = 23, fill = "white", size = 3
    ) +
    coord_flip() +
    theme_bw() +
    guides(fill = "none", colour = "none") +
    theme(
      legend.title = element_text(size = 9),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      axis.text.y = element_text(size = 12, colour = "black", hjust = 0.5, angle = 90),
      axis.text.x = element_text(size = 10, colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
    ) +
    scale_fill_manual(values = cbpl) +
    scale_colour_manual(values = cbpl) +
    scale_color_npg() +
    scale_fill_npg()
  
  return(p)
}


#--------plot function for moderator analysis--------#
# Function to run mod_results and prepare mod_table with sd and df
run_mod_results <- function(model, mod, group = "Strain", data) {
  library(dplyr)
  library(metafor)
  
  # Validate inputs
  if (!inherits(model, c("rma.mv", "rma", "robust.rma"))) {
    stop("Model must be of class rma.mv, rma, or robust.rma")
  }
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  if (!mod %in% names(data)) {
    stop("Moderator '", mod, "' not found in data")
  }
  
  # Run mod_results
  results <- mod_results(model, mod = mod, group = group, data = data)
  class(results) <- "list"
  
  # Add sd and df to mod_table
  results$mod_table <- results$mod_table %>%
    mutate(
      sd = sqrt(model$se[[1]]^2 + sum(model$sigma2)),
      df = model$ddf
    )
  
  return(results)
}

# Function to extract components from prepare_mod_data output
extract_mod_components <- function(results) {
  # Extract components
  mod_table <- results$mod_table
  data_trim <- results$data_trim
  legend <- results$legend
  group_no <- results$group_no
  color <- results$color
  mod_table$color2 <- results$color2
  fill <- color
  
  # Return as list
  return(list(
    mod_table = mod_table,
    data_trim = data_trim,
    legend = legend,
    group_no = group_no,
    color = color,
    fill = fill
  ))
}

# Function to create base ggplot with geom_quasirandom
plot_base_quasirandom <- function(data_trim, fill_color) {
  library(ggplot2)
  library(ggbeeswarm)
  
  p <- ggplot2::ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = data_trim,
      ggplot2::aes(y = yi, x = moderator, size = scale, color = color),
      alpha = 1, shape = 21, col = "#999999", fill = fill_color, width = 0.5
    )
  
  return(p)
}

#-------------------------------------------------------------------------------------------------
# modified from ihttps://github.com/davidsjoberg/ggsankey/blob/main/R/sankey.R
#-------------------------------------------------------------------------------------------------
utils::globalVariables(c(".", ".data", "x", "node", "next_node", "next_x", "..r"))
# importFrom(ggplot2, "%+replace%")
#' @importFrom ggplot2 %+replace%

# ** Support functions ----------
prepare_params <- function(...) {
  # Prepare aesthics for flow lines
  flow.aes <- list(...)
  removes <- names(flow.aes) %>%
    stringr::str_extract_all(., "(?<=flow.).*") %>% unlist()
  removes2 <- names(flow.aes) %>%
    stringr::str_subset(., "node") %>% unlist()
  flow.aes[c(removes, removes2)] <- NULL
  names(flow.aes) <- names(flow.aes) %>%
    stringr::str_replace_all("flow.", "")
  
  # Prepare aesthics for node boxes
  node.aes <- list(...)
  removes <- names(node.aes) %>%
    stringr::str_extract_all(., "(?<=node.).*") %>% unlist()
  removes2 <- names(node.aes) %>%
    stringr::str_subset(., "flow") %>% unlist()
  node.aes[c(removes, removes2)] <- NULL
  names(node.aes) <- names(node.aes) %>%
    stringr::str_replace_all(., "node.", "")
  
  return(list(flow.aes, node.aes))
}

find_default_space <- function(.df) {
  .df %>%
    dplyr::group_by(.data$n_x) %>%
    dplyr::summarise(n_groups = dplyr::n_distinct(.data$node),
                     freq = sum(.data$freq, na.rm = TRUE)) %>%
    dplyr::mutate(v = .data$freq / .data$n_groups / 4) %>%
    dplyr::pull(.data$v) %>%
    max()
}

sigmoid <- function(x_from, x_to, y_from, y_to, smooth = 5, n = 300) {
  x <- seq(-smooth, smooth, length = n)
  y <- exp(x) / (exp(x) + 1)
  out <- data.frame(x = (x + smooth) / (smooth * 2) * (x_to - x_from) + x_from,
                    y = y * (y_to - y_from) + y_from)
}


#long format -----------------------------------------------------------------
dlong <- function(.df, ..., value = NULL) {
  if("..r" %in% names(.df)) stop("The column name '..r' is not allowed")
  .vars <- dplyr::quos(...)
  
  if(!missing(value)) {
    value_var <- dplyr::enquo(value)
    out <- .df %>%
      dplyr::select(!!!.vars, value = !!value_var) %>%
      dplyr::mutate(..r = dplyr::row_number()) %>%
      tidyr::gather(x, node, -..r, -value) %>%
      dplyr::arrange(.data$..r) %>%
      dplyr::group_by(.data$..r) %>%
      dplyr::mutate(next_x = dplyr::lead(.data$x),
                    next_node = dplyr::lead(.data$node)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-..r) %>%
      dplyr::relocate(value, .after = dplyr::last_col())
  } else {
    out <- .df %>%
      dplyr::select(!!!.vars) %>%
      dplyr::mutate(..r = dplyr::row_number()) %>%
      tidyr::gather(x, node, -..r) %>%
      dplyr::arrange(.data$..r) %>%
      dplyr::group_by(.data$..r) %>%
      dplyr::mutate(next_x = dplyr::lead(.data$x),
                    next_node = dplyr::lead(.data$node)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-..r)
  }
  
  levels <- unique(out$x)
  
  out %>%
    dplyr::mutate(dplyr::across(c(x, next_x), ~factor(., levels = levels)))
}


#' @title sankey_themes
#' @name theme_sankey
#' @aliases theme_alluvial
#' @aliases theme_sankey_bump
#'
#' @description Minimal themes for sankey, alluvial and sankey bump plots
#'
#' @param base_size base font size, given in pts.
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect elements
#'
#' @export
theme_sankey <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black",
                                            size = ggplot2::rel(1)),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }

#' @rdname theme_sankey
#' @export
theme_alluvial <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }

#' @rdname theme_sankey
#' @export
theme_sankey_bump <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line("gray90")
        )
    }
  }


# FLOW LAYER ---------
StatSankeyFlow <- ggplot2::ggproto("StatSankeyFlow", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        flow_data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                          dplyr::summarise(flow_freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                        
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        flow_data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                          dplyr::summarise(flow_freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                        
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE),, .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      df <- data %>%
                                                        dplyr::left_join(flow_data, by = c("n_x", "node"))
                                                      
                                                      
                                                      
                                                      flows <- df %>%
                                                        dplyr::left_join(df %>%
                                                                           dplyr::select(n_x, node, ymin_end = ymin, ymax_end = ymax, xmin_end = xmin, xmax_end = xmax) %>%
                                                                           dplyr::distinct(),
                                                                         by = c("n_next_x" = "n_x", "next_node" = "node")) %>%
                                                        tidyr::drop_na(n_x, node, next_node, n_next_x, ymax_end, ymin_end, xmax_end, xmin_end) %>%
                                                        dplyr::mutate(r = dplyr::row_number()) %>%
                                                        dplyr::arrange(n_x, -r) %>%
                                                        dplyr::select(-r) %>%
                                                        dplyr::group_by(n_x, node) %>%
                                                        dplyr::mutate(cum_flow_freq = cumsum(flow_freq) - flow_freq) %>%
                                                        dplyr::ungroup() %>%
                                                        dplyr::group_by(n_x, n_next_x, node, next_node) %>%
                                                        dplyr::mutate(flow_start_ymax = ymax - cum_flow_freq,
                                                                      flow_start_ymin = flow_start_ymax - flow_freq)
                                                      
                                                      flows <- flows %>%
                                                        dplyr::arrange(n_x, n_next_x, next_node) %>%
                                                        dplyr::group_by(n_next_x, next_node) %>%
                                                        dplyr::mutate(cum_flow_freq_end = cumsum(flow_freq) - flow_freq) %>%
                                                        dplyr::mutate(flow_end_ymax = ymax_end - cum_flow_freq_end,
                                                                      flow_end_ymin = flow_end_ymax - flow_freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      flows <- flows %>%
                                                        dplyr::select(-n_x, -node, -freq, -ymax, -ymin, -xmin, -n_next_x, -next_node, -flow_freq, -ymin_end, -ymax_end, -xmax_end, -cum_flow_freq, -cum_flow_freq_end) %>%
                                                        dplyr::mutate(group = dplyr::row_number())
                                                      
                                                      flows %>%
                                                        dplyr::mutate(smooth = params$smooth) %>%
                                                        as.data.frame()
                                                    })
                                     
                                     
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     
                                     out1 <- sigmoid(data$xmax, data$xmin_end, data$flow_start_ymax, data$flow_end_ymax,
                                                     smooth = data$smooth)
                                     out2 <- sigmoid(data$xmin_end, data$xmax, data$flow_end_ymin, data$flow_start_ymin,
                                                     smooth = data$smooth)
                                     dplyr::bind_rows(out1, out2)
                                   }
)


# FLOW SANKEYBUMP LAYER ---------
StatSankeyBumpFlow <- ggplot2::ggproto("StatSankeyBumpFlow", ggplot2::Stat,
                                       extra_params = c("na.rm", "type", "space", "smooth"),
                                       
                                       setup_data = function(data, params) {
                                         
                                         purrr::map_dfr(unique(data$PANEL),
                                                        ~{
                                                          data <- data %>% dplyr::filter(PANEL == .x)
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(nodes = paste(node, x)) %>%
                                                            dplyr::arrange(x, -value) %>%
                                                            dplyr::mutate(bbb = dplyr::row_number()) %>%
                                                            dplyr::arrange(bbb) %>%
                                                            dplyr::mutate(nodes = fct_reorder(nodes, value, mean)) %>%
                                                            dplyr::arrange(node, x) %>%
                                                            dplyr::group_by(node) %>%
                                                            dplyr::mutate(next_x = dplyr::lead(x),
                                                                          node = nodes,
                                                                          next_node = dplyr::lead(nodes)) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::arrange(x, node)
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                          
                                                          if(!("value" %in% names(data))) {
                                                            flow_data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                              dplyr::summarise(flow_freq = dplyr::n(), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                            
                                                            data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::select(-n_next_x, -next_node) %>%
                                                              dplyr::group_by_all() %>%
                                                              dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                          } else {
                                                            flow_data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                              dplyr::summarise(flow_freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                            
                                                            data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::select(-n_next_x, -next_node) %>%
                                                              dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                              dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                          }
                                                          
                                                          if(is.null(params$space)) {
                                                            params$space <- find_default_space(data)
                                                          }
                                                          
                                                          data <- data %>%
                                                            dplyr::group_by(n_x) %>%
                                                            dplyr::arrange(node) %>%
                                                            dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                          ymin = ymax - freq) %>%
                                                            dplyr::ungroup()
                                                          
                                                          if(params$type == "sankey") {
                                                            data <- data %>%
                                                              dplyr::group_by(n_x) %>%
                                                              dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                            ymax = ymax - max(ymax)/2) %>%
                                                              dplyr::ungroup()
                                                          } else if (params$type == "alluvial"){
                                                            data <- data
                                                          }
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(xmin = n_x,
                                                                          xmax = n_x)
                                                          
                                                          df <- data %>%
                                                            dplyr::left_join(flow_data, by = c("n_x", "node"))
                                                          
                                                          flows <- df %>%
                                                            dplyr::left_join(df %>%
                                                                               dplyr::select(n_x, node, ymin_end = ymin, ymax_end = ymax, xmin_end = xmin, xmax_end = xmax, flow_freq_end = flow_freq) %>%
                                                                               dplyr::distinct(),
                                                                             by = c("n_next_x" = "n_x", "next_node" = "node")) %>%
                                                            tidyr::drop_na(n_x, node, next_node, n_next_x, ymax_end, ymin_end, xmax_end, xmin_end) %>%
                                                            dplyr::mutate(r = dplyr::row_number()) %>%
                                                            dplyr::arrange(n_x, -r) %>%
                                                            dplyr::select(-r) %>%
                                                            dplyr::group_by(n_x, node) %>%
                                                            dplyr::mutate(cum_flow_freq = cumsum(flow_freq) - flow_freq) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::group_by(n_x, n_next_x, node, next_node) %>%
                                                            dplyr::mutate(flow_start_ymax = ymax - cum_flow_freq,
                                                                          flow_start_ymin = flow_start_ymax - flow_freq)
                                                          
                                                          flows <- flows %>%
                                                            dplyr::arrange(n_x, n_next_x, next_node) %>%
                                                            dplyr::group_by(n_next_x, next_node) %>%
                                                            dplyr::mutate(cum_flow_freq_end = cumsum(flow_freq_end) - flow_freq_end) %>%
                                                            dplyr::mutate(flow_end_ymax = ymax_end - cum_flow_freq_end,
                                                                          flow_end_ymin = flow_end_ymax - flow_freq_end) %>%
                                                            dplyr::ungroup()
                                                          
                                                          flows <- flows %>%
                                                            dplyr::select(-n_x, -node, -freq, -ymax, -ymin, -xmin, -n_next_x, -next_node, -flow_freq, -ymin_end, -ymax_end, -xmax_end, -cum_flow_freq, -cum_flow_freq_end) %>%
                                                            dplyr::mutate(group = dplyr::row_number())
                                                          
                                                          flows %>%
                                                            rowwise() %>%
                                                            dplyr::mutate(..groupqq = stringr::str_remove(nodes, as.character(x))) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::group_by(..groupqq) %>%
                                                            dplyr::mutate(group = dplyr::cur_group_id()) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::select(-..groupqq) %>%
                                                            dplyr::mutate(smooth = params$smooth) %>%
                                                            as.data.frame()
                                                        })
                                       },
                                       
                                       compute_group = function(data, scales) {
                                         
                                         out1 <- purrr::map_dfr(1:nrow(data), ~{
                                           datat <- data %>% dplyr::slice(.x)
                                           sigmoid(datat$xmax, datat$xmin_end, datat$flow_start_ymax, datat$flow_end_ymax,
                                                   smooth = datat$smooth)
                                         }) %>%
                                           dplyr::arrange(x)
                                         out2 <- purrr::map_dfr(1:nrow(data), ~{
                                           datat <- data %>% dplyr::slice(.x)
                                           sigmoid(datat$xmin_end, datat$xmax, datat$flow_end_ymin, datat$flow_start_ymin,
                                                   smooth = datat$smooth)
                                         }) %>%
                                           dplyr::arrange(-x)
                                         
                                         dplyr::bind_rows(out1, out2)
                                       }
)

# TEXT LAYER -------
StatSankeyText <- ggplot2::ggproto("StatSankeyText", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(x = n_x,
                                                                      y = ymin + (ymax - ymin)/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      
                                                      return(as.data.frame(data))
                                                    })
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


# NODE LAYER -------
StatSankeyNode <- ggplot2::ggproto("StatSankeyNode", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      return(as.data.frame(data))
                                                    })
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


sankey_p <- function(mapping = NULL,
                     data = NULL,
                     position = "identity",
                     na.rm = FALSE,
                     show.legend = NA,
                     space = NULL,
                     type = "sankey",
                     width = .1,
                     smooth = 8,
                     inherit.aes = TRUE,
                     ...
) {
  params_list <- prepare_params(...)
  
  list(
    flow = ggplot2::layer(
      stat = StatSankeyFlow,
      data = data,
      mapping = mapping,
      geom = "polygon",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[1]]
        )
      )
    ),
    
    node = ggplot2::layer(
      stat = StatSankeyNode,
      data = data,
      mapping = mapping,
      geom = ggplot2::GeomRect,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[2]]
        )
      )
    )
  )
  
  
}


sankey_p_label <- function(mapping = NULL,
                           data = NULL,
                           position = "identity",
                           na.rm = FALSE,
                           show.legend = NA,
                           space = NULL,
                           type = "sankey",
                           width = .1,
                           inherit.aes = TRUE,
                           ...) {
  # Prepare aesthics for label
  label.aes <- list(...)
  
  list(
    label = ggplot2::layer(
      stat = StatSankeyText,
      data = data,
      mapping = mapping,
      geom = "label",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          type = type,
          label.aes
        )
      )
    )
  )
}

#-------------heterogeneity---------------#
#' @title cvh1_ml
#' @description CVH1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH1. TODO - we need to cite original CVH1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH1_eng_1 <- cvh1_ml(english_MA, boot = 10)
#' CVH1_eng_2 <- cvh1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH1_fish_1 <- cvh1_ml(model, boot = 10)
#' CVH1_fish_2 <- cvh1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH1_lim_1 <- cvh1_ml(lim_MR, boot = 10)
#' CVH1_lim_2 <- cvh1_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh1_ml <- function(model,
                    boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh1_ml cannot take models with heterogeneous variance.")}
  
  CVH1s <- ml_cvh1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CVH1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH1 <- ml_cvh1(tmp)
      })
    } else{
      CVH1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH1 <- ml_cvh1(tmp)
        return(CVH1) })
    }
    # Summarise the bootstrapped distribution.
    CVH1s_each_95 <- data.frame(t(apply(CVH1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH1s <-  round(CVH1s_each_95, digits = 3)
    colnames(CVH1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVH1s)
}


#' @title ml_cvh1
#' @description Calculated CVH1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh1 <- function(model){
  
  # total cvh1
  CVH1_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
  # cvh1 at different levels
  CVH1_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
  names(CVH1_each) <- paste0("CVH1_", model$s.names)
  names(CVH1_total) <- "CVH1_Total"
  
  CVH1s <- c(CVH1_total, CVH1_each)
  
  return(CVH1s)
}

#' @title cvh2_ml
#' @description CVH2 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH2. TODO - we need to cite original CVH2 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH2_eng_1 <- cvh2_ml(english_MA, boot = 10)
#' CVH2_eng_2 <- cvh2_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH2_fish_1 <- cvh2_ml(model, boot = 10)
#' CVH2_fish_2 <- cvh2_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH2_lim_1 <- cvh2_ml(lim_MR, boot = 10)
#' CVH2_lim_2 <- cvh2_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh2_ml <- function(model,
                    boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh2_ml cannot take models with heterogeneous variance.")}
  
  CVH2s <- ml_cvh2(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CVH2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH2 <- ml_cvh2(tmp)
      })
    } else{
      CVH2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH2 <- ml_cvh2(tmp)
        return(CVH2) })
    }
    # Summarise the bootstrapped distribution.
    CVH2s_each_95 <- data.frame(t(apply(CVH2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH2s <-  round(CVH2s_each_95, digits = 3)
    colnames(CVH2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVH2s)
}


#' @title ml_cvh2
#' @description Calculated CVH2 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh2 <- function(model){
  
  # total cvh2
  CVH2_total <- sum(model$sigma2) / (model$beta[[1]])^2
  # cvh2 at different levels
  CVH2_each <-  model$sigma2 / (model$beta[[1]])^2
  names(CVH2_each) <- paste0("CVH2_", model$s.names)
  names(CVH2_total) <- "CVH2_Total"
  
  CVH2s <- c(CVH2_total, CVH2_each)
  
  return(CVH2s)
}

#' @title m1_ml
#' @description M1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M1 - TODO - we need to cite original M1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M1. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M1_eng_1 <- m1_ml(english_MA, boot = 10)
#' M1_eng_2 <- m1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M1_fish_1 <- m1_ml(model, boot = 10)
#' M1_fish_2 <- m1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M1_lim_1 <- m1_ml(lim_MR, boot = 10)
#' M1_lim_2 <- m1_ml(lim_MR)
#' }
#' @references TODO
#' @export

m1_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m1_ml cannot take models with heterogeneous variance.")}
  
  M1s <- ml_m1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M1 <- ml_m1(tmp)
      })
    } else{
      M1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M1 <- ml_m1(tmp)
        return(M1) })
    }
    # Summarise the bootstrapped distribution.
    M1s_each_95 <- data.frame(t(apply(M1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M1s <-  round(M1s_each_95, digits = 3)
    colnames(M1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M1s)
}


#' @title ml_m1
#' @description Calculated M1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m1 <- function(model){
  
  # total m1
  M1_total <- sqrt(sum(model$sigma2)) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  # m1 at different levels
  M1_each <-  sqrt(model$sigma2) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  names(M1_each) <- paste0("M1_", model$s.names)
  names(M1_total) <- "M1_Total"
  
  M1s <- c(M1_total, M1_each)
  
  return(M1s)
}
#' @title m2_ml
#' @description M2 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M2 - TODO - we need to cite original M2 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M2_eng_1 <- m2_ml(english_MA, boot = 10)
#' M2_eng_2 <- m2_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M2_fish_1 <- m2_ml(model, boot = 10)
#' M2_fish_2 <- m2_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M2_lim_1 <- m2_ml(lim_MR, boot = 10)
#' M2_lim_2 <- m2_ml(lim_MR)
#' }
#' @references TODO
#' @export

m2_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m2_ml cannot take models with heterogeneous variance.")}
  
  M2s <- ml_m2(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M2 <- ml_m2(tmp)
      })
    } else{
      M2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M2 <- ml_m2(tmp)
        return(M2) })
    }
    # Summarise the bootstrapped distribution.
    M2s_each_95 <- data.frame(t(apply(M2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M2s <-  round(M2s_each_95, digits = 3)
    colnames(M2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M2s)
}


#' @title ml_m2
#' @description Calculated CV for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m2 <- function(model){
  
  # total m2
  M2_total <- sum(model$sigma2) / ( (model$beta[[1]])^2 + sum(model$sigma2) )
  # m2 at different levels
  M2_each <-  model$sigma2 / ( (model$beta[[1]])^2 + sum(model$sigma2) )
  names(M2_each) <- paste0("M2_", model$s.names)
  names(M2_total) <- "M2_Total"
  
  M2s <- c(M2_total, M2_each)
  
  return(M2s)
}

#-----------visualize heterogeneity-----------
# Get the results from the model
pred_dist_data <- function(mod_sim){
  # Get sigmas
  sigma2 <- mod_sim$sigma2
  
  # Get sampling error for mean
  se2_mu <- mod_sim$se^2
  
  # Get mean estimate
  mu <- mod_sim$b[1]
  
  # Names of random effect levels
  names <- c(attr(mod_sim$s.names, "names"), 'total')
  
  # Sigmas for prediction interval at each level
  cum_sigma2 <- c(sigma2 + se2_mu, sum(sigma2, se2_mu))
  
  # Create the dataframe
  data <- data.frame(group =  names, 
                     mean = mu,
                     sd = sqrt(cum_sigma2))
  return(data)
}

#### Calculate the proportion of effects beyond some threshold
prop_beyond <- function(data, threshold = 0.2, tail = c("above", "below")) {
  tail = match.arg(tail)
  
  # Get the mean and sd
  mean <- data$mean
  sd <- data$sd
  
  # Get the proportion beyond the threshold
  if (tail == "above") {
    proportion <- pnorm(threshold, mean, sd, lower.tail = FALSE)
  } else if (tail == "below") {
    proportion <- pnorm(threshold, mean, sd, lower.tail = TRUE)
  } else {
    stop("tail must be either 'above' or 'below'")
  }
  
  # Return the proportion
  return(round(proportion*100, 2))
}

# Plotting prediction distributions
pred_distri <- function(data) {
  mapply(function(mean, sd, group) {
    # new one
    list(stat_function(fun = dnorm, geom = "area", aes(color = group, fill = group), alpha = 0.1, args = list(mean = mean, sd = sd)), 
         stat_function(fun = dnorm, geom = "area", aes(color = group, fill = group), alpha = 0.1, args = list(mean = mean, sd = sd)))
    
    
  }, 
  mean = data$mean, 
  sd = data$sd, 
  group = data$group
  )
}

# plotting shaded prediction distribution
pred_distri_shaded <- function(x, m, sd, threshold) {
  y = dnorm(x, mean = m, sd = sd)
  y[x < threshold] <- NA
  return(y)
}


#--------------------------t distribution--------------------------#
library(extraDistr)

#### Calculate the proportion of effects beyond some threshold
propT_beyond <- function(data, df, threshold = 0.2, tail = c("above", "below")) {
  tail = match.arg(tail)
  
  # Get the mean and sd
  m <- data$mean
  sd <- data$sd
  
  # Get the proportion beyond the threshold
  if (tail == "above") {
    proportion <- extraDistr::plst(threshold, df = df, mu = m, sigma = sd, lower.tail = FALSE)
  } else if (tail == "below") {
    proportion <- extraDistr::plst(threshold, df= df, mu = m, sigma = sd, lower.tail = TRUE)
  } else {
    stop("tail must be either 'above' or 'below'")
  }
  
  # Return the proportion
  return(proportion*100)
}

# Plotting prediction distributions
pred_Tdistri <- function(data) {
  
  pred_Tdistri <- function(x, m, sd, df) {
    y = extraDistr::dlst(x, df = df, mu = m, sigma = sd)
    return(y)}
  
  mapply(function(df, mean, sd, group) {
    list(stat_function(fun = pred_Tdistri, geom = "area", aes(color = group, fill = group), alpha = 0.1, args = list(df = df, m = mean, sd = sd)), 
         stat_function(fun = pred_Tdistri, geom = "area", aes(color = group, fill = group), alpha = 0.1, args = list(df = df, m = mean, sd = sd)))
    
    
  }, 
  df = df,
  mean = data$mean, 
  sd = data$sd, 
  group = data$group
  )
}


# plotting prediction t distribution distribution
pred_Tdistri <- function(x, m, sd, df) {
  y = extraDistr::dlst(x, df = df, mu = m, sigma = sd)
  return(y)
}


# plotting shaded prediction distribution
pred_Tdistri_shaded <- function(x, m, sd, df, threshold) {
  y = extraDistr::dlst(x, df = df, mu = m, sigma = sd)
  y[x < threshold] <- NA
  return(y)
}

#--------extraction----------#
est <- function(mod) {
  df <- data.frame(est = round(mod$beta,3),
                   lb = round(mod$ci.lb,3),
                   ub = round(mod$ci.ub,3),
                   t = round(mod$zval,3),
                   df = mod$ddf,
                   p = round(mod$pval,4))
  df$levels <- rownames(df)
  return(df[c(7,1:6)])
}
