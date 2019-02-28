###############################################################################
##############        Spotted Hyena eRRBS DNA Methylation:       ##############
##############             Maternal care (FAS data)              ##############
##############                 By: Zach Laubach                  ##############
##############               created: 14 Feb 2019                ##############
##############             last updated: 28 Feb 2019             ##############
###############################################################################


  ### PURPOSE: This code is desingned to analyze ERRBS data and associations of 
  # maternal care from FAS data and genome wide DNA methylation in hyenas.
  
  
    # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Data management
    # 4: Univariate analyses
    # 5: Model 'close proximity' behaviors
    # 6: Model 'nursing' behaviors
    # 7: Model 'grooming' behaviors
    # 8: Format variables for MACUA 
    # 9: Export data files 



###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 clear global environment
    rm(list = ls())
  
  
  ### 1.2 Install and load packages 
    ## a) Data Manipulation and Descriptive Stats Packages
      # Check for tidyverse and install if not already installed
        if (!'tidyverse' %in% installed.packages()[,1]){
          install.packages ('tidyYverse')
        }
      # load tidyverse packages
        library ('tidyverse')
      
      # Check for sqldf and install if not already installed
        if (!'sqldf' %in% installed.packages()[,1]){
          install.packages ('sqldf')
        }
        options(gsubfn.engine = "R") #fixes tcltk bug; run before require sqldf
      # load tidyverse packages
        library ('sqldf')
      
      # Check for lubridate and install if not already installed
        if (!'lubridate' %in% installed.packages()[,1]){
          install.packages ('lubridate')
        }
      # load lubridate packages
        library ('lubridate') 
        
      # NOTE: used to make joins based on conditional statements
      # Check for fuzzyjoin and install if not already installed
#        if (!'fuzzyjoin' %in% installed.packages()[,1]){
#          install.packages ('fuzzyjoin')
#        }
      # load fuzzyjoin packages
#        library ('fuzzyjoin')
        
      # Check for here and install if not already installed
        if (!'here' %in% installed.packages()[,1]){
          install.packages ('here')
        }
      # load here packages
        library ('here')
        
    ## b) Modeling Packages
      # Check for broom and install if not already installed
        if (!'broom' %in% installed.packages()[,1]){
          install.packages ('broom')
        }
      # load broom packages
        library ('broom')
        
      # Check for glmmTMB and install if not already installed
        # for zero inlfated mixed models
        if (!'glmmTMB' %in% installed.packages()[,1]){
          install.packages ('glmmTMB')
        }
      # load glmmTMB packages
        library ('glmmTMB')
        
      # Check for bbmle and install if not already installed
        # for AICtab
        if (!'bbmle' %in% installed.packages()[,1]){
          install.packages ('bbmle')
        }
      # load bbmle packages
        library ('bbmle')

        
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 3.5.2 (2018-12-20)
    # Platform: x86_64-apple-darwin15.6.0 (64-bit)
    # Running under: macOS  Mojave 10.14.3
    
  
  ### 1.4 Set working directory 
    setwd(here())

    
  ### 1.5 Set file paths for data importing and exporting
    ## a) The path to errbs data
      rrbs_vars_data_path <- paste("~/R/R_wd/fisi/project/6_hy_RRBS/",
                              "data/", sep = '')
    
    ## b) The path to maternal care FAS data
      mat_fas_path <- paste("~/R/R_wd/fisi/project/fas_maternal_care/",
                                   "fas/output/", sep = '')    
      
    ## c) The path to Access Fisi backend data
      general_data_path <- paste("~/R/R_wd/fisi/project/0_data/",
                              sep = '')
    
    ## d) Source scripts path
      source_path <- paste("~/Git/source_code/")
      
  
  ### 1.6 Source functions
    ## a) all_char_to_lower function
      source(file = paste0(source_path, "all_char_to_lower.R"))
      
    ## b) format_var_names function
      source(file = paste0(source_path, "format_var_names.R"))
      
    ## c) format_var_names_dash function
      source(file = paste0(source_path, "format_var_names_dash.R"))  
      
    ## d) fix_dates_and_times function
    #NOTE: Not working
    #     source(file = paste0(source_path, "fix_dates_and_times.R")) 
      

      
###############################################################################
##############                  2. Import data                   ##############
###############################################################################    
      
  ### 2.1 Import data files (with readr)
      
    ## a) Import data from copy of sample seleciton file, which includes list 
      # of hyenas and variables
      rrbs_vars <- read_csv(paste(rrbs_vars_data_path,
                                  "hy_rrbs_variables.csv", sep = ''))
        
    ## b) check that all variables are of appropriate class    
      sapply(rrbs_vars, class)
      
    ## c) Import maternal care FAS data
      mat_fas <- read_csv(paste(mat_fas_path,
                                 "fas_counts.csv", sep = ''))  
      
    ## d) Import tblHyena data files
      tblHyena <- read_csv(paste(general_data_path, "1_output_tidy_tbls/",
                                  "tblHyenas.csv", sep = ''))
      
    ## e)  Import tblFemalerank data file
      tblFemalerank <- read_csv(paste(general_data_path, "1_output_tidy_tbls/",
                                      "tblFemalerank.csv", sep = ''))
      
    ## f)  Import tblReprostates
      repro_state <- read_csv(paste(general_data_path, 
                                      "reprostates_EDS.csv", sep = ''))
      
      

###############################################################################
##############                 3. Data management                ##############
###############################################################################

  ### 3.1 Tidy rrbs_vars  
    ## a) view the variable classes to determine which need to be modified
      spec(rrbs_vars)
      
    ## b) Convert all text to lower case
      rrbs_vars <- AllCharactersToLower(rrbs_vars)
      
    ## c) Format variable names (lowercase and separated by '.')
      rrbs_vars <- FormatVarNames(rrbs_vars)
      
    ## d) Convert dates stored as character (e.g. 8/3/05) to formatted dates
      rrbs_vars$darting.date <- as.Date(rrbs_vars$darting.date, 
                                        format = "%m/%d/%y")
      
      rrbs_vars$cub.dob <- as.Date(rrbs_vars$cub.dob, 
                                        format = "%m/%d/%y")
      
    ## e) Create an estimated age in months by subtracting birthdate from
      # darting date using lubridate and dividing by average # days in month 
      rrbs_vars <- rrbs_vars %>%
        mutate(dart.age.mon = round((interval(cub.dob, 
                                        darting.date) %/% days(1) / 30.44), 1))
   
    ## f) drop unecessary variables
      rrbs_vars <- rrbs_vars %>%
        select (- c(orig.box, orig.cell, sample.type, check.out.state,
                    check.out.vol.u.l., dna.conc.ng.u.l., initial.library.prep,
                    mat.care.sessions, cub.sex, age.at.darting, grooming..blup,
                    blup.tertile)) 
      
    ## g) convert hy.id to character class 
      rrbs_vars$hy.id <- as.character(rrbs_vars$hy.id)
      
      
  ### 3.2 Tidy mat_fas
    ## a) view the variable classes to determine which need to be modified
      spec(mat_fas)
      
    ## b) Convert all text to lower case
      mat_fas <- AllCharactersToLower(mat_fas)
      
    ## c) Format variable names (lowercase and separated by '.')
      mat_fas <- FormatVarNames(mat_fas)
      
    ## d) rename variables
      mat_fas <-  mat_fas %>%
        rename('hy.id' = 'cub')
      
    ## e) convert hy.id to character class  
      mat_fas$hy.id <- as.character(mat_fas$hy.id)  
      
    ## f) rename variables
      mat_fas <-  mat_fas %>%
        rename('fas.date' = 'date')
      
  
  ### 3.3 Tidy tblHyena
    ## a) view the variable classes to determine which need to be modified
      spec(tblHyena)
    
    ## b) Convert all text to lower case
      tblHyena <- AllCharactersToLower(tblHyena)
    
    ## c) Format variable names (lowercase and separated by '.')
      tblHyena <- FormatVarNames(tblHyena)
    
    ## d) Convert dates stored as character (e.g. 03-aug-05) to formatted dates
     tblHyena$first.seen <-  as.Date(tblHyena$first.seen, 
                                       format = "%d-%b-%y")
     tblHyena$den.grad <-  as.Date(tblHyena$den.grad, 
                                     format = "%d-%b-%y")
     tblHyena$disappeared <-  as.Date(tblHyena$disappeared, 
                                        format = "%d-%b-%y")
     tblHyena$birthdate <- as.Date(tblHyena$birthdate, 
                                     format = "%d-%b-%y")
     tblHyena$death.date <-  as.Date(tblHyena$death.date, 
                                       format = "%d-%b-%y")
     tblHyena$weaned <-  as.Date(tblHyena$weaned, 
                                   format = "%d-%b-%y")
    
    ## e) rename variables
       tblHyena <- tblHyena %>%
          rename('hy.id' = 'id')
    
    ## f) convert hy.id to character class  
       tblHyena$hy.id <- as.character(tblHyena$hy.id)
      
 
  ### 3.4 Tidy tblFemalerank
    ## a) Pattern recognize numbers from Year variable and copy; gets rid of
      # unwanted text characters
      tblFemalerank$rank.year <- as.numeric(regmatches
                                            (tblFemalerank$year,
                                              gregexpr("[[:digit:]]+",
                                                       tblFemalerank$year)))
           
    ## b) rename 'id' variable as 'mom' 
       tblFemalerank <- rename_(tblFemalerank, "mom" = "id")
       
       
    ### 3.5 Tidy repro_state       
       ## a) view the variable classes to determine which need to be modified
       spec(repro_state)
       
       ## b) Convert all text to lower case
       repro_state <- AllCharactersToLower(repro_state)
       
       ## c) Format variable names (lowercase and separated by '.')
       repro_state <- FormatVarNames(repro_state)
       
       ## d) rename variables
       repro_state <-  repro_state %>%
         rename('hy.id' = 'mom')
       
       ## e) convert hy.id to character class  
       repro_state$hy.id <- as.character( repro_state$hy.id)     
       
             
  ### 3.6 Join Data sets
    ## a) Semi join mat_fas by rrbs_vars, retains all rows in 
      # mat_fas that has a match id in rrbs_vars   
      fas_data <- semi_join(mat_fas,
                              rrbs_vars, by = "hy.id")
      
    ## b) A list of hyenas having both ERRBS and FAS data  
      hys <- unique(fas_data$hy.id)
      # the number of those hyenas
      length(hys)
  
    ## c) Join variables from rrbs_vars to fas_data
      fas_data <- fas_data  %>%
        left_join(select(rrbs_vars, c(hy.id, darting.date, cub.dob, 
                                      dart.age.mon)), 
                  by = c("hy.id" = "hy.id"))
      
    ## c) Join variables from tblHyenas to fas_data
      fas_data <- fas_data  %>%
        left_join(select(tblHyena, c(hy.id, number.littermates, litrank)),
                  by = c("hy.id" = "hy.id")) 
      
      
  ### 3.7 Tidy fas_data
    ## a) Create a varialbe, 'fas.vs.dart,' which indicates whether fas data 
      # were collected before vs after the darting date
      fas_data <- fas_data  %>%
        mutate(fas.vs.dart = case_when(fas_data$fas.date >= 
                                         fas_data$darting.date ~ c("after"),
                                       fas_data$fas.date < 
                                         fas_data$darting.date ~ 
                                         c("before")))
      
      
    ## b) Create a varialbe, 'fec.vs.dart.days,' which indicates how many days
      # between the fas.date and darting date
      fas_data <- fas_data  %>%
        mutate(fas.vs.dart.days = (fas_data$fas.date - 
                                     fas_data$darting.date)) 
     
    ## c) Create a varialbe, 'fas.age.mon' which is an estimated age in 
      # months by subtracting cub.dob from fas.date using lubridate
      fas_data <- fas_data %>%
        mutate(fas.age.mon = round((interval(cub.dob, 
                                         fas.date) %/% days(1)/ 30.44), 1))
    
    ## d) Extract month from fas.date
      # Use lubridate to extract the month during which a fas sample was 
      # collected and make a new variable  
      fas_data$fas.mon <- month(fas_data$fas.date)
      
    ## e) Create a varialbe, 'migratn.seas,' which indicates if a fas sample
      # was collected in migration (June 1 - Oct 31)
      fas_data <- fas_data  %>%
        mutate(migratn.seas = ifelse(fas.mon >= 6 & fas.mon<= 10, 'migration', 
                                     'none'))
      
    ## f) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of migratn.seas variable and sets the reference  
      # level to 'none'
      fas_data <- transform(fas_data,
                            migratn.seas = factor(migratn.seas,
                                              levels = c("none", 
                                                         "migration")))  
      
    ## g) Extract am vs. pm from stop.fas time
      # Use lubridate to extract the month during which a fas sample was 
      # collected and make a new variable
      fas_data <- fas_data  %>%
        mutate(fas.am.pm = ifelse(lubridate::am(stop.fas), 'am', 
                                     'pm'))
      
    ## h) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of fas.am.pm variable and sets the reference  
      # level to 'am' 
      fas_data <- fas_data %>%
        transform(fas.am.pm = factor(fas.am.pm,
                                     levels = c("am","pm"))) 
      
    ## i) Subset fas_data based on one 12 month (including mon. 12) cut-off  
      # for FAS samples collected prior to darting.data
      fas_data <- fas_data  %>%
        filter(fas.vs.dart == 'before' & fas.age.mon <= 13)
    #*** NOTE *** This is a data inclusion cut-off decision. 
          # ALL FAS ocurr before darting date (ERRBS sample) and when 
          # cub ages were less than 13 months 
      
    ## j) Subset fas_data based on minimum overlap time of 5 min 
      # for FAS samples when both mom and cub present together
      fas_data <- fas_data  %>%
        filter(overlap >= 5)
    #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL FAS in which there is a minimum overlap time of 5 min 
      # for FAS samples when both mom and cub present together 
      
    ## k) Litter size tidying
      fas_data <- fas_data  %>%
        replace_na(list(number.littermates = 0))
        
    #*** MANUAL DATA CHECK *** there is no lit size info for chco, reb, goof,
      # or maa in tblHyena (27 feb 2019). These animals do not share a 
      # dob with any other animals in tblHyenas, so they are manually 
      # assigne number.litter mates == 0
      
    ## l) Make a new 2 level nomial variable lit.size
      fas_data <- fas_data  %>%
        mutate(lit.size = case_when(number.littermates == 0 ~ c("single"),
                                    number.littermates == 1 ~ c("twin")))
                                    
      
      
  ### 3.8 Combine repro_state w fas_data
    ## a) Make an empty data frame
      # This is an empty data frame that can store the overlaping fas_data
      # and repro_states data
      fas_repro_data  <- c()   
     
    ## b) For loop to find overlap
      # Iterate over fecal data (hy.id and fas.date), to find interseciton
      # with repro_state data
      for (i in 1:nrow(fas_data)) { 
        
        # loop through 1:n IDs in fas_data
        mom = paste (fas_data$mom[i])  
        
        # loop through 1:n dates in fas_data
        fas.date <-(fas_data$fas.date[i])
      
        # create a dataframe to store the rows from repro_states where the
        # mom matches hy.id, and fas.date is in between cycle start and stop 
        overlap_fas_repro <- filter(repro_state, mom == hy.id & 
                                       fas.date >= cycle.start & 
                                       fas.date <= cycle.stop)
        
        # Control flow
          # if there is no id match and date overlap, 
          # then go to next loop iteration in fas_data
          if (nrow(overlap_fas_repro) < 1) {
            next
          }
      
        # add the fas.date onto the overlap_fas_repro data
        overlap_fas_repro <- cbind(overlap_fas_repro, fas.date, mom)
        
        # add the filtered overlap_poop_repro data to a new dataframe
        # over each iteration of the loop
        fas_repro_data <- rbind(fas_repro_data, 
                                  overlap_fas_repro)
      }
      
    ## c) Join repro state data to FAS data 
      fas_data <- fas_data %>%
        left_join(select(fas_repro_data, c(mom, fas.date, state, 
                                             cycle.start, cycle.stop, 
                                             trimester, parity)),
                  by = c("mom" = "mom",
                         "fas.date" = "fas.date")) 
      
      
  ### 3.9 Re-tidy FAS data
    ## a) Change state to character
      fas_data$state <- as.character(fas_data$state)
      
    ## b) Subset data to only include moms and cubs when lactating
      fas_data <- fas_data  %>%
        filter(state == "l")
      
    #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL FAS in which mom's are lactating 
      
    ## c) Re-code *nominal* factor (with ordered levels)  
      fas_data <- transform(fas_data, 
                             state = factor(state))
  
    ## d) Re-code parity as a two level variable
      # primiparous = prim v.s. multiparous = mult
      fas_data <- fas_data  %>%
        mutate(parity.binary = ifelse (fas_data$parity == 'a', 'prim', 'mult'))
       
    ## e) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of parity variable and sets the reference level 
      # to 'a'
      fas_data <- transform(fas_data, 
                             parity.binary = factor(parity.binary,
                                              levels = c("prim", "mult")))
    # *** NOTE *** levels may vary by data set
      
      

###############################################################################
##############              4. Univariate analyses               ##############
###############################################################################      
      
  ### 4.1 Overview
    # Generate summary / descriptive stats for variables in the data frame
      summary(fas_data)
      
      
  ### 4.2 Visualize (and transform as needed) raw data        
      ## a) Histogram Outcome (nos minutes in close proximity)
      ggplot(data=fas_data, aes(x=c)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 35, by = 2.5), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 35)) +
        labs(title= "Histogram for mother and cub close proximity") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Close proximity (num. mins.)", y="Frequency") 
      
      ## b) Histogram Outcome (ratio in close proximity : overlap time together)
      ggplot(data=fas_data, aes(x=(c/overlap))) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1, by = 0.01), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 1)) +
        labs(title= "Histogram for mother and cub close proximity
             vs overlap time together") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Ratio (close proximity : overlap", y="Frequency") 
  
  ### 4.3 Visualize (and transform as needed) raw data        
    ## a) Histogram Outcome (num minutes spent nursing)
      ggplot(data=fas_data, aes(x=n)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 40, by = 2.5), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 40)) +
        labs(title= "Histogram for amount of minutes spent nursing") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Nursing (num. mins.)", y="Frequency") 
      
    ## b) Histogram Outcome (ratio nursing : overlap time together)
      ggplot(data=fas_data, aes(x=(n/overlap))) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1, by = 0.01), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 1)) +
        labs(title= "Histogram for time spent nursing
             vs overlap time together") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Ratio (nursing : overlap", y="Frequency") 
      
      
    ### 4.4 Visualize (and transform as needed) raw data        
      ## a) Histogram Outcome (num minutes spent grooming)
      ggplot(data=fas_data, aes(x=g)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 15, by = 1), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 15)) +
        labs(title= "Histogram for amount minutes spent grooming") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Grooming (num. mins.)", y="Frequency")   
    
      
    ## b) Histogram Outcome (ratio grooming : overlap time together)
      ggplot(data=fas_data, aes(x=(g/overlap))) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1, by = 0.01), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 1)) +
        labs(title= "Histogram for time spent grooming
             vs overlap time together") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Ratio (grooming : overlap", y="Frequency") 
      
  
      
###############################################################################
##############       5. Model 'close proximity' behaviors        ##############
###############################################################################  
      
  ### 5.1 Close proximity model parameterization 
    ## a) Close proximity zero inflated, poisson distributed response
      close.zipoisson <- glmmTMB(c ~ fas.age.mon + fas.am.pm +  # sampling var
                                   #migratn.seas +              # mat care var
                                   #lit.size + parity.binary +  # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = fas_data,
                                 ziformula = ~1,
                                 family = list(family = 'poisson', 
                                               link = 'log'))
      # Model summary estimates
      summary(close.zipoisson)
    
    ## b) Close proximity zero inflated, negative bionomial (NB2) 
      # distributed response
      close.zinegbinom2 <- update(close.zipoisson, 
                                 family = list(family = 'nbinom2',
                                               link = 'log'))
      
      # Model summary estimates
      summary(close.zinegbinom2)
      
    ## c) Close proximity zero inflated, negative bionomial (NB1) 
      # distributed response
      close.zinegbinom1 <- update(close.zipoisson, 
                                 family = list(family = 'nbinom1',
                                               link = 'log'))
      
      # Model summary estimates
      summary(close.zinegbinom1)
      
    ## d) Close proximity poisson distributed response; NO zero inflation
      close.poisson <- glmmTMB(c ~ fas.age.mon + fas.am.pm +  # sampling var
                                 #migratn.seas +              # mat care var
                                 #lit.size + parity.binary +  # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = fas_data,
                                 ziformula = ~0,
                                 family = list(family = 'poisson', 
                                               link = 'log'))
      # Model summary estimates
      summary(close.poisson)
      
    ## e) Close proximity, negative bionomial (NB2); NO zero inflation 
      # distributed response
      close.negbinom2 <- update(close.poisson, 
                                  family = list(family = 'nbinom2',
                                                link = 'log'))
      
      # Model summary estimates
      summary(close.negbinom2)
      
    ## f) Close proximity negative bionomial (NB1); NO zero inflation 
      # distributed response
      close.negbinom1 <- update(close.poisson, 
                                family = list(family = 'nbinom1',
                                              link = 'log'))
      
      # Model summary estimates
      summary(close.negbinom1)
      
      
  ### 5.2 Close proximity model fit comparisons
    ## e) Compare model fit with AICtab (from bbmle)
      AICtab(close.zipoisson, close.zinegbinom2, close.zinegbinom1, 
             close.poisson, close.negbinom2, close.negbinom1)
    
      
  ### 5.3 Close proximity BLUP extractions

    # NOTE: Use when there are repeated measuresments for a variable
       # that is to be used as an explanatory variable in another analysis.
       # Can control for other variables that bias estimates of explanatory
       # variable
          
    # NOTE: BLUPs are conditional modes from a generalized linear model
      # (according to Doug Bates). Here we use the glmmTMB package 
      # to fit a zero-inflated, poisson/negative bionomial model
      # BLUP = fixef(intrcpt) + ranef 
       
      # Calucate the BLUPs, individual variation in close proximity between a
      # mom and cub from the best fitting model (above)
      
    ## a) Generate best fit model summary 
      ranef(close.zinegbinom1) # random effect
      fixef(close.zinegbinom1) # fixed effect
      coef(close.zinegbinom1) # fixed effect
      
    ## b) extract BLUPs from mixed model object
      close.blups <- as.data.frame(ranef(close.zinegbinom1)) # extract ranef as 
      # a dataframe, BLUPs = rand effects + intercept (from poiss/neg. binom)
      
    ## c) Rename variables in blups table
      close.blups <- close.blups %>%
        rename('hy.id' = 'grp') %>%
        rename('c.ranef' = 'condval') %>%
        rename('c.ranef.sd' = 'condsd') %>%
        select(c('hy.id', 'c.ranef', 'c.ranef.sd'))
      
    ## d) extract fixed effect (intercept) from poisson/neg. binomial model
      close.intrcpt <- (fixef(close.zinegbinom1)[[1]])[[1]] # fixed effect
    #**** INTERPETATION ****#  
      # Interpet as expected log count of behavior when controlling for /
      # holding constant the effects ofcovariates (or if exponentiated 
      # the intercept is the incident rate of the expected count 
      # of behavior) as a proporiton of time overlap (the offset)
      
    ## e) extract fixed effect (intercept) from zero inflation model
      zi.close.intrcpt <- (fixef(close.zinegbinom1)[[2]])[[1]] # fixed effect
      
    ## f) Create a new variable that is ranef plus both poisson/neg. binomial
      # model intercept. This provides estimates of individual level
      # variation in proportion of time spent doing a behavior vs time mom
      # and cub overlapped (present together)
      close.blups <-  close.blups  %>%
        mutate(c.blups = close.intrcpt + c.ranef) %>%
        mutate(c.expontd.blups = exp(c.blups))
    #**** INTERPETATION ****#  
      # Interpet as the conditional mode or the log counts / incident rate
      # of expected counts as a proportion of the overlap time (offset) for
      # each individual...while holding constant effect of other covariates
      
      
      
###############################################################################
##############            6. Model 'nursing' behaviors           ##############
###############################################################################
      
  ### 6.1 Nursing model parameterization 
    ## a) Nursing zero inflated, poisson distributed response
      nurse.zipoisson <- glmmTMB(n ~ fas.age.mon + fas.am.pm +  # sampling var
                                   #migratn.seas +              # mat care var
                                   #lit.size + parity.binary +  # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = fas_data,
                                 ziformula = ~1,
                                 family = list(family = 'poisson', 
                                               link = 'log'))
    # Model summary estimates
      summary(nurse.zipoisson)
      
    ## b) Nursing zero inflated, negative bionomial (NB2) 
      # distributed response
      nurse.zinegbinom2 <- update(nurse.zipoisson, 
                                  family = list(family = 'nbinom2',
                                                link = 'log'))
      
    # Model summary estimates
      summary(nurse.zinegbinom2)
      
    ## c) Nursing zero inflated, negative bionomial (NB1) 
      # distributed response
      nurse.zinegbinom1 <- update(nurse.zipoisson, 
                                  family = list(family = 'nbinom1',
                                                link = 'log'))
      
    # Model summary estimates
      summary(nurse.zinegbinom1)
      
    ## d) Nursing poisson distributed response; NO zero inflation
      nurse.poisson <- glmmTMB(n ~ fas.age.mon + fas.am.pm +  # sampling var
                                 #migratn.seas +              # mat care var
                                 #lit.size + parity.binary +  # mat care var
                                 offset(log(overlap)) + (1|hy.id),
                               data = fas_data,
                               ziformula = ~0,
                               family = list(family = 'poisson', 
                                             link = 'log'))
      # Model summary estimates
      summary(nurse.poisson)
      
    ## e) Nursing, negative bionomial (NB2); NO zero inflation 
      # distributed response
      nurse.negbinom2 <- update(nurse.poisson, 
                                family = list(family = 'nbinom2',
                                              link = 'log'))
      
    # Model summary estimates
      summary(nurse.negbinom2)
      
    ## f) Nursing negative bionomial (NB1); NO zero inflation 
      # distributed response
      nurse.negbinom1 <- update(nurse.poisson, 
                                family = list(family = 'nbinom1',
                                              link = 'log'))
      
    # Model summary estimates
      summary(nurse.negbinom1)
      
      
  ### 6.2 Nursing model fit comparisons
    ## e) Compare model fit with AICtab (from bbmle)
      AICtab(nurse.zipoisson, nurse.zinegbinom2, nurse.zinegbinom1, 
             nurse.poisson, nurse.negbinom2, nurse.negbinom1)
      
      
  ### 6.3 Nursing BLUP extractions
      
    # NOTE: Use when there are repeated measuresments for a variable
    # that is to be used as an explanatory variable in another analysis.
    # Can control for other variables that bias estimates of explanatory
    # variable
      
    # NOTE: BLUPs are conditional modes from a generalized linear model
    # (according to Doug Bates). Here we use the glmmTMB package 
    # to fit a zero-inflated, poisson/negative bionomial model
    # BLUP = fixef(intrcpt) + ranef 
      
    # Calucate the BLUPs, individual variation in nursing between a
    # mom and cub from the best fitting model (above)
      
    ## a) Generate best fit model summary 
      ranef(nurse.zipoisson) # random effect
      fixef(nurse.zipoisson) # fixed effect
      coef(nurse.zipoisson) # BLUPs
      
    ## b) extract BLUPs from mixed model object
      nurse.blups <- as.data.frame(ranef(nurse.zipoisson)) # extract ranef as 
      # a dataframe, BLUPs = rand effects + intercept (from poiss/neg. binom)
      
    ## c) Rename variables in blups table
      nurse.blups <- nurse.blups %>%
        rename('hy.id' = 'grp') %>%
        rename('n.ranef' = 'condval') %>%
        rename('n.ranef.sd' = 'condsd') %>%
        select(c('hy.id', 'n.ranef', 'n.ranef.sd'))
      
    ## d) extract fixed effect (intercept) from poisson/neg. binomial model
      nurse.intrcpt <- (fixef(nurse.zipoisson)[[1]])[[1]] # fixed effect
    #**** INTERPETATION ****#  
      # Interpet as expected log count of behavior when controlling for /
      # holding constant the effects ofcovariates (or if exponentiated 
      # the intercept is the incident rate of the expected count 
      # of behavior) as a proporiton of time overlap (the offset)
      
    ## e) extract fixed effect (intercept) from zero inflation model
      zi.nurse.intrcpt <- (fixef(nurse.zipoisson)[[2]])[[1]] # fixed effect
      
    ## f) Create a new variable that is ranef plus both poisson/neg. binomial
      # model intercept. This provides estimates of individual level
      # variation in proportion of time spent doing a behavior vs time mom
      # and cub overlapped (present together)
      nurse.blups <-  nurse.blups  %>%
        mutate(n.blups = nurse.intrcpt + n.ranef) %>%
        mutate(n.expontd.blups = exp(n.blups))    
    #**** INTERPETATION ****#  
      # Interpet as the conditional mode or the log counts / incident rate
      # of expected counts as a proportion of the overlap time (offset) for
      # each individual...while holding constant effect of other covariates

      
      
###############################################################################
##############           7. Model 'grooming' behaviors           ##############
###############################################################################
      
  ### 7.1 Grooming model parameterization 
    ## a) Grooming zero inflated, poisson distributed response
      groom.zipoisson <- glmmTMB(g ~ fas.age.mon + fas.am.pm +  # sampling var
                                   #migratn.seas +              # mat care var
                                   #lit.size + parity.binary +  # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = fas_data,
                                 ziformula = ~1,
                                 family = list(family = 'poisson', 
                                               link = 'log'))
      # Model summary estimates
      summary(groom.zipoisson)
      
    ## b) Grooming zero inflated, negative bionomial (NB2) 
      # distributed response
      groom.zinegbinom2 <- update(groom.zipoisson, 
                                  family = list(family = 'nbinom2',
                                                link = 'log'))
      
      # Model summary estimates
      summary(groom.zinegbinom2)
      
    ## c) Grooming zero inflated, negative bionomial (NB1) 
      # distributed response
      groom.zinegbinom1 <- update(groom.zipoisson, 
                                  family = list(family = 'nbinom1',
                                                link = 'log'))
      
      # Model summary estimates
      summary(groom.zinegbinom1)
      
    ## d) Grooming poisson distributed response; NO zero inflation
      groom.poisson <- glmmTMB(g ~ fas.age.mon + fas.am.pm +  # sampling var
                                 #migratn.seas +              # mat care var
                                 #lit.size + parity.binary +  # mat care var
                                 offset(log(overlap)) + (1|hy.id),
                               data = fas_data,
                               ziformula = ~0,
                               family = list(family = 'poisson', 
                                             link = 'log'))
      # Model summary estimates
      summary(groom.poisson)
      
    ## e) Grooming, negative bionomial (NB2); NO zero inflation 
      # distributed response
      groom.negbinom2 <- update(groom.poisson, 
                                family = list(family = 'nbinom2',
                                              link = 'log'))
      
      # Model summary estimates
      summary(groom.negbinom2)
      
    ## f) Grooming negative bionomial (NB1); NO zero inflation 
      # distributed response
      groom.negbinom1 <- update(groom.poisson, 
                                family = list(family = 'nbinom1',
                                              link = 'log'))
      
      # Model summary estimates
      summary(groom.negbinom1)
      
      
  ### 7.2 Grooming model fit comparisons
    ## e) Compare model fit with AICtab (from bbmle)
      AICtab(groom.zipoisson, groom.zinegbinom2, groom.zinegbinom1, 
             groom.poisson, groom.negbinom2, groom.negbinom1)
      
      
  ### 7.3 Grooming BLUP extractions
      
      # NOTE: Use when there are repeated measuresments for a variable
      # that is to be used as an explanatory variable in another analysis.
      # Can control for other variables that bias estimates of explanatory
      # variable
      
      # NOTE: BLUPs are conditional modes from a generalized linear model
      # (according to Doug Bates). Here we use the glmmTMB package 
      # to fit a zero-inflated, poisson/negative bionomial model
      # BLUP = fixef(intrcpt) + ranef 
      
      # Calucate the BLUPs, individual variation in grooming between a
      # mom and cub from the best fitting model (above)
      
    ## a) Generate best fit model summary 
      ranef(groom.zinegbinom1) # random effect
      fixef(groom.zinegbinom1) # fixed effect; coef in glmmTMB gives fixef
      coef(groom.zinegbinom1) # BLUPs
      
    ## b) extract BLUPs from mixed model object
      groom.blups <- as.data.frame(ranef(groom.zinegbinom1)) # extract ranef as 
      # a dataframe, BLUPs = rand effects + intercept (from poiss/neg. binom)               
      
      
    ## c) Rename variables in blups table
      groom.blups <- groom.blups %>%
        rename('hy.id' = 'grp') %>%
        rename('g.ranef' = 'condval') %>%
        rename('g.ranef.sd' = 'condsd') %>%
        select(c('hy.id', 'g.ranef', 'g.ranef.sd'))
      
    ## d) extract fixed effect (intercept) from poisson/neg. binomial model
      groom.intrcpt <- (fixef(groom.zinegbinom1)[[1]])[[1]] # fixed effect
    #**** INTERPETATION ****#  
      # Interpet as expected log count of behavior when controlling for /
      # holding constant the effects ofcovariates (or if exponentiated 
      # the intercept is the incident rate of the expected count 
      # of behavior) as a proporiton of time overlap (the offset)
      
    ## e) extract fixed effect (intercept) from zero inflation model
      zi.groom.intrcpt <- (fixef(groom.zinegbinom1)[[2]])[[1]] # fixed effect
      
    ## f) Create a new variable that is ranef plus both poisson/neg. binomial
      # model intercept. This provides estimates of individual level
      # variation in proportion of time spent doing a behavior vs time mom
      # and cub overlapped (present together)
      groom.blups <-  groom.blups  %>%
        mutate(g.blups = groom.intrcpt + g.ranef) %>%
        mutate(g.expontd.blups = exp(g.blups))  
    #**** INTERPETATION ****#  
      # Interpet as the conditional mode or the log counts / incident rate
      # of expected counts as a proportion of the overlap time (offset) for
      # each individual...while holding constant effect of other covariates
      
      
      
###############################################################################
##############           8. Format variables for MACUA           ##############
###############################################################################      
  
  ### 8.1 Join the BLUPs for each behavior into a single data frame
    ## a) Left join close.blups to nurse.blups  
      fas_blups <- left_join(close.blups,
                             nurse.blups, by = 'hy.id')
      
    ## b) Left join close.blups to nurse.blups    
      fas_blups <- left_join(fas_blups,
                             groom.blups, by = 'hy.id')
      
      
  ### 8.2 Make dataframe for txt. file export to MACUA (for EWAS)
      ## a) Manual data clean up drop animals that failed ERRBS library QAQC
      rrbs_vars <- rrbs_vars %>%
        filter (!grepl('tato', hy.id))
    
      ## b) Left join fas_blups to rrbs_vars
      rrbs_vars <- rrbs_vars %>%
        left_join(fas_blups, by = 'hy.id')
    
        
  ### 8.3 Arrange data set by hy.id abc order
      ## This is necessary because the order of rows matches data to a 
      ## specific animal
      rrbs_vars <- arrange(rrbs_vars, hy.id)
    
        
  ### 8.4 Subset explanatory/predictor variables and format for use in MACAU
    ## a) Maternal rank formatted
     mat_rank_predictor_hy_n29 <- rrbs_vars %>%
        select(mom.rank) 
    
    ## b) Close proximity BLUPs formatted
      close_prox_predictor_hy_n29 <- rrbs_vars %>%
        select(c.expontd.blups) 
      
    ## c) Nursing BLUPs formatted
      nurse_predictor_hy_n29 <- rrbs_vars %>%
        select(n.expontd.blups) 
      
    ## d) Grooing BLUPs formatted
      groom_predictor_hy_n29 <- rrbs_vars %>%
        select(g.expontd.blups)   

            
  ### 8.5 Subset explanatory/predictor variables and format for use in MACAU
    ## a) Covariates (variables to control analyses)
      mat_rank_care_covars_hy_n29 <- rrbs_vars %>%
        mutate(intercept = 1) %>%
        select(intercept, dart.age.mon) 
       
 
           
###############################################################################
##############                9. Export data files               ##############
###############################################################################
      
  ### 9.1 Export FAS BLUPs to csv     
      # Save and export tables as a .cvs spreadsheet and named with today's
      # date. Files are saved in the 'output' folder in the working directory.
      
    ## a) Generate File Names
      # For each table that will be saved as a .csv file, first generate a 
      # file name to save each table
      # here, to paste the folder path, followed by the file name 
      csv.file.name.FAS.BLUPs <- paste0(here("data", 
                                                    "fas_blups.csv")) 
     
    ## b) Save Tables 
      # Save data frame as a .csv file (a spreadsheet/table) into the 
      # output data folder in the working directory.
      write.csv (fas_blups, file = csv.file.name.FAS.BLUPs)


  ### 9.2 Export explanatory/predictor variables to txt     
      # Save and export part of table as a .txt file. Files are formatted 
      # for use in MACAU as explanatory/predictor variables
            
    ## a) Generate maternal rank file name
      # Use here, to paste the folder path, followed by the file name 
      mat.rank.pred <- paste0(here("data", 
                                   "mat_rank_predictor_hy_n29.txt")) 
      
    ## b) Save maternal rank file 
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(mat_rank_predictor_hy_n29, file = mat.rank.pred, 
                  row.names = F, col.names = F)    
      
    ## c) Generate close proximity file name
      # Use here, to paste the folder path, followed by the file name 
      close.prox.pred <- paste0(here("data", 
                                     "close_prox_predictor_hy_n29.txt")) 
      
    ## d) Save close proximity file
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(close_prox_predictor_hy_n29, file = close.prox.pred, 
                  row.names = F, col.names = F)    
      
    ## e) Generate nursing file name
      # Use here, to paste the folder path, followed by the file name 
      nurse.pred <- paste0(here("data", 
                                "nurse_predictor_hy_n29.txt")) 
      
    ## f) Save nursing file
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(nurse_predictor_hy_n29, file = nurse.pred, 
                  row.names = F, col.names = F)    
      
    ## g) Generate grooming file name
      # Use here, to paste the folder path, followed by the file name 
      groom.pred <- paste0(here("data", "groom_predictor_hy_n29.txt")) 
      
    ## h) Save grooming file
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(groom_predictor_hy_n29, file = groom.pred, 
                  row.names = F, col.names = F)     
      
      
  ### 9.3 Export covariate variables to txt     
    # Save and export part of table as a .txt file. Files are formatted 
    # for use in MACAU as covariate variables
      
    ## a) Generate maternal rank/care covariate file name
      # Use here, to paste the folder path, followed by the file name 
      mat.rank.care.covar <- paste0(here('data', 
                                         'mat_rank_care_covars_hy_n29.txt')) 
      
    ## b) Save maternal rank file 
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(mat_rank_care_covars_hy_n29, file = mat.rank.care.covar, 
                  row.names = F, col.names = F)  

      
      
      
      
      
      
      
      