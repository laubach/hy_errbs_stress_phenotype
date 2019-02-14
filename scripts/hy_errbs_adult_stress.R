###############################################################################
##############        Spotted Hyena eRRBS DNA Methylation:       ##############
##############               Adult stress phenotype              ##############
##############                 By: Zach Laubach                  ##############
##############               created: 20 Nov 2018                ##############
##############             last updated: 8 Feb 2019              ##############
###############################################################################


  ### PURPOSE: This code is desingned to analyze ERRBS data and associations of 
  # genome wide DNA methylation and adult stress phenotype in hyenas.
  
  
    # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Data management
    # 4: Model fecal corticosterone 
    # 5: Export data files




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
        
      # Check for lme4 and install if not already installed
        if (!'lme4' %in% installed.packages()[,1]){
          install.packages ('lme4')
        }
      # load lme4 packages
        library ('lme4')
        

  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 3.5.1 (2018-07-02)
    # Platform: x86_64-apple-darwin15.6.0 (64-bit)
    # Running under: macOS  10.14
    
  
  ### 1.4 Set working directory 
    setwd(here())
  
  
  ### 1.5 Set file paths for data importing and exporting
    ## a) The path to other variables to be modeled with rrbs data
      rrbs_vars_data_path <- paste("~/R/R_wd/fisi/project/6_hy_RRBS/",
                              "data/", sep = '')
      
    ## b) The path to prey data
      general_data_path <- paste("~/R/R_wd/fisi/project/0_data/",
                              sep = '')
    
    ## c) Source scripts path
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
      
    ## c) Import fecal horomens data files, including fecal repository 
      fecal_horm <- read_csv(paste(rrbs_vars_data_path,
                                  "tblFecalHormones.csv", sep = ''))
    ## d) Fecal repository
      fecal_repos <- read_csv(paste(rrbs_vars_data_path,
                                    "tblFecalRepository.csv", sep = ''))
    ## e) tblReprostats
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
      
      
  ### 3.2 Tidy fecal_horm
    ## a) view the variable classes to determine which need to be modified
      spec(fecal_horm)
      
    ## b) Convert all text to lower case
      fecal_horm <- AllCharactersToLower(fecal_horm)
      
    ## c) Format variable names (lowercase and separated by '.')
      fecal_horm <- FormatVarNames(fecal_horm)
      
  
  ### 3.3 Tidy fecal_repos
    ## a) view the variable classes to determine which need to be modified
      spec(fecal_repos)
    
    ## b) Convert all text to lower case
      fecal_repos <- AllCharactersToLower(fecal_repos)
    
    ## c) Format variable names (lowercase and separated by '.')
      fecal_repos <- FormatVarNames(fecal_repos)
    
    ## d) Convert dates stored as character (e.g. 8/3/05) to formatted dates
      fecal_repos$poop.date <- as.Date(fecal_repos$poop.date, 
                                        format = "%m/%d/%y")
    
    ## e) rename variables
      fecal_repos <- fecal_repos %>%
        rename('hy.id' = 'hyena.id')
    
    ## f) convert hy.id to character class  
      fecal_repos$hy.id <- as.character(fecal_repos$hy.id)
      
 
  ### 3.4 Join Data sets
    ## a) Semi join fecal_repos by rrbs_vars, retains all rows in 
      # fecal_repos that has a match id in rrbs_vars   
      fecal_repos <- semi_join(fecal_repos,
                              rrbs_vars, by = "hy.id")
      
    
    ## b) Left join fecal_repos fecal_horm, retains all columns from both
      # fecal_repos that has a match id in rrbs_vars   
      fecal_data <- left_join(fecal_repos,
                               fecal_horm, by = "fecal.sample.id")
      
      
    ## c) add cub.dob and darting.date to fecal_data
      dates_df <- select(rrbs_vars, c(hy.id, darting.date, cub.dob)) 
      
      fecal_data <- left_join(fecal_data,
                              dates_df, by = "hy.id")
      
    ## d) remove fecal_horm, fecal_repos, and dates_df from workspace
      rm(fecal_horm)
      rm(fecal_repos)
      rm(dates_df)

      
  ### 3.5 Tidy fecal_data
    ## a) Create a varialbe, 'fec.vs.dart,' which indicates whether a fecal 
      # sample was collected before vs after the darting date
      fecal_data <- fecal_data  %>%
        mutate(fec.vs.dart = case_when(fecal_data$poop.date >= 
                                         fecal_data$darting.date ~ c("after"),
                                       fecal_data$poop.date < 
                                         fecal_data$darting.date ~ 
                                         c("before")))
      
    ## b) Identify which hyenas don't have an 'after' fecal sample
      hy_list <- fecal_data %>% 
        filter(grepl("after", fec.vs.dart))%>% 
        distinct(hy.id)
      
     # no_after_fec <- anti_join(fecal_data, hy_list, by = "hy.id")
     
    ## c) Create a varialbe, 'fec.vs.dart.days,' which indicates how many days
      # between the poop.date and darting date
      fecal_data <- fecal_data  %>%
        mutate(fec.vs.dart.days = (fecal_data$poop.date - 
                                     fecal_data$darting.date)) 
      
    ## d) Create a varialbe, 'fecal.age.days,' which indicates how old hyena was 
      # at the poop.date in days
      fecal_data <- fecal_data  %>%
        mutate(fecal.age.days = (fecal_data$poop.date - 
                                     fecal_data$cub.dob))   
    
    ## e) Extract month from poop.date
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable  
        fecal_data$poop.mon <- month(fecal_data$poop.date)
      
    ## f) Create a varialbe, 'migratn.seas,' which indicates if a poop sample
      # was collected in migration (June 1 - Oct 31)
      fecal_data <- fecal_data  %>%
        mutate(migratn.seas = ifelse(poop.mon >= 6 & poop.mon<= 10, 'migration', 
                                     'none'))
      
    ## g) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of migratn.seas variable and sets the reference  
      # level to 'none'
      fecal_data <- transform( fecal_data, 
                               migratn.seas = factor(migratn.seas,
                                              levels = c("none", 
                                                         "migration")))  
               
    ## h) Convert poop.time to a datetime class
      fecal_data$poop.time <- as.POSIXct(fecal_data$poop.time, 
                                       format = "%m/%d/%y %H:%M")
      
    ## i) Extract am vs. pm from poop.time
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable
      fecal_data <- fecal_data  %>%
        mutate(poop.am.pm = ifelse(lubridate::am(poop.time), 'am', 
                                     'pm'))
      
    ## j) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of poop.am.pm variable and sets the reference  
      # level to 'am' 
      fecal_data <- transform( fecal_data, 
                               poop.am.pm = factor(poop.am.pm,
                                                     levels = c("am", 
                                                                "pm"))) 
      
    ## k) Subset fecal_data based on one year cut-off for fecal samples 
      # collected prior to darting.data
      fecal_data <- fecal_data  %>%
        filter(fecal.age.days >= 365)
    #*** NOTE *** This is a data inclusion cut-off decision. 
          # ALL FECAL CORT. when HYENAS >1yr


  ### 3.6 Tidy repro_state       
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
      
      
  ### 3.7 Combine repro_state w fecal_data
    ## a) Make an empty data frame
      # This is an empty data frame that can store the overlaping fecal_data
      # and repro_states data
      fecal_repro_data  <- c()   
     
    ## b) For loop to find overlap
      # Iterate over fecal data (hy.id and poop.date), to find interseciton
      # with repro_state data
      for (i in 1:nrow(fecal_data)) { 
        
        # loop through 1:n IDs in fecal_data
        id = paste (fecal_data$hy.id[i])  
        
        # loop through 1:n dates in fecal_data
        poop.date <-(fecal_data$poop.date[i])
      
        # create a dataframe to store the rows from repro_states where the
        # id matches hy.id, and poop.date is in between cycle start and stop 
        overlap_poop_repro <- filter(repro_state, id == hy.id & 
                                       poop.date >= cycle.start & 
                                       poop.date <= cycle.stop)
        
        # Control flow
          # if there is no id match and date overlap, 
          # then go to next loop iteration in fecal_data
          if (nrow(overlap_poop_repro) < 1) {
            next
          }
      
        # add the poop.date onto the overlap_poop_repro data
        overlap_poop_repro <- cbind(overlap_poop_repro, poop.date, id)
        
        # add the filtered overlap_poop_repro data to a new dataframe
        # over each iteration of the loop
        fecal_repro_data <- rbind(fecal_repro_data, 
                                  overlap_poop_repro)
      }
      
    ## c) Join repro state data to fecal data 
      fecal_data <- fecal_data %>%
        left_join(select(fecal_repro_data, c(hy.id, poop.date, state, 
                                             cycle.start, cycle.stop, 
                                             trimester, parity)),
                  by = c("hy.id" = "hy.id",
                         "poop.date" = "poop.date")) 
      
      
  ### 3.8 Re-tidy fecal data
    ## a) Change state to character
      fecal_data$state <- as.character(fecal_data$state)
      
    ## b) Replaces NA with repro state
      # *** NOTE *** Hyena's less than ~750 days old are not in tblReprostates,  
          # but are by default n = nulliparous. Animals older than ~750 days
          # sometimes have missing data on repro state, possibly because
          # cub goes missing - HERE WE MADE DECISION to classify these
          # animals' rerpro state as o = other
      fecal_data <- fecal_data  %>%
        mutate(state = case_when(!is.na(fecal_data$state)
                                 ~ state,
                                 is.na(fecal_data$state) &
                                   fecal.age.days < 750
                                 ~ c("n"),
                                 is.na(fecal_data$state) &
                                   fecal.age.days > 750
                                 ~ c("o")))
      
    ## c) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of state variable and sets the reference level 
      # to 'n' makes this
      fecal_data <- transform( fecal_data, 
                             state = factor(state,
                                              levels = c("n", "p", "l", "o")))
      
      

###############################################################################
##############           4. Model fecal corticosterone           ##############
###############################################################################      
      
  ### 4.1 Overview
    # Generate summary / estimates for repeated measures fecal corticosterone 
    # to be used as independent variable in models.
  
  ### 4.2 Visualize (and transform as needed) raw data        
      ## a) Histogram Outcome (fecal corticosterone)
      ggplot(data=fecal_data, aes(x=corticosterone.ng.g)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 750, by = 10), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0,750)) +
        labs(title= "Histogram for fecal corticosterone") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Corticosterone (ng/g)", y="Frequency") 
   class(adult_fec_cort$cort)
    ## b) Natural log transformation
      fecal_data$corticosterone.ng.g.log <- log(fecal_data$corticosterone.ng.g)
      
    ## c) Histogram Outcome (log fecal corticosterone)  
      ggplot(data=fecal_data, aes(x=corticosterone.ng.g.log)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(-0.5, 7, by = 0.05), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(-0.5,7)) +
        labs(title= "Histogram for log fecal corticosterone") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Log Corticosterone (ng/g)", y="Frequency") 
      
      
  ### 4.3 Mean summary
    ## a) Calucate the average corticosterone measure for each hyena 
      # Start from one year of age and including all additonal fecal samples.
      fecal_data_avg <- fecal_data  %>%
        group_by(hy.id) %>%
        summarize(n.fec.corticost = sum(!is.na(corticosterone.ng.g)),
                  mean.fec.corticost = round(mean(corticosterone.ng.g, 
                                            na.rm = T), 2),
                  sd.fec.corticost = round(sd(corticosterone.ng.g, 
                                              na.rm = T), 2))

    ## b) Calucate the average of log corticosterone measure for each hyena 
      # Start from one year of age and including all additonal fecal samples.
      log_fecal_data_avg <- fecal_data  %>%
        group_by(hy.id) %>%
        summarize(n.log.fec.corticost = sum(!is.na(corticosterone.ng.g.log)),
                  mean.log.fec.corticost = round(mean(corticosterone.ng.g.log, 
                                                  na.rm = T), 2),
                  sd.log.fec.corticost = round(sd(corticosterone.ng.g.log, 
                                              na.rm = T), 2))
  
      
  ### 4.4 BLUPs from mixed model linear regression
    
    # NOTE: Use when there are repeated measuresments for a variable
       # that is to be used as an explanatory variable in another analysis.
       # Can control for other variables that bias estimates of explanatory
       # variable
          
    # NOTE: BLUPs are conditional means from linear model with a
      # Gaussian distribution (according to Doug Bates)
      # BLUP = fixef(intrcpt) + ranef 
       
    ## a) Calucate the BLUPs, individual variation in fecal adult cort.    
      fec.cort.lmm <- lme4::lmer(corticosterone.ng.g.log ~ fecal.age.days + 
                                   state + poop.am.pm + migratn.seas + 
                                   (1|hy.id ), data = fecal_data)
      
    ## b) Generate mixed model summary 
      summary(fec.cort.lmm)
      ranef(fec.cort.lmm) # random effect
      fixef(fec.cort.lmm) # fixed effect
      
    ## c) extract BLUPs from mixed model object
      blups <- coef(fec.cort.lmm)[[1]] # extract BLUPs as a dataframe
                                       # BLUPs = rand ef. + fix ef. (intercept)
      
    ## d) Use tibble to add row names as their own column
      blups <- rownames_to_column(blups, "id") 
      
    ## e) Rename variables in blups table
      blups <- rename(blups, 'log_cort' = '(Intercept)') 
      
    ## f) Create a new dataframe that includes hyean id, the log cort BLUPs,
      # and the exponeniated (biological scale) cort BLUPs
      adult_fec_cort <- as.tibble(cbind(id = blups$id,
                              log.cort = blups$log_cort, 
                              cort = exp(blups$log_cort)), round = 4)
      
    ## g) Coerce from character to numeric class
      adult_fec_cort$log.cort <- as.numeric(adult_fec_cort$log.cort)
      adult_fec_cort$cort <- as.numeric(adult_fec_cort$cort)
   
      
  ### 4.5 Graph adult fecal cort conditional averages (BLUPs)
      ## a) Histogram Outcome (fecal corticosterone BLUPs)  
      ggplot(data=adult_fec_cort, aes(x=cort)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(20, 100, by = 7), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(20, 100)) +
        labs(title= "Histogram of adult fecal corticosterone BLUPs") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Corticosterone BLUPs (ng/g)", y="Frequency") 
      
      ## b) Histogram Outcome (log fecal corticosterone)  
      ggplot(data=adult_fec_cort, aes(x=log.cort)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(2.5, 5, by = 0.15), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(2.5, 5)) +
        labs(title= "Histogram of log adult fecal corticosterone BLUPs") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Log Corticosterone BLUPs (ng/g)", y="Frequency") 
      

      
###############################################################################
##############                5. Export data files               ##############
###############################################################################
      
  ### 10.1 Export fecal corticosterone BLUPs to csv     
      # Save and export tables as a .cvs spreadsheet and named with today's
      # date. Files are saved in the 'output' folder in the working directory.
      
    ## a) Generate File Names
      # For each table that will be saved as a .csv file, first generate a 
      # file name to save each table
      # here, we paste the folder path, followed by the file name 
      csv.file.name.fecal.cort.BLUPs <- paste0(here("output", 
                                                    "fecal_cort_BLUPs.csv")) 
      
    ## b) Save Tables 
      # Save data frame as a .csv file (a spreadsheet/table) into the 
      # output data folder in the working directory.
      write.csv (adult_fec_cort, file = csv.file.name.fecal.cort.BLUPs)
      
      
      
      
      
      
      
      
      
      
      
      
      