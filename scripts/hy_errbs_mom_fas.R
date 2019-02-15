###############################################################################
##############        Spotted Hyena eRRBS DNA Methylation:       ##############
##############             Maternal care (FAS data)              ##############
##############                 By: Zach Laubach                  ##############
##############               created: 14 Feb 2018                ##############
##############             last updated: 14 Feb 2019            ##############
###############################################################################


  ### PURPOSE: This code is desingned to analyze ERRBS data and associations of 
  # maternal care from FAS data and genome wide DNA methylation in hyenas.
  
  
    # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Data management



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
  
    ## c) add cub.dob and darting.date to fas_data
      dates_df <- select(rrbs_vars, c(hy.id, darting.date, cub.dob, 
                                      dart.age.mon)) 
      
      fas_data <- left_join(fas_data,
                              dates_df, by = "hy.id")
      
      
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
      fas_data <- transform( fas_data, 
                              fas.am.pm = factor(fas.am.pm,
                                                     levels = c("am", 
                                                                "pm"))) 
      
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
##############                10. Export data files               ##############
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
      
      
      
      
      
      
      
      
      
      
      
      
      