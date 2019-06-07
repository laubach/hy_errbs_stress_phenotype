###############################################################################
##############        Spotted Hyena eRRBS DNA Methylation:       ##############
##############               Adult stress phenotype              ##############
##############                 By: Zach Laubach                  ##############
##############               created: 20 Nov 2018                ##############
##############             last updated: 7 June 2019             ##############
###############################################################################


  ### PURPOSE: This code is desingned to analyze ERRBS data and associations of 
  # genome wide DNA methylation and adult stress phenotype in hyenas.
  
  
    # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Data management
    # 4: Join tables and re-tidy fecal luma data 
    # 5: Univariate analyses
    # 6: Bivariate data exploration  
    # 7: Bivariate analyses
    # 8: Model fecal corticosterone
    # 9: Format variables for MACUA 
    # 10: Export data files



###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 clear global environment
    rm(list = ls())
  

  ### 1.2 Install and load Mara Hyena Project packages 
    ## a) Load the Mara Hyena Project data files from github
    # Check for devtools and install if not already installed
      if(!'devtools' %in% row.names(installed.packages())){
        install.packages('devtools')
        }

    # load hyenadata package
      library('hyenadata')


  ### 1.3 Install and load packages 
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
        
      # Check for here and install if not already installed
        if (!'here' %in% installed.packages()[,1]){
          install.packages ('here')
        }
      # load here packages
        library ('here')
        
    ## b) Graph Plotting and Visualization Packages
      # Check for ggplot2 and install if not already installed
        if (!'ggplot2' %in% installed.packages()[,1]){
          install.packages ('ggplot2')
        }
      # load ggplot2 packages
        library ('ggplot2')
        
      # Check for gridExtra and install if not already installed
        if (!'gridExtra' %in% installed.packages()[,1]){
          install.packages ('gridExtra')
        }
      # load gridExtra packages
        library ('gridExtra')
        
      # Check for dotwhisker and install if not already installed
        # used with broom to graph beta estimates
        if (!'dotwhisker' %in% installed.packages()[,1]){
          install.packages ('dotwhisker')
        }
      # load dotwhisker packages
        library ('dotwhisker')
        
    ## c) Modeling Packages
      # Check for broom and install if not already installed
        if (!'broom' %in% installed.packages()[,1]){
          install.packages ('broom')
        }
      # load broom packages
        library ('broom')
      
      # Check for nlme and install if not already installed
        if (!'nlme' %in% installed.packages()[,1]){
          install.packages ('nlme')
        }
      # load nlme packages
        library ('nlme')  
        
      # Check for lme4 and install if not already installed
        if (!'lme4' %in% installed.packages()[,1]){
          install.packages ('lme4')
        }
      # load lme4 packages
        library ('lme4')
        

  ### 1.4 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 3.5.1 (2018-07-02)
    # Platform: x86_64-apple-darwin15.6.0 (64-bit)
    # Running under: macOS  10.14
    
  
  ### 1.5 Set working directory 
    setwd(here())
  
  
  ### 1.6 Set file paths for data importing and exporting
    ## a) The path to other variables to be modeled with rrbs data
      rrbs_vars_data_path <- paste("~/R/R_wd/fisi/project/6_hy_RRBS/",
                              "data/", sep = '')
      
    ## b) The path to prey data
      general_data_path <- paste("~/R/R_wd/fisi/project/0_data/",
                              sep = '')
    
    ## c) Source scripts path
      source_path <- paste("~/Git/source_code/")
      
  
  ### 1.7 Source functions
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
      
  ### 2.1 Import RRBS data files (with readr)
    ## a) Import data from copy of sample seleciton file, which includes list 
      # of hyenas and variables
      rrbs_vars <- read_csv(paste(rrbs_vars_data_path,
                                  "hy_rrbs_variables.csv", sep = ''))
        
    ## b) check that all variables are of appropriate class    
      sapply(rrbs_vars, class)
      
      
  ### 2.2 Import Access Fisi data files
    ## a) Import Access data backend from Mara Hyena Project data package
      # Use hyenadata package 'load_all_tables' function to load all tables
      hyenadata::load_all_tables()
      
    ## b) Import working version of tblReprostats
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
   
    ## d) drop unecessary variables
      rrbs_vars <- rrbs_vars %>%
        select (- c(orig.box, orig.cell, sample.type, check.out.state,
                    check.out.vol.u.l., dna.conc.ng.u.l., initial.library.prep,
                    mat.care.sessions, cub.sex, age.at.darting, grooming..blup,
                    blup.tertile, cub.dob)) 
      
    ## e) convert hy.id to character class 
      rrbs_vars$hy.id <- as.character(rrbs_vars$hy.id)
      
      
  ### 3.2 Tidy tblFecalHormones
    ## a) rename data frame
      fecal_horm <- tblFecalHormones
      
    ## b) Convert all text to lower case
      fecal_horm <- AllCharactersToLower(fecal_horm)
      
    ## c) Format variable names (lowercase and separated by '.')
      fecal_horm <- FormatVarNames(fecal_horm)
      
    ## d) Convert hormone concentrations to numeric
      # make list of variable namges that contain 'ng.g'
      horm_columns <- fecal_horm %>%
        select(contains('ng.g')) %>%
        colnames()
      # convert variables (from horm_columns list) to character first 
      fecal_horm <- fecal_horm %>%
        mutate_at(horm_columns, funs(as.character)) 
      
      # convert variables (from horm_columns list) to numeric 
      fecal_horm <- fecal_horm %>%
        mutate_at(horm_columns, funs(as.numeric))
      
      
  ### 3.3 Tidy fecal_repos
    ## a) rename data frame
      fecal_repos <- tblFecalRepository
      
    ## b) Convert all text to lower case
      fecal_repos <- AllCharactersToLower(fecal_repos)
      
    ## c) Format variable names (lowercase and separated by '.')
      fecal_repos <- FormatVarNames(fecal_repos)
      
    ## d) Rename hyena.id as hy.id
      fecal_repos <- fecal_repos %>%
        rename('hy.id' = 'hyena.id')
      
      
  ### 3.4 Tidy repro_state
    ## a) Convert all text to lower case
      repro_state <- AllCharactersToLower(repro_state)
      
    ## b) Format variable names (lowercase and separated by '.')
      repro_state <- FormatVarNames(repro_state)  
      
    ## c) Rename mom as hy.id 
      repro_state <- repro_state %>%
        rename('hy.id' = 'mom')
      
    ## d) Remove first numbering column
      repro_state <- repro_state %>%
        subset(select = -c(x1))
      
      
  ### 3.5 Tidy hyenas
    ## a) rename data frame
      hyenas <- tblHyenas
      
    ## b) Convert all text to lower case
      hyenas <- AllCharactersToLower(hyenas)
      
    ## c) Format variable names (lowercase and separated by '.')
      hyenas <- FormatVarNames(hyenas)    
    
    ## d) Rename id as hy.id 
      hyenas <- hyenas %>%
        rename('hy.id' = 'id') 
      
      
  ### 3.6 Tidy life_hist
    ## a) rename data frame
      life_hist <- tblLifeHistory.wide
      
    ## b) Convert all text to lower case
      life_hist <- AllCharactersToLower(life_hist)
      
    ## c) Format variable names (lowercase and separated by '.')
      life_hist <- FormatVarNames(life_hist) 
      
    ## d) Rename id as hy.id 
      life_hist <- life_hist %>%
        rename('hy.id' = 'id') 
      
    ## e) convert hy.id to character class  
      life_hist$hy.id <- as.character(life_hist$hy.id) 
      
      
  ### 3.7 tidy NA values  
    ## a) Create a vector of values representing NA
      na_strings <- c('NaN', 'NA', 'na')
      
    ## b) Use naniar to replace all values representing NA with actual NA
      life_hist <- life_hist %>%
        replace_with_na_all(condition = ~.x %in% na_strings)
      
    ## c) Convert all '...error' variables from character into numeric
      life_hist <- life_hist %>%
        mutate_at(vars(contains('error')), funs(as.numeric))
      
      
  ### 3.8 Clean global environment
    ## a) Remove extra tables/dataframes
      rm(list = ls(pattern = 'tbl'))
      
      
      
###############################################################################
##############     4. Join tables and re-tidy fecal luma data    ##############
############################################################################### 
      
      
  ### 4.1 Make the fecal_cort (fecal cort) dataframe
    ## a) Left join fecal_repos to fecal_horm, retains all columns from both
      fecal_cort <- fecal_horm %>%
        left_join(select(fecal_repos, c(fecal.sample.id, hy.id, kaycode,
                                        poop.date, poop.time)),
                  by = "fecal.sample.id")
      
    ## b) add cub.dob and darting.date to fecal_cort
      dates_df <- select(rrbs_vars, c(hy.id, darting.date)) 
      
    ## c) Convert dates stored as character (e.g. 8/3/05) to formatted dates
      dates_df$darting.date <- as.Date(dates_df$darting.date, 
                                        format = "%m/%d/%y")
      
    ## d) Left join dates_df to fecal_horm, retains all columns from both 
      fecal_cort <- left_join(fecal_cort,
                              dates_df, by = "hy.id")
      
    ## e) remove fecal_horm, fecal_repos, and dates_df from workspace
      rm(fecal_horm)
      rm(fecal_repos)
      rm(dates_df)
      
    ## f) Left join life_hist to fecal_cort,
      fecal_cort <- fecal_cort %>% 
        left_join(life_hist, by = c('hy.id' = 'hy.id'))
      
    ## g) Left join hyenas to fecal_cort,
      fecal_cort <- fecal_cort %>% 
        left_join(select(hyenas, c(hy.id, sex, status, mom, dad, 
                                   number.littermates, litrank)),
                  by = c('hy.id' = 'hy.id'))  
      
  
  ### 4.2 Tidy fecal_cort
    ## a) Create an estimated age in months by subtracting birthdate from
      # darting date using lubridate and dividing by average # days in month 
      fecal_cort <- fecal_cort  %>%
        mutate(dart.age.mon = round((interval(dob, 
                                              darting.date) %/% days(1) 
                                     / 30.44), 1))  
      
    ## b) Create a varialbe, 'fec.vs.dart,' which indicates whether a fecal 
      # sample was collected before vs after the darting date
      fecal_cort <- fecal_cort  %>%
        mutate(fec.vs.dart = case_when(fecal_cort$poop.date >= 
                                         fecal_cort$darting.date ~ c("after"),
                                       fecal_cort$poop.date < 
                                         fecal_cort$darting.date ~ 
                                         c("before")))
      
    ## c) Identify which hyenas have an 'after' fecal sample
      hy_list <- fecal_cort %>% 
        filter(grepl("after", fec.vs.dart))%>% 
        distinct(hy.id)
      
     # no_after_fec <- anti_join(fecal_cort, hy_list, by = "hy.id")
     
    ## d) Create a varialbe, 'poop.from.dart.mon,' which indicates how many days
      # between the poop.date and darting date
      fecal_cort <- fecal_cort  %>%
        mutate(poop.from.dart.mon = round((interval(darting.date, 
                                              poop.date) %/% days(1) 
                                     / 30.44), 1))
      
    ## e) Create a varialbe, 'poop.age.mon,' which indicates how old hyena was 
      # at the poop.date in days
      fecal_cort <- fecal_cort  %>%
        mutate(poop.age.mon = round((interval(dob, 
                                              poop.date) %/% days(1) 
                                     / 30.44), 1))
      
      #*** HACK *** removed negate poop.age.mon... dob looks wrong
      fecal_cort <- fecal_cort  %>%
        filter(poop.age.mon > 1) 
    
    ## f) Extract month from poop.date
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable  
        fecal_cort$poop.mon <- month(fecal_cort$poop.date)
      
    ## g) Create a varialbe, 'migratn.seas,' which indicates if a poop sample
      # was collected in migration (June 1 - Oct 31)
      fecal_cort <- fecal_cort  %>%
        mutate(migratn.seas = ifelse(poop.mon >= 6 & poop.mon<= 10, 'migration', 
                                     'none'))
      
    ## h) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of migratn.seas variable and sets the reference  
      # level to 'none'
      fecal_cort <- transform(fecal_cort, 
                               migratn.seas = factor(migratn.seas,
                                              levels = c("none", 
                                                         "migration")))  
               
    ## i) Convert poop.time to a datetime class
      fecal_cort$poop.time <- as.POSIXct(paste(fecal_cort$poop.date,
                                                    fecal_cort$poop.time), 
                                              format = '%Y-%m-%d %H:%M:%S')
      
    ## j) Extract am vs. pm from poop.time
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable
      fecal_cort <- fecal_cort  %>%
        mutate(poop.am.pm = ifelse(lubridate::am(poop.time), 'am', 
                                     'pm'))
      
    ## k) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of poop.am.pm variable and sets the reference  
      # level to 'am' 
      fecal_cort <- transform( fecal_cort, 
                               poop.am.pm = factor(poop.am.pm,
                                                     levels = c("am", 
                                                                "pm")))     


  ### 4.3 Combine repro_state w fecal_cort
    ## a) Make an empty data frame
      # This is an empty data frame that can store the overlaping fecal_cort
      # and repro_states data
      fecal_repro_data  <- c()   
     
    ## b) For loop to find overlap
      # Iterate over fecal data (hy.id and poop.date), to find interseciton
      # with repro_state data
      for (i in 1:nrow(fecal_cort)) { 
        
        # loop through 1:n IDs in fecal_cort
        id = paste (fecal_cort$hy.id[i])  
        
        # loop through 1:n dates in fecal_cort
        poop.date <-(fecal_cort$poop.date[i])
      
        # create a dataframe to store the rows from repro_states where the
        # id matches hy.id, and poop.date is in between cycle start and stop 
        overlap_poop_repro <- filter(repro_state, id == hy.id & 
                                       poop.date >= cycle.start & 
                                       poop.date <= cycle.stop)
        
        # Control flow
          # if there is no id match and date overlap, 
          # then go to next loop iteration in fecal_cort
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
      fecal_cort <- fecal_cort %>%
        left_join(select(fecal_repro_data, c(hy.id, poop.date, state, 
                                             cycle.start, cycle.stop, 
                                             trimester, parity)),
                  by = c("hy.id" = "hy.id",
                         "poop.date" = "poop.date")) 
      
      
  ### 4.4 Re-tidy fecal data
    ## a) Change state to character
      fecal_cort$state <- as.character(fecal_cort$state)
      
    ## b) Replaces NA with repro state
      # *** NOTE *** Hyena's less than ~750 days (~24 mon) are not in   
          # tblReprostates, but are by default n = nulliparous. Animals older 
          # than ~24 mon some times have missing data on repro state, 
          # possibly becausecub goes missing - HERE WE MADE DECISION to 
          # classify these animals' rerpro state as o = other
      fecal_cort <- fecal_cort  %>%
        mutate(state = case_when(!is.na(fecal_cort$state)
                                 ~ state,
                                 is.na(fecal_cort$state) &
                                   poop.age.mon < 24
                                 ~ c("n"),
                                 is.na(fecal_cort$state) &
                                   poop.age.mon > 24
                                 ~ c("o")))
      
    ## c) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of state variable and sets the reference level 
      # to 'n' makes this
      fecal_cort <- transform(fecal_cort, 
                             state = factor(state,
                                              levels = c("n", "p", "l", "o")))
      
    ## d) Extract year from fas.date
      # Use lubridate to extract the year during which an FAS was done
      fecal_cort$poop.yr <- year(fecal_cort$poop.date)
      
    ## e) Restrict the fecal_cort data set to match rrbs data
      # Subset the fecal_cort data based on overlapping date ranges starting
      # with the earliest dob of animals in the RRBS data set...in this case 
      # it is 2007 
      fecal_cort <- fecal_cort  %>%
        filter(poop.yr > 2006)
      
      #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL poop which occur in an overlapping period with RRBS data
      
    ## f) Subset fecal_cort based on 13 month cut-off for fecal samples
      # from females only
      fecal_cort <- fecal_cort  %>%
        filter(grepl('f', sex) & poop.age.mon >= 13)
      #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL FECAL CORT. when fecal sample collected when HYENAS >= 13 months
      # (subadult and adult animals)
      
      
###############################################################################
##############              5. Univariate analyses               ##############
###############################################################################      
      
  ### 5.1 Overview
      # Visualize data and generate summary statistics for fecal corticosterone 
      # as the outcome and eRRBS methylation as the explanatory variable.
 
  ### 5.2 Visualize (and transform as needed) raw data        
      ## a) Histogram Outcome (fecal corticosterone)
      ggplot(data=fecal_cort, aes(x=corticosterone.ng.g)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 750, by = 10), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0,750)) +
        labs(title= "Histogram for fecal corticosterone") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Corticosterone (ng/g)", y="Frequency") 

    ## b) Natural log transformation
      fecal_cort$corticosterone.ng.g.log <- log(fecal_cort$corticosterone.ng.g)
      
    ## c) Histogram Outcome (log fecal corticosterone)  
      ggplot(data=fecal_cort, aes(x=corticosterone.ng.g.log)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(-0.5, 7, by = 0.05), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(-0.5,7)) +
        labs(title= "Histogram for log fecal corticosterone") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Log Corticosterone (ng/g)", y="Frequency") 
      
      
  ### 5.3 Mean summary
    ## a) Calucate the average corticosterone measure for each hyena 
      # Start from one year of age and including all additonal fecal samples.
      fecal_cort_avg <- fecal_cort  %>%
        group_by(hy.id) %>%
        summarize(n.fec.corticost = sum(!is.na(corticosterone.ng.g)),
                  mean.fec.corticost = round(mean(corticosterone.ng.g, 
                                            na.rm = T), 2),
                  sd.fec.corticost = round(sd(corticosterone.ng.g, 
                                              na.rm = T), 2))

    ## b) Calucate the average of log corticosterone measure for each hyena 
      # Start from one year of age and including all additonal fecal samples.
      log_fecal_cort_avg <- fecal_cort  %>%
        group_by(hy.id) %>%
        summarize(n.log.fec.corticost = sum(!is.na(corticosterone.ng.g.log)),
                  mean.log.fec.corticost = round(mean(corticosterone.ng.g.log, 
                                                  na.rm = T), 2),
                  sd.log.fec.corticost = round(sd(corticosterone.ng.g.log, 
                                              na.rm = T), 2))
  
      
      
###############################################################################
##############           6. Bivariate data exploration           ##############
###############################################################################       
    
  ### 6.1 Bivariate statistics fecal cort. by age at fecal collection
    ## a) Plot fecal corticosterone by age at fecal collection
      # NOTE: These summaries contain non-independent measures of fecal cort
      ggplot(data = subset(fecal_cort, !is.na(x = poop.age.mon)),
#   ************** Scatterplot START **************
             aes(x = poop.age.mon, y = corticosterone.ng.g)) +
        geom_point(shape = 1) +
        geom_smooth(method = loess, se = F) + # Add smooth curve best fit lines
#   *************** Scatterplot END ***************
        
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Fecal corticosterone by 
             age on fecal collection date") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "darkgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("Fecal Corticosterone (ng/g)") +
        xlab("Age (months)")
      
    ## b) Save Plot
      # use ggsave to save the plot
      ggsave("cort_by_age_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
      
  ### 6.2  Bivariate statistics fecal cort. by reproductive state
    ## a) Summary stats corticosterone.ng.g by reproductive state
      # NOTE: These summaries contain non-independent measures of fecal cort
      cort_by_repro_state <- fecal_cort %>%
        group_by (state) %>%
        summarise (n.samp.id = n(),
                   n.hy.id =  n_distinct(hy.id),
                   avg = round (mean(corticosterone.ng.g, na.rm = T), 2),
                   median =  round (quantile(corticosterone.ng.g, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(corticosterone.ng.g, na.rm = T), 2))  
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/cort_by_repro_state.pdf"),
          height = 4, width = 7)
      grid.table(cort_by_repro_state)
      dev.off()
      
      ## c) Plot fecal corticosterone by reproductive state
      # NOTE: These summaries contain non-independent measures of fecal cort
      ggplot(data = subset(fecal_cort, !is.na(x = state)), 
             aes(x = state, y = corticosterone.ng.g, color = state)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Fecal corticosterone by 
             reproductive state") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "darkgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("Fecal Corticosterone (ng/g)") +
        xlab("Reproductive state")
      
      ## d) Save Plot
      # use ggsave to save the plot
      ggsave("cort_by_repro_state_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)    
      
      
  ### 6.3  Bivariate statistics fecal cort. by fecal collection time of day 
      # (am/pm)
    ## a) Summary stats corticosterone.ng.g by time of day 
      # NOTE: These summaries contain non-independent measures of fecal cort
      cort_by_am_pm <- fecal_cort %>%
        group_by (poop.am.pm) %>%
        summarise (n.samp.id = n(),
                   n.hy.id =  n_distinct(hy.id),
                   avg = round (mean(corticosterone.ng.g, na.rm = T), 2),
                   median =  round (quantile(corticosterone.ng.g, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(corticosterone.ng.g, na.rm = T), 2))  
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/cort_by_am_pm.pdf"),
          height = 3, width = 5)
      grid.table(cort_by_am_pm)
      dev.off()
      
    ## c) Plot fecal corticosterone by time of day
      # NOTE: These summaries contain non-independent measures of fecal cort
      ggplot(data = subset(fecal_cort, !is.na(x = poop.am.pm)), 
             aes(x = poop.am.pm, y = corticosterone.ng.g, color = poop.am.pm)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Fecal corticosterone by fecal
             collection time of day") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "darkgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("Fecal Corticosterone (ng/g)") +
        xlab("Time of day (am/pm)")
      
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("cort_by_am_pm_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)    
      
      
    ### 6.4  Bivariate statistics fecal cort. by migration season (on date the
      # sample was collected)
      ## a) Summary stats corticosterone.ng.g by migratn.seas.fec
      # NOTE: These summaries contain non-independent measures of fecal cort
      cort_by_migratn <- fecal_cort %>%
        group_by (migratn.seas) %>%
        summarise (n.samp.id = n(),
                   n.hy.id =  n_distinct(hy.id),
                   avg = round (mean(corticosterone.ng.g, na.rm = T), 2),
                   median =  round (quantile(corticosterone.ng.g, 
                                             c(.5), na.rm = T), 2),
                   sd = round (sd(corticosterone.ng.g, na.rm = T), 2))  
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(paste0(here(),"/output/cort_by_migratn.pdf"),
          height = 3, width = 5)
      grid.table(cort_by_migratn)
      dev.off()
      
    ## c) Plot fecal corticosterone by migration season
      # NOTE: These summaries contain non-independent measures of fecal cort
      ggplot(data = subset(fecal_cort, !is.na(x = migratn.seas)), 
             aes(x = migratn.seas, y = corticosterone.ng.g, 
                 color = migratn.seas)) + 
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs(title = "Fecal corticosterone by 
             migration season") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(legend.position = "none") + # remove legend
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = "white")) +
        # add major axes
        theme(axis.line = element_line(colour = "darkgrey", 
                                       size = 1, linetype = "solid")) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face="bold", color="black", 
                                         size=18, angle=0,
                                         margin = margin(t = 0, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face="bold", color="black", 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        ylab("Fecal Corticosterone (ng/g)") +
        xlab("Migration season")
      
    ## d) Save Plot
      # use ggsave to save the plot
      ggsave("cort_by_migratn_plot.pdf", plot = last_plot(), device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 7, height = 5,
             units = c("in"), dpi = 300, limitsize = TRUE)    
      
       
         
###############################################################################
##############               7. Bivariate analyses               ##############
###############################################################################  
      
  ### 7.1 Overview
    ## a) Asssess variable relationship to formally identify potential 
      # confounding variables, precision variables, and effect modification 
      # (significant interactions) so that models can be appropriately and 
      # efficiently parameterized. 
      
      
  ### 7.2 Assess potential precision variables
    ## a) Bivariate Regression: Fecal cort. by age when fecal sample collected
      # uses 'nmle' package, which will provided p-value estimates
      fec.age.lme <- lme(corticosterone.ng.g.log ~ poop.age.mon, 
                         random = ~1|hy.id, 
                         subset(fecal_cort,!is.na(x = poop.age.mon)))
      
      # Summary and parameter estimates
      summary(fec.age.lme) 
      intervals(fec.age.lme, which = "fixed")
      
    ## b) Bivariate Regression: Fecal cort. by reproductive state when
      # fecal sample collected
      # uses 'nmle' package, which will provided p-value estimates
      repro.state.lme <- lme(corticosterone.ng.g.log ~ state, 
                             random = ~1|hy.id, 
                             subset(fecal_cort,!is.na(x = state)))
      
      # Summary and parameter estimates
      summary(repro.state.lme) 
      intervals(repro.state.lme, which = "fixed")
      
    ## c) Bivariate Regression: Fecal cort. by time of day (am/pm) when
      # fecal sample collected
      # uses 'nmle' package, which will provided p-value estimates
      am.pm.lme <- lme(corticosterone.ng.g.log ~ poop.am.pm, 
                       random = ~1|hy.id, 
                       subset(fecal_cort,!is.na(x = poop.am.pm)))
      
      # Summary and parameter estimates
      summary(am.pm.lme) 
      intervals(am.pm.lme, which = "fixed")  
      
    ## d) Bivariate Regression: Fecal cort. by migration status when
      # fecal sample collected
      # uses 'nmle' package, which will provided p-value estimates
      migratn.lme <- lme(corticosterone.ng.g.log ~ migratn.seas, 
                         random = ~1|hy.id, 
                         subset(fecal_cort,!is.na(x = migratn.seas)))
      
      # Summary and parameter estimates
      summary(migratn.lme) 
      intervals(migratn.lme, which = "fixed")  
      
      
            
###############################################################################
##############           8. Model fecal corticosterone           ##############
###############################################################################      
      
      
  ### 8.1 BLUPs from mixed model linear regression
    
    # NOTE: Use when there are repeated measuresments for a variable
       # that is to be used as an explanatory variable in another analysis.
       # Can control for other variables that bias estimates of explanatory
       # variable
          
    # NOTE: BLUPs are conditional means from linear model with a
      # Gaussian distribution (according to Doug Bates)
      # BLUP = fixef(intrcpt) + ranef 
       
    ## a) Calucate the BLUPs, individual variation in fecal adult cort.    
      fec.cort.lmm <- lme4::lmer(corticosterone.ng.g.log ~ poop.age.mon + 
                                   state + poop.am.pm + #migratn.seas + 
                                   (1|hy.id ), data = fecal_cort)
      
    ## b) Generate mixed model summary 
      summary(fec.cort.lmm)
      ranef(fec.cort.lmm) # random effect
      fixef(fec.cort.lmm) # fixed effect
      
    ## c) extract BLUPs from mixed model object
      blups <- coef(fec.cort.lmm)[[1]] # extract BLUPs as a dataframe
                                       # BLUPs = rand ef. + fix ef. (intercept)
      
    ## d) Use tibble to add row names as their own column
      blups <- rownames_to_column(blups, "hy.id") 
  
    ## e) Rename variables in blups table
      blups <- rename(blups, 'log_cort' = '(Intercept)') 
      
    ## f) Create a new dataframe that includes hyena id, the log cort BLUPs,
      # and the exponeniated (biological scale) cort BLUPs
      fec_cort_blups <- as_tibble(cbind(hy.id = blups$hy.id,
                              log.cort = blups$log_cort, 
                              cort = exp(blups$log_cort)), round = 4)
      
    ## g) Coerce from character to numeric class
      fec_cort_blups$log.cort <- as.numeric(fec_cort_blups$log.cort)
      fec_cort_blups$cort <- as.numeric(fec_cort_blups$cort)
   
      
  ### 8.2 Graph adult fecal cort conditional averages (BLUPs)
      ## a) Histogram Outcome (fecal corticosterone BLUPs)  
      ggplot(data=fec_cort_blups, aes(x=cort)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(20, 100, by = 7), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(20, 100)) +
        labs(title= "Histogram of adult fecal corticosterone BLUPs") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Corticosterone BLUPs (ng/g)", y="Frequency") 
      
      ## b) Histogram Outcome (log fecal corticosterone)  
      ggplot(data=fec_cort_blups, aes(x=log.cort)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(2.5, 5, by = 0.15), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(2.5, 5)) +
        labs(title= "Histogram of log adult fecal corticosterone BLUPs") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Log Corticosterone BLUPs (ng/g)", y="Frequency") 
      
      

###############################################################################
##############           9. Format variables for MACUA           ##############
###############################################################################      
      
  ### 9.1 Make dataframe for txt. file export to MACUA (for EWAS)
    ## a) Manual data clean up drop animals that failed ERRBS library QAQC
      cort_rrbs_vars <- rrbs_vars %>%
        filter (!grepl('tato', hy.id))
      
    ## b) Left join fec_cort_blups to cort_rrbs_vars
      cort_rrbs_vars <- cort_rrbs_vars %>%
        left_join(fec_cort_blups, by = 'hy.id')
      
    ## c) Group by hy.id and save age in months at time of darting
      age_errbs <- fecal_cort %>%
        group_by(hy.id) %>%
        summarise(age.dart.mon = first(dart.age.mon))
      
    ## d) Left join fecal_cort$dart.age to cort_rrbs_vars
      cort_rrbs_vars <- cort_rrbs_vars %>%
        left_join(age_errbs, by = 'hy.id') 
      
    ## e) Left join fecal_cort_avg to cort_rrbs_vars
      cort_rrbs_vars <- cort_rrbs_vars %>%
        left_join(fecal_cort_avg, by = 'hy.id')
      
    ## f) Left join fecal_cort_avg to fec_cort_blups
      fec_cort_blups_avg <- fec_cort_blups %>%
        left_join(fecal_cort_avg, by = 'hy.id')
    
      
  ### 9.2 Arrange data set by hy.id abc order
      ## This is necessary because the order of rows matches data to a 
      ## specific animal
      cort_rrbs_vars <- arrange(cort_rrbs_vars, hy.id)
      
      
  ### 9.3 Subset explanatory/predictor variables and format for use in MACAU
    ## a) Fecal corticosterone formatted
      cort_predictor_hy_n29 <- cort_rrbs_vars %>%
        select(cort) 
      
    ## b) Nat. log fecal corticosterone formatted
      log_cort_predictor_hy_n29 <- cort_rrbs_vars %>%
        select(log.cort) 
      

  ### 9.4 Subset explanatory/predictor variables and format for use in MACAU
    ## a) Covariates (variables to control analyses)
      fec_cort_covars_hy_n29 <- cort_rrbs_vars %>%
        mutate(intercept = 1) %>%
        select(intercept, age.dart.mon) 

      
      
###############################################################################
##############                10. Export data files              ##############
###############################################################################
     
  ### 10.1 Export fecal and rrbs data to csv     
    # Save and export tables as a .cvs spreadsheet. 
    # Files are saved in the 'output' folder in the working directory.
      
    ## a) Generate file name for fecal_cort
      # this the full fas data set with repeated fas samples
      csv.file.name.fec.cort <- paste0(here("data",
                                               "fecal_cort.csv")) 
    ## b) Save mat_car_fas table as csv
      write.csv (fecal_cort, file = csv.file.name.fec.cort)
      
    ## c) Generate file name for fec_cort_blups_avg
      # this is the data set summarized using blups and with covariates
      csv.file.name.fec.cort.blups.avg <- paste0(here("data",
                                             "fec_cort_blups_avg.csv")) 
    ## d) Save mat_car_fas_blups table as csv
      write.csv (fec_cort_blups_avg, file = csv.file.name.fec.cort.blups.avg) 
      
    ## e) Generate file name for mat_care_rrbs
      # includes fas data and covariates only for animals with rrbs data
      csv.file.name.cort.rrbs <- paste0(here("data",
                                                "cort_rrbs_vars.csv"))
      
    ## f) Save mat_care_rrbs table as csv
      write.csv (cort_rrbs_vars, file = csv.file.name.cort.rrbs)
      
      
  ### 10.2 Export explanatory/predictor variables to txt     
      # Save and export part of table as a .txt file. Files are formatted 
      # for use in MACAU as explanatory/predictor variables
      
    ## a) Generate fecal cort file name
      # Use here, to paste the folder path, followed by the file name 
      cort.pred <- paste0(here("data", 
                                   "cort_predictor_hy_n29.txt")) 
      
    ## b) Save fecal cort file 
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(cort_predictor_hy_n29, file = cort.pred, 
                  row.names = F, col.names = F)    
      
    ## c) Generate fecal cort file name
      # Use here, to paste the folder path, followed by the file name 
      log.cort.pred <- paste0(here("data", 
                               "log_cort_predictor_hy_n29.txt")) 
      
    ## d) Save fecal cort file 
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(log_cort_predictor_hy_n29, file = log.cort.pred, 
                  row.names = F, col.names = F)  
          
      
  ### 10.3 Export covariate variables to txt     
      # Save and export part of table as a .txt file. Files are formatted 
      # for use in MACAU as covariate variables
      
    ## a) Generate fecal cort covariate file name
      # Use here, to paste the folder path, followed by the file name 
      cort.covar <- paste0(here('data',
                                'fec_cort_covars_hy_n29.txt')) 
      
    ## b) Save fecal cort file 
      # Save part of dataframe as a .txt file into the 
      # output data folder in the working directory.
      write.table(fec_cort_covars_hy_n29, file = cort.covar, 
                  row.names = F, col.names = F)  
      
      
      
      