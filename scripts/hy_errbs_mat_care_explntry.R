###############################################################################
##############        Spotted Hyena eRRBS DNA Methylation:       ##############
##############             Maternal care (FAS data)              ##############
##############                 By: Zach Laubach                  ##############
##############               created: 14 Feb 2019                ##############
##############             last updated: 6 June 2019             ##############
###############################################################################


### PURPOSE: This code is desingned to analyze ERRBS data and associations of 
# maternal care from FAS data and genome wide DNA methylation in hyenas.

    # Code Blocks
    # 1: Configure workspace
    # 2: Import data
    # 3: Tidy individual tables
    # 4: Join tables and re-tidy fas data
    # 5: Model 'close proximity' behaviors
    # 6: Model 'nursing' behaviors
    # 7: Model 'grooming' behaviors
    # 8: Format variables for downstream analyses 
    # 9: Export data files 



###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 Set global environment
    ## a) clear global environment
      rm(list = ls())

    ## b) prevent R from automatically reading charater strins as factors
      options(stringsAsFactors = FALSE)
  
  ### 1.2 Install and load Mara Hyena Project packages 
    ## a) Load the Mara Hyena Project data files from github
      # Check for devtools and install if not already installed
        if(!'devtools' %in% row.names(installed.packages())){
          install.packages('devtools')
        }
    ## b) Use devtools to install hyeanadata package from the MaraHyenaProject
      # on github
        devtools::install_github('MaraHyenaProject/hyenadata',
                      auth_token = "e31091e1f6277964e3a52cf3e62ffdc59ea1cadc")
      # load hyenadata package
        library('hyenadata')

  ### 1.2 Install and load packages CRAN packages  
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
        
      # Check for naniar and install if not already installed
        # used to replace values with NA
        if (!'naniar' %in% installed.packages()[,1]){
          install.packages ('naniar')
        }
      # load naniar package
        library ('naniar') 
        
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

        
  ### 1.4 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 3.5.2 (2018-12-20)
    # Platform: x86_64-apple-darwin15.6.0 (64-bit)
    # Running under: macOS  Mojave 10.14.3
    
  
  ### 1.5 Set working directory 
    setwd(here())

    
  ### 1.6 Set file paths for data importing and exporting
    ## a) The path to errbs data
    rrbs_vars_data_path <- paste("~/R/R_wd/fisi/project/6_hy_RRBS/",
                                 "data/", sep = '')
    
    ## b) The path to maternal care FAS data
      mat_fas_path <- paste("~/R/R_wd/fisi/project/fas_maternal_care/",
                                   "fas/output/", sep = '')    
      
    ## c) Source scripts path
      source_path <- paste("~/Git/source_code/")
      
    ## d) The path to Access Fisi backend data
      general_data_path <- paste("~/R/R_wd/fisi/project/0_data/",
                                 sep = '')
      
  
  ### 1.7 Source functions
    ## a) all_char_to_lower function
      source(file = paste0(source_path, "all_char_to_lower.R"))
      
    ## b) format_var_names function
      source(file = paste0(source_path, "format_var_names.R"))
      
    ## c) format_var_names_dash function
      source(file = paste0(source_path, "format_var_names_dash.R"))  
      
      

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
      
      
  ### 2.2 Import maternal care FAS data
    ## a) Import maternal care FAS data which has undergone QAQC (see R 
      # script fas_mat_care_calc_annotated.R in the directory 
      # R...fas_matneral_care)
      mat_fas <- read_csv(paste(mat_fas_path,
                                 "fas_counts.csv", sep = ''))  
      
      
  ### 2.3 Import Access Fisi data files
    ## a) Import Access data backend from Mara Hyena Project data package
      # Use hyenadata package 'load_all_tables' function to load all tables
      hyenadata::load_all_tables()
      
    ## b) Import working version of tblReprostats
      repro_state <- read_csv(paste(general_data_path,
                                    "reprostates_EDS.csv", sep = ''))
      
      

###############################################################################
##############              3. Tidy individual tables            ##############
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
                                              darting.date) %/% days(1) / 
                                       30.44), 1))
     
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
      
    ## b) remove first numbering column
      mat_fas <- mat_fas %>%
        subset(select = -c(X1))
      
    ## c) rename variables
      mat_fas <-  mat_fas %>%
        rename('hy.id' = 'cub')
      
    ## d) Convert all text to lower case
      mat_fas <- AllCharactersToLower(mat_fas)
      
    ## d) convert hy.id to character class  
  #    mat_fas$hy.id <- as.character(mat_fas$hy.id)  
      
    ## e) rename variables
      mat_fas <-  mat_fas %>%
        rename('fas.date' = 'date') %>%
        rename('fas.mom' = 'mom')
      
 
  ### 3.3 Tidy repro_state       
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
      repro_state$hy.id <- as.character(repro_state$hy.id)    
      
    ## f) remove first numbering column
      repro_state <-  repro_state %>%
        subset(select = -c(x1))
      
      
  ### 3.4 Tidy hyenas
    ## a) rename data frame
      hyenas <- tblHyenas
      
    ## b) Convert all text to lower case
      hyenas <- AllCharactersToLower(hyenas)
      
    ## c) Format variable names (lowercase and separated by '.')
      hyenas <- FormatVarNames(hyenas)    
      
    ## d) Rename id as hy.id 
      hyenas <- hyenas %>%
        rename('hy.id' = 'id') 
      
      
  ### 3.5 Tidy life_hist
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
      
      
  ### 3.6 tidy NA values  
    ## a) Create a vector of values representing NA
      na_strings <- c('NaN', 'NA', 'na')
      
    ## b) Use naniar to replace all values representing NA with actual NA
      life_hist <- life_hist %>%
        replace_with_na_all(condition = ~.x %in% na_strings)
      
    ## c) Convert all '...error' variables from character into numeric
      life_hist <- life_hist %>%
        mutate_at(vars(contains('error')), funs(as.numeric))
      
      
  ### 3.7 Clean global environment
    ## a) Remove extra tables/dataframes
      rm(list = ls(pattern = 'tbl'))
      
      
      
###############################################################################
##############     4. Join tables and re-tidy fas / luma data    ##############
###############################################################################    
             
  ### 4.1 Join Data sets
    ## a) Left join life_hist to mat_care_fas,
      mat_care_fas <- mat_fas %>% 
        left_join(life_hist, by = c('hy.id' = 'hy.id'))
  
    ## b) Left join hyenas to mat_care_fas,
      mat_care_fas <- mat_care_fas %>% 
        left_join(select(hyenas, c(hy.id, sex, status, mom, dad, 
                                   number.littermates, litrank)),
                  by = c('hy.id' = 'hy.id'))  
      
    ## c) A list of hyenas having FAS data  
      hys <- unique(mat_care_fas$hy.id)
      # the number of those hyenas
      length(hys)

      
  ### 4.2 Tidy mat_care_fas
    ## a) Create a varialbe, 'fas.age.mon' which is an estimated age in 
      # months by subtracting dob.date from fas.date using lubridate
      mat_care_fas <- mat_care_fas  %>%
        mutate(fas.age.mon = round((interval(dob, 
                                         fas.date) %/% days(1)/ 30.44), 1))
  
    ## b) Extract month from fas.date
      # Use lubridate to extract the month during which a fas sample was 
      # collected and make a new variable  
      mat_care_fas$fas.mon <- month(mat_care_fas$fas.date)
      
    ## c) Create a varialbe, 'migratn.seas.fas,' which indicates if a fas 
      # sample was collected in migration (June 1 - Oct 31)
      mat_care_fas <- mat_care_fas  %>%
        mutate(migratn.seas.fas = ifelse(fas.mon >= 6 & fas.mon<= 10, 
                                         'migration', 'none'))
      
    ## d) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of migratn.seas.fas variable and sets the reference  
      # level to 'none'
      mat_care_fas <- transform(mat_care_fas,
                            migratn.seas.fas = factor(migratn.seas.fas,
                                              levels = c("none", 
                                                         "migration")))  
      
    ## e) Extract am vs. pm from stop.fas time
      # Use lubridate to extract the month during which a fas sample was 
      # collected and make a new variable
      mat_care_fas <- mat_care_fas  %>%
        mutate(fas.am.pm = ifelse(lubridate::am(stop.fas), 'am', 
                                     'pm'))
      
    ## f) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of fas.am.pm variable and sets the reference  
      # level to 'am' 
      mat_care_fas <- mat_care_fas %>%
        transform(fas.am.pm = factor(fas.am.pm,
                                     levels = c("am","pm"))) 
      
    ## g) Extract year from dob
      # Use lubridate to extract the year during which a hyena was born 
      # and make a new variable  
      mat_care_fas$dob.yr <- year(mat_care_fas$dob)
      
    ## h) Subset mat_care_fas based on 12 month (including mon. 12) cut-off  
      # for FAS samples collected prior to darting.data
      mat_care_fas <- mat_care_fas  %>%
        filter(fas.age.mon < 13)
    #*** NOTE *** This is a data inclusion cut-off decision. 
          # ALL FAS ocurr before darting date (LUMA sample) and when 
          # cub ages were less than 13 months 
      
    ## i) Subset mat_care_fas based on minimum overlap time of 5 min 
      # for FAS samples when both mom and cub present together
      mat_care_fas <- mat_care_fas  %>%
        filter(overlap >= 5)
    #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL FAS in which there is a minimum overlap time of 5 min 
      # for FAS samples when both mom and cub present together 
         
    ## j) Update number.littermates (a new 2 level nomial variable)
      mat_care_fas <- mat_care_fas  %>%
        mutate(number.littermates = case_when(number.littermates == 0 
                                              ~ c("single"),
                                              number.littermates == 1 
                                              ~ c("twin")))
      
      
  ### 4.3 Combine repro_state w/ mat_care_fas
    ## a) Make an empty data frame
      # This is an empty data frame that can store the overlaping mat_care_fas
      # and repro_states data
      fas_repro_data  <- c()   
     
    ## b) For loop to find overlap
      # Iterate over luma_fas_data (hy.id and fas.date), to find interseciton
      # with repro_state data
      for (i in 1:nrow(mat_care_fas)) { 
        
        # loop through 1:n IDs in mat_care_fas
        mom = paste (mat_care_fas$mom[i])  
        
        # loop through 1:n dates in mat_care_fas
        fas.date <-(mat_care_fas$fas.date[i])
      
        # create a dataframe to store the rows from repro_states where the
        # mom matches hy.id, and fas.date is in between cycle start and stop 
        overlap_fas_repro <- filter(repro_state, mom == hy.id & 
                                       fas.date >= cycle.start & 
                                       fas.date <= cycle.stop)
        
        # Control flow
          # if there is no id match and date overlap, 
          # then go to next loop iteration in mat_care_fas
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
      mat_care_fas <- mat_care_fas  %>%
        left_join(select(fas_repro_data, c(mom, fas.date, state, 
                                             cycle.start, cycle.stop, 
                                             trimester, parity, cub1, cub2)),
                  by = c("mom" = "mom",
                         "fas.date" = "fas.date")) %>%
        distinct(.keep_all = F) # add this to remove duplicate rows caused by 
                                # joining data where key is duplicated in 
                                # left table

      
  ### 4.4 Re-tidy FAS data
    ## a) Change state to character
      mat_care_fas$state <- as.character(mat_care_fas$state)
      
    ## b) Subset data to only include moms and cubs when possibly lactating
      mat_care_fas <- mat_care_fas  %>%
        filter(!grepl('o', state))
    #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL FAS in which mom's are not 'other' (...when no chance of lactating)
      
    ## c) Re-code *nominal* factor (with ordered levels)  
      mat_care_fas <- transform(mat_care_fas, 
                             state = factor(state))
  
    ## d) Re-code parity as a two level variable
      # primiparous = prim v.s. multiparous = mult
      mat_care_fas <- mat_care_fas  %>%
        mutate(parity.binary = ifelse(mat_care_fas$parity == 'a', 
                                      'prim', 'mult'))
       
    ## e) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of parity variable and sets the reference level 
      # to 'a'/ 'primiparous'
      mat_care_fas <- transform(mat_care_fas, 
                             parity.binary = factor(parity.binary,
                                              levels = c("prim", "mult")))
    # *** NOTE ***  this will (April 4, 2019) exclude serena hyenas,
      # which are not currently in tblReprostates
      
    ## f) Rename all variables from repro_states to reflect that these
      # variable pertain to the mom...not the cub / hy.id
      mat_care_fas <- mat_care_fas %>%
        rename('mom.state' = 'state') %>%
        rename('mom.cycle.start' = 'cycle.start') %>%
        rename('mom.cycle.stop' = 'cycle.stop') %>%
        rename('mom.trimester' = 'trimester') %>%
        rename('mom.parity' = 'parity') %>%
        rename('mom.cub1' = 'cub1') %>%
        rename('mom.cub2' = 'cub2') %>%
        rename('mom.parity.binary' = 'parity.binary')
      
    ## g) A list of hyenas having FAS data  
      hys <- unique(mat_care_fas$hy.id)
      # the number of those hyenas
      length(hys) 
      
    ## h) Extract year from fas.date
      # Use lubridate to extract the year during which an FAS was done
      mat_care_fas$fas.yr <- year(mat_care_fas$fas.date)
      
    ## i) Restrict the mat_care_fas data set to match rrbs data
      # Subset the mat_care_fas data based on overlapping date ranges starting
      # with the earliest dob of animals in the RRBS data set...in this case 
      # it is 2007 but FAS data are 2012-2013
      mat_care_fas <- mat_care_fas  %>%
        filter(fas.yr > 2006)
      
      #*** NOTE *** This is a data inclusion cut-off decision. 
      # ALL FAS which occur in an overlapping period with RRBS data
      
      
      
###############################################################################
##############              5. Univariate analyses               ##############
###############################################################################      
      
  ### 5.1 Overview
    ## a) Generate summary / descriptive stats for all variables in the data frame
      summary(mat_care_fas)
      
    ## b) Across all fas trials during which a hyena was <= 1yr, summarize data
      # grouped by hy.id (number of fas trials, first and last fas date etc.) 
      mat_care_fas_summry <- mat_care_fas %>%
        arrange(stop.fas) %>%
        group_by(hy.id) %>%
        summarize(n.fas = sum(!is.na(hy.id)),
                  sex = first(sex),
                  avg.cub.mins = round(mean(cub.mins, 
                                         na.rm = T), 2),
                  sd.cub.mins = round(sd(cub.mins, 
                                     na.rm = T), 2),
                  avg.mom.mins = round(mean(mom.mins, 
                                            na.rm = T), 2),
                  sd.mom.mins = round(sd(mom.mins, 
                                         na.rm = T), 2),
                  avg.overlap = round(mean(overlap, 
                                            na.rm = T), 2),
                  sd.overlap = round(sd(overlap, 
                                         na.rm = T), 2),
                  avg.c.mins = round(mean(c, na.rm = T), 2),
                  sd.c = round(sd(c, na.rm = T), 2),
                  prop.c = round((avg.c.mins/avg.overlap), 2),
                  avg.n.mins = round(mean(n, na.rm = T), 2),
                  sd.n = round(sd(n, na.rm = T), 2),
                  prop.n = round((avg.n.mins/avg.overlap), 2),
                  avg.g.mins = round(mean(g, na.rm = T), 2),
                  sd.g = round(sd(g, na.rm = T), 2),
                  prop.g = round((avg.g.mins/avg.overlap), 2),
                  first.fas.date = first(fas.date),
                  last.fas.date = last(fas.date),
                  fas.format = first(fas.format),
                  first.fas.age.mon = first(fas.age.mon),
                  last.fas.age.mon = last(fas.age.mon),
                  first.fas.yr = first(fas.yr),
                  # sample.id = first(sample.id),
                  # darting.date = first(darting.date),
                  # dart.age.mon = first(dart.age.mon),
                  dob = first(dob),
                  dob.yr = first(dob.yr),
                  dob.event.status = first(dob.event.status),
                  dob.event.data = first(dob.event.data),
                  dob.error = first(dob.error),
                  dob.notes = first(dob.notes),
                  dfs = first(dfs),
                  dfs.event.status = first(dfs.event.status),
                  dfs.event.data = first(dfs.event.data),
                  dfs.error = first(dfs.error),
                  dfs.notes = first(dfs.notes),
                  dengrad = first(dengrad),
                  dengrad.event.status = first(dengrad.event.status),
                  dengrad.event.data = first(dengrad.event.data),
                  dengrad.error = first(dengrad.error),
                  dengrad.notes = first(dengrad.notes),
                  weaned = first(weaned),
                  weaned.event.status = first(weaned.event.status),
                  weaned.event.data = first(weaned.event.data),
                  weaned.error = first(weaned.error),
                  weaned.notes = first(weaned.notes),
                  change.clan = first(change.clan),
                  change.clan.event.status = first(change.clan.event.status),
                  change.clan.event.data = first(change.clan.event.data),
                  change.clan.error = first(change.clan.error),
                  change.clan.notes = first(change.clan.notes),
                  change.clan2 = first(change.clan2),
                  change.clan2.event.status = first(change.clan2.event.status),
                  change.clan2.event.data = first(change.clan2.event.data),
                  change.clan2.error = first(change.clan2.error),
                  change.clan2.notes = first(change.clan2.notes),
                  disappeared = first(disappeared),
                  disappeared.event.status = first(disappeared.event.status),
                  disappeared.event.data = first(disappeared.event.data),
                  disappeared.error = first(disappeared.error),
                  disappeared.notes = first(disappeared.notes),
                  mom = first(mom),
                  dad = first(dad),
                  number.littermates = first(number.littermates),
                  litrank = first(litrank),
                  mom.parity = first(mom.parity),
                  mom.parity.binary = first(mom.parity.binary),
                  mom.cub1 = first(mom.cub1),
                  mom.cub2 = first(mom.cub2))
      
    ## c) Summarize proportions of time spent in close proximity, nursing,
      # and grooming
      mat_care_prop_summry <- mat_care_fas_summry %>%
        summarize(n.prop = sum(!is.na(hy.id)),
                  avg.prop.c = round(mean(prop.c, 
                                            na.rm = T), 2),
                  sd.prop.c = round(sd(prop.c, 
                                         na.rm = T), 2),
                  avg.prop.n = round(mean(prop.n, 
                                          na.rm = T), 2),
                  sd.prop.n = round(sd(prop.n, 
                                       na.rm = T), 2),
                  avg.prop.g = round(mean(prop.g, 
                                          na.rm = T), 2),
                  sd.prop.g = round(sd(prop.g, 
                                       na.rm = T), 2))
      
    ## d) Create an estimated duration in months over which period of time
      # fas data were collected (months between first and last fas date)
      mat_care_fas_summry <-  mat_care_fas_summry %>% 
        mutate(fas.durtn.mon = round(
          (interval(first.fas.date, last.fas.date) %/% 
             days(1)/ 30.44), 2))
                  
                  
  ### 5.2 Visualize (and transform as needed) raw data        
      ## a) Histogram Outcome (nos minutes in close proximity)
      ggplot(data=mat_care_fas, aes(x=c)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 35, by = 2.5), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 35)) +
        labs(title= "Histogram for mother and cub close proximity") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Close proximity (num. mins.)", y="Frequency") 
      
    ## b) Histogram Outcome (ratio in close proximity : overlap time together)
      ggplot(data=mat_care_fas, aes(x=(c/overlap))) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1, by = 0.01), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 1)) +
        labs(title= "Histogram for mother and cub close proximity
             vs overlap time together") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Ratio (close proximity : overlap", y="Frequency") 
  
  ### 5.3 Visualize (and transform as needed) raw data        
    ## a) Histogram Outcome (num minutes spent nursing)
      ggplot(data=mat_care_fas, aes(x=n)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 40, by = 2.5), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 40)) +
        labs(title= "Histogram for amount of minutes spent nursing") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Nursing (num. mins.)", y="Frequency") 
      
    ## b) Histogram Outcome (ratio nursing : overlap time together)
      ggplot(data=mat_care_fas, aes(x=(n/overlap))) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1, by = 0.01), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 1)) +
        labs(title= "Histogram for time spent nursing
             vs overlap time together") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Ratio (nursing : overlap", y="Frequency") 
      
      
  ### 5.4 Visualize (and transform as needed) raw data        
    ## a) Histogram Outcome (num minutes spent grooming)
      ggplot(data=mat_care_fas, aes(x=g)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 15, by = 1), 
                       col="black",
                       fill = "dark grey") +
        xlim(c(0, 15)) +
        labs(title= "Histogram for amount minutes spent grooming") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x="Grooming (num. mins.)", y="Frequency")   
    
      
    ## b) Histogram Outcome (ratio grooming : overlap time together)
      ggplot(data=mat_care_fas, aes(x=(g/overlap))) + 
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
                                   migratn.seas.fas +           # sampling var
                                   # sex +                        # mat care var
                                   # number.littermates +         # mat care var  
                                   # mom.parity.binary +          # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                   # below is hy.id nested within fas.mom
                                   # offset(log(overlap)) + (1|fas.mom/hy.id),
                                 data = mat_care_fas,
                                 ziformula = ~1,
                                 family = list(family = 'poisson', 
                                               link = 'log'))
      # Model summary estimates
      summary(close.zipoisson)
        # When generating BLUPs, control for FAS sampling variables, i.e. 
        # variable that can affect mom cub interction but not later life 
        # DNA methylation. 
        # After exponentiating, the beta estimates are differences in 
        # proportions of mins of behavior/total mins of overlap between 
        # groups (cat. variable) or per 1 unit change (cont. variable)?
         
      
      
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
                                 migratn.seas.fas +           # sampling var
                                 # sex +                        # mat care var
                                 # number.littermates +         # mat care var  
                                 # mom.parity.binary +          # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = mat_care_fas,
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
      close.intrcpt <- (fixef(close.zinegbinom2)[[1]])[[1]] # fixed effect
    #**** INTERPETATION ****#  
      # Interpet as expected log count of behavior when controlling for /
      # holding constant the effects ofcovariates (or if exponentiated 
      # the intercept is the incident rate of the expected count 
      # of behavior) as a proporiton of time overlap (the offset)
      
    ## e) extract fixed effect (intercept) from zero inflation model
      zi.close.intrcpt <- (fixef(close.zinegbinom2)[[2]])[[1]] # fixed effect
      
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
                                   migratn.seas.fas +           # sampling var
                                   # sex +                        # mat care var
                                   # number.littermates +         # mat care var  
                                   # mom.parity.binary +          # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = mat_care_fas,
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
                                 migratn.seas.fas +           # sampling var
                                 # sex +                        # mat care var
                                 # number.littermates +         # mat care var  
                                 # mom.parity.binary +          # mat care var
                                 offset(log(overlap)) + (1|hy.id),
                               data = mat_care_fas,
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
      ranef(nurse.zinegbinom1) # random effect
      fixef(nurse.zinegbinom1) # fixed effect
      coef(nurse.zinegbinom1) # BLUPs
      
    ## b) extract BLUPs from mixed model object
      nurse.blups <- as.data.frame(ranef(nurse.zinegbinom1)) # extract ranef as 
      # a dataframe, BLUPs = rand effects + intercept (from poiss/neg. binom)
      
    ## c) Rename variables in blups table
      nurse.blups <- nurse.blups %>%
        rename('hy.id' = 'grp') %>%
        rename('n.ranef' = 'condval') %>%
        rename('n.ranef.sd' = 'condsd') %>%
        select(c('hy.id', 'n.ranef', 'n.ranef.sd'))
      
    ## d) extract fixed effect (intercept) from poisson/neg. binomial model
      nurse.intrcpt <- (fixef(nurse.zinegbinom1)[[1]])[[1]] # fixed effect
    #**** INTERPETATION ****#  
      # Interpet as expected log count of behavior when controlling for /
      # holding constant the effects ofcovariates (or if exponentiated 
      # the intercept is the incident rate of the expected count 
      # of behavior) as a proporiton of time overlap (the offset)
      
    ## e) extract fixed effect (intercept) from zero inflation model
      zi.nurse.intrcpt <- (fixef(nurse.zinegbinom1)[[2]])[[1]] # fixed effect
      
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
                                   migratn.seas.fas +           # sampling var
                                   # sex +                        # mat care var
                                   # number.littermates +         # mat care var  
                                   # mom.parity.binary +          # mat care var
                                   offset(log(overlap)) + (1|hy.id),
                                 data = mat_care_fas,
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
                                 migratn.seas.fas +           # sampling var
                                 # sex +                        # mat care var
                                 # number.littermates +         # mat care var  
                                 # mom.parity.binary +          # mat care var
                                 offset(log(overlap)) + (1|hy.id),
                               data = mat_care_fas,
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
##############    8. Format variables for downstream analyses    ##############
###############################################################################      
  
  ### 8.1 Join the BLUPs for each behavior into a single data frame
    ## a) Left join close.blups to nurse.blups  
      mat_care_fas_blups <- left_join(close.blups,
                             nurse.blups, by = 'hy.id')
      
    ## b) Left join groom.blups to mat_care_fas_blups    
      mat_care_fas_blups <- left_join(mat_care_fas_blups,
                             groom.blups, by = 'hy.id')
      
    ## c) Left join mat_care_fas_summry to mat_care_fas_blups  
      mat_care_fas_blups <- left_join(mat_care_fas_blups,
                                      mat_care_fas_summry, by = 'hy.id')
      
      
  ### 8.2 Left join mat_care_fas by rrbs_vars
    ## a) Select variabels of importance from rrbs data
      rrbs_vars <- rrbs_vars %>%
        select(c(hy.id, sample.id, darting.date, 
                           dart.age.mon, mom.rank)) 
      
    ## b) Left join
      # retains all rows in mat_care_fas that has a match id in rrbs_vars   
      mat_care_rrbs <- rrbs_vars %>%
        left_join(mat_care_fas_blups, by = c("hy.id" = "hy.id"))
      
      
  ### 8.3 Tidy maternal care
    ## a) Create a variable, 'dart.from.fas,' which indicates whether fas data 
      # were collected before vs after the darting date
      mat_care_rrbs <- mat_care_rrbs  %>%
        mutate(dart.from.fas = case_when(mat_care_rrbs$first.fas.date >= 
                                           mat_care_rrbs$darting.date ~ 
                                         c("after"),
                                         mat_care_rrbs$first.fas.date < 
                                           mat_care_rrbs$darting.date ~ 
                                         c("before")))
      
    ## b) Create a varialbe, 'fec.vs.dart.days,' which indicates how many days
      # between the fas.date and darting date
      mat_care_rrbs <- mat_care_rrbs  %>%
        mutate(dart.from.fas.days = (mat_care_rrbs$darting.date - 
                                     mat_care_rrbs$first.fas.date)) 
      

  ### 8.4 Make dataframes for txt. file export to MACUA (for EWAS)
    ## a) Manual data clean up drop animals that failed ERRBS library QAQC
      mat_care_rrbs <- mat_care_rrbs %>%
        filter (!grepl('tato', hy.id))
      
    ## b) Arrange data set by hy.id abc order
      ## This is necessary because the order of rows matches data to a 
      ## specific animal
      mat_care_rrbs <- arrange(mat_care_rrbs, hy.id)
    
    ## c) Maternal rank formatted
      mat_rank_predictor_hy_n29 <- mat_care_rrbs %>%
        select(mom.rank) 
      
    ## d) Close proximity BLUPs formatted
      close_prox_predictor_hy_n29 <- mat_care_rrbs %>%
        select(c.expontd.blups) 
      
    ## e) Nursing BLUPs formatted
      nurse_predictor_hy_n29 <- mat_care_rrbs %>%
        select(n.expontd.blups) 
      
    ## e) Grooing BLUPs formatted
      groom_predictor_hy_n29 <- mat_care_rrbs %>%
        select(g.expontd.blups)   
      
      
  ### 8.5 Subset explanatory/predictor variables and format for use in MACAU
    ## a) Covariates (variables to control analyses)
      mat_rank_care_covars_hy_n29 <- mat_care_rrbs %>%
        mutate(intercept = 1) %>%
        #mutate(intercept = ifelse(is.na(c.expontd.blups), NA, 1)) %>%
        select(intercept, dart.age.mon, number.littermates, mom.parity.binary) 
      
    ## b) Fill in missing number.littermates
      mat_rank_care_covars_hy_n29$number.littermates <-
        ifelse((is.na(mat_rank_care_covars_hy_n29$number.littermates) &
                  !is.na(mat_rank_care_covars_hy_n29$mom.parity.binary)), 
               'single', mat_rank_care_covars_hy_n29$number.littermates)
      
      
      
      
      
###############################################################################
##############                9. Export data files               ##############
###############################################################################
      
  ### 9.1 Export FAS BLUPs to csv     
      # Save and export tables as a .cvs spreadsheet and named with today's
      # date. Files are saved in the 'output' folder in the working directory.
      
    ## a) Generate file name for mat_car_fas
      # this the full fas data set with repeated fas samples
      csv.file.name.FAS.mat.car <- paste0(here("data",
                                               "mat_care_fas.csv")) 
    ## b) Save mat_car_fas table as csv
      write.csv (mat_care_fas, file = csv.file.name.FAS.mat.car)
      
    ## c) Generate file name for mat_car_fas_blups
      # this is the data set summarized using blups and with covariates
      csv.file.name.FAS.BLUPs <- paste0(here("data",
                                             "mat_care_fas_blups.csv")) 
    ## d) Save mat_car_fas_blups table as csv
      write.csv (mat_care_fas_blups, file = csv.file.name.FAS.BLUPs) 
    
    ## e) Generate file name for mat_care_rrbs
      # includes fas data and covariates only for animals with rrbs data
      csv.file.name.mat.car.rrbs <- paste0(here("data",
                                               "mat_care_rrbs.csv"))
      
    ## f) Save mat_care_rrbs table as csv
      write.csv (mat_care_rrbs, file = csv.file.name.mat.car.rrbs)
      

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
      
      
    
