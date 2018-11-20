###############################################################################
##############        Spotted Hyena eRRBS DNA Methylation:       ##############
##############               Adult stress phenotype              ##############
##############                 By: Zach Laubach                  ##############
##############             last updated: 20 Nov 2018             ##############
###############################################################################


### PURPOSE: This code is desingned to analyze ERRBS data and associations of 
# genome wide DNA methylation and adult stress phenotype in hyenas.


# Code Blocks
# 1: Configure workspace
# 2: Import data
# 3: Data management
# 4: Univariate analyses
# 5: Data transformations 
# 6: Bi-variate analyses 
# 7: Re-tidy data for analyses
# 8: Cub models
# 9: Subadult models
# 10: Adult models



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

# Check for here and install if not already installed
if (!'here' %in% installed.packages()[,1]){
  install.packages ('here')
}
# load here packages
library ('here')