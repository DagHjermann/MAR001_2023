# MAR001_2023

Scripts for analysis of levels and trends of marine contaminants at the European level for ETC, 2023  
  
Note: the scripts are based on the 2021 version  
- 
  

## Script overview    

* The order of scripts is as below (note that it doesn't strictly follow the script numbers):    
    - 01 - read raw data from the two sources    
    - 03 - based on data from 01, remove duplicated data (data that are in both sources) and make 'cleaned' data   
    - 04 - based on data from 03, make per-year medians  
    - 05 - based on data from 04 and 07, perform time series analysis  
    - 07 - based on data from 04, perform status classification    
    - 06 - based on results from 05 and 07, make output tables and graphs  
    
* Function scripts (contains functions used by other scripts - should not normally need to be opened)   
    - 02_Fix_station_duplicates_functions.R (used by script 03)   
    - 05_Time_trend_regression_functions.R (used by script 05)  
    
* Input data  
    - Folder Input_data
    
* Other files
    - Thresholds/MGR_05_df_thresholds.xlsx - used by scr 07, contains thesholds (EQS etc) used for classification (low/medium/high concentration)    
    
* Major changes from the 2021 scripts:  
    - script 01: There should now be Portuguese data after 2010 -> Script just check that there are PT data after 2010  
    - script 03: TO DO  

## Scripts, details      

* 01: Read raw data from 'Input_data'   

    - two data sources, ICES and EMODnet  
    - note that we read from two different EMODnet versions, a and b, see 'readme' in 'Input_data'. 
    We save both at the end of 01, but end up using a.
    - columns fixed: Country and UNIT (converts all concentrations to ug/kg)   
    - keeps only records with parameters that are in both data sets
    - removes TBTIN when TBSN+ is availiable (part 3, 94 records)
    - a few records (ca 100) are deleted for having strange units (part 5)  

* 03: Reads data from 01, checks duplication, and produces data for script 04

    - Small fixes in the start: cleans one huge Hg value, fixes unit of fat weight
    - Select data from bivalves of 8 species (part 3 + 4a): 304546 rows out of a total of 1553527 (dat_2)  
    - dat_3: remove ICES duplicates  
    - dat_4: EMODnet data summarized per station/year  
    - dat_5: Add PCB7 + PCB6 to data 
    - dat_6: Delete French spring data   
    - dat_7: Add species-specific DW and fat weight to data   
    - Part 7: For duplicates, the data are split in two, those with common station name (in both ICES and Emodnet) and those not  
    - For the ones with duplicates, there is a job to pick data from the different years and in the right basis  
    - Afterwards, data with and without duplicates are put together again
    
    *) Species used: Cerastoderma edule (common cockle), Mya arenaria (softshell clam), Ruditapes philippinarum (manila clam), 
    Mytilus edulis (blue mussel), Mytilus galloprovincialis (Mediteranean mussel), Crassostrea gigas (Pacific oyster), 
    Ostrea edulis (native oyster), Macoma balthica (Baltic clam)
      
* 04: Reads data from 03, makes per-year medians, and saves it for script 05 and 07  

* 05: Reads data from 04 and 07, perform regression, and saves it for script 06
    
    - Most of the script is for the regressions, then dat_status_2 is added in order to get 
    combined status + trend for each station  
    - Series with no data < LOQ: Ordinary regression
    - Series with some data (annual medians) < LOQ: Bayesian modelling using JAGS  
    - Criteria for trend analysis (see part 3): 1. Only data after 2005 used, 2. The end of the time series should be in 2015 or later
    3. The length of the time series should be minimum 5 and maximum 10 years   

    - As the script is now, it seems to need 'dat_status_2' in memory... (from script 07)

* 06: Reads data from 05 and 07 (especially) and makes all output tables and graphs used to write assessment text  
    
    - The calculation of 'relative class' is done in this script   
    - There is also a version '06_MAR001_results_2021_BEFORE_UPDATE_2022-06-17.Rmd' using the old, erronous data from scr 07

* 07: Reads data from 04, perform status classification, combines the result and produces data for script 05

    - Criteria for status calculation: needs at least 3 years of data since 2010  
    - Note: for the first version of the assessment, this criteria was not implemented properly   
    - This was corrected on 2022-06-17 (commit 1426b98)   
    - For comparison, the old saved file from this script was saved as '07_dat_status_trend_BEFORE_UPDATE_2022-06-17.rds'  
    - For comparison, this is used in the result file '06_MAR001_results_2021_BEFORE_UPDATE_2022-06-17.Rmd'  


