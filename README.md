# Bur oak Forecasting
Scripts for Ecological forecasting model for phenology at the arboretum

# Script Folders

Each of these folders has their own README inside

## Mortonarb_Buroak_Forecast

Purpose: A folder for running our shiny app

# Scripts

## 1_NPN_download.R 

Purpose: This script downloads bur oak NPN data, cleans them, and subsets by proximity to source populations

Inputs: Mother Tree Locations csv which contains the location of all mother trees

Outputs: Filtered_NPN.csv which contains NPN bur oak observations with a column to flag if they are in the source area or not

         Filtered_NPN_Sum_Stats.csv which is a summary table of the NPN bur oak observations
         
         Bur_Oak_Source_Area_NPN.pdf which contains an image of all NPN bur oak observations and their proximity to source
         
Notes: Creating the spatial hull for OK requires a manual adjustment. It currently works but is an oddity and could be an issue down the line

## 2_Source_mapping.R

Purpose: This script creates figures showing the mother tree source populations locations and NPN observations

Inputs: Filtered_NPN.csv which contains NPN bur oak observations with a column to flag if they are in the source area or not

         Mother Tree Locations csv which contains the location of all mother trees
         
Outputs: NPN_Bur_Oak_Observation_Locations.png which shows all bur oak observations in NPN

          Bur_Oak_Source_Area_NPN_points.pdf which shows the convex hull around mother trees and nearby NPN points
          
          Bur_Oak_Source_Area_Expanded_NPN_Range.pdf whcih shows the extended range of all the mother tree hulls
          
Notes: This script only produces figures and isn't required in the overall workflow

## 3_NPN_weather_matching.R

Purpose: This script takes our filter NPN data and matches them with weather metrics

Inputs: Filtered_NPN.csv which contains NPN bur oak observations with a column to flag if they are in the source area or not

Outputs: Full_Bur_Obs.csv which contains all npn observations nearby source populations matched with weather metrics

         Daymet_clean_data.csv which contains yearly weather data for each site of observation
         
## 4_Model_Fitting.R

Purpose: This script takes our filtered NPN data and fits our different model paramters.

Inputs: Full_Bur_Obs.csv which contains all npn observations nearby source populations matched with weather metrics

         Daymet_clean_data.csv which contains yearly weather data for each site of observation
         
Outputs: Parameter distributions for models

Notes: Current models are YDAY null model, GDD5, PTTGDD5, linear GTmean

## Output_Visualization.R

Purpose: This script makes intial visuals for the model output.

Inputs:
         
Outputs: 

Notes: No longer used now that shiny app is working

## Test_models.R

Purpose: This script takes our filtered NPN data and tries to fit more complex models

Inputs: Full_Bur_Obs.csv which contains all npn observations nearby source populations matched with weather metrics

        Daymet_clean_data.csv which contains yearly weather data for each site of observation
        
Outputs: Parameter distributions for models

# Functions

## weather_calc.R

Purpose: This function serves to calculate weather statistics of interest. 

Currently Growing degree days at 5C, Growing degree days at 0C, Number of chill days, and Growing season mean temperature
         
Inputs: Lubridate package

Outputs:This function will take a data frame of daily weather data and produce the following summary statistics

GDD5 = Growing degree days at 5 degrees C 

GDD0 = Growing degree days at 0 degrees C

CDD0 = Chilling degree days at 0 degrees C

NCD = Number of chilling days 

GTmean = Growing season mean temperature from Jan 1st to March 31st

PTTGDD = Growing degree day at 5 degrees C * the amount of light per day
  
Notes: The defaults for this funcion are

Starting year of interest                       y_start = 1975

Ending year of interest                         y_end = 2019

Julian yday for start of growing season         g_start = 1

Julian yday for end of growing season           g_end = 90

## met_download_CFS.R

## met_download_GHCN.R

## met_gapfill.R
