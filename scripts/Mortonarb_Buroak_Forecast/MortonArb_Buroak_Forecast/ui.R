#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Phenology forecasting app
# Purpose: This script is the ui side of the shiny app
# Inputs: Weather_ArbCOOP_historical_latest.csv created by M1_Meterology_download.R
#         Weather_Arb_forecast_ensemble_latest.csv created by M1_Meterology_download.R
#         Species_Name_Catalogue.csv created by SP_Names.R
# Outputs: 
# Notes: The first half of this script is taken from M2_Meteorology_Graphing_Bayes.R
#-----------------------------------------------------------------------------------------------------------------------------------#
library(shiny)
library(ggplot2)
library(plotly)
library(stringr)
library(shinyWidgets)
# -------------------------------------
# Load in the data and calculate ensemble spread
# -------------------------------------
path.in <- "data/"

#Reading in the list of species names
mod.catalogue<- read.csv(file.path( path.in, "Model_List.csv"))

#Reading in the list of available past forecasts: This system will only work for one year this time. I should split it by year into different folders down the line
fc.df <- read.csv(file.path(path.in, "Old_Forecast_List.csv"))

#This puts everything on a single page where I can specify the organizaiton easily
fluidPage(
  
  titlePanel("The Morton Arboretum: Spring Oak Budburst Forecast"),
  
  #This defines what will be in the first row of the page
  #Currently that is our explanation of what the app is
  fluidRow(
    p(h2("     This app provides a visualization of the predicted dates of the budburst phenophase for Bur oak species in the common garden experiment at The Morton Arboretum utilizing data from", 
         tags$a(href="https://www.usanpn.org/usa-national-phenology-network", "the National Phenology Network"))),
    p(h3(strong("YDAY"), " is fit using the day of year of previous observations. It functions as a null model")),
    p(h3(strong("GDD5"), " is fit using the cumulative growing degree days (daily mean temp - 5C) of previous observations. This is the most common and simple phenology model")),
    p(h3(strong("PTTGDD5")," is fit using the cumulative photothermal growing degree days (GDD5 * daylength) of previous observations")),
    p(h3(strong("GTmean"), " is a linear model that fits an intercept and slope using the mean daily temperature from Jan 1st to March 31st as the independent variable")),
  ),
  

  fluidRow(br()), #This line is exclusively to give white space between the description and the options
  
  #The second row is the naming convention selection and previous forecast slider
  fluidRow(
    #Allowing the choice between scientific and common: Maybe remove given the audience and functionality quirks
    column(width = 4, pickerInput('Model','Choose a Model: ', choices = c(mod.catalogue$Model.name), selected= "GDD5", options = list('live-search' = FALSE), multiple = T)),
    
    #Allowing for a slider between the forecasts
    column(width = 4, sliderTextInput("Forecast date", "Previous forecasts", choices=fc.df$Date, selected = as.character(max(fc.df$Date))))),

  
  #Fifth row
  fluidRow(
    #The graphs themselves. The only show up when the submit button is pressed
    plotOutput("weather")),
  
  #Fifth row
  fluidRow(
    #The graphs themselves. The only show up when the submit button is pressed
    uiOutput("plot.thresh.ui")),
  
)

