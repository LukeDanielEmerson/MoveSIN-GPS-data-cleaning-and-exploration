## ########################################################################## ##
## 01. Directory creation                                                     ##                                 
##                                                                            ##
## Project: MoveSIN GPS data cleaning and exploration                         ##
##                                                                            ##
## Created by: Luke Emerson                                                   ##
## Created: 25th April 2026                                                   ##
##                                                                            ##
## Edited by: Luke Emerson                                                    ##
## Edited: 19th June 2026                                                     ##
## ########################################################################## ##

# Check working directory is correct
getwd()

# Define the folders to create -------------------------------------------------
folders <- c("Scripts", 
             "Telemetry_data", 
             "RDS", 
             "Spatial_layers", 
             "Figures", 
             "Outputs", 
             "Literature", 
             "Documentation")

# Create folders if they don't already exist
lapply(folders, function(x) if(!dir.exists(x)) dir.create(x))
