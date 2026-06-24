## ########################################################################## ##
## Extra. GPS data cleaning considerations                                    ##
##                                                                            ##
## Project: MoveSIN GPS data cleaning and exploration                         ##
##                                                                            ##
## Created by: Luke Emerson                                                   ##
## Created: 23rd June 2026                                                    ##
##                                                                            ##
## Edited by: Luke Emerson                                                    ##
## Edited: 23rd June 2026                                                     ##
## ########################################################################## ##


# Load packages ----------------------------------------------------------------
library(gridExtra)
library(grid)

# GPS data cleaning considerations --------------------------------------------- 
#Although standard cleaning procedures apply to most workflows using telemetry
#data.Telemetry data cleaning should be tailored to the biological question, the
#movement metric being estimated, and the assumptions of the analytical method.

grid.newpage()

grid.text(
  "GPS data cleaning considerations",
  x = 0.5, y = 0.85,
  gp = gpar(fontsize = 24, fontface = "bold")
)

grid.text(
  paste(
    "Although standard cleaning procedures apply to most workflows using telemetry data,",
    "cleaning decisions should be tailored to:",
    "",
    "• the biological question",
    "• the movement metric being estimated",
    "• the assumptions of the analytical method",
    sep = "\n"
  ),
  x = 0.5, y = 0.5,
  just = "center",
  gp = gpar(fontsize = 16)
)


# Cleaning steps required regardless of analysis -------------------------------
# 1. Remove observations with missing coordinates
# 2. Account for duplicate timestamps
# 3. Remove impossible timestamps
# 4. Remove observations with obvious coordinate errors

grid.newpage()

grid.text(
  "Cleaning steps required regardless of analysis",
  x = 0.5, y = 0.85,
  gp = gpar(fontsize = 24, fontface = "bold")
)

grid.text(
  paste(
    "1. Remove observations with missing coordinates",
    "2. Account for duplicate timestamps",
    "3. Remove impossible timestamps",
    "4. Remove observations with obvious coordinate errors",
    sep = "\n"
  ),
  x = 0.5, y = 0.5,
  just = "center",
  gp = gpar(fontsize = 16)
)


# Common spatial ecology workflows using telemetry data ------------------------
# Potential other analytical approaches to determine "what animals do", "how they
# move through landscape", "how they interact", "how much energy they expend".

title <- "Common spatial ecology workflows using telemetry data to explore space use, habitat selection, movement modelling, behavioural state, activity patterns"

# Create table data
cleaning_table <- data.frame(
  "Cleaning or preprocessing step" = c(
    "Remove missing coordinates (NA lat/long)",
    "Convert lat/long to projected coordinates (UTM)",
    "Sort locations chronologically by individual",
    "Remove duplicate timestamps within individual",
    "Check sampling interval / fix rate",
    "Remove large temporal gaps",
    "Remove unrealistic step lengths",
    "Remove unrealistic speeds",
    "Account for GPS error / fix quality",
    "Remove stationary collars / mortality periods",
    "Check GPS drift / jumps",
    "Bearing / turning angle filtering",
    "Regularise sampling interval",
    "Thin locations to independence",
    "Align CRS and extract environmental covariates",
    "Generate random locations",
    "Generate random steps",
    "Model movement process",
    "Assess model diagnostics"
  ),
  
  MCP = c(
    "Required","Recommended","Recommended","Required","Useful",
    "Sometimes","Recommended","Recommended","Useful",
    "Often","Recommended","No","Optional","Sometimes",
    "No","No","No","No","Minimal"
  ),
  
  KDE = c(
    "Required","Recommended","Recommended","Required","Useful",
    "Sometimes","Recommended","Recommended","Useful",
    "Recommended","Recommended","No","Optional","Sometimes",
    "No","No","No","No","Minimal"
  ),
  
  wAKDE = c(
    "Required","Recommended","Essential","Required","Essential",
    "Sometimes","Important","Important","Important",
    "Essential","Important","No","No","No",
    "No","No","No","Essential","Essential"
  ),
  
  RSF = c(
    "Required","Recommended","Essential","Required","Important",
    "Sometimes","Important","Important","Useful",
    "Essential","Important","Sometimes","Sometimes","Sometimes",
    "Essential","Essential","No","No","Essential"
  ),
  
  SSF = c(
    "Required","Recommended","Essential","Required","Essential",
    "Essential","Essential","Essential","Important",
    "Essential","Important","Essential","Essential","No",
    "Essential","No","Essential","Sometimes","Essential"
  ),
  
  HMM = c(
    "Required","Recommended","Essential","Required","Essential",
    "Essential","Recommended","Recommended","Important",
    "Essential","Important","Important","Essential","No",
    "Sometimes","No","No","Essential","Essential"
  ),
  
  "Activity pattern" = c(
    "Required","Recommended","Essential","Required","Essential",
    "Important","Useful","Useful","Sometimes",
    "Essential","Important","No","Recommended","No",
    "No","No","No","No","Important"
  )
)

# Create table
table_plot <- tableGrob(
  cleaning_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 12,
    padding = unit(c(5,5), "mm")
  )
)

# Create title
title_grob <- textGrob(
  title,
  gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  just = "center"
)

# Display
grid.newpage()

grid.arrange(
  title_grob,
  table_plot,
  heights = unit(c(0.08,0.92), "npc")
)


