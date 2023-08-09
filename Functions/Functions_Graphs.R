#---------------------------------------------------------------------------------------------------------------
# Theme for the graphs 
# Last modified on 27-06-2023
#---------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(ggthemes)
library(extrafont)

## Creation of a theme for the graphs

theme_graphs <- function() {
  theme(plot.title = element_text(family="serif", size = 24, face = "bold", 
                                  hjust=0.5, margin = margin(0, 0, 10, 0)),
        plot.subtitle = element_text(family="serif", colour = "#000000", size = 14, 
                                     hjust=0.5, margin = margin(0, 0, 10, 0)),
        plot.caption = element_text(family="serif", colour = "#000000", size = 14, 
                                    hjust=1, margin = margin(10, 0, 10, 0)),
        plot.background = element_rect(fill = "#FFFFFF"),
        plot.margin = margin(0.2, 0.5, 0.2, 0.5, "cm"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.background = element_rect(fill = "#FFFFFF", colour = "#000000", linetype = "solid"), 
        panel.grid.major.x = element_line(colour = "#C9C9C9", linetype = "dotted"),
        panel.grid.major.y = element_line(colour = "#C9C9C9", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(family="serif", face="bold", size=26),
        strip.background = element_rect(fill="#FFFFFF", colour="#000000", linewidth = 0.5),
        axis.title.x = element_text(family="serif", size = 24, colour = "#000000", 
                                    hjust=0.5, face = "bold", margin = margin(10, 0, 0, 0)), 
        axis.title.y = element_text(family="serif", size = 24, colour = "#000000", 
                                    face = "bold", margin = margin(0, 10, 0, 0)), 
        axis.text = element_text(family="serif", size = 20, colour = "#000000"),
        axis.line.y = element_line(colour = "#000000"),
        axis.line.x = element_line(colour = "#000000"),
        axis.ticks = element_line(colour = "#000000", linewidth = 1),    
        legend.title = element_text(family="serif", size = 20, colour = "#000000", face = "bold"),
        legend.text = element_text(family="serif", size = 20, colour = "#000000"),
        legend.background = element_rect(fill = "#FFFFFF", colour = "#000000", linewidth = 0.3, linetype = "solid"), 
        legend.key = element_rect(fill = NA), 
        legend.position = "bottom", 
        legend.direction = "horizontal"
  )
}

# Convert SOCSIM months to calendar years. 
asYr <- function(month, last_month, final_sim_year) {
  return(final_sim_year - trunc((last_month - month)/12))
}