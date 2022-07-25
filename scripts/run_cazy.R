library(tidyverse)
library(lubridate)

setwd(paste0("~"))
##Read in data
#DRAM function heatmap form
function_heatmap_form_cazy <- read_tsv("data/cazy_rules.tsv") 
#dbCAN ids in MQHQ bins
dbcan_mqhq <- read_tsv("data/dbcan_mqhq.tsv") 

#apply DRAM SOM rules for function presence absence on dbcan_mqhq
#using function heatmap form
#must include an oligo cleaving and a backbone cleaving for TRUE

#get unique list of polymers
substrate_names <- function_heatmap_form_cazy %>%
  pull(function_name) %>%
  str_remove(" \\(.+\\)") %>%
  str_remove(" Cleavage| cleavage") %>%
  str_remove(" Oligo| oligo") %>%
  str_remove(" Backbone| backbone") %>%
  unique()



dram_product <- buildProduct(substrate_names,
                             function_heatmap_form_cazy,
                             dbcan_mqhq)


