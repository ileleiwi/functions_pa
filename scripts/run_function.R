library(tidyverse)


setwd(paste0("~"))

#Data 
test <- read_tsv("data/test_data.tsv") 
#rules
rules <- read_tsv("data/sugar_scfa_resp_rules.tsv") %>%
  as.data.frame()

#functions
source("scripts/rule_functions.R")

#build checkRule function
checkRule <- makecheckRule(rules)

#generate dataframe
x <- evaluateCountsDf(test, rules)
