library(tidyverse)


setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/functions_pa/"))


#Data 
picr_ko <- read_tsv("data/test_data.tsv") 
#rules
rules <- read_tsv("Clean_Data/picr_rules.tsv") %>%
  as.data.frame()

#functions
source("Scripts/rule_functions.R")

#build checkRule function
checkRule <- makecheckRule(rules)

#generate dataframe
x <- evaluateCountsDf(picr, rules)







##################################################################



##################################################################








#Rules

# Sulfate Reduction 
#   must have Dsr
#   Dsr - K11180|EC1.8.99.5|K11181|EC1.8.99.5
#   Asr - K00380|EC1.8.1.2|K00381|EC1.8.1.2|K00392|EC1.8.7.1
# 
# Fumarate Reduction
#   must have one of the following KOs
#   K00244,K00245,K00246,K00247
# 
# Aerobic
#   must have [50% CytochromeCOxidase OR 50% CytochromeAa3600MenaquinolOxidase OR 50% CytochromeOUbiquinolOxidase] AND ETC50
# 
#   ETC50
#     must have 50% ComplexIA OR 50% ComplexIB OR 50% ComplexIC
#     ComplexIA - K00330,[K00331&K00332&[K00333|K00331]&[K13378|K13380]],K00334,K00335,K00336,K00337,K00338,K00339,K00340,[K00341&K00342|K15863],K00343
#     ComplexIB - K05574,K05582,K05581,K05579,K05572,K05580,K05578,K05576,K05577,K05575,K05573,K05583,K05584,K05585
#     ComplexIC - K03945,K03946,K03947,K03948,K03949,K03950,K03951,K03952,K03953,K03954,K03955,K03956,K11352,K11353
#   
#   CytochromeCOxidase - K02275,K02274,K02276|K15408,K02277
#   CytochromeAa3600MenaquinolOxidase - K02827,K02826,K02828,K02829
#   CytochromeOUbiquinolOxidase - K02297,K02298,K02299,K02300
# 
# 
# Tetrathionate Reduction
#   must have [1 of ttrrs AND 1 of ttrabc] OR 2 of ttrabc
#   ttrabc - K08359,K08358,K08357
#   ttrrs - K13040|K13041
# 
# Microaerophillic
#   must have [50% CytochromeBDUbiquinolOxidase OR 50% CytochromeCOxidaseCBB3] AND ETC50
#     
#     CytochromeBDUbiquinolOxidase - K00425,K00426,K00424|K22501
#     CytochromeCOxidaseCBB3 - K00404|K00405|K15862,K00407,K00406
#     
# Denitrification (Not ETC)
#   must have [NitrateReductase OR NitriteReductase OR NitricOxideReductase OR NitrousOxideReductase] AND does not have ETC50
#     NitrateReductase - K00370|K00371|K00374|K02567|K02568|EC1.7.5.1|EC1.9.6.1
#     NitriteReductase - K00368|K15864|EC1.7.2.1
#     NitricOxideReductase - K04561|K02305|EC1.7.2.5
#     NitrousOxideReductase - K00376|EC1.7.2.4