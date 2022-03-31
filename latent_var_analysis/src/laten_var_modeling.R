################################################################################
############### Tsinghua Computional Social Science R Workshop 2022 Spring
############### Latent Variable Modeling
############### Wenquan Wu
############### 2022-04-01 
################################################################################
# load ----
library(tidyverse)
load('latent_var_analysis/data/wvs_wave7.rdata')
data <- `WVS_Cross-National_Wave_7_R_v3_0`

dat_eco <- data %>% 
  select(
    Q106,
    Q107,
    Q108,
    Q109,
    Q110
  )

# Principal components analysis (PCA) ----
  


# Structural Equation Modeling (SEM) ----

# Item Response Theory (IRT) ----