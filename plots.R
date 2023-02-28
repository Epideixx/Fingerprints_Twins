##############################
###       PLOTS         ######
##############################

library(ggseg)
library(ggplot2)
library(patchwork)
library(ggsegDesterieux)
library(dplyr)
library(colorspace)

########## ICC ##########

ICC <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/ICC_and_Heritability/ICC_avg.csv')
ICC_final <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Destrieux_atlas.csv')

ICC_final_DELTA <- ICC_final
ICC_final_DELTA$icc = ICC$DELTA
ICC_final_DELTA$freqband = "Delta"

ICC_final_THETA <- ICC_final
ICC_final_THETA$icc = ICC$THETA
ICC_final_THETA$freqband = "Theta"

ICC_final_ALPHA <- ICC_final
ICC_final_ALPHA$icc = ICC$ALPHA
ICC_final_ALPHA$freqband = "Alpha"

ICC_final_BETA <- ICC_final
ICC_final_BETA$icc = ICC$BETA
ICC_final_BETA$freqband = "Beta"

ICC_final_GAMMA <- ICC_final
ICC_final_GAMMA$icc = ICC$GAMMA
ICC_final_GAMMA$freqband = "Gamma"

ICC_final_HIGH_GAMMA <- ICC_final
ICC_final_HIGH_GAMMA$icc = ICC$HIGH.GAMMA
ICC_final_HIGH_GAMMA$freqband = "High Gamma"

ICC = rbind(ICC_final_DELTA, ICC_final_THETA, ICC_final_ALPHA, ICC_final_BETA, ICC_final_GAMMA, ICC_final_HIGH_GAMMA)

ICC$freqband= factor(ICC$freqband, levels = c("Delta", "Theta","Alpha", "Beta", "Gamma", "High Gamma"))


###### Heritability #######

Heritability <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/ICC_and_Heritability/Heritability_avg.csv')
Heritability_final <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Destrieux_atlas.csv')

Heritability_final_DELTA <- Heritability_final
Heritability_final_DELTA$h = Heritability$DELTA
Heritability_final_DELTA$freqband = "Delta"

Heritability_final_THETA <- Heritability_final
Heritability_final_THETA$h = Heritability$THETA
Heritability_final_THETA$freqband = "Theta"

Heritability_final_ALPHA <- Heritability_final
Heritability_final_ALPHA$h = Heritability$ALPHA
Heritability_final_ALPHA$freqband = "Alpha"

Heritability_final_BETA <- Heritability_final
Heritability_final_BETA$h = Heritability$BETA
Heritability_final_BETA$freqband = "Beta"

Heritability_final_GAMMA <- Heritability_final
Heritability_final_GAMMA$h = Heritability$GAMMA
Heritability_final_GAMMA$freqband = "Gamma"

Heritability_final_HIGH_GAMMA <- Heritability_final
Heritability_final_HIGH_GAMMA$h = Heritability$HIGH.GAMMA
Heritability_final_HIGH_GAMMA$freqband = "High Gamma"

Heritability = rbind(Heritability_final_DELTA, Heritability_final_THETA, Heritability_final_ALPHA, Heritability_final_BETA, Heritability_final_GAMMA, Heritability_final_HIGH_GAMMA)

Heritability$freqband= factor(Heritability$freqband, levels = c("Delta", "Theta","Alpha", "Beta", "Gamma", "High Gamma"))


Heritability %>%
  group_by(freqband) %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = h)) +
  geom_sf(show.legend = TRUE) +
  facet_wrap( ~ freqband) + scale_fill_binned_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1) +
  theme_void() + scale_color_manual('white')






