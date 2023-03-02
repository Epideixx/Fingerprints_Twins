##############################
###       PLOTS         ######
##############################

library(ggseg)
library(ggplot2)
library(ggridges)
library(patchwork)
library(ggsegDesterieux)
library(dplyr)
library(colorspace)
library(cowplot)


####### DISTRIBUTIONS CORR PER BAND ########

corr <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/PSD_correlations/All_correlation_every_freq.csv')
corr <- corr[,2:length(corr)]

autocorr <- corr[,1:3]
autocorr <- stack(autocorr)
autocorr <- na.omit(autocorr)
autocorr$ind <- recode_factor(autocorr$ind, 'Autocorr_MZ_BROADBAND' = 'MZ', 'Autocorr_DZ_BROADBAND' = 'DZ', 'Autocorr_NT_BROADBAND' = 'NT' )


crossocorr <- corr[,4:6]
crossocorr <- stack(crossocorr)
crossocorr <- na.omit(crossocorr)
crossocorr$ind <- recode_factor(crossocorr$ind, 'Crosscorr.MZ_BROADBAND' = 'MZ', 'Crosscorr.DZ_BROADBAND' = 'DZ', 'Crosscorr.NT_BROADBAND' = 'NT' )

# For Broadband
p1 <- ggplot(autocorr, aes(x = values)) +
  geom_histogram() +
  facet_grid(ind ~ .) +
  xlim(0.85, 1.00)

p2 <- ggplot(crossocorr, aes(x = values)) +
  geom_histogram() +
  facet_grid(ind ~ ., scales = "free")+
  xlim(0.7, 0.95)

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

###
corr <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/PSD_correlations/All_correlation_every_freq.csv')
corr <- corr[,2:length(corr)]

corr <- corr[,1:6]
corr <- stack(corr)
corr <- na.omit(corr)

corr$ind = as.character( corr$ind )
corr$corrtype = sapply(strsplit(corr$ind, "_|\\."), '[', 1)
corr$Group = sapply(strsplit(corr$ind, "_|\\."), '[', 2)

corr$corrtype <- recode_factor(corr$corrtype, 'Autocorr' = 'Autocorrelation', 'Crosscorr' = 'Cross-Correlation')
corr$Group = factor(corr$Group, levels = c("NT", "DZ", "MZ"))

corr %>%
  ggplot(aes(x=values,y=Group, fill=Group)) +
  geom_density_ridges()+
  theme(legend.position = "none")+
  facet_wrap(~corrtype, scale = "free_x")+
  xlab("Correlation")+
  ylab("Twin Category")


####### PLOT CORRELATIONS PER ROI ########

corr <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/PSD_correlations/Dataset_1_VS_Dataset2/ROI_correlations.csv')
corr <- corr[2:149, ]

corr_final <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Destrieux_atlas.csv')

corr_MZ <- corr_final
corr_MZ$corr <- corr$Crosscorr.MZ
corr_MZ$typeTwin <- "Crosscorrelation MZ"

corr_DZ <- corr_final
corr_DZ$corr <- corr$Crosscorr.DZ
corr_DZ$typeTwin <- "Crosscorrelation DZ"

corr_final <- rbind(corr_MZ, corr_DZ)

corr_final$typeTwin <- factor(corr_final$typeTwin, levels = c("Crosscorrelation MZ", "Crosscorrelation DZ"))

corr_final %>%
  group_by(typeTwin) %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = corr)) +
  geom_sf(show.legend = TRUE) +
  facet_wrap( ~ typeTwin) + scale_fill_continuous_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = "Correlation") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/Corr_ROI_dataset_1_2_continuous.pdf', device = "pdf")

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

ICC %>%
  group_by(freqband) %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = icc)) +
  geom_sf(show.legend = TRUE) +
  facet_wrap( ~ freqband) + scale_fill_continuous_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = "ICC") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/ICC_per_band_continuous.pdf', device = "pdf")

ICC %>%
  group_by(freqband) %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = icc)) +
  geom_sf(show.legend = TRUE) +
  facet_wrap( ~ freqband) + scale_fill_binned_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = "ICC") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/ICC_per_band_binned.pdf', device = "pdf")

ICC <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/ICC_and_Heritability/ICC_avg.csv')
ICC_final <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Destrieux_atlas.csv')

ICC_final_avg_broad <- ICC_final
ICC_final_avg_broad$icc = ICC$total_avg

ICC_final_avg_broad %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = icc)) +
  geom_sf(show.legend = TRUE) +
  scale_fill_continuous_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = "ICC") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/ICC_mean_band_continuous.pdf', device = "pdf")


###### Heritability #######

Heritability <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/ICC_and_Heritability/Heritability_mean_avg_per_band.csv')
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
  facet_wrap( ~ freqband) + scale_fill_continuous_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = " h²") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/Heritability_per_band_continuous.pdf', device = "pdf")

Heritability %>%
  group_by(freqband) %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = h)) +
  geom_sf(show.legend = TRUE) +
  facet_wrap( ~ freqband) + scale_fill_binned_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = " h²") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/Heritability_per_band_binned.pdf', device = "pdf")



Heritability <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Results/ICC_and_Heritability/Heritability_avg.csv')
Heritability_final <- read.csv(file = 'Neuro_McGill/Fingerprints_Twins/Destrieux_atlas.csv')
Heritability_final_Broadband <- Heritability_final
Heritability_final_Broadband$h <- Heritability$total_avg

Heritability_final_Broadband %>%
  brain_join(desterieux) %>%
  reposition_brain(hemi ~ side) %>%
  ggplot(aes(fill = h)) +
  geom_sf(show.legend = TRUE) +
  scale_fill_continuous_sequential(palette= 'Rocket', rev= FALSE, begin = 0, end = 1, name = " h²") +
  theme_void() + scale_color_manual('white')

ggsave('Neuro_McGill/Fingerprints_Twins/Results/Figures/Heritability_mean_band_continuous.pdf', device = "pdf")


