##############################
###       PLOTS         ######
##############################

library(viridis)
library(ggseg)
library(ggplot2)
library(ggridges)
library(patchwork)
library(dplyr)
library(colorspace)
library(cowplot)
library(ggsegDK)
library(ggpubr)

magma_3 = c("#51127c", "#b73779", "#fc8961")

# --------------------------------------
#         PLOT ACCURACIES
# --------------------------------------

# Import the data
data <- read.csv("Results_Log_DK/PSD_accuracy/All_accuracies_every_freq_stacked.csv")

# Define the order of the frequency bands
freq_order <- c("BROADBAND", "DELTA", "THETA", "ALPHA", "BETA", "GAMMA", "HIGH GAMMA")

# Define the order of the account types
acc_order <- c("Autocorr", "Crosscorr MZ", "Crosscorr DZ")

# Create the aesthetic boxplot
ggplot(data, aes(x=factor(FreqBand, levels = freq_order), y=values, fill=factor(AccType, levels = acc_order))) + 
  geom_boxplot()+
  scale_fill_manual(values = magma_3, labels = c("Self-matching accuracy", "MZ matching accuracy", "DZ matching accuracy")) +
  labs(x = "Frequency Band", y = "Accuracy", color = "Accuracy Type", fill = "Accuracy Type") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 10), panel.background = element_rect(fill='white', colour='black'))

# Define the output file path and name
out_file <- "Results_Log_DK/Figures/accuracy_box_plot.pdf"

# Save the plot as PDF
ggsave(filename = out_file, width = 10, height = 6, units = "in")


# ------ Bar plot -------

# Calculate means and standard errors
means <- aggregate(values ~ FreqBand + AccType, data, mean)
se <- aggregate(values ~ FreqBand + AccType, data, function(x) sd(x)/sqrt(length(x)))
sd <- aggregate(values ~ FreqBand + AccType, data, function(x) sd(x))

means$se = se$values
means$sd = sd$values


# Create the bar plot with error bars
ggplot(data = means, aes(x = factor(FreqBand, levels = freq_order), y = values, fill = factor(AccType, levels = acc_order))) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = values - sd, ymax = values + sd), position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_manual(values = magma_3, labels = c("Self-matching accuracy", "MZ matching accuracy", "DZ matching accuracy")) +
  labs(x = "Frequency Band", y = "Accuracy", color = "Accuracy Type", fill = "Accuracy Type") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 10), panel.background = element_rect(fill='white', colour='black'))

# Define the output file path and name
out_file <- "Results_Log_DK/Figures/accuracy_bar_plot.pdf"

# Save the plot as PDF
ggsave(filename = out_file, width = 10, height = 6, units = "in")


# --------------------------------------
#         PLOT CORRELATIONS
# --------------------------------------

all_corr = read.csv("Results_Log_DK/PSD_correlations/All_correlation_every_freq_stacked.csv")

ggplot(all_corr, aes(x=factor(FreqBand, levels = freq_order), y=values, fill=paste(TypeCorr, TwinType))) + 
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) + 
  scale_fill_discrete(name="Correlation Type") +
  labs(x = "Frequency Band", y = "Pearson Correlation")

ggsave('Results_Log_DK/Figures/Corr_per_Band_Twin_Auto_and_Cross_DK.pdf', device = "pdf")



bands <- unique(all_corr$FreqBand)
theme_set(theme_bw())

for (band in bands) {
  
  auto_corr <- all_corr %>% 
    filter(FreqBand == band & TypeCorr == "Autocorr")
  
  p1 <- ggplot(auto_corr, aes(x = values, fill = TwinType, group = TwinType, col = TwinType)) +
    geom_density(alpha = 0.55) +
    scale_fill_manual(values = magma_3) +
    scale_color_manual(values = magma_3) + 
    xlab("Correlation") +
    ylab("Density") +
    ggtitle("Correlation with oneself") +
    theme(plot.title = element_text(face = "bold", size = 14), panel.background = element_rect(fill='white', colour='black'))
  
  cross_corr <- all_corr %>% 
    filter(FreqBand == band & TypeCorr == "Crosscorr")
  
  p2 <- ggplot(cross_corr, aes(x = values, fill = TwinType, group = TwinType, col = TwinType)) +
    geom_density(alpha = 0.55) +
    scale_fill_manual(values = magma_3) +
    scale_color_manual(values = magma_3) + 
    xlab("Correlation") +
    ylab("Density") +
    ggtitle("Correlation with one's twin") +
    theme(plot.title = element_text(face = "bold", size = 14), panel.background = element_rect(fill='white', colour='black'))
  
  p1 + p2 + plot_layout(guides = "collect")
  
  ggsave(filename = paste('Results_Log_DK/Figures/Corr_density', band ,'DK.pdf'), device = "pdf", width = 14, height = 6, units = "in")
  
}


# --------------------------------------
#             PLOT ICC
# --------------------------------------

# Read in atlas and ICC data
atlas <- read.csv('new_Data/dka_HCP_test.csv')
atlas <- atlas[, c(1, 3, 4, 5)]

ICC <- read.csv('Results_Log_DK/ICC_and_Heritability/ICC.csv', header = TRUE)
ICC <- ICC[-1]

# Define X data as a data frame with the row means for each frequency band
X <- data.frame(delta = rowMeans(ICC[,1:8]), theta = rowMeans(ICC[,9:16]), alpha = rowMeans(ICC[,17:26]),
                beta = rowMeans(ICC[,27:60]), gamma = rowMeans(ICC[,61:100]), hgamma = rowMeans(ICC[,101:301]))

# Add row means to the atlas data frame for each frequency band and broadband
atlas$Delta <- rowMeans((ICC[,1:8]))
atlas$Theta <- rowMeans((ICC[,9:16]))
atlas$Alpha <- rowMeans((ICC[,17:26]))
atlas$Beta <- rowMeans((ICC[,27:60]))
atlas$Gamma <- rowMeans((ICC[,61:100]))
atlas$High_Gamma <- rowMeans((ICC[,101:301]))
atlas$Broadband <- rowMeans((X))

# Replace any values in the data that are less than 0.4 with 0.4
data4plot <- cbind(stack(atlas[,5:10]), atlas[,1:4])
data4plot$values[data4plot$values < 0.4] <- 0.4

# Plot ICC narrowbands for each network using the DK7_200 atlas
data4plot %>%
  group_by(ind) %>% 
  brain_join(dk) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Inferno', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Inferno', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

# Save the narrowbands ICC plot as a pdf
ggsave('Results_Log_DK/Figures/ICC_narrowbands_DK.pdf', device = "pdf")

# Replace any values in the broadband column that are less than 0.6 with 0.6
atlas$Broadband[atlas$Broadband < 0.6] <- 0.6

# Plot ICC broadband using the DK7_200 atlas
atlas %>%
  brain_join(dk) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill = Broadband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0.6,0.8)) +  
  scale_color_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0.6,1)) +
  theme_void() + 
  scale_color_manual('white')

# Save the broadband ICC plot as a pdf
ggsave('Results_Log_DK/Figures/ICC_broadbands_DK.pdf', device = "pdf")



# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')


ggplot(atlas, aes(network, Broadband, fill=network)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('fingerprinting features')

ggsave('Results_Log_DK/Figures/Boxplot_ICC_vs_Networks_DK.pdf', device = "pdf")



# --------------------------------------
#         PLOT HERITABILITY
# --------------------------------------

# Read in atlas and Heritability data
atlas <- read.csv('new_Data/dka_HCP_test.csv')
atlas <- atlas[, c(1, 3, 4, 5)]

H <- read.csv('Results_Log_DK/ICC_and_Heritability/heritability_mean.csv', header = TRUE)
H <- H[-1]

# Define h data as a data frame with the row means for each frequency band
h <- data.frame(delta = rowMeans(H[,1:8]), theta = rowMeans(H[,9:16]), alpha = rowMeans(H[,17:26]),
                beta = rowMeans(H[,27:60]), gamma = rowMeans(H[,61:100]), hgamma = rowMeans(H[,101:301]))

# Add row means to the atlas data frame for each frequency band and broadband
atlas$Delta <- rowMeans((H[,1:8]))
atlas$Theta <- rowMeans((H[,9:16]))
atlas$Alpha <- rowMeans((H[,17:26]))
atlas$Beta <- rowMeans((H[,27:60]))
atlas$Gamma <- rowMeans((H[,61:100]))
atlas$High_Gamma <- rowMeans((H[,101:301]))
atlas$Broadband <- rowMeans((h))

# Replace any values in the data that are less than 0.4 with 0.4
data4plot <- cbind(stack(atlas[,5:10]), atlas[,1:4])
data4plot$values[data4plot$values < 0.0] <- 0.0
data4plot$values[data4plot$values > 1.0] <- 1.0

# Plot H narrowbands for each network using the DK7_200 atlas
data4plot %>%
  group_by(ind) %>% 
  brain_join(dk) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill = values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  
  scale_fill_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0,1)) +  
  scale_color_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0,1)) +
  theme_void() + 
  scale_color_manual('white')

# Save the narrowbands H plot as a pdf
ggsave('Results_Log_DK/Figures/Heritability_narrowbands_DK.pdf', device = "pdf")

# Replace any values in the broadband column that are less than 0.6 with 0.6
atlas$Broadband[atlas$Broadband < 0.0] <- 0.0
atlas$Broadband[atlas$Broadband > 1.0] <- 1.0

# Plot H broadband using the DK7_200 atlas
atlas %>%
  brain_join(dk) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill = Broadband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0.0,1.0)) +  
  scale_color_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0.0,1)) +
  theme_void() + 
  scale_color_manual('white')

# Save the broadband H plot as a pdf
ggsave('Results_Log_DK/Figures/Heritability_broadbands_DK.pdf', device = "pdf")



# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')


ggplot(atlas, aes(network, Broadband, fill=network)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('Heritability features')

ggsave('Results_Log_DK/Figures/Boxplot_Heritability_vs_Networks_DK.pdf', device = "pdf")

# --------------------------------------
#         CORR ICC vs HERITABILITY
# --------------------------------------

X$broadband = rowMeans(X[1:6])
h$broadband = rowMeans(h[1:6])

cor.test(X$broadband, h$broadband)
cor.test(X$delta, h$delta)
cor.test(X$theta, h$theta)
cor.test(X$alpha, h$alpha)
cor.test(X$beta, h$beta)
cor.test(X$gamma, h$gamma)
cor.test(X$hgamma, h$hgamma)


lm3=lm(scale(X$alpha) ~ scale(h$alpha))
summary(lm3)

someData_scatter= data.frame(ICC= X$alpha, h2=h$alpha)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_Alpha.pdf", device = "pdf", width = 8, height = 6, units = "in")

### P-value test ###
shuffle_idx = read.csv("new_Data/permuted_indexes_of_DK_atlas_SPINs&Twirl.csv", header = TRUE)
shuffle_idx = shuffle_idx[-1]
orig=cor.test(X$alpha, h$alpha)
permuted_corr=c()
for (i in 1000:2001){
  cor_temp=cor.test(X$alpha, h$alpha[shuffle_idx[,i] + 1])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
}
sum(orig$estimate < permuted_corr)/1000


lm3=lm(scale(X$beta) ~ scale(h$beta))
summary(lm3)

someData_scatter= data.frame(ICC= X$beta, h2=h$beta)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_Beta.pdf", device = "pdf", width = 8, height = 6, units = "in")



lm3=lm(scale(X$delta) ~ scale(h$delta))
summary(lm3)

someData_scatter= data.frame(ICC= X$delta, h2=h$delta)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_Delta.pdf", device = "pdf", width = 8, height = 6, units = "in")



lm3=lm(scale(X$theta) ~ scale(h$theta))
summary(lm3)

someData_scatter= data.frame(ICC= X$theta, h2=h$theta)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_Theta.pdf", device = "pdf", width = 8, height = 6, units = "in")



lm3=lm(scale(X$gamma) ~ scale(h$gamma))
summary(lm3)

someData_scatter= data.frame(ICC= X$gamma, h2=h$gamma)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_Gamma.pdf", device = "pdf", width = 8, height = 6, units = "in")



lm3=lm(scale(X$hgamma) ~ scale(h$hgamma))
summary(lm3)

someData_scatter= data.frame(ICC= X$hgamma, h2=h$hgamma)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_High_Gamma.pdf", device = "pdf", width = 8, height = 6, units = "in")


lm3=lm(scale(X$broadband) ~ scale(h$broadband))
summary(lm3)

someData_scatter= data.frame(ICC= X$broadband, h2=h$broadband)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Salient Fingerprinting Features (ICC)") + ylab("Heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none")

ggsave("Results_Log_DK/Figures/Corr_ICC_vs_Heritability_Broadband.pdf", device = "pdf", width = 8, height = 6, units = "in")


### P-value test ###
shuffle_idx = read.csv("new_Data/permuted_indexes_of_DK_atlas_SPINs&Twirl.csv", header = TRUE)
shuffle_idx = shuffle_idx[-1]
orig=cor.test(X$broadband, h$broadband)
permuted_corr=c()
for (i in 1000:2001){
  cor_temp=cor.test(X$broadband, h$broadband[shuffle_idx[,i] + 1])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
}
sum(orig$estimate < permuted_corr)/1000

# --------------------------------------
#  CORR NEURORECEPTORS vs HERITABILITY
# --------------------------------------

neuro_data = read.csv("new_Data/DK2018_200Parcels_7Networks_Neuromaps.csv")

neuro_receptors = neuro_data[,16:34]

H <- read.csv('Results_Log_DK/ICC_and_Heritability/heritability_mean.csv', header = TRUE)
H <- H[-1]

neuro_data$new_region <- paste(neuro_data$region, neuro_data$Yeo, sep= '_')


# Define h data as a data frame with the row means for each frequency band
h <- data.frame(delta = rowMeans(H[,1:8]), theta = rowMeans(H[,9:16]), alpha = rowMeans(H[,17:26]),
                beta = rowMeans(H[,27:60]), gamma = rowMeans(H[,61:100]), hgamma = rowMeans(H[,101:301]))

h$broadband <- rowMeans((h[1:6]))

### Heatmap correlations done on Python ###


lm3=lm(scale(h$broadband) ~ scale(neuro_receptors$D1))
summary(lm3)

someData_scatter= data.frame(h2= h$broadband, receptor = neuro_receptors$D1)

ggplot(someData_scatter, aes(x=h2 , y = receptor, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Heritability") + ylab("D1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave("Results_Log_DK/Figures/Corr_Heritability_D1_Broadband_DK.pdf", device = "pdf", width = 8, height = 6, units = "in")


lm3=lm(scale(h$broadband) ~ scale(neuro_receptors$D2))
summary(lm3)

someData_scatter= data.frame(h2= h$broadband, receptor = neuro_receptors$D2)

ggplot(someData_scatter, aes(x=h2 , y = receptor, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Heritability") + ylab("D2")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave("Results_Log_DK/Figures/Corr_Heritability_D2_Broadband_DK.pdf", device = "pdf", width = 8, height = 6, units = "in")


lm3=lm(scale(h$broadband) ~ scale(neuro_receptors$MOR))
summary(lm3)

someData_scatter= data.frame(h2= h$broadband, receptor = neuro_receptors$MOR)

ggplot(someData_scatter, aes(x=h2 , y = receptor, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Heritability") + ylab("MOR")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave("Results_Log_DK/Figures/Corr_Heritability_MOR_Broadband_DK.pdf", device = "pdf", width = 8, height = 6, units = "in")


lm3=lm(scale(h$broadband) ~ scale(neuro_receptors$NET))
summary(lm3)

someData_scatter= data.frame(h2= h$broadband, receptor = neuro_receptors$NET)

ggplot(someData_scatter, aes(x=h2 , y = receptor, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("Heritability") + ylab("NET")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave("Results_Log_DK/Figures/Corr_Heritability_NET_Broadband_DK.pdf", device = "pdf", width = 8, height = 6, units = "in")


# --------------------------------------
#     FUNC GRADIENT vs HERITABILITY
# --------------------------------------


neuro_data = read.csv("new_Data/DK2018_200Parcels_7Networks_Neuromaps.csv")

func_grad = neuro_data[,14:15]

H <- read.csv('Results_Log_DK/ICC_and_Heritability/heritability_mean.csv', header = TRUE)
H <- H[-1]

neuro_data$new_region <- paste(neuro_data$region, neuro_data$Yeo, sep= '_')


# Define h data as a data frame with the row means for each frequency band
h <- data.frame(delta = rowMeans(H[,1:8]), theta = rowMeans(H[,9:16]), alpha = rowMeans(H[,17:26]),
                beta = rowMeans(H[,27:60]), gamma = rowMeans(H[,61:100]), hgamma = rowMeans(H[,101:301]))

h$broadband <- rowMeans((h[1:6]))

func_grad$heritability[func_grad$heritability < 0.2] <- 0.2
func_grad$heritability[func_grad$heritability > 0.75] <- 0.75

ggplot(func_grad, aes(fcgradient01, fcgradient02)) +
  geom_point(aes(colour = heritability), size = 4.5, alpha = 1) +
  scale_fill_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0.2,0.75)) +  
  scale_color_continuous_sequential(palette = 'Inferno', rev = FALSE, limits = c(0.2,0.75)) +
  theme(text = element_text(size = 16), axis.text = element_text(size = 10), panel.background = element_rect(fill='white', colour='black'))

ggsave("Results_Log_DK/Figures/Func_Grad_and_Heritability_Broadband_DK.pdf", device = "pdf", width = 8, height = 6, units = "in")



