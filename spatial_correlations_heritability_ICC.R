
library(ggplot2)
library(ggsegDesterieux)
library(ggseg)
library(dplyr)
library(colorspace)
library(ggsegSchaefer)
library(ggpubr)



# # # # # # # # # # # # # # # # # # # # # # #
## plot the schaefer atlas
## plot the Schaefer

########################## Schaefer atlas PSD #######################################

library(ggsegSchaefer)
atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')

psd_1= read.csv('~/Documents/HCP_twin_projec/outputs/Scaefer_PSD_m1_participantmean.csv', header = FALSE)
atlas$Ldelta= rowMeans(log10(psd_1[,1:8]))
atlas$Ltheta= rowMeans(log10(psd_1[,9:16]))
atlas$Lalpha= rowMeans(log10(psd_1[,17:26]))
atlas$Lbeta= rowMeans(log10(psd_1[,27:60]))
atlas$Lgamma= rowMeans(log10(psd_1[,61:100]))
atlas$LHgamma= rowMeans(log10(psd_1[,101:301]))


atlas %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Ldelta)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') 


atlas[,35:40]= lapply(atlas[,35:40], scales::rescale)
data4plot=cbind(stack(atlas[,35:40]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Brain_map_power.pdf', device = "pdf")



# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')
atlas[,35:40]= lapply(atlas[,35:40], scales::rescale)

data4plot=cbind(stack(atlas[,35:40]), atlas[,5:7])
colnames(data4plot)[4]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_power.pdf', device = "pdf")


atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')

psd_2= read.csv('~/Documents/HCP_twin_projec/outputs/Scaefer_PSD_m2_participantmean.csv', header = FALSE)
atlas$Ldelta= rowMeans(log10(psd_2[,1:8]))
atlas$Ltheta= rowMeans(log10(psd_2[,9:16]))
atlas$Lalpha= rowMeans(log10(psd_2[,17:26]))
atlas$Lbeta= rowMeans(log10(psd_2[,27:60]))
atlas$Lgamma= rowMeans(log10(psd_2[,61:100]))
atlas$LHgamma= rowMeans(log10(psd_2[,101:301]))


atlas %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Lbeta)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') 


atlas[,35:40]= lapply(atlas[,35:40], scales::rescale)
data4plot=cbind(stack(atlas[,35:40]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Brain_map_power_2.pdf', device = "pdf")



# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')
atlas[,35:40]= lapply(atlas[,35:40], scales::rescale)

data4plot=cbind(stack(atlas[,35:40]), atlas[,5:7])
colnames(data4plot)[4]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_power_2.pdf', device = "pdf")




atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')

psd_3= read.csv('~/Documents/HCP_twin_projec/outputs/Scaefer_PSD_m3_participantmean.csv', header = FALSE)
atlas$Ldelta= rowMeans(log10(psd_3[,1:8]))
atlas$Ltheta= rowMeans(log10(psd_3[,9:16]))
atlas$Lalpha= rowMeans(log10(psd_3[,17:26]))
atlas$Lbeta= rowMeans(log10(psd_3[,27:60]))
atlas$Lgamma= rowMeans(log10(psd_3[,61:100]))
atlas$LHgamma= rowMeans(log10(psd_3[,101:301]))


atlas %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Lbeta)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') 


atlas[,35:40]= lapply(atlas[,35:40], scales::rescale)
data4plot=cbind(stack(atlas[,35:40]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Brain_map_power_3.pdf', device = "pdf")


# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')
atlas[,35:40]= lapply(atlas[,35:40], scales::rescale)

data4plot=cbind(stack(atlas[,35:40]), atlas[,5:7])
colnames(data4plot)[4]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_power_3.pdf', device = "pdf")



########################## Schaefer atlas ICC ##########################

atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')
atlas$new_region=paste(atlas$region, atlas$Yeo, sep= '_')

ICC= read.csv('~/Documents/HCP_twin_projec/ICC_Scahefer.csv', header = TRUE)
ICC= ICC[-1]
#atlas=atlas[order(atlas$new_region),]

# Define X and Y data
X=data.frame(delta= rowMeans(ICC[,1:8]),theta= rowMeans(ICC[,9:16]),alpha= rowMeans(ICC[,17:26]),
             beta= rowMeans(ICC[,27:60]), gamma= rowMeans(ICC[,61:100]), hgamma= rowMeans(ICC[,101:301]))

atlas$delta= rowMeans((ICC[,1:8]))
atlas$theta= rowMeans((ICC[,9:16]))
atlas$alpha= rowMeans((ICC[,17:26]))
atlas$beta= rowMeans((ICC[,27:60]))
atlas$gamma= rowMeans((ICC[,61:100]))
atlas$high_gamma= rowMeans((ICC[,101:301]))
atlas$Boradband= rowMeans((X))

data4plot=cbind(stack(atlas[,36:41]), atlas[,c(1:3, 6)])

data4plot$values[data4plot$values<0.4] =0.4
data4plot %>%
  group_by(ind) %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) + viridis::scale_fill_viridis(option="magma", limits=c(0.4, 0.9)) +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_ICC_narrowband.pdf', device = "pdf")


atlas$Boradband[atlas$Boradband<0.6] =0.6
atlas %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", limits=c(0.6, 0.9)) +
  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,0.8)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_ICC_broadband.pdf', device = "pdf")

# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')


ggplot(atlas, aes(Yeo, Boradband, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('fingerprinting features')

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_ICC.pdf', device = "pdf")

########################## Schaefer atlas bootstrapped ICC ##########################

atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')
atlas$new_region=paste(atlas$region, atlas$Yeo, sep= '_')

ICCb= read.csv('~/Documents/HCP_twin_projec/ICC_bootstrap_avg_total.csv', header = TRUE)
ICCb= ICCb[-1]
#atlas=atlas[order(atlas$new_region),]

# Define X and Y data
Xb=data.frame(delta= rowMeans(ICCb[,1:8]),theta= rowMeans(ICCb[,9:16]),alpha= rowMeans(ICCb[,17:26]),
             beta= rowMeans(ICCb[,27:60]), gamma= rowMeans(ICCb[,61:100]), hgamma= rowMeans(ICCb[,101:301]))


ICCbb= read.csv('~/Documents/HCP_twin_projec/ICC_without_twins.csv', header = TRUE)
ICCbb= ICCbb[-1]
#atlas=atlas[order(atlas$new_region),]

# Define X and Y data
Xbb=data.frame(delta= rowMeans(ICCbb[,1:8]),theta= rowMeans(ICCbb[,9:16]),alpha= rowMeans(ICCbb[,17:26]),
              beta= rowMeans(ICCbb[,27:60]), gamma= rowMeans(ICCbb[,61:100]), hgamma= rowMeans(ICCbb[,101:301]))



atlas$delta= rowMeans((ICCb[,1:8]))
atlas$theta= rowMeans((ICCb[,9:16]))
atlas$alpha= rowMeans((ICCb[,17:26]))
atlas$beta= rowMeans((ICCb[,27:60]))
atlas$gamma= rowMeans((ICCb[,61:100]))
atlas$high_gamma= rowMeans((ICCb[,101:301]))
atlas$Boradband= rowMeans((Xb))

data4plot=cbind(stack(atlas[,36:41]), atlas[,c(1:3, 6)])

data4plot$values[data4plot$values<0.4] =0.4
data4plot %>%
  group_by(ind) %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) + viridis::scale_fill_viridis(option="magma", limits=c(0.4, 0.9)) +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_ICC_narrowband_boostrapped.pdf', device = "pdf")


atlas$Boradband[atlas$Boradband<0.6] =0.6
atlas %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", limits=c(0.6, 0.9)) +
  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,0.8)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_ICC_broadband_boostrapped.pdf', device = "pdf")

# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')


ggplot(atlas, aes(Yeo, Boradband, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('fingerprinting features')

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_ICC_boostrapped.pdf', device = "pdf")


cor.test(rowMeans(X), rowMeans(Xb)) # .995


########################## Schaefer atlas heritability ##########################

atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')
atlas$new_region=paste(atlas$region, atlas$Yeo, sep= '_')

H= read.csv('~/Documents/HCP_twin_projec/heritability_final.csv', header = TRUE)
H= H[-1]
#atlas=atlas[order(atlas$new_region),]


# Define X and Y data
h=data.frame(delta= rowMeans(H[,1:8]),theta= rowMeans(H[,9:16]),alpha= rowMeans(H[,17:26]),
             beta= rowMeans(H[,27:60]), gamma= rowMeans(H[,61:100]), hgamma= rowMeans(H[,101:301]))

atlas$delta= rowMeans((H[,1:8]))
atlas$theta= rowMeans((H[,9:16]))
atlas$alpha= rowMeans((H[,17:26]))
atlas$beta= rowMeans((H[,27:60]))
atlas$gamma= rowMeans((H[,61:100]))
atlas$high_gamma= rowMeans((H[,101:301]))
atlas$Boradband= rowMeans((h))

data4plot=cbind(stack(atlas[,36:41]), atlas[,c(1:3, 6)])

data4plot$values[data4plot$values< -0.4] = -0.4
data4plot$values[data4plot$values> 1.5] = 1.5
data4plot %>%
  group_by(ind) %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  viridis::scale_fill_viridis(option="magma", limits=c(-0.4, 1.5)) + #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_H_narrowband_June2023.pdf', device = "pdf")

data4plot$values[data4plot$values< -0.4] = -0.4
data4plot$values[data4plot$values> 1.2] = 1.2
atlas %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", limits=c(-0.4, 1.2)) +
  #scale_fill_continuous_sequential(palette= 'Magma', rev= FALSE, limits= c(-0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_H_broadband_June2023.pdf', device = "pdf")


# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')

ggplot(atlas, aes(Yeo, Boradband, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('heritability')

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_H_June2023.pdf', device = "pdf")




########################## correlate heritability & ICC ##########################

cor.test(rowMeans(Xb), rowMeans(h))

lm3=lm(scale(rowMeans(X)) ~ scale(rowMeans(h)))
summary(lm3)
sjPlot::tab_model(lm3)

permuted_index= read.csv('/Users/jason/Desktop/Schaefer_neuromaps/permuted_indexes_of_schaefer_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1

#gradient offset
orig=cor.test(rowMeans(Xb), rowMeans(h))
permuted_corr=c()
for (i in 1001:2000){
  
  cor_temp=cor.test(rowMeans(Xb), rowMeans(h)[permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}

sum(orig$estimate < permuted_corr)/1000
# permuted p value = 0.026

someData_scatter= data.frame(ICC= rowMeans(Xb), h2=rowMeans(h))

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic2() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_broadband_ICC_vs_h_June2023.pdf', device = "pdf", width = 7, height = 7, dpi = 800)

cor.test(Xb$delta, h$delta)
cor.test(Xb$theta, h$theta)
cor.test(Xb$alpha, h$alpha)
cor.test(Xb$beta, h$beta)
cor.test(Xb$gamma, h$gamma)
cor.test(Xb$hgamma, h$hgamma)


lm3=lm(scale(Xb$alpha) ~ scale(h$alpha))
summary(lm3)

someData_scatter= data.frame(ICC= X$alpha, h2=h$alpha)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic2() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_alpha_ICC_vs_h_June2023.pdf', device = "pdf", width = 7, height = 7, dpi = 800)



lm3=lm(scale(Xb$beta) ~ scale(h$beta))
summary(lm3)

someData_scatter= data.frame(ICC= X$beta, h2=h$beta)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic2() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_beta_ICC_vs_h_June2023.pdf', device = "pdf", width = 7, height = 7, dpi = 800)


delta=cor.test(Xb$delta, h$delta)
theta=cor.test(Xb$theta, h$theta)
alpha=cor.test(Xb$alpha, h$alpha)
beta=cor.test(Xb$beta, h$beta)
gamma=cor.test(Xb$gamma, h$gamma)
hgamma=cor.test(Xb$hgamma, h$hgamma)


data4plot= data.frame(corrs= c(delta$estimate,theta$estimate,alpha$estimate,
                    beta$estimate, gamma$estimate, hgamma$estimate), 
           CIlower= c(delta$conf.int[1], theta$conf.int[1], alpha$conf.int[1], 
           beta$conf.int[1], gamma$conf.int[1], hgamma$conf.int[1]),
CIlupper= c(delta$conf.int[2], theta$conf.int[2], alpha$conf.int[2], 
beta$conf.int[2], gamma$conf.int[2], hgamma$conf.int[2]), band= c('delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))

data4plot$band = factor( data4plot$band, levels=c('delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))

ggplot(data4plot, aes(x=band , y = corrs, colour= band, fill=band)) + geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_errorbar(colour="#000000", aes(ymin=CIlower, ymax=CIlupper),  width=.2, position=position_dodge(.9)) + geom_point(stat='identity',position=position_dodge(.9), size= 10)   + 
  scale_fill_manual(values=cbbPalette)  + scale_color_manual(values=cbbPalette) + ggpubr::theme_classic2() + xlab("") + ylab("correlation between icc & heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/ICC_heritability_corr_June2023.pdf', device = "pdf", height = 5, width = 7)
ggsave('~/Documents/HCP_twin_projec/figures/ICC_heritability_corrJune2023.jpeg', device = "jpeg", height = 5, width = 7)





#gradient offset
orig=cor.test(Xb$delta, h$delta)
permuted_corr_delta=c()
for (i in 1001:2000){
  
  cor_temp=cor.test(Xb$delta, h$delta[permuted_index[,i]])
  permuted_corr_delta= c(permuted_corr_delta, cor_temp$estimate)
  
}

sum(orig$estimate > permuted_corr_delta)/1000 #0.101


orig=cor.test(Xb$theta, h$theta)
permuted_corr_theta=c()
for (i in 2001:3000){
  
  cor_temp=cor.test(Xb$theta, h$theta[permuted_index[,i]])
  permuted_corr_theta= c(permuted_corr_theta, cor_temp$estimate)
  
}

sum(orig$estimate < permuted_corr_theta)/1000 #0.477


orig=cor.test(Xb$alpha, h$alpha)
permuted_corr_alpha=c()
for (i in 3001:4000){
  
  cor_temp=cor.test(Xb$alpha, h$alpha[permuted_index[,i]])
  permuted_corr_alpha= c(permuted_corr_alpha, cor_temp$estimate)
  
}

sum(orig$estimate < permuted_corr_alpha)/1000 #0


orig=cor.test(Xb$beta, h$beta)
permuted_corr_beta=c()
for (i in 4001:5000){
  
  cor_temp=cor.test(Xb$beta, h$beta[permuted_index[,i]])
  permuted_corr_beta= c(permuted_corr_beta, cor_temp$estimate)
  
}

sum(orig$estimate < permuted_corr_beta)/1000 #0

orig=cor.test(Xb$gamma, h$gamma)
permuted_corr_gamma=c()
for (i in 5001:6000){
  
  cor_temp=cor.test(Xb$gamma, h$gamma[permuted_index[,i]])
  permuted_corr_gamma= c(permuted_corr_gamma, cor_temp$estimate)
  
}

sum(orig$estimate < permuted_corr_gamma)/1000 # 0.069


orig=cor.test(Xb$hgamma, h$hgamma)
permuted_corr_hgamma=c()
for (i in 6001:7000){
  
  cor_temp=cor.test(Xb$hgamma, h$hgamma[permuted_index[,i]])
  permuted_corr_hgamma= c(permuted_corr_hgamma, cor_temp$estimate)
  
}

sum(orig$estimate < permuted_corr_hgamma)/1000 # 0.371

p.adjust(c(0.101,0.477,0,0,0.069,0.371), 'fdr')


dataperm= data.frame(null= c(permuted_corr_delta, permuted_corr_theta, permuted_corr_alpha,
                               permuted_corr_beta, permuted_corr_gamma, permuted_corr_hgamma), band= rep(c('delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'), each=1000))

dataperm$band = factor( data4plot$band, levels=c('delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))



ggplot(data4plot, aes(x=band , y = corrs, colour= band, fill='null')) + geom_hline(yintercept=0, linetype="dashed", color = "grey") +
 geom_boxplot(data= dataperm, aes(x=band , y = null, alpha=0.5 )) + geom_point(stat='identity',position=position_dodge(.9), size= 10)   + 
  scale_fill_manual(values='#E8E8E8')  + scale_color_manual(values=cbbPalette) + ggpubr::theme_classic2() + xlab("") + ylab("correlation between icc & heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/ICC_heritability_corr_June2023_with_nulls.pdf', device = "pdf", height = 5, width = 7)
ggsave('~/Documents/HCP_twin_projec/figures/ICC_heritability_corrJune2023_with_nulls.jpeg', device = "jpeg", height = 5, width = 7)




########################## correlate heritability & ICC bootstrapped ##########################

cor.test(rowMeans(Xb), rowMeans(h))

lm3=lm(scale(rowMeans(Xb)) ~ scale(rowMeans(h)))
summary(lm3)
sjPlot::tab_model(lm3)

someData_scatter= data.frame(ICC= rowMeans(Xb), h2=rowMeans(h))

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic2() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_broadband_ICC_bootstrapped_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)

cor.test(Xb$delta, h$delta)
cor.test(Xb$theta, h$theta)
cor.test(Xb$alpha, h$alpha)
cor.test(Xb$beta, h$beta)
cor.test(Xb$gamma, h$gamma)
cor.test(Xb$hgamma, h$hgamma)


lm3=lm(scale(Xb$alpha) ~ scale(h$alpha))
summary(lm3)

someData_scatter= data.frame(ICC= Xb$alpha, h2=h$alpha)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic2() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_alpha_ICC_bootstrapped_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)



lm3=lm(scale(Xb$beta) ~ scale(h$beta))
summary(lm3)

someData_scatter= data.frame(ICC= Xb$beta, h2=h$beta)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic2() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_beta_ICC_bootstrapped_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)


delta=cor.test(Xb$delta, h$delta)
theta=cor.test(Xb$theta, h$theta)
alpha=cor.test(Xb$alpha, h$alpha)
beta=cor.test(Xb$beta, h$beta)
gamma=cor.test(Xb$gamma, h$gamma)
hgamma=cor.test(Xb$hgamma, h$hgamma)


data4plot= data.frame(corrs= c(delta$estimate,theta$estimate,alpha$estimate,
                               beta$estimate, gamma$estimate, hgamma$estimate), 
                      CIlower= c(delta$conf.int[1], theta$conf.int[1], alpha$conf.int[1], 
                                 beta$conf.int[1], gamma$conf.int[1], hgamma$conf.int[1]),
                      CIlupper= c(delta$conf.int[2], theta$conf.int[2], alpha$conf.int[2], 
                                  beta$conf.int[2], gamma$conf.int[2], hgamma$conf.int[2]), band= c('delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))

data4plot$band = factor( data4plot$band, levels=c('delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))

ggplot(data4plot, aes(x=band , y = corrs, colour= band, fill=band)) + geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_errorbar(colour="#000000", aes(ymin=CIlower, ymax=CIlupper),  width=.2, position=position_dodge(.9)) + geom_point(stat='identity',position=position_dodge(.9), size= 10)   + 
  scale_fill_manual(values=cbbPalette)  + scale_color_manual(values=cbbPalette) + ggpubr::theme_classic2() + xlab("") + ylab("correlation between icc & heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/ICC_bootstrapped_heritability_corr.pdf', device = "pdf", height = 5, width = 7)
ggsave('~/Documents/HCP_twin_projec/figures/ICC_bootstrapped_heritability_corr.jpeg', device = "jpeg", height = 5, width = 7)


