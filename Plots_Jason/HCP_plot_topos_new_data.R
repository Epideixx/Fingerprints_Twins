
library(ggplot2)
library(ggsegDesterieux)
library(ggseg)
library(dplyr)
library(colorspace)
library(ggsegSchaefer)
library(ggpubr)


# # # # # # # # # # # # # # # # # # # # # # #
## plot the Destrieux atlas
#################################################################

atlas= read.csv('~/Desktop/Destrieux_atlas_test_HCP.csv')

psd_1= read.csv('~/Documents/HCP_twin_projec/outputs/New_PSD_m1_participantmean.csv', header = FALSE)

atlas$Ldelta= rowMeans(log10(psd_1[,1:8]))
atlas$Ltheta= rowMeans(log10(psd_1[,9:16]))
atlas$Lalpha= rowMeans(log10(psd_1[,17:26]))
atlas$Lbeta= rowMeans(log10(psd_1[,27:60]))
atlas$Lgamma= rowMeans(log10(psd_1[,61:100]))
atlas$LHgamma= rowMeans(log10(psd_1[,101:301]))

data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Topo_rest_1.pdf', device = "pdf")

atlas[,22:27]= lapply(atlas[,22:27], scales::rescale)
data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind ) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Topo_rest_1_free_scale.pdf', device = "pdf")



atlas= read.csv('~/Desktop/Destrieux_atlas_test_HCP.csv')

psd_2= read.csv('~/Documents/HCP_twin_projec/outputs/New_PSD_m2_participantmean.csv', header = FALSE)

atlas$Ldelta= rowMeans(log10(psd_2[,1:8]))
atlas$Ltheta= rowMeans(log10(psd_2[,9:16]))
atlas$Lalpha= rowMeans(log10(psd_2[,17:26]))
atlas$Lbeta= rowMeans(log10(psd_2[,27:60]))
atlas$Lgamma= rowMeans(log10(psd_2[,61:100]))
atlas$LHgamma= rowMeans(log10(psd_2[,101:301]))



data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Topo_rest_2.pdf', device = "pdf")

atlas[,22:27]= lapply(atlas[,22:27], scales::rescale)
data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind ) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Topo_rest_2_free_scale.pdf', device = "pdf")

atlas= read.csv('~/Desktop/Destrieux_atlas_test_HCP.csv')

psd_3= read.csv('~/Documents/HCP_twin_projec/outputs/New_PSD_m3_participantmean.csv', header = FALSE)

atlas$Ldelta= rowMeans(log10(psd_3[,1:8]))
atlas$Ltheta= rowMeans(log10(psd_3[,9:16]))
atlas$Lalpha= rowMeans(log10(psd_3[,17:26]))
atlas$Lbeta= rowMeans(log10(psd_3[,27:60]))
atlas$Lgamma= rowMeans(log10(psd_3[,61:100]))
atlas$LHgamma= rowMeans(log10(psd_3[,101:301]))



data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Topo_rest_3.pdf', device = "pdf")

atlas[,22:27]= lapply(atlas[,22:27], scales::rescale)
data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind ) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Topo_rest_3_free_scale.pdf', device = "pdf")

# # # # # # # # # # # # # # # # # # # # # # #
## plot the schaefer atlas
## plot the Schaefer
#################################################################

library(ggsegSchaefer)
atlas= read.csv('~/Documents/HCP_twin_projec/Schaefer_atlas.csv')

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


atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)
data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
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
atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)

data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_power.pdf', device = "pdf")


atlas= read.csv('~/Documents/HCP_twin_projec/Schaefer_atlas.csv')

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


atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)
data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
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
atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)

data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_power_2.pdf', device = "pdf")




atlas= read.csv('~/Documents/HCP_twin_projec/Schaefer_atlas.csv')

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


atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)
data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Brain_map_power_3.pdf', device = "pdf")





## relative 
#################################################################

atlas= read.csv('~/Documents/HCP_twin_projec/Schaefer_atlas.csv')

psd_1= read.csv('~/Documents/HCP_twin_projec/outputs/Scaefer_Relative_PSD_m1_participantmean.csv', header = FALSE)
atlas$Ldelta= rowMeans((psd_1[,1:8]))
atlas$Ltheta= rowMeans((psd_1[,9:16]))
atlas$Lalpha= rowMeans((psd_1[,17:26]))
atlas$Lbeta= rowMeans((psd_1[,27:60]))
atlas$Lgamma= rowMeans((psd_1[,61:100]))
atlas$LHgamma= rowMeans((psd_1[,101:301]))


atlas %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Lalpha)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') 


atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)
data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Relative_Brain_map_power.pdf', device = "pdf")

# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')
atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)

data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Relative_box_plot_power.pdf', device = "pdf")




atlas= read.csv('~/Documents/HCP_twin_projec/Schaefer_atlas.csv')

psd_2= read.csv('~/Documents/HCP_twin_projec/outputs/Scaefer_Relative_PSD_m2_participantmean.csv', header = FALSE)
atlas$Ldelta= rowMeans((psd_2[,1:8]))
atlas$Ltheta= rowMeans((psd_2[,9:16]))
atlas$Lalpha= rowMeans((psd_2[,17:26]))
atlas$Lbeta= rowMeans((psd_2[,27:60]))
atlas$Lgamma= rowMeans((psd_2[,61:100]))
atlas$LHgamma= rowMeans((psd_2[,101:301]))


atlas %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Lalpha)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') 


atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)
data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Relative_Brain_map_power_2.pdf', device = "pdf")


# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')
atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)

data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Relative_box_plot_power_2.pdf', device = "pdf")


psd_3= read.csv('~/Documents/HCP_twin_projec/outputs/Scaefer_Relative_PSD_m3_participantmean.csv', header = FALSE)
atlas$Ldelta= rowMeans((psd_3[,1:8]))
atlas$Ltheta= rowMeans((psd_3[,9:16]))
atlas$Lalpha= rowMeans((psd_3[,17:26]))
atlas$Lbeta= rowMeans((psd_3[,27:60]))
atlas$Lgamma= rowMeans((psd_3[,61:100]))
atlas$LHgamma= rowMeans((psd_3[,101:301]))


atlas %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Lalpha)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') 


atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)
data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Relative_Brain_map_power_3.pdf', device = "pdf")

# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')
atlas[,8:13]= lapply(atlas[,8:13], scales::rescale)

data4plot=cbind(stack(atlas[,8:13]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_Relative_box_plot_power_3.pdf', device = "pdf")

##################################################################


library(ggplot2)
library(ggsegDesterieux)
library(ggseg)
library(dplyr)
library(colorspace)
library(ggsegSchaefer)


# # # # # # # # # # # # # # # # # # # # # # #
## plot the Destrieux atlas
#################################################################

atlas= read.csv('Destrieux_atlas.csv')

ICC= read.csv('Results_Log_Destrieux/ICC_and_Heritability/ICC.csv', header = TRUE)
ICC= ICC[-1]

# Define X and Y data
X=data.frame(delta= rowMeans(ICC[,1:8]),theta= rowMeans(ICC[,9:16]),alpha= rowMeans(ICC[,17:26]),
             beta= rowMeans(ICC[,27:60]), gamma= rowMeans(ICC[,61:100]), hgamma= rowMeans(ICC[,101:301]))

atlas$Ldelta= rowMeans((ICC[,1:8]))
atlas$Ltheta= rowMeans((ICC[,9:16]))
atlas$Lalpha= rowMeans((ICC[,17:26]))
atlas$Lbeta= rowMeans((ICC[,27:60]))
atlas$Lgamma= rowMeans((ICC[,61:100]))
atlas$LHgamma= rowMeans((ICC[,101:301]))
atlas$Boradband= rowMeans((X))

data4plot=cbind(stack(atlas[,10:15]), atlas[,1:3])

data4plot$values[data4plot$values<0.4] =0.4
data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/ICC_narrowband.pdf', device = "pdf")


atlas$Boradband[atlas$Boradband<0.6] =0.6
atlas %>%
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,0.8)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,1)) +
  theme_void() + scale_color_manual('white')


ggsave('~/Documents/HCP_twin_projec/figures/ICC_broadband.pdf', device = "pdf")

# heritability
atlas= read.csv('~/Desktop/Destrieux_atlas_test_HCP.csv')

H= read.csv('~/Documents/HCP_twin_projec/heritability.csv', header = TRUE)
H= H[-1]

# Define X and Y data
h=data.frame(delta= rowMeans(H[,1:8]),theta= rowMeans(H[,9:16]),alpha= rowMeans(H[,17:26]),
             beta= rowMeans(H[,27:60]), gamma= rowMeans(H[,61:100]), hgamma= rowMeans(H[,101:301]))

atlas$Ldelta= rowMeans((H[,1:8]))
atlas$Ltheta= rowMeans((H[,9:16]))
atlas$Lalpha= rowMeans((H[,17:26]))
atlas$Lbeta= rowMeans((H[,27:60]))
atlas$Lgamma= rowMeans((H[,61:100]))
atlas$LHgamma= rowMeans((H[,101:301]))
atlas$Boradband= rowMeans((h))

data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])
data4plot$values[data4plot$values< -0.4] = -0.4
data4plot$values[data4plot$values> 1] = 1
data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/H_narrowband.pdf', device = "pdf")


atlas$Boradband[atlas$Boradband< -0.4] = -0.4
atlas$Boradband[atlas$Boradband> 1] = 1
atlas %>%
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/H_broadband.pdf', device = "pdf")



cor.test(rowMeans(X), rowMeans(h))

lm3=lm(scale(rowMeans(X)) ~ scale(rowMeans(h)))
summary(lm3)
sjPlot::tab_model(lm3)

someData_scatter= data.frame(ICC= rowMeans(X), h2=rowMeans(h))

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

ggsave('~/Documents/HCP_twin_projec/figures/scatter_plot_broadband_ICC_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)



lm3=lm(scale(X$alpha) ~ scale(h$alpha))
summary(lm3)
#sjPlot::tab_model(lm3)

someData_scatter= data.frame(ICC= X$alpha, h2=h$alpha)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/scatter_plot_alpha_ICC_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)



lm3=lm(scale(X$beta) ~ scale(h$beta))
summary(lm3)
#sjPlot::tab_model(lm3)

someData_scatter= data.frame(ICC= X$alpha, h2=h$alpha)

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/scatter_plot_beta_ICC_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)


# TO DO: plot delta ICC MZ vs DZ, look at ICC of CamCan


atlas= read.csv('~/Desktop/Destrieux_atlas_test_HCP.csv')

ICCMZ= read.csv('~/Documents/HCP_twin_projec/ICC_only_MZ.csv', header = TRUE)
ICCMZ= ICCMZ[-1]

ICCDZ= read.csv('~/Documents/HCP_twin_projec/ICC_only_DZ.csv', header = TRUE)
ICCDZ= ICCDZ[-1]

ICC= ICCMZ-ICCDZ

# Define X and Y data
X=data.frame(delta= rowMeans(ICC[,1:8]),theta= rowMeans(ICC[,9:16]),alpha= rowMeans(ICC[,17:26]),
             beta= rowMeans(ICC[,27:60]), gamma= rowMeans(ICC[,61:100]), hgamma= rowMeans(ICC[,101:301]))

atlas$Ldelta= rowMeans((ICC[,1:8]))
atlas$Ltheta= rowMeans((ICC[,9:16]))
atlas$Lalpha= rowMeans((ICC[,17:26]))
atlas$Lbeta= rowMeans((ICC[,27:60]))
atlas$Lgamma= rowMeans((ICC[,61:100]))
atlas$LHgamma= rowMeans((ICC[,101:301]))
atlas$Boradband= rowMeans((X))

data4plot=cbind(stack(atlas[,22:27]), atlas[,1:3])

data4plot$values[data4plot$values< -0.5] =0.5
data4plot %>%
  group_by(ind) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_diverging(palette= 'Red-Green', rev= FALSE, limits= c(-0.5,0.5)) +  scale_fill_continuous_diverging(palette= 'Red-Green', rev= FALSE, limits= c(-0.5,0.5)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/ICC_narrowband_MZ_minus_DZ.pdf', device = "pdf")


atlas %>%
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_diverging(palette= 'Red-Green', rev= FALSE, limits= c(-0.5,0.5)) +  scale_fill_continuous_diverging(palette= 'Red-Green', rev= FALSE, limits= c(-0.5,0.5)) +
  theme_void() + scale_color_manual('white')


ggsave('~/Documents/HCP_twin_projec/figures/ICC_broadband_MZ_minus_DZ.pdf', device = "pdf")



# # # # # # # # # # # # # # # # # # # # # # #
## plot the Scahefer atlas
#################################################################

atlas= read.csv('new_Data/Schaefer_atlas.csv')
atlas$new_region=paste(atlas$region, atlas$Yeo, sep= '_')

ICC= read.csv('Results_Log_Schaefer/ICC_and_Heritability/ICC.csv', header = TRUE)
ICC= ICC[-1]
#atlas=atlas[order(atlas$new_region),]

# Define X and Y data
X=data.frame(delta= rowMeans(ICC[,1:8]),theta= rowMeans(ICC[,9:16]),alpha= rowMeans(ICC[,17:26]),
             beta= rowMeans(ICC[,27:60]), gamma= rowMeans(ICC[,61:100]), hgamma= rowMeans(ICC[,101:301]))

atlas$Ldelta= rowMeans((ICC[,1:8]))
atlas$Ltheta= rowMeans((ICC[,9:16]))
atlas$Lalpha= rowMeans((ICC[,17:26]))
atlas$Lbeta= rowMeans((ICC[,27:60]))
atlas$Lgamma= rowMeans((ICC[,61:100]))
atlas$LHgamma= rowMeans((ICC[,101:301]))
atlas$Boradband= rowMeans((X))

data4plot=cbind(stack(atlas[,9:14]), atlas[,c(1:3, 6)])

data4plot$values[data4plot$values<0.4] =0.4
data4plot %>%
  group_by(ind) %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_ICC_narrowband.pdf', device = "pdf")


atlas$Boradband[atlas$Boradband<0.6] =0.6
atlas %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,0.8)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.6,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_ICC_broadband.pdf', device = "pdf")

# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')


ggplot(atlas, aes(Yeo, Boradband, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('fingerprinting features')

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_ICC.pdf', device = "pdf")



# heritability
atlas= read.csv('new_Data/Schaefer_atlas.csv')
atlas$new_region=paste(atlas$region, atlas$Yeo, sep= '_')

H= read.csv('Results_Log_Schaefer/ICC_and_Heritability/heritability_mean.csv', header = TRUE)
H= H[-1]
#atlas=atlas[order(atlas$new_region),]


# Define X and Y data
h=data.frame(delta= rowMeans(H[,1:8]),theta= rowMeans(H[,9:16]),alpha= rowMeans(H[,17:26]),
             beta= rowMeans(H[,27:60]), gamma= rowMeans(H[,61:100]), hgamma= rowMeans(H[,101:301]))

atlas$Ldelta= rowMeans((H[,1:8]))
atlas$Ltheta= rowMeans((H[,9:16]))
atlas$Lalpha= rowMeans((H[,17:26]))
atlas$Lbeta= rowMeans((H[,27:60]))
atlas$Lgamma= rowMeans((H[,61:100]))
atlas$LHgamma= rowMeans((H[,101:301]))
atlas$Boradband= rowMeans((h))

data4plot=cbind(stack(atlas[,9:14]), atlas[,c(1:3, 6)])

data4plot$values[data4plot$values< -0.4] = -0.4
data4plot$values[data4plot$values> 1] = 1
data4plot %>%
  group_by(ind) %>% 
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ ind) +  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_H_narrowband.pdf', device = "pdf")

atlas$Boradband[atlas$Boradband< -0.4] = -0.4
atlas$Boradband[atlas$Boradband> 1] = 1
atlas %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =Boradband)) + 
  geom_sf(show.legend = TRUE) + 
  scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(-0.4,1)) +
  theme_void() + scale_color_manual('white')

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_H_broadband.pdf', device = "pdf")


# set colour palette
cbbPalette= c("#DBF132", "#15A705",'#4F89E0',"#D81B99",'#4d3d87','#115d80', '#FF0000')

ggplot(atlas, aes(Yeo, Boradband, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0.9, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + ylab('heritability')

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_ICC.pdf', device = "pdf")



data4plot=cbind(stack(atlas[,9:14]), atlas[,6:7])
colnames(data4plot)[3]= 'region'

ggplot(data4plot, aes(values, Yeo, fill=Yeo)) + 
  ggdist::stat_halfeye(adjust = .5, width = .1, .width = 0, justification = -.5, alpha=0.5, point_alpha= 0) + 
  geom_boxplot(width = .3, outlier.shape = NA, colour= '#888888') + ggpubr::theme_classic2() + scale_fill_manual(values=cbbPalette) + facet_wrap(~ ind)

ggsave('~/Documents/HCP_twin_projec/figures/Scahefer_box_plot_power_narrowband.pdf', device = "pdf", width = 10, height = 10)




cor.test(rowMeans(X), rowMeans(h))

lm3=lm(scale(rowMeans(X)) ~ scale(rowMeans(h)))
summary(lm3)
sjPlot::tab_model(lm3)

someData_scatter= data.frame(ICC= rowMeans(X), h2=rowMeans(h))

ggplot(someData_scatter, aes(x=ICC , y = h2, fill=cbbPalette[6])) + 
  geom_jitter(colour = cbbPalette[6]) + 
  stat_smooth(method = "lm", colour = 'black',fullrange = T) +
  labs(title = paste("Adj R2 = ", signif(summary(lm3)$adj.r.squared, 5),
                     "Intercept =", signif(lm3$coef[[1]],5 ),
                     " Slope =", signif(lm3$coef[[2]], 5),
                     " P =", signif(summary(lm3)$coef[2,4], 5) )) + scale_fill_manual(values=cbbPalette[6])  +
  theme_classic() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_broadband_ICC_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)

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
  theme_classic() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_alpha_ICC_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)




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
  theme_classic() +   xlab("salient fingerprinting features (ICC)") + ylab("heritability")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_scatter_plot_beta_ICC_vs_h.pdf', device = "pdf", width = 7, height = 7, dpi = 800)




