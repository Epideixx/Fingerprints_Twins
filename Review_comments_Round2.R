library(ggplot2)
library(ggseg)
library(dplyr)
library(colorspace)
library(ggsegSchaefer)
library(ggpubr)

# load in libraries and set up colours
cbbPalette= c("#F0792C", "#f5ec6c",'#4d3d87',"#76D7C4","#268236", '#4d3d87','#593d80', '#3b49e3')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Check how sample size etc biases heritability estimates
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')
atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')
atlas$new_region=paste(atlas$region, atlas$Yeo, sep= '_')

H= read.csv('~/Documents/HCP_twin_projec/heritability_final.csv', header = TRUE)
H= H[-1]
#atlas=atlas[order(atlas$new_region),]

# group data together into 6 cannonical bands
h=data.frame(delta= rowMeans(H[,1:8]),theta= rowMeans(H[,9:16]),alpha= rowMeans(H[,17:26]),
             beta= rowMeans(H[,27:60]), gamma= rowMeans(H[,61:100]), hgamma= rowMeans(H[,101:301]))

atlas$delta= rowMeans((H[,1:8]))
atlas$theta= rowMeans((H[,9:16]))
atlas$alpha= rowMeans((H[,17:26]))
atlas$beta= rowMeans((H[,27:60]))
atlas$gamma= rowMeans((H[,61:100]))
atlas$high_gamma= rowMeans((H[,101:301]))
atlas$Boradband= rowMeans((h))

data4plot=cbind(stack(atlas[,38:43]), atlas[,c(1:3,6)])
data4plot$values= as.numeric(data4plot$values)
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

#ggsave('~/Documents/HCP_twin_projec/figures/Schaefer_H_narrowband_June2023.pdf', device = "pdf")


H_pearson= read.csv('/Users/jason/Downloads/Fingerprints_Twins-main-2/Results_Log_Schaefer_test/Only_GT/PSD_Heritability/pearson/heritability_mean.csv', header = TRUE)
H_pearson= H_pearson[-1]
#atlas=atlas[order(atlas$new_region),]

# Define X and Y data
h_p=data.frame(delta= rowMeans(H_pearson[,1:8]),theta= rowMeans(H_pearson[,9:16]),alpha= rowMeans(H_pearson[,17:26]),
               beta= rowMeans(H_pearson[,27:60]), gamma= rowMeans(H_pearson[,61:100]), hgamma= rowMeans(H_pearson[,101:301]))

diag(cor(h_p, h))

# Read in demographics
subj_HC= read.csv('~/Documents/HCP_twin_projec/outputs/subject_codes.csv')
demo_twin= read.csv('~/Documents/HCP_twin_projec/ID_to_Name.csv', header = TRUE)

demo_hcp= read.csv('~/Documents/HCP_twin_projec/HCP_zygocity_just_key_measures.csv', header = TRUE)
demo2= read.csv('~/Documents/HCP_twin_projec/HCP_demo_test.csv')
demo_hcp$Gender= plyr::mapvalues(demo_hcp$Subject,demo2$Subject, demo2$Gender )

demo_hcp= demo_hcp[demo_hcp$Subject %in%  subj_HC$Ids,]
order=  demo_hcp$Subject %in%  subj_HC$Ids
demo_hcp= demo_hcp[ order(subj_HC$Ids),]

demo_hcp$twin_ID=demo_twin$New_Name[order(demo_twin$ID)]
demo_hcp$Group= substr(demo_hcp$twin_ID, start = 6, stop=7)

demo_hcp$Group = plyr::mapvalues(demo_hcp$Group, c('MZ', 'DZ', 'in'), c('MZ', 'DZ', 'NotTwin'))

demo_hcp_MZ= demo_hcp[demo_hcp$Group=='MZ',]
demo_hcp_MZ= demo_hcp_MZ[order(demo_hcp_MZ$twin_ID),]

demo_hcp_DZ= demo_hcp[demo_hcp$Group=='DZ',]
demo_hcp_DZ= demo_hcp_DZ[order(demo_hcp_DZ$twin_ID),]

demo_hcp_Non= demo_hcp[demo_hcp$Group=='NotTwin',]

# read data and compute data for demographic table
psd_1=read.csv('~/Documents/HCP_twin_projec/outputs/record_1.csv', header = TRUE)
psd_2=read.csv('~/Documents/HCP_twin_projec/outputs/record_2.csv', header = TRUE)
psd_3=read.csv('~/Documents/HCP_twin_projec/outputs/record_3.csv', header = TRUE)

psd_1=log(psd_1)
psd_2=log(psd_2)
psd_3=log(psd_3)

psd_1=psd_1[,-1]
psd_2=psd_2[,-1]
psd_3=psd_3[,-1]

# define function to compute ICC
icc <- function(df, n = 89, k = 2) {
  df_b = n-1
  df_w = n*(k-1)
  
  x_w_mean =  rowMeans(df)
  x_g_mean = mean(unlist(df))
  ss_t = sum(((df - x_g_mean)^2))
  ss_w = sum((df - (x_w_mean))^2)
  ss_b = ss_t - ss_w
  ms_b = ss_b / df_b
  ms_w = ss_w / df_w
  icc = (ms_b - ms_w) / (ms_b + ((k-1)*ms_w))
  
}

# create vector that you can iterate over
demo_hcp$uniqueTwinID=gsub("A", "_A", demo_hcp$twin_ID)
demo_hcp$uniqueTwinID=gsub("B", "_B", demo_hcp$uniqueTwinID)

order_of_filesMZ= c("Twin_MZ_1_A",  "Twin_MZ_1_B", "Twin_MZ_2_A",  "Twin_MZ_2_B" , "Twin_MZ_3_A",  "Twin_MZ_3_B",  "Twin_MZ_4_A",  "Twin_MZ_4_B" ,
                    "Twin_MZ_5_A",  "Twin_MZ_5_B",  "Twin_MZ_6_A" , "Twin_MZ_6_B",  "Twin_MZ_7_A",  "Twin_MZ_7_B" ,
                    "Twin_MZ_8_A",  "Twin_MZ_8_B" , "Twin_MZ_9_A" , "Twin_MZ_9_B",  "Twin_MZ_10_A", "Twin_MZ_10_B", "Twin_MZ_11_A", "Twin_MZ_11_B",
                    "Twin_MZ_12_A", "Twin_MZ_12_B", "Twin_MZ_13_A" ,"Twin_MZ_13_B" ,"Twin_MZ_14_A", "Twin_MZ_14_B",
                    "Twin_MZ_15_A" ,"Twin_MZ_15_B" ,"Twin_MZ_16_A", "Twin_MZ_16_B", "Twin_MZ_18_A", "Twin_MZ_18_B"
)

order_of_filesDZ= c("Twin_DZ_1_A",  "Twin_DZ_1_B",  "Twin_DZ_2_A",  "Twin_DZ_2_B",  "Twin_DZ_3_A",  "Twin_DZ_3_B", 
                    "Twin_DZ_5_A",  "Twin_DZ_5_B" , "Twin_DZ_6_A",  "Twin_DZ_6_B",  "Twin_DZ_7_A",  "Twin_DZ_7_B",
                    "Twin_DZ_8_A",  "Twin_DZ_8_B",  "Twin_DZ_9_A",  "Twin_DZ_9_B", "Twin_DZ_11_A",
                    "Twin_DZ_11_B", "Twin_DZ_12_A", "Twin_DZ_12_B","Twin_DZ_13_A", "Twin_DZ_13_B" )

# make an index of MZ twin A and Bs 
temp=match(order_of_filesMZ,demo_hcp$uniqueTwinID)
MZ_order=temp
twin_MZ_A= temp[seq(1,34,2)]
twin_MZ_B= temp[seq(2,34,2)]

temp=match(order_of_filesDZ,demo_hcp$uniqueTwinID)
DZ_order=temp
twin_DZ_A= temp[seq(1,22,2)]
twin_DZ_B= temp[seq(2,22,2)]

output_MZ= matrix(, nrow = 60200, ncol=1)
output_DZ= matrix(, nrow = 60200, ncol=1)

# zscore all features before computing ICC
broadband_1_Mz= t(scale(t(psd_1[MZ_order,])))
broadband_2_Mz= t(scale(t(psd_2[MZ_order,])))
broadband_3_Mz= t(scale(t(psd_3[MZ_order,])))

# zscore all features before computing ICC
broadband_1_Dz= t(scale(t(psd_1[DZ_order,])))
broadband_2_Dz= t(scale(t(psd_2[DZ_order,])))
broadband_3_Dz= t(scale(t(psd_3[DZ_order,])))

# iterate acorss all features to compute heritability 
for( i in 1:60200){
  
  #. Monozygotic twins
  x1= c(broadband_1_Mz[seq(1,34,2),i],broadband_2_Mz[seq(1,34,2),i], broadband_3_Mz[seq(1,34,2),i],
        broadband_1_Mz[seq(1,34,2),i],broadband_2_Mz[seq(1,34,2),i], broadband_1_Mz[seq(1,34,2),i],
        broadband_3_Mz[seq(1,34,2),i],broadband_2_Mz[seq(1,34,2),i], broadband_3_Mz[seq(1,34,2),i])
  x2= c(broadband_1_Mz[seq(2,34,2),i],broadband_2_Mz[seq(2,34,2),i], broadband_3_Mz[seq(2,34,2),i],
        broadband_2_Mz[seq(2,34,2),i],broadband_1_Mz[seq(2,34,2),i], broadband_3_Mz[seq(2,34,2),i],
        broadband_1_Mz[seq(2,34,2),i],broadband_3_Mz[seq(2,34,2),i], broadband_2_Mz[seq(2,34,2),i])
  df_temp= data.frame(x1= x1, x2=x2)
  output_MZ[i,1]=icc(df_temp, 9*length(twin_MZ_B), dim(df_temp)[2])
  
  #. dizygotic twins
  x1= c(broadband_1_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_3_Dz[seq(1,22,2),i],
        broadband_1_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_1_Dz[seq(1,22,2),i],
        broadband_3_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_3_Dz[seq(1,22,2),i])
  x2= c(broadband_1_Dz[seq(2,22,2),i],broadband_2_Dz[seq(2,22,2),i], broadband_3_Dz[seq(2,22,2),i],
        broadband_2_Dz[seq(2,22,2),i],broadband_1_Dz[seq(2,22,2),i], broadband_3_Dz[seq(2,22,2),i],
        broadband_1_Dz[seq(2,22,2),i],broadband_3_Dz[seq(2,22,2),i], broadband_2_Dz[seq(2,22,2),i])
  df_temp= data.frame(x1= x1, x2=x2)
  output_DZ[i,1]=icc(df_temp, 9*length(twin_DZ_B), dim(df_temp)[2])
  
  
}

heritability = rowMeans(2*(output_MZ-output_DZ))
heritability_mat<- matrix(heritability, nrow = 200, byrow = TRUE)

write.csv(heritability_mat, '/Users/jason/Documents/HCP_twin_projec/heritability_Round2.csv')


# Now let us redo this process with a subsample of MZ to match the size of DZ
temp=match(order_of_filesMZ,demo_hcp$uniqueTwinID)
MZ_order=temp[1:22]
twin_MZ_A= temp[seq(1,22,2)]
twin_MZ_B= temp[seq(2,22,2)]

temp=match(order_of_filesDZ,demo_hcp$uniqueTwinID)
DZ_order=temp[1:22]
twin_DZ_A= temp[seq(1,22,2)]
twin_DZ_B= temp[seq(2,22,2)]

output_MZ= matrix(, nrow = 60200, ncol=1)
output_DZ= matrix(, nrow = 60200, ncol=1)

# zscore all features before computing ICC
broadband_1_Mz= t(scale(t(psd_1[MZ_order,])))
broadband_2_Mz= t(scale(t(psd_2[MZ_order,])))
broadband_3_Mz= t(scale(t(psd_3[MZ_order,])))

# zscore all features before computing ICC
broadband_1_Dz= t(scale(t(psd_1[DZ_order,])))
broadband_2_Dz= t(scale(t(psd_2[DZ_order,])))
broadband_3_Dz= t(scale(t(psd_3[DZ_order,])))

# compute heritability
for( i in 1:60200){
  
  #. Monozygotic twins
  x1= c(broadband_1_Mz[seq(1,22,2),i],broadband_2_Mz[seq(1,22,2),i], broadband_3_Mz[seq(1,22,2),i],
        broadband_1_Mz[seq(1,22,2),i],broadband_2_Mz[seq(1,22,2),i], broadband_1_Mz[seq(1,22,2),i],
        broadband_3_Mz[seq(1,22,2),i],broadband_2_Mz[seq(1,22,2),i], broadband_3_Mz[seq(1,22,2),i])
  x2= c(broadband_1_Mz[seq(2,22,2),i],broadband_2_Mz[seq(2,22,2),i], broadband_3_Mz[seq(2,22,2),i],
        broadband_2_Mz[seq(2,22,2),i],broadband_1_Mz[seq(2,22,2),i], broadband_3_Mz[seq(2,22,2),i],
        broadband_1_Mz[seq(2,22,2),i],broadband_3_Mz[seq(2,22,2),i], broadband_2_Mz[seq(2,22,2),i])
  df_temp= data.frame(x1= x1, x2=x2)
  output_MZ[i,1]=icc(df_temp, 9*length(twin_MZ_B), dim(df_temp)[2])
  
  #. dizygotic twins
  x1= c(broadband_1_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_3_Dz[seq(1,22,2),i],
        broadband_1_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_1_Dz[seq(1,22,2),i],
        broadband_3_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_3_Dz[seq(1,22,2),i])
  x2= c(broadband_1_Dz[seq(2,22,2),i],broadband_2_Dz[seq(2,22,2),i], broadband_3_Dz[seq(2,22,2),i],
        broadband_2_Dz[seq(2,22,2),i],broadband_1_Dz[seq(2,22,2),i], broadband_3_Dz[seq(2,22,2),i],
        broadband_1_Dz[seq(2,22,2),i],broadband_3_Dz[seq(2,22,2),i], broadband_2_Dz[seq(2,22,2),i])
  df_temp= data.frame(x1= x1, x2=x2)
  output_DZ[i,1]=icc(df_temp, 9*length(twin_DZ_B), dim(df_temp)[2])
  
  
}

heritability_sub = rowMeans(2*(output_MZ-output_DZ))
heritability_sub<- matrix(heritability_sub, nrow = 200, byrow = TRUE)

# group together data within 6 cannonical bands
h=data.frame(delta= rowMeans(heritability[,1:8]),theta= rowMeans(heritability[,9:16]),alpha= rowMeans(heritability[,17:26]),
             beta= rowMeans(heritability[,27:60]), gamma= rowMeans(heritability[,61:100]), hgamma= rowMeans(heritability[,101:301]))

hsub=data.frame(delta= rowMeans(heritability_sub[,1:8]),theta= rowMeans(heritability_sub[,9:16]),alpha= rowMeans(heritability_sub[,17:26]),
                beta= rowMeans(heritability_sub[,27:60]), gamma= rowMeans(heritability_sub[,61:100]), hgamma= rowMeans(heritability_sub[,101:301]))
# look at correlation
diag(cor(h, hsub))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# permute data to get heritability  p values
# make an index of MZ twin A and Bs 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
heritability_mat= read.csv('~/Documents/HCP_twin_projec/heritability_Round2.csv', header = TRUE)
heritability_mat= heritability_mat[-1]

heritability_band=as.matrix(data.frame(delta= rowMeans(heritability_mat[,1:8]),theta= rowMeans(heritability_mat[,9:16]),alpha= rowMeans(heritability_mat[,17:26]),
                                       beta= rowMeans(heritability_mat[,27:60]), gamma= rowMeans(heritability_mat[,61:100]), hgamma= rowMeans(heritability_mat[,101:301])))

temp=match(order_of_filesMZ,demo_hcp$uniqueTwinID)
MZ_order=temp
twin_MZ_A= temp[seq(1,34,2)]
twin_MZ_B= temp[seq(2,34,2)]

temp=match(order_of_filesDZ,demo_hcp$uniqueTwinID)
DZ_order=temp
twin_DZ_A= temp[seq(1,22,2)]
twin_DZ_B= temp[seq(2,22,2)]

heritability_perm= matrix(, nrow = 500, ncol=60200)
heritability_perm_bands= array(0, dim=c(500,200,6))

# run 500 permutations of the heritability estimate 
for (iperm in 1:500){
  
  output_MZ= matrix(, nrow = 60200, ncol=1)
  output_DZ= matrix(, nrow = 60200, ncol=1)
  
  temp=sample(1:17,11)
  DZ_order2=DZ_order
  MZ_order2=DZ_order2
  MZ_order2[seq(1,22,2)]=MZ_order[c(2*(temp)-1)]
  MZ_order2[seq(2,22,2)]=MZ_order[c(2*(temp))]
  
  # get random perm of data
  # permute labels of MZ and DZ PAIRS! 
  Dz_rand= sample(1:11,6)
  Mz_rand= sample(1:11,6)
  
  Dz_rand=c(2*(Dz_rand)-1, 2*Dz_rand)
  Mz_rand=c(2*(Mz_rand)-1, 2*Mz_rand)
  
  MZ_order_perm= MZ_order2
  MZ_order_perm[Mz_rand] = DZ_order2[Dz_rand]
  
  DZ_order_perm= DZ_order2
  DZ_order_perm[Dz_rand] = MZ_order2[Mz_rand]
  
  # zscore all features before computing ICC
  broadband_1_Mz= t(scale(t(psd_1[MZ_order_perm,])))
  broadband_2_Mz= t(scale(t(psd_2[MZ_order_perm,])))
  broadband_3_Mz= t(scale(t(psd_3[MZ_order_perm,])))
  
  # zscore all features before computing ICC
  broadband_1_Dz= t(scale(t(psd_1[DZ_order_perm,])))
  broadband_2_Dz= t(scale(t(psd_2[DZ_order_perm,])))
  broadband_3_Dz= t(scale(t(psd_3[DZ_order_perm,])))
  
  # compute heritability for two new permuted labels 
  for( i in 1:60200){
    
    #. Monozygotic twins
    x1= c(broadband_1_Mz[seq(1,22,2),i],broadband_2_Mz[seq(1,22,2),i], broadband_3_Mz[seq(1,22,2),i],
          broadband_1_Mz[seq(1,22,2),i],broadband_2_Mz[seq(1,22,2),i], broadband_1_Mz[seq(1,22,2),i],
          broadband_3_Mz[seq(1,22,2),i],broadband_2_Mz[seq(1,22,2),i], broadband_3_Mz[seq(1,22,2),i])
    x2= c(broadband_1_Mz[seq(2,22,2),i],broadband_2_Mz[seq(2,22,2),i], broadband_3_Mz[seq(2,22,2),i],
          broadband_2_Mz[seq(2,22,2),i],broadband_1_Mz[seq(2,22,2),i], broadband_3_Mz[seq(2,22,2),i],
          broadband_1_Mz[seq(2,22,2),i],broadband_3_Mz[seq(2,22,2),i], broadband_2_Mz[seq(2,22,2),i])
    df_temp= data.frame(x1= x1, x2=x2)
    output_MZ[i,1]=icc(df_temp, 9*length(twin_DZ_B), dim(df_temp)[2])
    
    #. dizygotic twins
    x1= c(broadband_1_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_3_Dz[seq(1,22,2),i],
          broadband_1_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_1_Dz[seq(1,22,2),i],
          broadband_3_Dz[seq(1,22,2),i],broadband_2_Dz[seq(1,22,2),i], broadband_3_Dz[seq(1,22,2),i])
    x2= c(broadband_1_Dz[seq(2,22,2),i],broadband_2_Dz[seq(2,22,2),i], broadband_3_Dz[seq(2,22,2),i],
          broadband_2_Dz[seq(2,22,2),i],broadband_1_Dz[seq(2,22,2),i], broadband_3_Dz[seq(2,22,2),i],
          broadband_1_Dz[seq(2,22,2),i],broadband_3_Dz[seq(2,22,2),i], broadband_2_Dz[seq(2,22,2),i])
    df_temp= data.frame(x1= x1, x2=x2)
    output_DZ[i,1]=icc(df_temp, 9*length(twin_DZ_B), dim(df_temp)[2])
    
    
  }
  
  heritability_perm[iperm,] = rowMeans(2*(output_MZ-output_DZ))
  temph=matrix(heritability_perm[iperm,], nrow = 200, byrow = TRUE)
  
  # for computational efficency save only the averages within a band
  heritability_perm_bands[iperm,,]=as.matrix(data.frame(delta= rowMeans(temph[,1:8]),theta= rowMeans(temph[,9:16]),alpha= rowMeans(temph[,17:26]),
                                                        beta= rowMeans(temph[,27:60]), gamma= rowMeans(temph[,61:100]), hgamma= rowMeans(temph[,101:301])))
  
  print(iperm)
}

# save data as csv because we worked so hard!
write.csv(heritability_perm_bands, '/Users/jason/Documents/HCP_twin_projec/permutations_of_heritability_bands_Round2.csv')

heritability_band=as.matrix(data.frame(delta= rowMeans(heritability_mat[,1:8]),theta= rowMeans(heritability_mat[,9:16]),alpha= rowMeans(heritability_mat[,17:26]),
                                       beta= rowMeans(heritability_mat[,27:60]), gamma= rowMeans(heritability_mat[,61:100]), hgamma= rowMeans(heritability_mat[,101:301])))

# now compute p values by assessing if heritability is larger (one tailed test) 
p_values= array(0, dim=c(200,6))
for ( b in 1:6){
  p_values[,b]=apply(heritability_perm_bands[,,b] > heritability_band[,b],2,sum)/500
}

write.csv(p_values, '/Users/jason/Documents/HCP_twin_projec/permutations_of_heritability_bands_pvalues_Round2.csv')

# now let us plot result over the brain
atlas= read.csv('~/Desktop/Schaefer_neuromaps/Schaefer2018_200Parcels_7Networks_Neuromaps.csv')
atlas$delta_pval= p_values[,1]
atlas$theta_pval= p_values[,2]
atlas$alpha_pval= p_values[,3]
atlas$beta_pval= p_values[,4]
atlas$gamma_pval= p_values[,5]
atlas$Hgamma_pval= p_values[,6]

data4plot=cbind(stack(atlas[,37:42]), atlas[,6:7])
colnames(data4plot)[3]= 'region'


data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", direction=-1) +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Review_heritability_pval_Round2.pdf', device = "pdf")


data4plot$values[data4plot$values>0.05] =0.07
data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", limits=c(0, 0.07), direction=-1) +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Review_heritability_pval_threshold_Round2.pdf', device = "pdf")


data4plot$values[data4plot$values>0.05] =1
data4plot$values[data4plot$values<=0.05] =0.0
data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =values)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", limits=c(0, 1), direction=-1) +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Review_heritability_pval_threshold_2_Round2.pdf', device = "pdf")



# let us plot the orginal heritability estimates but thresholded by the p values 
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

data4plot=cbind(stack(atlas[,37:42]), stack(atlas[,43:48]),atlas[,6:7])
colnames(data4plot)= c('pvalues', 'ind', 'heritability', 'band', 'region', "x.mni" )

data4plot$heritability[data4plot$pvalues>0.05] =NaN

data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =heritability)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma") +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Review_heritability_pval_threshold_3__Round2.pdf', device = "pdf")



data4plot %>% group_by(ind) %>%
  brain_join(schaefer7_200) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =heritability)) + 
  geom_sf(show.legend = TRUE) + viridis::scale_fill_viridis(option="magma", limits= c(0, 1.5)) +  #scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +  scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE, limits= c(0.4,1)) +
  theme_void() + scale_color_manual('white') + facet_wrap( ~ ind )

ggsave('~/Documents/HCP_twin_projec/figures/Review_heritability_pval_threshold_4__Round2.pdf', device = "pdf")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#   Match two random unrelated individuals and compute accuracy 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# get data for non twins
psd_1_nontwin=psd_1[demo_hcp$Group == 'NotTwin',]
psd_2_nontwin=psd_2[demo_hcp$Group == 'NotTwin',]
psd_3_nontwin=psd_3[demo_hcp$Group == 'NotTwin',]


# randomly pair up with different reocrding and compute matching Acc
Acc_1= c()
Acc_2= c()
for (i in 1:100){
corr_psd_NonTwin=cor(t(psd_1_nontwin[,]),t(psd_2_nontwin[sample(1:25,25),]))
tt=apply(corr_psd_NonTwin, 1, which.max)
Acc_1[i]=sum(seq(1:25)==tt)/25

tt=apply(corr_psd_NonTwin, 2, which.max)
Acc_2[i]=sum(seq(1:25)==tt)/25

}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,]),t(psd_3_nontwin[sample(1:25,25),]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,]),t(psd_1_nontwin[sample(1:25,25),]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2[i]=sum(seq(1:25)==tt)/25
  
}

index_delta=unlist(lapply(1:200, function(x) ((x-1)*301)+(1:8)))
index_theta=unlist(lapply(1:200, function(x) ((x-1)*301)+(9:16))) # theta
index_alpha=unlist(lapply(1:200, function(x) ((x-1)*301)+(17:26))) # alpha
index_beta=unlist(lapply(1:200, function(x) ((x-1)*301)+(27:60))) # beta
index_gamma=unlist(lapply(1:200, function(x) ((x-1)*301)+(61:100))) # gamma
index_hgamma=unlist(lapply(1:200, function(x) ((x-1)*301)+(101:301))) # Hgamma


Acc_1_delta= c()
Acc_2_delta= c()
for (i in 1:100){
  corr_psd_NonTwin=cor(t(psd_1_nontwin[,index_delta]),t(psd_2_nontwin[sample(1:25,25),index_delta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_delta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_delta[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,index_delta]),t(psd_3_nontwin[sample(1:25,25),index_delta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_delta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_delta[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,index_delta]),t(psd_1_nontwin[sample(1:25,25),index_delta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_delta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_delta[i]=sum(seq(1:25)==tt)/25
  
}


Acc_1_theta= c()
Acc_2_theta= c()
for (i in 1:100){
  corr_psd_NonTwin=cor(t(psd_1_nontwin[,index_theta]),t(psd_2_nontwin[sample(1:25,25),index_theta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_theta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_theta[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,index_theta]),t(psd_3_nontwin[sample(1:25,25),index_theta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_theta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_theta[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,index_theta]),t(psd_1_nontwin[sample(1:25,25),index_theta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_theta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_theta[i]=sum(seq(1:25)==tt)/25
  
}


Acc_1_alpha= c()
Acc_2_alpha= c()
for (i in 1:100){
  corr_psd_NonTwin=cor(t(psd_1_nontwin[,index_alpha]),t(psd_2_nontwin[sample(1:25,25),index_alpha]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_alpha[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_alpha[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,index_alpha]),t(psd_3_nontwin[sample(1:25,25),index_alpha]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_alpha[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_alpha[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,index_alpha]),t(psd_1_nontwin[sample(1:25,25),index_alpha]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_alpha[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_alpha[i]=sum(seq(1:25)==tt)/25
  
}



Acc_1_beta= c()
Acc_2_beta= c()
for (i in 1:100){
  corr_psd_NonTwin=cor(t(psd_1_nontwin[,index_beta]),t(psd_2_nontwin[sample(1:25,25),index_beta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_beta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_beta[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,index_beta]),t(psd_3_nontwin[sample(1:25,25),index_beta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_beta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_beta[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,index_beta]),t(psd_1_nontwin[sample(1:25,25),index_beta]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_beta[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_beta[i]=sum(seq(1:25)==tt)/25
  
}


Acc_1_gamma= c()
Acc_2_gamma= c()
for (i in 1:100){
  corr_psd_NonTwin=cor(t(psd_1_nontwin[,index_gamma]),t(psd_2_nontwin[sample(1:25,25),index_gamma]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_gamma[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_gamma[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,index_gamma]),t(psd_3_nontwin[sample(1:25,25),index_gamma]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_gamma[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_gamma[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,index_gamma]),t(psd_1_nontwin[sample(1:25,25),index_gamma]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_gamma[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_gamma[i]=sum(seq(1:25)==tt)/25
  
}


Acc_1_hgamma= c()
Acc_2_hgamma= c()
for (i in 1:100){
  corr_psd_NonTwin=cor(t(psd_1_nontwin[,index_hgamma]),t(psd_2_nontwin[sample(1:25,25),index_hgamma]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_hgamma[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_hgamma[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 101:200){
  corr_psd_NonTwin=cor(t(psd_2_nontwin[,index_hgamma]),t(psd_3_nontwin[sample(1:25,25),index_hgamma]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_hgamma[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_hgamma[i]=sum(seq(1:25)==tt)/25
  
}

for (i in 201:300){
  corr_psd_NonTwin=cor(t(psd_3_nontwin[,index_hgamma]),t(psd_1_nontwin[sample(1:25,25),index_hgamma]))
  tt=apply(corr_psd_NonTwin, 1, which.max)
  Acc_1_hgamma[i]=sum(seq(1:25)==tt)/25
  
  tt=apply(corr_psd_NonTwin, 2, which.max)
  Acc_2_hgamma[i]=sum(seq(1:25)==tt)/25
  
}

# --------------------------------------
#         PLOT ACCURACIES
# --------------------------------------

data4plot_acc= data.frame( acc = c(mean(c(Acc_1, Acc_2)),
                                   mean(c(Acc_1_delta, Acc_2_delta)),
                                   mean(c(Acc_1_theta, Acc_2_theta)),
                                   mean(c(Acc_1_alpha, Acc_2_alpha)),
                                   mean(c(Acc_1_beta, Acc_2_beta)),
                                   mean(c(Acc_1_gamma, Acc_2_gamma)),
                                   mean(c(Acc_1_hgamma, Acc_2_hgamma)) ),
                           lower= c(quantile(c(Acc_1, Acc_2), c(0.025, 0.975))[1],
                                 quantile(c(Acc_1_delta, Acc_2_delta), c(0.025, 0.975))[1],
                                 quantile(c(Acc_1_theta, Acc_2_theta), c(0.025, 0.975))[1],
                                 quantile(c(Acc_1_alpha, Acc_2_alpha), c(0.025, 0.975))[1],
                                 quantile(c(Acc_1_beta, Acc_2_beta), c(0.025, 0.975))[1],
                                 quantile(c(Acc_1_gamma, Acc_2_gamma), c(0.025, 0.975))[1],
                                 quantile(c(Acc_1_hgamma, Acc_2_hgamma), c(0.025, 0.975))[1] ),
                           upper= c(quantile(c(Acc_1, Acc_2), c(0.025, 0.975))[2],
                                 quantile(c(Acc_1_delta, Acc_2_delta), c(0.025, 0.975))[2],
                                 quantile(c(Acc_1_theta, Acc_2_theta), c(0.025, 0.975))[2],
                                 quantile(c(Acc_1_alpha, Acc_2_alpha), c(0.025, 0.975))[2],
                                 quantile(c(Acc_1_beta, Acc_2_beta), c(0.025, 0.975))[2],
                                 quantile(c(Acc_1_gamma, Acc_2_gamma), c(0.025, 0.975))[2],
                                 quantile(c(Acc_1_hgamma, Acc_2_hgamma), c(0.025, 0.975))[2] ),
                           band = c('broadband', 'delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))

data4plot_acc$band = factor(data4plot_acc$band, levels = c('broadband', 'delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'))

ggplot(data=data4plot_acc, aes(x=band, y=acc, fill='1',  width=.5)) +  coord_cartesian(ylim=c(0,1))+
  geom_bar(stat = "identity", position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2 ,position=position_dodge(.5))  +
  theme_classic2() + scale_fill_manual(values=cbbPalette) + labs(y="Accuracy (%)", x="",title="") + 
  theme(text = element_text(size=40), legend.title=element_blank())

ggsave('~/Documents/HCP_twin_projec/figures/Review_AccuracyNonTwin_Round2.pdf', device = "pdf", width =10, height =5)


setwd('/Users/jason/Downloads/Fingerprints_Twins-main-2/')
# Import the data
data <- read.csv("Results_Log_Schaefer_test/Only_GT/PSD_Accuracy/All_accuracies_every_freq_stacked.csv")

# Define the order of the frequency bands
freq_order <- c("BROADBAND", "DELTA", "THETA", "ALPHA", "BETA", "GAMMA", "HIGH GAMMA")

# Define the order of the account types
acc_order <- c("Autocorr", "Crosscorr MZ", "Crosscorr DZ")


# ------ Bar plot -------

# Calculate means and standard errors
means <- aggregate(values ~ FreqBand + AccType, data, mean)
se <- aggregate(values ~ FreqBand + AccType, data, function(x) sd(x)/sqrt(length(x)))
sd <- aggregate(values ~ FreqBand + AccType, data, function(x) sd(x))
conf_int_02_5 <- aggregate(values ~ FreqBand + AccType, data, function(x) quantile(x, c(.025)))
conf_int_97_5 <- aggregate(values ~ FreqBand + AccType, data, function(x) quantile(x, c(.975)))

means$se = se$values
means$sd = sd$values 
means$conf_int_02_5 = conf_int_02_5$values
means$conf_int_97_5 = conf_int_97_5$values

means$values[ means$AccType == 'Autocorr'] = data4plot_acc$acc
means$conf_int_02_5[ means$AccType == 'Autocorr'] = data4plot_acc$lower
means$conf_int_97_5[ means$AccType == 'Autocorr'] = data4plot_acc$upper

# Create the bar plot with error bars
ggplot(data = means, aes(x = factor(FreqBand, levels = freq_order), y = values, fill = factor(AccType, levels = acc_order))) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = conf_int_02_5, ymax = conf_int_97_5), position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_manual(values = magma_3, labels = c("Individual differentiation", "MZ pair differentiation", "DZ pair differentiation")) +
  ylim(-0.01, 1) +
  labs(x = "Frequency Band", y = "Accuracy", color = "", fill = "") +
  theme_classic() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 10), panel.background = element_rect(fill='white', colour='black'))

# Define the output file path and name
out_file <- "Results_Log_Schaefer/Only_GT/Figures/accuracy_bar_plot_Review.pdf"

# Save the plot as PDF
ggsave(filename = out_file, width = 10, height = 6, units = "in")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# as per suggestion, run permutations over the alignemnt between GO categories and ICC
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# read in ICC
ICCb= read.csv('~/Documents/HCP_twin_projec/ICC_bootstrap_avg_total.csv', header = TRUE)
ICCb= ICCb[-1]

# average within cannonical bands
Xb=data.frame(delta= rowMeans(ICCb[,1:8]),theta= rowMeans(ICCb[,9:16]),alpha= rowMeans(ICCb[,17:26]),
              beta= rowMeans(ICCb[,27:60]), gamma= rowMeans(ICCb[,61:100]), hgamma= rowMeans(ICCb[,101:301]))



## Run analysis where we permute the ICC and corr with each GO category
# read in genetic data obtained from abagen
gene_expression= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/All_genes_expression_genes.csv')

index= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/row_ids.csv', header = FALSE)

gene_expression=gene_expression[index$V1,]
gene_names= read.csv('/Users/jason/Documents/HCP_twin_projec/top50_per_gene_loadings.csv')

pos_gene_expression= rowMeans(gene_expression[,which(colnames(gene_expression) %in% gene_names$pos_genes)])

neg_gene_expression= rowMeans(gene_expression[,which(colnames(gene_expression) %in% gene_names$neg_genes)])

permuted_index= read.csv('/Users/jason/Desktop/Schaefer_neuromaps/permuted_indexes_of_schaefer_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1

GO_negative= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/negative_top_GO_biological_process_terms.csv')

GO_negative$pvals_permuted=1

# fot the top negative categories 
for (g in 1:35){
  # get the avg gene expression across cortex (spatial map) for these genes GO
  list_temp=unique(unlist(strsplit(GO_negative$Genes[g], split = " "))) 
  colnums2Avg=which(colnames(gene_expression) %in% list_temp)
  gene_express_obs= rowMeans(gene_expression[,colnums2Avg])
  
  # get the obsered correlation
  orig=cor(gene_express_obs, rowMeans(Xb[,-2]))
  # run permutations based on the SPIN tests
  permuted_corr=c()
  for (i in 1:1001){
    
    cor_temp=cor.test(rowMeans(Xb[permuted_index[,i],-2]), gene_express_obs)
    permuted_corr= c(permuted_corr, cor_temp$estimate)
    
  }
  
  # see if the observed alignemnt between the GO category is greater than the SPINS
  # note one tailed test here since we know the direction of the effect already (i.e., negative)
  GO_negative$pvals_permuted[g]=(sum(orig > permuted_corr)+1)/1001
}

# rinse and repeat for the positive genes 
GO_positive= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/positive_top_GO_biological_process_terms.csv')
GO_positive$pvals_permuted=1

for (g in 1:35){
  list_temp=unique(unlist(strsplit(GO_positive$Genes[g], split = " ")))
  colnums2Avg=which(colnames(gene_expression) %in% list_temp)
  gene_express_obs= rowMeans(gene_expression[,colnums2Avg])
  
  orig=cor(gene_express_obs, rowMeans(Xb[,-2]))
  permuted_corr=c()
  for (i in 1:1001){
    
    cor_temp=cor.test(rowMeans(Xb[permuted_index[,i],-2]), gene_express_obs)
    permuted_corr= c(permuted_corr, cor_temp$estimate)
    
  }
  
  GO_positive$pvals_permuted[g]=(sum(orig < permuted_corr)+1)/1001
}

# now that we have the pvalues let us plot the effects and make it pretty 

GO_negative$Pathway = factor(GO_negative$Pathway, levels = GO_negative$Pathway[order(GO_negative$FoldEnrichment)])

ggplot(GO_negative, aes(x= log10(FoldEnrichment), y=Pathway, colour = log10(pvals_permuted), fill = log10(pvals_permuted), size = abs(log10(EnrichmentFDR)))) +
  geom_point() + ggpubr::theme_classic2() + scale_fill_gradient(low= '#40E0D0', high = '#37706a') +
  scale_color_gradient(low= '#40E0D0', high = '#37706a')

ggsave('~/Documents/HCP_twin_projec/figures/negativegenesTop_reviewer.pdf', device = "pdf", width= 20, height=9)


GO_positive$Pathway = factor(GO_positive$Pathway, levels = GO_positive$Pathway[order(GO_positive$FoldEnrichment)])

ggplot(GO_positive, aes(x= log10(FoldEnrichment), y=Pathway, colour = log10(pvals_permuted), fill = log10(pvals_permuted), size = abs(log10(EnrichmentFDR)))) +
  geom_point() + ggpubr::theme_classic2() + scale_fill_gradient(low= '#F315AD', high = '#61214d') +
  scale_color_gradient(low= '#F315AD', high = '#61214d')

ggsave('~/Documents/HCP_twin_projec/figures/positivegenesTop_reviewer.pdf', device = "pdf", width= 20, height=9)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# plots of self- other- and twin- similarity of brain profiles 
# do this for the broadband data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# need to compute histogram of within and between corrs of subject 

corr_psd=cor(t(psd_1),t(psd_2)) 
demo_hcp$psd_self_corr12=diag(corr_psd)
demo_hcp$psd_identifiability12=diag(apply(corr_psd, 2, scale))

corr_psd_13=cor(t(psd_1),t(psd_3)) 
demo_hcp$psd_self_corr13=diag(corr_psd_13)
demo_hcp$psd_identifiability13=diag(apply(corr_psd_13, 2, scale))

corr_psd_23=cor(t(psd_2),t(psd_3)) 
demo_hcp$psd_self_corr23=diag(corr_psd_23)
demo_hcp$psd_identifiability23=diag(apply(corr_psd_23, 2, scale))

data4plot=cbind(stack(demo_hcp[, c(41, 43, 45)]), demo_hcp$Group)
colnames(data4plot)[3]= 'group'

ggplot(data4plot, aes(DescTools::FisherZ(values), color= group, fill=group, alpha=0.5)) +
  geom_density(aes( color= group, fill=group), alpha= 0.5) +
  #geom_histogram(aes(y=..density.., color= group, fill=group), alpha= 0.5, bins = 20) +
  theme_minimal() +
  labs(x = "self correlation", y = "denstity")  + theme_classic2() +
  ggtitle("") + scale_fill_manual(values=cbbPalette[c(1,2,3)]) +
  scale_color_manual(values=cbbPalette[c(1,2,3)]) 

ggsave('~/Documents/HCP_twin_projec/figures/Review_hist_self_corr_broadband.pdf', device = "pdf")


temp=match(order_of_filesMZ,demo_hcp$uniqueTwinID)
corr_psd_MZ=cor(t(psd_1[temp,]),t(psd_2[temp,])) 
corss_corrs_MZ_12= c(diag(corr_psd_MZ[seq(1,34,2), seq(2,34,2)]), diag(corr_psd_MZ[seq(2,34,2), seq(1,34,2)]))
corr_psd_MZ[seq(1,34,2), seq(2,34,2)]= NaN
corr_psd_MZ[seq(2,34,2), seq(1,34,2)]= NaN
diag(corr_psd_MZ)= NaN
other_corrs_MZ_12= c(corr_psd_MZ)

corr_psd_MZ=cor(t(psd_2[temp,]),t(psd_3[temp,])) 
corss_corrs_MZ_23= c(diag(corr_psd_MZ[seq(1,34,2), seq(2,34,2)]), diag(corr_psd_MZ[seq(2,34,2), seq(1,34,2)]))
corr_psd_MZ[seq(1,34,2), seq(2,34,2)]= NaN
corr_psd_MZ[seq(2,34,2), seq(1,34,2)]= NaN
diag(corr_psd_MZ)= NaN
other_corrs_MZ_23= c(corr_psd_MZ)

corr_psd_MZ=cor(t(psd_1[temp,]),t(psd_3[temp,])) 
corss_corrs_MZ_13= c(diag(corr_psd_MZ[seq(1,34,2), seq(2,34,2)]), diag(corr_psd_MZ[seq(2,34,2), seq(1,34,2)]))
corr_psd_MZ[seq(1,34,2), seq(2,34,2)]= NaN
corr_psd_MZ[seq(2,34,2), seq(1,34,2)]= NaN
diag(corr_psd_MZ)= NaN
other_corrs_MZ_13= c(corr_psd_MZ)


temp=match(order_of_filesDZ,demo_hcp$uniqueTwinID)
corr_psd_DZ=cor(t(psd_1[temp,]),t(psd_2[temp,])) 
corss_corrs_DZ_12= c(diag(corr_psd_DZ[seq(1,22,2), seq(2,22,2)]), diag(corr_psd_DZ[seq(2,22,2), seq(1,22,2)]))
corr_psd_DZ[seq(1,22,2), seq(2,22,2)]= NaN
corr_psd_DZ[seq(2,22,2), seq(1,22,2)]= NaN
diag(corr_psd_DZ)= NaN
other_corrs_DZ_12= c(corr_psd_DZ)

corr_psd_DZ=cor(t(psd_2[temp,]),t(psd_3[temp,])) 
corss_corrs_DZ_23= c(diag(corr_psd_DZ[seq(1,22,2), seq(2,22,2)]), diag(corr_psd_DZ[seq(2,22,2), seq(1,22,2)]))
corr_psd_DZ[seq(1,22,2), seq(2,22,2)]= NaN
corr_psd_DZ[seq(2,22,2), seq(1,22,2)]= NaN
diag(corr_psd_DZ)= NaN
other_corrs_DZ_23= c(corr_psd_DZ)


corr_psd_DZ=cor(t(psd_1[temp,]),t(psd_3[temp,])) 
corss_corrs_DZ_13= c(diag(corr_psd_DZ[seq(1,22,2), seq(2,22,2)]), diag(corr_psd_DZ[seq(2,22,2), seq(1,22,2)]))
corr_psd_DZ[seq(1,22,2), seq(2,22,2)]= NaN
corr_psd_DZ[seq(2,22,2), seq(1,22,2)]= NaN
diag(corr_psd_DZ)= NaN
other_corrs_DZ_13= c(corr_psd_DZ)


data4plot= data.frame(values= c(corss_corrs_MZ_12, corss_corrs_MZ_23, corss_corrs_MZ_13, 
                                corss_corrs_DZ_12, corss_corrs_DZ_23, corss_corrs_DZ_13),
                      group= c(rep('MZ', 17*3), rep('DZ', 11*3)))

ggplot(data4plot, aes(DescTools::FisherZ(values), color= group, fill=group, alpha=0.5)) +
  geom_density(aes( color= group, fill=group), alpha= 0.5) +
  #geom_histogram(aes(y=..density.., color= group, fill=group), alpha= 0.5, bins = 20) +
  theme_minimal() +
  labs(x = "twin pair correlation", y = "denstity")  + theme_classic2() +
  ggtitle("") + scale_fill_manual(values=cbbPalette[c(1,2,3)]) + 
  scale_color_manual(values=cbbPalette[c(1,2,3)]) 

ggsave('~/Documents/HCP_twin_projec/figures/Review_hist_twinpair_corr_broadband.pdf', device = "pdf")


data4plot2= data.frame(values= c(other_corrs_MZ_12, other_corrs_MZ_23, other_corrs_MZ_13, 
                                 other_corrs_DZ_12, other_corrs_DZ_23, other_corrs_DZ_13),
                       group= c(rep('MZ other', 34*34),rep('MZ other', 34*34),rep('MZ other', 34*34),
                                rep('DZ other', 22*22),rep('DZ other', 22*22),rep('DZ other', 22*22)))

ggplot(data4plot2, aes(DescTools::FisherZ(values), color= group, fill=group, alpha=0.5)) +
  geom_density(aes( color= group, fill=group), alpha= 0.5) +
  #geom_histogram(aes(y=..density.., color= group, fill=group), alpha= 0.5, bins = 20) +
  theme_minimal() +
  labs(x = "twin pair correlation", y = "denstity")  + theme_classic2() +
  ggtitle("") + scale_fill_manual(values=cbbPalette[c(1,2,3)]) + 
  scale_color_manual(values=cbbPalette[c(1,2,3)]) 

ggsave('~/Documents/HCP_twin_projec/figures/Review_hist_other_corr_broadband.pdf', device = "pdf")


data4plot=rbind(data4plot2, data4plot)

ggplot(data4plot, aes(DescTools::FisherZ(values), color= group, fill=group, alpha=0.5)) +
  geom_density(aes( color= group, fill=group), alpha= 0.5) +
  #geom_histogram(aes(y=..density.., color= group, fill=group), alpha= 0.5, bins = 20) +
  theme_minimal() +
  labs(x = "twin pair correlation", y = "denstity")  + theme_classic2() +
  ggtitle("") + scale_fill_manual(values=cbbPalette[c(1,2,3,4)]) + 
  scale_color_manual(values=cbbPalette[c(1,2,3,4)]) 

ggsave('~/Documents/HCP_twin_projec/figures/Review_hist_other_corr_broadband_all_data.pdf', device = "pdf")



temp=demo_hcp$Group == "NotTwin"
corr_psd_NT=cor(t(psd_1[temp,]),t(psd_2[temp,])) 
diag(corr_psd_NT)= NaN
other_corrs_NT_12= c(corr_psd_NT)

temp=demo_hcp$Group == "NotTwin"
corr_psd_NT=cor(t(psd_1[temp,]),t(psd_3[temp,])) 
diag(corr_psd_NT)= NaN
other_corrs_NT_13= c(corr_psd_NT)

temp=demo_hcp$Group == "NotTwin"
corr_psd_NT=cor(t(psd_2[temp,]),t(psd_3[temp,])) 
diag(corr_psd_NT)= NaN
other_corrs_NT_23= c(corr_psd_NT)


data4plot3= data.frame(values= c(other_corrs_NT_12, other_corrs_NT_23, other_corrs_NT_13),
                       group= c(rep('NonTwin', 25*25)))

ggplot(data4plot3, aes(DescTools::FisherZ(values), color= group, fill=group, alpha=0.5)) +
  geom_density(aes( color= group, fill=group), alpha= 0.5) +
  #geom_histogram(aes(y=..density.., color= group, fill=group), alpha= 0.5, bins = 20) +
  theme_minimal() +
  labs(x = "twin pair correlation", y = "denstity")  + theme_classic2() +
  ggtitle("") + scale_fill_manual(values=cbbPalette[c(1,2,3)]) + 
  scale_color_manual(values=cbbPalette[c(1,2,3)]) 

ggsave('~/Documents/HCP_twin_projec/figures/Review_hist_other_corr_broadband_nonTwin.pdf', device = "pdf")

data4plot3=rbind(data4plot2, data4plot3)

ggplot(data4plot3, aes(DescTools::FisherZ(values), color= group, fill=group, alpha=0.5)) +
  geom_density(aes( color= group, fill=group), alpha= 0.5) +
  #geom_histogram(aes(y=..density.., color= group, fill=group), alpha= 0.5, bins = 20) +
  theme_minimal() +
  labs(x = "twin pair correlation", y = "denstity")  + theme_classic2() +
  ggtitle("") + scale_fill_manual(values=cbbPalette[c(1,2,3,4)]) + 
  scale_color_manual(values=cbbPalette[c(1,2,3,4)]) 

ggsave('~/Documents/HCP_twin_projec/figures/Review_hist_other_corr_broadband_all_data.pdf', device = "pdf")



