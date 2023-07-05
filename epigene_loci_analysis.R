
# read list of genes
epigen= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/CoRSIVs_epigenome_areas_variation.csv')

# save promoter, body and three prime seperately 
promoter_genes= unique(unlist(strsplit(epigen$promoter....3kb.TSS..overlapping.genes, ',')))

body_genes= unique(unlist(strsplit(epigen$gene.body.overlapping.genes, ',')))

threepimre_genes= unique(unlist(strsplit(epigen$three.prime.region.....3kb.TES..overlapping.genes, ',')))

# get gene loadings for positive, negative, and all genes used in PLS
positive_genes= read.csv('~/Documents/HCP_twin_projec/abagen_analysis/positive_genes.csv')

negative_genes= read.csv('~/Documents/HCP_twin_projec/abagen_analysis/negative_genes.csv')

all_genes= read.csv('~/Documents/HCP_twin_projec/abagen_analysis/gene_names.csv', header = FALSE)



############### over representation of promoter region #########
intersect(positive_genes$genenames, promoter_genes)
# 43 genes in promoter region for positive set
intersect(negative_genes$genenames, promoter_genes)
# 30 genes in promoter region for negative set
intersect(all_genes$V1, promoter_genes)
# 164 genes in promoter region for all genes

# permutations of all gene set positive genes
rand_values=c()
for (i in 1:10000){
  
  rand_genes=all_genes$V1[sample(9104)[1:2272]]
  rand_values[i]=length(intersect(rand_genes, promoter_genes))

}

sum(rand_values < 43)/10000


# permutations of all gene set negative genes
rand_values=c()
for (i in 1:10000){
  
  rand_genes=all_genes$V1[sample(9104)[1:2280]]
  rand_values[i]=length(intersect(rand_genes, promoter_genes))
  
}

sum(rand_values < 30)/10000



# plot hsitogram 
library(ggplot2)
ggplot() + geom_density(aes(x= rand_values, color= 'grey', fill='grey'), alpha= 0.5) +
  ggpubr::theme_classic2() + geom_vline(xintercept = 30) +
  labs(x = "self correlation", y = "denstity")  + 
  ggtitle("") + scale_fill_manual(values= '#D3D3D3') + scale_color_manual(values='#D3D3D3') 

ggsave('~/Documents/HCP_twin_projec/figures/epigenetics_rand_perm_promoter_negative_gene.pdf', device = "pdf")


############### over representation of gene body methyl #########
intersect(positive_genes$genenames, body_genes)
# 368 genes in promoter region for positive set
intersect(negative_genes$genenames, body_genes)
# 357 genes in promoter region for negative set
intersect(all_genes$V1, body_genes)
# 1449 genes in promoter region for all genes

# permutations of all gene set positive genes
rand_values=c()
for (i in 1:10000){
  
  rand_genes=all_genes$V1[sample(9104)[1:2272]]
  rand_values[i]=length(intersect(rand_genes, body_genes))
  
}

hist(rand_values)
sum(rand_values < 368)/10000

# permutations of all gene set negative genes
rand_values=c()
for (i in 1:10000){
  
  rand_genes=all_genes$V1[sample(9104)[1:2280]]
  rand_values[i]=length(intersect(rand_genes, body_genes))
  
}

sum(rand_values < 357)/10000




############### over representation of three prime region #########
intersect(positive_genes$genenames, threepimre_genes)
# 38 genes in promoter region for positive set
intersect(negative_genes$genenames, threepimre_genes)
# 45 genes in promoter region for negative set
intersect(all_genes$V1, threepimre_genes)
# 194 genes in promoter region for all genes

# permutations
rand_values=c()
for (i in 1:10000){
  
  rand_genes=all_genes$V1[sample(9104)[1:2272]]
  rand_values[i]=length(intersect(rand_genes, threepimre_genes))
  
}

hist(rand_values)
sum(rand_values < 38)/10000


# plot hsitogram 
library(ggplot2)
ggplot() + geom_density(aes(x= rand_values, color= 'grey', fill='grey'), alpha= 0.5) +
  ggpubr::theme_classic2() + geom_vline(xintercept = 30) +
  labs(x = "self correlation", y = "denstity")  + 
  ggtitle("") + scale_fill_manual(values= '#D3D3D3') + scale_color_manual(values='#D3D3D3') 

ggsave('~/Documents/HCP_twin_projec/figures/epigenetics_rand_three_promoter_positive_genes.pdf', device = "pdf")


# permutations
rand_values=c()
for (i in 1:10000){
  
  rand_genes=all_genes$V1[sample(9104)[1:2280]]
  rand_values[i]=length(intersect(rand_genes, threepimre_genes))
  
}

sum(rand_values < 45)/10000

# none survive FDR correction 
p.adjust(c(0.6178,0.0155,0.6419,0.3284,0.0312, 0.2462), 'fdr')




