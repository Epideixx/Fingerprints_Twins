
library(dplyr)
library(ggplot2)

negativegenes= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/negative_top_genes.csv')
positivegenes= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/positive_top_genes.csv')

negativegenes$Pathway = factor(negativegenes$Pathway, levels = negativegenes$Pathway[order(positivegenes$FoldEnrichment)])


ggplot(negativegenes, aes(x= log10(FoldEnrichment), y=Pathway, colour = Loadings, fill = Loadings, size = abs(log10(EnrichmentFDR)))) +
  geom_point() + ggpubr::theme_classic2() + scale_fill_gradient(low= '#40E0D0', high = '#2A796D') +
  scale_color_gradient(low= '#40E0D0', high = '#2A796D')

ggsave('~/Documents/HCP_twin_projec/figures/negativegenesTop.pdf', device = "pdf", width= 20, height=9)


positivegenes$Pathway = factor(positivegenes$Pathway, levels = positivegenes$Pathway[order(positivegenes$FoldEnrichment)])

ggplot(positivegenes, aes(x= log10(FoldEnrichment), y=Pathway, colour = Loadings, fill = Loadings, size = abs(log10(EnrichmentFDR)))) +
  geom_point() + ggpubr::theme_classic2() + scale_fill_gradient(low= '#F315AD', high = '#84075D') +
  scale_color_gradient(low= '#F315AD', high = '#84075D')

ggsave('~/Documents/HCP_twin_projec/figures/positivegenesTop.pdf', device = "pdf", width= 20, height=9)



data= read.csv('~/Documents/HCP_twin_projec/Behav_Neurosynth_loadings.csv')
dataplot=data[c(1:10, 114:123), ]
dataplot$Sign= factor(dataplot$Sign)
dataplot$terms = factor(dataplot$terms, levels = dataplot$terms)

ggplot(dataplot, aes(x=terms, y= loadings, colour = Sign, fill = Sign)) + 
  geom_bar(stat = "identity", position=position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.0 ,position=position_dodge(.5), color= 'grey') +
  ggpubr::theme_classic2() + scale_color_manual(values=c('#40E0D0','#F315AD')) + scale_fill_manual(values=c('#40E0D0','#F315AD')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('')

ggsave('~/Documents/HCP_twin_projec/figures/TopNeuroSynthTerms.pdf', device = "pdf", width= 10, height=4)

