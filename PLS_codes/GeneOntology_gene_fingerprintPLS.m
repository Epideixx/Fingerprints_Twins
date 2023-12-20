% This script finds the biological processes with which two gene sets are
% most involved. For each biological process (category), the mean loading
% of genes is computed (called the category score). Significance is
% assessed against a null model of category scores. This null model is
% constructed using the loadings of the same genes of interest,
% where null loadings come from PLS analysis on a
% spatial-autocorrelation preserving permutation of one matrix.

% see https://github.com/benfulcher/GeneSetEnrichmentAnalysis for more
% details on GO annotations

% cross corr usc with Y and vsc with X

%% load
clear all
clc
cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('ICC_schaefer_fullcohort.mat')      % node by ICC matrix
load('gene_expression_abagen.mat')           % relevant node indices
load('SPINStwirls.mat')           % spin test indices
load('coordinates.mat')          % (x,y,z) coordinates for brain regions
load('gene_names.mat')          


expressiongenes= expressiongenes(rowids,:); % fix row ids align by region
spins= permutedindexesofschaeferatlasSPINsTwirl+1;

ICCSchaeferband= [mean(ICCScahefer(:,1:8), 2), mean(ICCScahefer(:,9:16), 2),...
  mean(ICCScahefer(:,17:26), 2), mean(ICCScahefer(:,27:60), 2), mean(ICCScahefer(:,61:100), 2),...
  mean(ICCScahefer(:,101:301), 2)];

load('./GOAnnotationDirect-human.mat')
load('./GOTerms_BP.mat')

%% get genes with entrezID

T = table2cell(readtable('gene_entrez_ids')); % load entrezID of genes

gene_name = genenames.genenames;     % get relevant gene names
entrezIDs = zeros(size(gene_name));

idx = [];
for k = 1:length(gene_name)                                                % for each gene
    if ismember(gene_name(k), T(:,1))                                      % if the gene has an entrezID
        entrezIDs(k) = cell2mat(T(find(strcmp(gene_name(k), T(:,1))),2));  % store the entrezID
        idx = [idx;k];                                                     % also store the index of the gene
    end
end
%entrezIDs = entrezIDs(entrezIDs ~= 0);                                     % remove all genes without entrezID
entrezIDsNONID = entrezIDs(entrezIDs ~= 0); % this will be our background genes to compare to in the enrichment analysis 

%% get category scores

% set up PLS analysis
addpath(genpath('./Pls/'));

% set up PLS analysis

X = zscore(expressiongenes);
Y = zscore(ICCSchaeferband);

nnodes = 200;
nterms = 6;
ngenes = length(expressiongenes);

% behav pls
option.method = 3;
option.num_boot = 0;
option.num_perm = 0;    
option.stacked_behavdata = Y;

exp{1} = X(:,:);

result = pls_analysis(exp, nnodes, 1, option); % empirical category scores come from here

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(expressiongenes(:,k),result.vsc(:,1));
end

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading

gpos_ID = entrezIDs(ipos(gpos_idx));   % these are the entrezIDs of the genes in the positive set
gneg_ID = entrezIDs(ineg(gneg_idx));  % these are the entrezIDs of the genes in the negative set

gpos_ID_nonzero= gpos_ID(gpos_ID~=0);
gneg_ID_nonzero= gneg_ID(gneg_ID~=0);

% find genes with no ids
% compute category score for them too 
ids_pos=ipos(gpos_idx);
ids_neg=ineg(gneg_idx);



%% get null category scores (old GO analysis) 
% 
% nperm = 1000;
% categoryScores_null = zeros(length(gload),2,nperm);
% 
% for m = 1:nperm % for each permutation
%     
%     % PLS analysis
%     option.stacked_behavdata = zscore(Y(spins(:,m),:)); % rows of neurosynth matrix has been permuted while preserving spatial autocorrelation
%     result = pls_analysis(exp, nnodes, 1, option);              % this is the null PLS result
%     
% 
%         gloadRAN = zeros(ngenes,1);
%     for k = 1:ngenes
%         gloadRAN(k) = corr(expressiongenes(:,k),result.vsc(:,1));
%     end
%     
%     categoryScores_null(:,1,m) = gloadRAN();
% 
% end
% 
% %% get biological processes
% 
% 
% % set up p-value template
% pvals = zeros(length(gene_name),1);
% 
% % calculate permuted p-value
% for k = 1:length(pvals)
%     
%     if gload(k) >0
%         pvals(k) = length(find(categoryScores_null(k,1,:) > gload(k)))/(nperm);
%     end
%     
%     if gload(k) <0
%         pvals(k) = length(find(categoryScores_null(k,1,:) < gload(k)))/(nperm);
%     end
%     
% end
% 
% 
% %% make a nice table with all info 
% 
% varNames = {'gene', 'EntrezID','p_value','loading'};
% NEGATIVE_tab=table(gene_name(ids_neg), entrezIDs(ids_neg),pvals(ids_neg), gload(ids_neg), 'VariableNames',varNames);
% 
% POSITIVE_tab=table(gene_name(ids_pos),entrezIDs(ids_pos), pvals(ids_pos), gload(ids_pos), 'VariableNames',varNames);
% 
% POSITIVE_tab_cleaned = POSITIVE_tab(POSITIVE_tab.p_value < 0.05 & POSITIVE_tab.EntrezID ~=0, :);
% 
% NEGATIVE_tab_cleaned = NEGATIVE_tab(NEGATIVE_tab.p_value < 0.05 & NEGATIVE_tab.EntrezID ~=0, :);

% save('./matlab_GO_workspace_new.mat')



%%
enrichmentallGOpositive=importfile('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/SupplementalData/GO_biological_processes_positive.csv');


% load gene enrichment files too 
load('./matlab_GO_workspace_new.mat')
enrichmentallGOpositiveClean=enrichmentallGOpositive(enrichmentallGOpositive.EnrichmentFDR <0.05, :);

Loadings=zeros(length(enrichmentallGOpositiveClean.Genes),1);
for i =1:length(enrichmentallGOpositiveClean.Genes)
    
    genes_temp=split(enrichmentallGOpositiveClean.Genes(i));
    temp_gload= nan(length(genes_temp),1);
    for k =1:length(genes_temp)

        if ~isempty(find(strcmp(genes_temp(k), gene_name)))
        
            temp_gload(k,1)=gload(find(strcmp(genes_temp(k), gene_name)));
        end
    end
   Loadings(i)= mean(temp_gload, 'omitnan' );
end

enrichmentallGOpositiveClean.Loadings=Loadings;

cutoff=sort(enrichmentallGOpositiveClean.EnrichmentFDR, "ascend");
enrichmentallGOpositiveClean2= enrichmentallGOpositiveClean(enrichmentallGOpositiveClean.EnrichmentFDR <=cutoff(35), : );
[B ,I]= sort(enrichmentallGOpositiveClean2.FoldEnrichment, 'descend');
enrichmentallGOpositiveClean2= enrichmentallGOpositiveClean2(I,:);


m1 = floor(70*0.5);
r = (0:m1-1)'/max(m1,1);
g = r;
r = [r; ones(m1+1,1)];
g = [g; 1; flipud(g)];
b = flipud(r);
c = [r g b]; 

enrichmentallGOpositiveClean2.Pathway=lower(enrichmentallGOpositiveClean2.Pathway);


figure
wordcloud(enrichmentallGOpositiveClean2,'Pathway','Loadings'); %, 'Color', c(71:-1:37,:));
title("gene process (positive)")

%saveas(gcf,'./wordclous_positive_GO_DEC2023.png')
%saveas(gcf,'./wordclous_positive_GO_DEC2023.pdf')
%saveas(gcf,'./wordclous_positive_GO_DEC2023.fig')


% negative 

enrichmentallGOnegative=importfile('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/SupplementalData/GO_biological_processes_negative.csv');


enrichmentallGOnegativeClean=enrichmentallGOnegative(enrichmentallGOnegative.EnrichmentFDR <0.05, :);

Loadings=zeros(length(enrichmentallGOnegativeClean.Genes),1);
for i =1:length(enrichmentallGOnegativeClean.Genes)
    
    genes_temp=split(enrichmentallGOnegativeClean.Genes(i));
    temp_gload= nan(length(genes_temp),1);
    for k =1:length(genes_temp)

        if ~isempty(find(strcmp(genes_temp(k), gene_name)))
        
            temp_gload(k,1)=gload(find(strcmp(genes_temp(k), gene_name)));
        end
    end
   Loadings(i)= mean(temp_gload, 'omitnan' );
end

enrichmentallGOnegativeClean.Loadings=Loadings;

cutoff=sort(enrichmentallGOnegativeClean.EnrichmentFDR, "ascend");
enrichmentallGOnegativeClean2= enrichmentallGOnegativeClean(enrichmentallGOnegativeClean.EnrichmentFDR <=cutoff(35), : );
[B ,I]= sort(enrichmentallGOnegativeClean2.FoldEnrichment, 'descend');
enrichmentallGOnegativeClean2= enrichmentallGOnegativeClean2(I,:);
enrichmentallGOnegativeClean2.Loadings = abs(enrichmentallGOnegativeClean2.Loadings);



m1 = 70*0.5;
r = (0:m1-1)'/max(m1-1,1);
g = r;
r = [r; ones(m1,1)];
g = [g; flipud(g)];
b = flipud(r);
c = [r g b]; 

enrichmentallGOnegativeClean2.Pathway=lower(enrichmentallGOnegativeClean2.Pathway);

figure
wordcloud(enrichmentallGOnegativeClean2,'Pathway','FoldEnrichment' );%, 'Color', c(35:-1:1,:));
title("gene process (negative)")

%saveas(gcf,'./wordclous_negative_GO_DEC2023.png')
%saveas(gcf,'./wordclous_negative_GO_DEC2023.pdf')
%saveas(gcf,'./wordclous_negative_GO_DEC2023.fig')


