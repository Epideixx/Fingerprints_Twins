
%% load
clear all
clc

HARgenes = table2cell(readtable('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/HAR_gene_list.csv')); % load specific cell type expression

%% load
cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('gene_expression_abagen.mat')           % relevant node indices
load('gene_names.mat')          
load('result_bandlimited.mat')
load('coordinates.mat')          % (x,y,z) coordinates for brain regions

label=genenames.genenames; % genes in AHBA

expressiongenes= expressiongenes(rowids,:); % fix row ids

genenames = cellstr(HARgenes(:,1));

%% get index of genes with specific cell type expression

for k = 1:length(genenames)                                           % for each gene with specific cell type expression 
    if ismember(genenames(k),label)                                   % if gene overlaps with AHBA genes
        HARgenes(k,4) = num2cell(find(strcmp(label,genenames(k))));  % add index of gene
    else 
        HARgenes(k,4) = {0};                                         % otherwise, add 0
    end
end
bad_idx = cell2mat(HARgenes(:,4))==0; % find indices of genes not in AHBA
HARgenes(bad_idx,:) = [];             % remove these genes


%% compute ratio of HAR genes in pos and neg loadings 

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
ngenes = length(expressiongenes);
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(expressiongenes(:,k),result.vsc(:,1));
end

%gload2 = corr(expressiongenes,result.vsc(:,1)); Does the same thing 

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading

ctd_ratios(1,1) = length(intersect(ipos(gpos_idx),cell2mat(HARgenes(:,4))))/length(gpos_idx);
ctd_ratios(1,2) = length(intersect(ineg(gneg_idx),cell2mat(HARgenes(:,4))))/length(gneg_idx);
ctd_ratios(1,3) = length(intersect([ipos(gpos_idx); ineg(gneg_idx)],cell2mat(HARgenes(:,4))))/length([ipos(gpos_idx); ineg(gneg_idx)]);

%% get null dist of genes
n = 1000;
ctd_null = zeros(2,n);
for k = 1:n % for each repetition
    
    % positive nulls
    y = datasample([1:length(gload)],length(gpos_idx),'Replace',false); % get random gene set the size of the positive gene set
    ctd_null(1,k) = length(intersect(y,cell2mat(HARgenes(:,4)) ))/length(gpos_idx); % find ratio of genes expressed 
    
    % negative nulls
    y = datasample([1:length(gload)],length(gneg_idx),'Replace',false);                      % get random gene set the size of the negative gene set
    ctd_null(2,k) = length(intersect(y, cell2mat(HARgenes(:,4)) ))/length(gneg_idx); % find ratio of genes expressed 

    
    % full gene nulls (top 50 neg and pos) 
    y = datasample([1:length(gload)],length([ipos(gpos_idx); ineg(gneg_idx)]),'Replace',false);                      % get random gene set the size of the negative gene set
    ctd_null(3,k) = length(intersect(y, cell2mat(HARgenes(:,4)) ))/length([ipos(gpos_idx); ineg(gneg_idx)]); % find ratio of genes expressed 
    
end

% get p-values with two-tailed significance test
p_ctd = ctd_ratios - mean(ctd_null,2); % mean centre
p_null = ctd_null - mean(ctd_null,2);
pval = zeros(1,2);
pval(1) = (1+(nnz(find(abs(p_null(1,:)) >= abs(p_ctd(1))))))/(n+1); % pval for positive gene set
pval(2) = (1+(nnz(find(abs(p_null(2,:)) >= abs(p_ctd(2))))))/(n+1); % pval for negative gene set
pval(3) = (1+(nnz(find(abs(p_null(3,:)) >= abs(p_ctd(3))))))/(n+1); % pval for negative gene set


pvalFDR2= mafdr(pval(:),'BHFDR',true);

% both negative and poistive loadings are enriched for Human acelerated
% genes HARs 


figure
scatter(1:3,ctd_ratios(:,2),200,'filled')
hold on
boxplot(squeeze(ctd_null)')
ylim([0.05 0.15])
xtickangle(90)


%% get coritical expression of HAR genes 

% find empirical cell type ratio
HARmap= mean(expressiongenes(:,cell2mat(HARgenes(:,4))),2);


%% let us look at the subset only associated to the brain 

clear all
clc

HARgenes = table2cell(readtable('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/HAR_gene_list.csv')); % load specific cell type expression

% load
cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('gene_expression_abagen.mat')           % relevant node indices
load('gene_names.mat')          
load('result_bandlimited.mat')
load('coordinates.mat')          % (x,y,z) coordinates for brain regions

label=genenames.genenames; % genes in AHBA

expressiongenes= expressiongenes(rowids,:); % fix row ids

HARgenes=HARgenes(cell2mat(HARgenes(:,3))==1,:);

genenames = cellstr(HARgenes(:,1));

%% get index of genes with specific cell type expression

for k = 1:length(genenames)                                           % for each gene with specific cell type expression 
    if ismember(genenames(k),label)                                   % if gene overlaps with AHBA genes
        HARgenes(k,4) = num2cell(find(strcmp(label,genenames(k))));  % add index of gene
    else 
        HARgenes(k,4) = {0};                                         % otherwise, add 0
    end
end
bad_idx = cell2mat(HARgenes(:,4))==0; % find indices of genes not in AHBA
HARgenes(bad_idx,:) = [];             % remove these genes



%% compute ratio of HAR genes in pos and neg loadings 

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
ngenes = length(expressiongenes);
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(expressiongenes(:,k),result.vsc(:,1));
end

%gload2 = corr(expressiongenes,result.vsc(:,1)); Does the same thing 

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading

ctd_ratios(1,1) = length(intersect(ipos(gpos_idx),cell2mat(HARgenes(:,4))))/length(gpos_idx);
ctd_ratios(1,2) = length(intersect(ineg(gneg_idx),cell2mat(HARgenes(:,4))))/length(gneg_idx);
ctd_ratios(1,3) = length(intersect([ipos(gpos_idx); ineg(gneg_idx)],cell2mat(HARgenes(:,4))))/length([ipos(gpos_idx); ineg(gneg_idx)]);


%% get null dist of genes
n = 1000;
ctd_null = zeros(2,n);
for k = 1:n % for each repetition
    
    % positive nulls
    y = datasample([1:length(gload)],length(gpos_idx),'Replace',false); % get random gene set the size of the positive gene set
    ctd_null(1,k) = length(intersect(y,cell2mat(HARgenes(:,4)) ))/length(gpos_idx); % find ratio of genes expressed 
    
    % negative nulls
    y = datasample([1:length(gload)],length(gneg_idx),'Replace',false);                      % get random gene set the size of the negative gene set
    ctd_null(2,k) = length(intersect(y, cell2mat(HARgenes(:,4)) ))/length(gneg_idx); % find ratio of genes expressed 

    
    % full gene nulls (top 50 neg and pos) 
    y = datasample([1:length(gload)],length([ipos(gpos_idx); ineg(gneg_idx)]),'Replace',false);                      % get random gene set the size of the negative gene set
    ctd_null(3,k) = length(intersect(y, cell2mat(HARgenes(:,4)) ))/length([ipos(gpos_idx); ineg(gneg_idx)]); % find ratio of genes expressed 
    
end

% get p-values with two-tailed significance test
p_ctd = ctd_ratios - mean(ctd_null,2); % mean centre
p_null = ctd_null - mean(ctd_null,2);
pval = zeros(1,2);
pval(1) = (1+(nnz(find(abs(p_null(1,:)) >= abs(p_ctd(1))))))/(n+1); % pval for positive gene set
pval(2) = (1+(nnz(find(abs(p_null(2,:)) >= abs(p_ctd(2))))))/(n+1); % pval for negative gene set
pval(3) = (1+(nnz(find(abs(p_null(3,:)) >= abs(p_ctd(3))))))/(n+1); % pval for negative gene set


pvalFDR2= mafdr(pval(:),'BHFDR',true);
% does not remain true when retsirected to HAR brain genes (BUT VERY VERY
% CLOSE) 


figure
scatter(1:3,ctd_ratios(:,2),200,'filled')
hold on
boxplot(squeeze(ctd_null)')
ylim([0.0 0.07])
xtickangle(90)


% find empirical cell type ratio
HARmap= mean(expressiongenes(:,cell2mat(HARgenes(:,4))),2);

