% This script organizes and parcellates data from BrainSpan, then computes
% estimated gene scores based on the gene weights computed in the original
% analyses (see scpt_genes_cog_pls.m). Note that the figure generation will
% run into an error if you don't have cbrewer added to your path
% (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)

% BrainSpan data can be downloaded in its original form from
% https://www.brainspan.org/static/download.html

%% load
clear all
close all

cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('gene_expression_abagen.mat')           % relevant node indices
load('gene_names.mat')          
load('result_bandlimited.mat')
load('SPINStwirls.mat')           % spin test indices

label=genenames.genenames; % genes in AHBA

spins= permutedindexesofschaeferatlasSPINsTwirl+1;


expressiongenes= expressiongenes(rowids,:); % fix row ids

cd('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/hansen_genescognition-master');
load('./BrainSpan/mapping.mat')  % mapping from 34 node parcellation to the 16 unique cortical regions included in BrainSpan

% AHBA harmonized files have been organized to include only genes included
% in AHBA
brainspan = readtable('./BrainSpan/gene_expression_AHBA_harmonized.csv'); % gene by sample matrix of gene expression
gene_info = readtable('./BrainSpan/gene_metadata_AHBA_harmonized.csv');   % metadata on genes
sample_info = readtable('./BrainSpan/samples_metadata.csv');              % metadata on tissue samples
brainspan = table2array(brainspan(2:end,2:end));              % remove column of gene names and row of sample IDs

%% remove non-cortical samples
% remove samples labeled 'not cortex' and amygdala samples

notcortex_idx = [];
sensorifugal = table2cell(sample_info(:,15)); % this column labels noncortical regions as 'Not_Cortex' 
sensorifugal = string(sensorifugal);

% find noncortical defined by sensorifugal
for k = 1:length(sensorifugal)
    if strcmp(sensorifugal(k),'Not_Cortex')
        notcortex_idx = [notcortex_idx; k];
    end
end

% find amygdala indices
amygdala = find(contains(table2cell(sample_info(:,8)),'amygdaloid complex'));
notcortex_idx = [notcortex_idx; amygdala];
notcortex_idx = sort(notcortex_idx);

% get all relevant (cortical) indices
cortex_idx = setdiff([1:length(sensorifugal)],notcortex_idx); 

% remove noncortical indices
brainspan(:,notcortex_idx) = [];


%% find pos neg loadings

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


%% get index of pos genes

% genes included in original analyses are stable (differential stability >
% 0.1). For comparability, we only include (available) stable genes as defined on the
% our schaefer parcellation

pos_label = label(ipos(gpos_idx)); % names of positive loaded genes
gene_label = table2cell(gene_info(:,4));     % names of genes in BrainSpan

notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)                            % for each gene in BrainSpan
    if ~ismember(gene_label(k), pos_label)            % if gene isn't stable
        notstable_idx = [notstable_idx; k];             % keep its index
    else                                                % if gene is stable
        i = find(ismember(pos_label,gene_label(k))); % get index of stable gene
        bspan_gidx033 = [bspan_gidx033; i];             % store it
    end
end

% remove nonstable genes
pos_brainspan= brainspan;
pos_brainspan(notstable_idx,:) = [];


neg_label = label(ineg(gneg_idx)); % names of positive loaded genes
gene_label = table2cell(gene_info(:,4));     % names of genes in BrainSpan

notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)                            % for each gene in BrainSpan
    if ~ismember(gene_label(k), neg_label)            % if gene isn't stable
        notstable_idx = [notstable_idx; k];             % keep its index
    else                                                % if gene is stable
        i = find(ismember(neg_label,gene_label(k))); % get index of stable gene
        bspan_gidx033 = [bspan_gidx033; i];             % store it
    end
end

% remove nonstable genes
neg_brainspan= brainspan;
neg_brainspan(notstable_idx,:) = [];


%% organize data by life stage
% samples are organized into 5 life stages: fetal, infant, child,
% adolescent, and adult

% get indices
fetal_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'fetal'));
infant_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'infant'));
child_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'child'));
adolescent_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adolescent'));
adult_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adult'));

% make gene expression matrices for each life stage
fetal = neg_brainspan(:,fetal_idx);
infant = neg_brainspan(:,infant_idx);
child = neg_brainspan(:,child_idx);
adolescent = neg_brainspan(:,adolescent_idx);
adult = neg_brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
M_NEG = {fetal, infant, child, adolescent, adult};

fetal = pos_brainspan(:,fetal_idx);
infant = pos_brainspan(:,infant_idx);
child = pos_brainspan(:,child_idx);
adolescent = pos_brainspan(:,adolescent_idx);
adult = pos_brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
M_POS = {fetal, infant, child, adolescent, adult};

%% parcellate brain regions

sample_regions = table2cell(sample_info(cortex_idx,8)); % get region name of each sample
regions = unique(string(sample_regions));               % get unique regions - this is the maximum number of brain regions in each gene expression matrix

% make value-based region mapping instead of string-based (for
% simplicity)
region_mapping = zeros(length(sample_regions),1);

for k = 1:length(regions)                               % for each unique region
    i = find(contains(sample_regions,regions(k)));      % find all samples from that region
    region_mapping(i) = k;                              % make index-based map of regions
end

% average expression of each gene from identical regions
% done separately for each life stage, where the number of unique regions
% with gene expression estimates varies across life stage
regionsIncluded = cell(5,1);
for k = 1:length(M_NEG)                                          % for each life stage
    mat = M_NEG{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M_NEG{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end

regionsIncluded = cell(5,1);
for k = 1:length(M_POS)                                          % for each life stage
    mat = M_POS{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M_POS{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end


%% track gene scores across development
% see how gene scores at the brain regions with gene expression estimates
% across all 5 life stages change with development

% get regions with available gene expression estimates across all life stages
reg = intersect(regionsIncluded{1},[regionsIncluded{2};regionsIncluded{3};regionsIncluded{4};regionsIncluded{5}]);
M_POS{1} = M_POS{1}(reg,:); % only need to change first life stage (16 regions)
M_NEG{1} = M_NEG{1}(reg,:); % only need to change first life stage (16 regions)


cellfun(@(x)mean(x, 2, 'omitnan'), M_NEG, 'UniformOutput', false)
temp=cellfun(@(x)mean(x, 2, 'omitnan'), M_POS, 'UniformOutput', false);

POS = horzcat(temp{:});

temp=cellfun(@(x)mean(x, 2, 'omitnan'), M_NEG, 'UniformOutput', false);

NEG = horzcat(temp{:});

%% repeat for random sample of genes 

% FIRST remove unstable genes

% genes included in original analyses are stable (differential stability >
% 0.1). For comparability, we only include (available) stable genes as defined on the
% our schaefer parcellation

stable_label = label; % names of stable genes from AHBA
gene_label = table2cell(gene_info(:,4));     % names of genes in BrainSpan

notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)                            % for each gene in BrainSpan
    if ~ismember(gene_label(k), stable_label)            % if gene isn't stable
        notstable_idx = [notstable_idx; k];             % keep its index
    else                                                % if gene is stable
        i = find(ismember(stable_label,gene_label(k))); % get index of stable gene
        bspan_gidx033 = [bspan_gidx033; i];             % store it
    end
end

% remove nonstable genes
brainspan(notstable_idx,:) = [];


nperm=1000;
ran_M= zeros(12,5,nperm);
for p= 1:nperm

    ranindex = randi([1 length(brainspan)],1, length(neg_brainspan));

    ran_brainspan=brainspan;
    ran_brainspan(ranindex,:) = [];

    fetal = ran_brainspan(:,fetal_idx);
    infant = ran_brainspan(:,infant_idx);
    child = ran_brainspan(:,child_idx);
    adolescent = ran_brainspan(:,adolescent_idx);
    adult = ran_brainspan(:,adult_idx);

    % organize matrices and indices
    M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
    M_ran = {fetal, infant, child, adolescent, adult};

    % average expression of each gene from identical regions
    % done separately for each life stage, where the number of unique regions
    % with gene expression estimates varies across life stage
    regionsIncluded = cell(5,1);
    for k = 1:length(M_ran)                                          % for each life stage
        mat = M_ran{k};                                              % get original gene expression matrix (genes x samples)
        r = region_mapping(M_idx{k});                            % get region mapping of this life stage
        regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
        mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
        i = 1;
        for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
            aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
            mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
            i = i+1;
        end
        M_ran{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
    end


    M_ran{1} = M_ran{1}(reg,:); % only need to change first life stage (16 regions)
    temp=cellfun(@(x)mean(x, 2, 'omitnan'), M_ran, 'UniformOutput', false);

    ran_M(:,:,p) = horzcat(temp{:});

end


%% make plots of gene expression over age groups

newcolors = [255, 52, 153
             255, 100, 153
             255, 52, 200
             255, 52, 100];
         

colororder(newcolors/255)
y = mean(mean(ran_M, 3)',2);
x = 1:numel(y);
std_dev = mean(std(ran_M, 0,3)',2);
curve1 = y + std_dev/ sqrt(1000);
curve2 = y - std_dev/ sqrt(1000);
x2 = [x, fliplr(x)];
inBetween = [curve1; fliplr(curve2')'];
fill(x2, inBetween, [0.8, 0.8, 0.8], 'LineWidth', 0.0001, 'FaceAlpha',0.3);
hold on;
plot(x, y, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);

plot(POS', '-.')
hold on;
plot(mean(POS)', 'LineWidth', 2, 'Color', [255, 52, 153]/255);
xticks(1:5)
xticklabels({'fetus','infant','child','adolescent','adult'})
ylabel('estimated gene expression')

newcolors = [0, 255, 255
             0, 220, 255
             0, 255, 220
             0, 230, 230];
         

colororder(newcolors/255)
y = mean(mean(ran_M, 3)',2);
x = 1:numel(y);
std_dev = mean(std(ran_M, 0,3)',2);
curve1 = y + std_dev/ sqrt(1000);
curve2 = y - std_dev/ sqrt(1000);
x2 = [x, fliplr(x)];
inBetween = [curve1; fliplr(curve2')'];
fill(x2, inBetween, [0.8, 0.8, 0.8], 'LineWidth', 0.0001, 'FaceAlpha',0.3);
hold on;
plot(x, y, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);

plot(NEG', '-.')
hold on;
plot(mean(NEG)', 'LineWidth', 2, 'Color', [0, 200, 200]/255);
xticks(1:5)
xticklabels({'fetus','infant','child','adolescent','adult'})
ylabel('estimated gene expression')



%% organize data by life stage
% samples are organized into 5 life stages: fetal, infant, child,
% adolescent, and adult

% get indices
fetal_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'fetal'));
infant_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'infant'));
child_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'child'));
adolescent_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adolescent'));
adult_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adult'));

% make gene expression matrices for each life stage
fetal = brainspan(:,fetal_idx);
infant = brainspan(:,infant_idx);
child = brainspan(:,child_idx);
adolescent = brainspan(:,adolescent_idx);
adult = brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
M = {fetal, infant, child, adolescent, adult};

%% parcellate brain regions

sample_regions = table2cell(sample_info(cortex_idx,8)); % get region name of each sample
regions = unique(string(sample_regions));               % get unique regions - this is the maximum number of brain regions in each gene expression matrix

% make value-based region mapping instead of string-based (for
% simplicity)
region_mapping = zeros(length(sample_regions),1);

for k = 1:length(regions)                               % for each unique region
    i = find(contains(sample_regions,regions(k)));      % find all samples from that region
    region_mapping(i) = k;                              % make index-based map of regions
end

% average expression of each gene from identical regions
% done separately for each life stage, where the number of unique regions
% with gene expression estimates varies across life stage
regionsIncluded = cell(5,1);
for k = 1:length(M)                                          % for each life stage
    mat = M{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end

%% get gene and term scores

% estimate gene scores by multiplying PLS-derived gene weights with gene
% expression matrices

gscore = cell(length(M),1);

% as per reviewer's request, use only genes with complete data instead of
% imputing with median expression (commented above, method for the preprint)
missing_genes = []; 
for k = 1:length(M) % for each life stage
    m = M{k}; % get the gene expression matrix (genes x brain regions)
    [~,i] = find(isnan(m)); % get 
    missing_genes = union(missing_genes, i);
end

bspan_gidx033(missing_genes) = [];

for k = 1:length(M)
    m = M{k};
    m(:,missing_genes) = [];
    M{k} = m;
    gscore{k} = m * result.u(bspan_gidx033);
end


%% track gene scores across development
% see how gene scores at the brain regions with gene expression estimates
% across all 5 life stages change with development

% get regions with available gene expression estimates across all life stages
reg = intersect(regionsIncluded{1},[regionsIncluded{2};regionsIncluded{3};regionsIncluded{4};regionsIncluded{5}]);
tmp = gscore;
tmp{1} = tmp{1}(reg); % only need to change first life stage (16 regions)
gscore_mat = [tmp{1} tmp{2} tmp{3} tmp{4} tmp{5}]; % organize gene scores

%% visualize

addpath(genpath('/Users/jason/Documents/cbrewer2'));

% change the colourmap if you don't want to download cbrewer
cm=cbrewer2('qual', 'Paired', 16, 'PiYG');

% gene scores across development

figure;
for k = 1:length(gscore_mat)                                    % for each brain region
    hold on
    plot(gscore_mat(k,:),'LineWidth',1.3,'Color',cm(reg(k),:))  % plot a curve of gene score in that brain region across all five life stages
end
legend(regions(reg),'Location','northwest');
ylim([-500, 100])
xticks(1:5)
xticklabels({'fetus','infant','child','adolescent','adult'})
ylabel('estimated gene scores')

%% RNAD PERMUTATIONS

cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('ICC_schaefer_fullcohort.mat')      % node by ICC
addpath(genpath('./Pls/'));


ICCSchaeferband= [mean(ICCScahefer(:,1:8), 2), mean(ICCScahefer(:,9:16), 2),...
  mean(ICCScahefer(:,17:26), 2), mean(ICCScahefer(:,27:60), 2), mean(ICCScahefer(:,61:100), 2),...
  mean(ICCScahefer(:,101:301), 2)];



% set up PLS analysis

X = zscore(expressiongenes);
Y = zscore(ICCSchaeferband);

nnodes = 200; % number of nodes/ ROIs 
ngenes = length(expressiongenes);
nterms= 6;
nperm=1000;
gscoreperm = cell(length(M),1);

gscore_permutations=zeros(12,5,nperm);

for p=1:nperm
% behav pls
option.method = 3;
option.num_boot = 0;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata = Y(spins(:,p),:);

exp{1} = X;
result_perm = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses

for k = 1:length(M)
    m = M{k};
    %m(:,missing_genes) = [];
    %M{k} = m;
    gscoreperm{k} = m * result_perm.u(bspan_gidx033);
end


gscoreperm{1} = gscoreperm{1}(reg); % only need to change first life stage (16 regions)
gscore_permutations(:,:,p)= horzcat(gscoreperm{:});

end


figure
y = mean(mean(gscore_permutations, 3)',2);
x = 1:numel(y);
std_dev = mean(std(gscore_permutations, 0,3)',2);
curve1 = y + std_dev/ sqrt(1000);
curve2 = y - std_dev/sqrt(1000);
x2 = [x, fliplr(x)];
inBetween = [curve1; fliplr(curve2')'];
fill(x2, inBetween, [0.8, 0.8, 0.8], 'LineWidth', 0.0001, 'FaceAlpha',0.3);
hold on;
plot(x, y, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);

for k = 1:length(gscore_mat)                                    % for each brain region
    hold on
    plot(gscore_mat(k,:),'LineWidth',1.3,'Color',cm(reg(k),:))  % plot a curve of gene score in that brain region across all five life stages
end
legend([''; ''; regions(reg)],'Location','northwest');
ylim([-500, 100])
xticks(1:5)
xticklabels({'fetus','infant','child','adolescent','adult'})
ylabel('estimated gene scores')

%% fit a bunch of linear models 

C = squeeze(num2cell(gscore_permutations, 2));
polyfits_permutations=cellfun(@(x)polyfit(1:5, x, 1), C, 'UniformOutput', false);

permuted_slopes= cellfun(@(x)x(1), polyfits_permutations, 'UniformOutput', false);

C = squeeze(num2cell(gscore_mat, 2));
polyfits_gscore=cellfun(@(x)polyfit(1:5, x, 1), C, 'UniformOutput', false);
polyfits_gscore= cellfun(@(x)x(1), polyfits_gscore, 'UniformOutput', false);

pval=sum(cell2mat(polyfits_gscore)< cell2mat(permuted_slopes),2)/1000

pvalFDR2= mafdr(pval,'BHFDR',true);

%% PLOT HISTOGRAMS OF PERMUTED SLOPES PER REGION

permuted_slopes=cell2mat(permuted_slopes);

figure
for l=1:12
subplot(12,1,l)
[f,xi] = ksdensity(permuted_slopes(l,:)); 
%f= [f; zeros(1,length(xi))]
%xi = [xi; xi];
fill(xi,f,cm(reg(l),:))
hold on;
xline(polyfits_gscore{l}, 'LineWidth', 2);
ylim([0,0.025])
xlim([-80, 150])
end


