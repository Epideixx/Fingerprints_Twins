
% PLS code (pls_analysis.m) can be downloaded at
% http://pls.rotman-baycrest.on.ca/source/ ("Latest PLS Applications")

%% load
clear all
clc
cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('neurosynth_terms.mat')      % node by term probability matrix
load('gene_expression_abagen.mat')           % relevant node indices
load('SPINStwirls.mat')           % spin test indices
load('coordinates.mat')          % (x,y,z) coordinates for brain regions
load('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/hansen_genescognition-master/terms.mat')
load('terms.mat') 

terms= terms.names;
neurosynth= atlschaefer2018res200neurosynth(rowids,:);
expressiongenes= expressiongenes(rowids,:); % fix row ids
spins= permutedindexesofschaeferatlasSPINsTwirl+1;


%% PLS analysis

addpath(genpath('./Pls/'));

% set up PLS analysis

X = zscore(expressiongenes);
Y = zscore(neurosynth);

nnodes = 200; % number of nodes/ ROIs 
ngenes = length(expressiongenes);
nterms= 123;

% behav pls
option.method = 3;
option.num_boot = 1000;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata = Y;

exp{1} = X;

result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
save('./result_neurosynth.mat','result', '-v7.3')

load('./result_neurosynth.mat')


%% spin test
% this code comes from pls_analysis.m and is modified to account for a
% spatial autocorrelation-preserving permutation test

nspins = 1000;                 % number of permutations ("spins")
s_spins = zeros(nterms,nspins); % singular values
option.method = 3;              % set up PLS
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X;
for k = 1:nspins    
    option.stacked_behavdata = Y(spins(:,k),:);  % permute neurosynth matrix
    
    datamatsvd=rri_xcor(option.stacked_behavdata,exp{1},0); % refer to pls_analysis.m
    [r,c] = size(datamatsvd);
    if r <= c
        [pu, sperm, pv] = svd(datamatsvd',0);
    else
        [pv, sperm, pu] = svd(datamatsvd,0);
    end
    
    %  rotate pv to align with the original v
    rotatemat = rri_bootprocrust(result.v,pv);
 
    %  rescale the vectors
    pv = pv * sperm * rotatemat;

    sperm = sqrt(sum(pv.^2));
    
    s_spins(:,k) = sperm;
end

sprob = zeros(nterms,1); % p-value for each latent variable

for k = 1:nterms % get permuted (via spin test) p-values
    sprob(k) = (1+(nnz(find(s_spins(k,:)>=result.s(k)))))/(1+nspins);
end  


%% plot variance explained

cb=(result.s.^2)/(sum(result.s.^2)); % calculate percent varaince explained 

figure
plot(1:10,cb(1:10)*100,'b.','MarkerSize', 30); % plot percent var explained
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
plot(1:10,cb(1:10)*100,'b-','LineWidth', 1.5); % plot percent var explained
xlabel("Component Number")
ylabel("Percent Covariance Explained (%)")

% first significant component explaines 67% of variance 
% corr b/w usc and vsc 1dst comp is .66


%% bootstrap the data to set CI for the var explained 

nnodes = 160; % number of nodes/ ROIs 
ngenes = length(expressiongenes);
nterms= 6;


var_explained_CI=[];
for b=1:1000
    
    bootind= randi(200,[1,160]);
    % set up PLS analysis

    X = zscore(expressiongenes(bootind,:));
    Y = zscore(neurosynth(bootind,:));
    % behav pls
    option.method = 3;
    option.num_boot = 0;
    option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
    option.stacked_behavdata = Y;

    exp{1} = X;

    result_CI = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
    cb_CI=(result_CI.s.^2)/(sum(result_CI.s.^2)); % calculate percent varaince explained 
    var_explained_CI(b)=cb_CI(1);

end

quantile(var_explained_CI, [0.025, 0.975])

%% distribute PLS-derived gene and term scores

% intrinsic (resting-state) networks from Yeo et al., 2011
rsn = Schaefer2018200Parcels7NetworksNeuromaps.Yeo

figure;               % distribute scores (as boxplots) in 7 networks
subplot(2,1,1)
boxplot(result.usc(:,1),rsn) % gene scores
title('gene scores')
subplot(2,1,2)
boxplot(result.vsc(:,1),rsn) % term scores
title('neurophy scores')

Schaefer2018200Parcels7NetworksNeuromaps.gene_score= result.usc(:,1); 
Schaefer2018200Parcels7NetworksNeuromaps.NeurPhy_score= result.vsc(:,1); 

writetable(Schaefer2018200Parcels7NetworksNeuromaps, './BEHAVIOUR_outputs_of_PLS_scores.csv');

%% get terms

[B,I] = sort(result.boot_result.orig_corr(:,1),'descend');  % sort loadings
npos = length(find(result.boot_result.orig_corr(:,1) > 0)); % number of positive loadings
nneg = length(find(result.boot_result.orig_corr(:,1) < 0)); % number of negative loadings
tpos = I(1:floor(0.25*npos));                                % get top 25% of positive loadings
tneg = I(end-floor(0.25*nneg):end);                          % get top 25% of negative loadings
pos_terms = terms(tpos);                                         % these are the positive terms contributing most
neg_terms = terms(tneg);                                         % these are the negative terms contributing most


varNames = {'term', 'loading'};
NEGATIVE_tab=table(terms(tneg)', abs(result.boot_result.orig_corr(tneg,1)), 'VariableNames',varNames);

colorMap = [linspace(0.5,1,30)', zeros(30,2)];

m1 = 40*0.5;
r = (0:m1-1)'/max(m1-1,1);
g = r;
r = [r; ones(m1,1)];
g = [g; flipud(g)];
b = flipud(r);
c = [r g b]; 

figure
wordcloud(NEGATIVE_tab,'term','loading', 'Color', c(length(tneg):-1:1,:));
title("behavioural terms (negative)")

saveas(gcf,'./wordclous_terms_negative_GO.png')
saveas(gcf,'./wordclous_terms_negative_GO.pdf')
saveas(gcf,'./wordclous_terms_negative_GO.fig')

varNames = {'term', 'loading'};
POSITIVE_tab=table(terms(tpos)', abs(result.boot_result.orig_corr(tpos,1)), 'VariableNames',varNames);

colorMap = [linspace(0.5,1,30)', zeros(30,2)];

m1 = floor(40*0.5);
r = (0:m1-1)'/max(m1,1);
g = r;
r = [r; ones(m1+1,1)];
g = [g; 1; flipud(g)];
b = flipud(r);
c = [r g b]; 

figure
wordcloud(POSITIVE_tab,'term','loading', 'Color', c(40:-1:40-length(tpos)+1, :));
title("behavioural terms (positive)")

saveas(gcf,'./wordclous_terms_positive_GO.png')
saveas(gcf,'./wordclous_terms_positive_GO.pdf')
saveas(gcf,'./wordclous_terms_positive_GO.fig')
