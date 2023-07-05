
% PLS code (pls_analysis.m) can be downloaded at
% http://pls.rotman-baycrest.on.ca/source/ ("Latest PLS Applications")

%% load
clear all
clc
cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('ICC_schaefer_fullcohort.mat')      % node by ICC
load('gene_expression_abagen.mat')           % relevant node indices
load('SPINStwirls.mat')           % spin test indices
load('coordinates.mat')          % (x,y,z) coordinates for brain regions

expressiongenes= expressiongenes(rowids,:); % fix row ids
spins= permutedindexesofschaeferatlasSPINsTwirl+1;

ICCSchaeferband= [mean(ICCScahefer(:,1:8), 2), mean(ICCScahefer(:,9:16), 2),...
  mean(ICCScahefer(:,17:26), 2), mean(ICCScahefer(:,27:60), 2), mean(ICCScahefer(:,61:100), 2),...
  mean(ICCScahefer(:,101:301), 2)];

%% PLS analysis

addpath(genpath('./Pls/'));

% set up PLS analysis

X = zscore(expressiongenes);
Y = zscore(ICCSchaeferband);

nnodes = 200; % number of nodes/ ROIs 
ngenes = length(expressiongenes);
nterms= 6;

% behav pls
option.method = 3;
option.num_boot = 1000;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata = Y;

exp{1} = X;

result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
save('./result_bandlimited.mat','result', '-v7.3')


load('./result_bandlimited.mat')


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
cb_spin=(s_spins.^2)./repmat((sum(s_spins.^2, 1)),[6,1]); % calculate percent varaince explained 


figure
hold on
boxplot(cb_spin'*100)
plot(1:6,cb*100,'b.','MarkerSize', 30); % plot percent var explained
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
plot(1:6,cb*100,'b-','LineWidth', 1.5); % plot percent var explained
xlabel("Component Number")
ylabel("Percent Covariance Explained (%)")
ylim([0 100])


saveas(gcf,'./NeuroPhys_percVarEx.png')
saveas(gcf,'./NeuroPhys_percVarEx.pdf')
%saveas(gcf,'./NeuroPhys_percVarEx.fig')

% first significant component explaines 85% of variance 


%% bootstrap the data to set CI for the var explained 

nnodes = 160; % number of nodes/ ROIs 
ngenes = length(expressiongenes);
nterms= 6;


var_explained_CI=[];
for b=1:1000
    
    bootind= randi(200,[1,160]);
    % set up PLS analysis

    X = zscore(expressiongenes(bootind,:));
    Y = zscore(ICCSchaeferband(bootind,:));
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


writetable(Schaefer2018200Parcels7NetworksNeuromaps, './outputs_ofPLS_scores.csv');

%% loadings of bands

l1=corr(ICCSchaeferband,result.usc(:,1)); 

bands= categorical({'1.delta', '2.theta', '3.alpha', '4.beta', '5.gamma', '6.high gamma'});
bands = reordercats(bands,cellstr(bands)');

figure
clear g
g(1,1)=gramm('x',bands,'y',l1, 'color', bands);
g(1,1).stat_summary('geom','bar','setylim',true);
g(1,1).set_title('Neurophysiology loadings ''geom'',''bar''');

g.draw();

saveas(gcf,'./NeuroPhys_loadings.png')
saveas(gcf,'./NeuroPhys_loadings.pdf')
saveas(gcf,'./NeuroPhys_loadings.fig')

