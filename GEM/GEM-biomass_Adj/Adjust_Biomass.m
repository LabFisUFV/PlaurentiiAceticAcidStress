%%Adjust biomass

%Load model

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE.mat');
model = ecModel;
%save('ecpaplaGEM_stage3_BRENDA_wostd_ALE.mat','model')
% Load biomass information

fid           = fopen('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\biomassCuration_Parwo.csv');
loadedData    = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);
BM.name       = loadedData{1};    BM.mets     = loadedData{2};
BM.pseudorxn  = loadedData{3};    BM.coeff    = loadedData{4};


% Nucleotides (DNA)
% Find out which rows contain the relevant information
indexes = find(contains(BM.pseudorxn, 'DNA'));
% Define new stoichiometries
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
% Change reaction
model = changeRxns(model, 'r_4050', equations, 1);

% Ribonucleotides (RNA)
indexes = find(contains(BM.pseudorxn, 'RNA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4049', equations, 1);

% Amino acids (protein)
indexes = find(contains(BM.pseudorxn, 'AA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4047', equations, 1);

% Carbohydrates
indexes = find(contains(BM.pseudorxn, 'carbohydrate'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4048', equations, 1);

% Lipid backbones
indexes = find(contains(BM.pseudorxn, 'backbone'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4063', equations, 1);


% Total protein content in the cell [g protein/gDw]
Ptot  = 0.4227;

% Fraction of enzymes in the model [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;

model.lb(strcmp(model.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);


%Test model
% model   = setParam(model,'lb','r_1714',0);
%model   = setParam(model,'eq','r_1718',-0.634);
%model   = setParam(model,'eq','r_1634',-0.861);
% model   = setParam(model,'lb','r_4046',0);
% 
% model   = setParam(model,'lb','r_2111',0.055);
% model   = setParam(model,'obj','r_2111',0);
%model   = setParam(model,'lb','r_4046',0);

sol  = solveLP(model,1)
printFluxes(model, sol.x)


%save('ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Ac.mat','model')
