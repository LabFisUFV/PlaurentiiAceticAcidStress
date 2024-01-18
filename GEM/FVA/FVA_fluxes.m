%FVA 

%%
% initCobraToolbox(false)

%% Cleaning the workspace and the command window
clear;clc

%% Par wo

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Wo.mat');
model_Par_wo = model; 

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.4227;

% Fraction of enzymes in the model [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;



%Collect exchange reactions 

exchangeRxns = model_Par_wo.rxns(endsWith(model_Par_wo.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


model_Par_wo = setParam(model_Par_wo, 'lb', exchangeRxns, 0);
model_Par_wo = setParam(model_Par_wo, 'lb', requiredRxns, -1000);

model_Par_wo = setParam(model_Par_wo, 'var', {'r_1992'}, -1.5756, 10);    % O2
model_Par_wo = setParam(model_Par_wo, 'var', {'r_1672'}, 1.7964, 10);    % CO2

model_Par_wo.lb(strcmp(model_Par_wo.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model_Par_wo = setParam(model_Par_wo, 'eq', {'r_1714'}, 0); % close glucose
model_Par_wo = setParam(model_Par_wo, 'lb', {'r_4046'}, 0.05) %change lower bound NGAM

model_Par_wo = setParam(model_Par_wo, 'var', 'prot_pool_exchange', -105.75, 10);

model_Par_wo   = setParam(model_Par_wo,'var','r_2111',0.043,5);
%model_Par_wo   = setParam(model_Par_wo,'obj','r_2111',1);
%model_Par_wo   = setParam(model_Par_wo,'obj','r_2111',0);

model_Par_wo = setParam(model_Par_wo, 'lb', {'r_1718'}, (-0.657-0.002)); % open xylose
model_Par_wo = setParam(model_Par_wo, 'ub', {'r_1718'}, (-0.657+0.002)); % open xylose


% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model_Par_wo.rxns,'prot_'));
toRemove1 = model_Par_wo.rxnNames(toRemove1);

rxnNamesList = model_Par_wo.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);


[minFlux, maxFlux] = fluxVariability(model_Par_wo, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_Par_wo.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear maxFlux minFlux ranges

model_conv = readCbModel('papla-GEM.xml');
[minFlux, maxFlux] = ecFVA(model_Par_wo, model_conv)
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(model_conv.rxns,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_Par_wo_mapped.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear
%% Par Ac

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Ac.mat');

model_Par_Ac = model;

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.3608;

% Fraction of enzymes in the model [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;

%Collect exchange reactions 

exchangeRxns = model_Par_Ac.rxns(endsWith(model_Par_Ac.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


model_Par_Ac = setParam(model_Par_Ac, 'lb', exchangeRxns, 0);
model_Par_Ac = setParam(model_Par_Ac, 'lb', requiredRxns, -1000);

model_Par_Ac = setParam(model_Par_Ac, 'var', {'r_1992'}, -2.2179, 10);    % O2
model_Par_Ac = setParam(model_Par_Ac, 'var', {'r_1672'}, 2.4992, 10);    % CO2

model_Par_Ac.lb(strcmp(model_Par_Ac.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model_Par_Ac = setParam(model_Par_Ac, 'eq', {'r_1714'}, 0); % close glucose
model_Par_Ac = setParam(model_Par_Ac, 'lb', {'r_4046'}, 0.05); %change lower bound NGAM

model_Par_Ac = setParam(model_Par_Ac, 'var', 'prot_pool_exchange', -90.2, 10);

model_Par_Ac   = setParam(model_Par_Ac,'var','r_2111',0.055,5);
%model_Par_Ac   = setParam(model_Par_Ac,'obj','r_2111',1);
%model_Par_Ac   = setParam(model_Par_Ac,'obj','r_2111',0);

model_Par_Ac = setParam(model_Par_Ac, 'lb', {'r_1718'}, (-0.554-0.009)); % open xylose
model_Par_Ac = setParam(model_Par_Ac, 'ub', {'r_1718'}, (-0.554+0.009)); % open xylose

model_Par_Ac = setParam(model_Par_Ac, 'lb', {'r_1634'}, (-0.757-0.007)); % open acetate
model_Par_Ac = setParam(model_Par_Ac, 'ub', {'r_1634'}, (-0.757+0.007)); % open acetate

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model_Par_Ac.rxns,'prot_'));
toRemove1 = model_Par_Ac.rxnNames(toRemove1);

rxnNamesList = model_Par_Ac.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);


[minFlux, maxFlux] = fluxVariability(model_Par_Ac, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_Par_Ac.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear maxFlux minFlux ranges

model_conv = readCbModel('papla-GEM.xml');
[minFlux, maxFlux] = ecFVA(model_Par_Ac, model_conv);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(model_conv.rxns,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_Par_Ac_mapped.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear




%% ATS wo

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_ATS_Wo.mat');

model_ATS_wo = model;


% Total protein content in the cell [g protein/gDw]
Ptot  = 0.3967;

% Fraction of enzymes in the model_ATS_wo [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;



%Collect exchange reactions 

exchangeRxns = model_ATS_wo.rxns(endsWith(model_ATS_wo.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


model_ATS_wo = setParam(model_ATS_wo, 'lb', exchangeRxns, 0);
model_ATS_wo = setParam(model_ATS_wo, 'lb', requiredRxns, -1000);

model_ATS_wo = setParam(model_ATS_wo, 'var', {'r_1992'}, -1.6266, 10);    % O2
model_ATS_wo = setParam(model_ATS_wo, 'var', {'r_1672'}, 1.8585, 10);    % CO2

model_ATS_wo.lb(strcmp(model_ATS_wo.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model_ATS_wo = setParam(model_ATS_wo, 'eq', {'r_1714'}, 0); % close glucose
model_ATS_wo = setParam(model_ATS_wo, 'lb', {'r_4046'}, 0.05) %change lower bound NGAM

model_ATS_wo = setParam(model_ATS_wo, 'var', 'prot_pool_exchange', -99.175, 10);

model_ATS_wo   = setParam(model_ATS_wo,'var','r_2111',0.045,5);
%model_ATS_wo   = setParam(model_ATS_wo,'obj','r_2111',1);
%model_ATS_wo   = setParam(model_ATS_wo,'obj','r_2111',0);

model_ATS_wo = setParam(model_ATS_wo, 'lb', {'r_1718'}, (-0.678-0.010)); % open xylose
model_ATS_wo = setParam(model_ATS_wo, 'ub', {'r_1718'}, (-0.678+0.010)); % open xylose

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model_ATS_wo.rxns,'prot_'));
toRemove1 = model_ATS_wo.rxnNames(toRemove1);

rxnNamesList = model_ATS_wo.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);


[minFlux, maxFlux] = fluxVariability(model_ATS_wo, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(rxnNamesList,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_ATS_wo.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear maxFlux minFlux ranges

model_conv = readCbModel('papla-GEM.xml');
[minFlux, maxFlux] = ecFVA(model_ATS_wo, model_conv)
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(model_conv.rxns,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_ATS_wo_mapped.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear
%% ATS Ac

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_ATS_Ac.mat');

model_ATS_Ac = model; 

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.3881;

% Fraction of enzymes in the model [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;



%Collect exchange reactions 

exchangeRxns = model_ATS_Ac.rxns(endsWith(model_ATS_Ac.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


model_ATS_Ac = setParam(model_ATS_Ac, 'lb', exchangeRxns, 0);
model_ATS_Ac = setParam(model_ATS_Ac, 'lb', requiredRxns, -1000);

model_ATS_Ac = setParam(model_ATS_Ac, 'var', {'r_1992'}, -2.5135, 10);    % O2
model_ATS_Ac = setParam(model_ATS_Ac, 'var', {'r_1672'}, 2.828, 10);    % CO2

model_ATS_Ac.lb(strcmp(model_ATS_Ac.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model_ATS_Ac = setParam(model_ATS_Ac, 'eq', {'r_1714'}, 0); % close glucose
model_ATS_Ac = setParam(model_ATS_Ac, 'lb', {'r_4046'}, 0.05); %change lower bound NGAM

model_ATS_Ac = setParam(model_ATS_Ac, 'var', 'prot_pool_exchange', -97.025, 10);

model_ATS_Ac   = setParam(model_ATS_Ac,'var','r_2111',0.055,5);
%model_ATS_Ac   = setParam(model_ATS_Ac,'obj','r_2111',1);
%model_ATS_Ac   = setParam(model_ATS_Ac,'obj','r_2111',0);

model_ATS_Ac = setParam(model_ATS_Ac, 'lb', {'r_1718'}, (-0.634-0.002)); % open xylose
model_ATS_Ac = setParam(model_ATS_Ac, 'ub', {'r_1718'}, (-0.634+0.002)); % open xylose

model_ATS_Ac = setParam(model_ATS_Ac, 'lb', {'r_1634'}, (-0.861-0.017)); % open acetate
model_ATS_Ac = setParam(model_ATS_Ac, 'ub', {'r_1634'}, (-0.861+0.017)); % open acetate

% Perform Flux Variability Analysis on enzyme usage
toRemove1 = find(contains(model_ATS_Ac.rxns,'prot_'));
toRemove1 = model_ATS_Ac.rxnNames(toRemove1);

rxnNamesList = model_ATS_Ac.rxns;
rxnNamesList = setdiff(rxnNamesList, toRemove1);


%[minFlux, maxFlux] = fluxVariability(model_ATS_Ac, 99, 'max', rxnNamesList);
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(model_conv.rxns,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_ATS_Ac_mapped.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear maxFlux minFlux ranges

model_conv = readCbModel('papla-GEM.xml');
[minFlux, maxFlux] = ecFVA(model_ATS_Ac, model_conv)
ranges = maxFlux - minFlux;

% Export results table
varNamesT = {'Rxns' 'ranges' 'minFlux' 'maxFlux'};
FVAtable  = table(model_conv.rxns,ranges,minFlux,maxFlux,'VariableNames', varNamesT);
FVA_filename = "FVA_fluxes_ATS_Ac_mapped.csv";
writetable(FVAtable, FVA_filename, 'Delimiter','\t')

clear

%% Plot enzyme usage FVA

%%
clc;clear

%%
Par_wo = readtable("FVA_fluxes_Par_wo.csv");
Par_Ac = readtable("FVA_fluxes_Par_Ac.csv");
ATS_wo = readtable("FVA_fluxes_ATS_wo.csv");
ATS_Ac = readtable("FVA_fluxes_ATS_Ac.csv");

%%
distributions = {Par_wo.ranges, Par_Ac.ranges, ATS_wo.ranges, ATS_Ac.ranges};
legends       = {'Parental strain w/o acetic acid', 'Parental strain with acetic acid','ATS I strain w/o acetic acid', 'ATS I with acetic acid'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);

%%
Par_wo = readtable("FVA_fluxes_Par_wo_mapped.csv");
Par_Ac = readtable("FVA_fluxes_Par_Ac_mapped.csv");
ATS_wo = readtable("FVA_fluxes_ATS_wo_mapped.csv");
ATS_Ac = readtable("FVA_fluxes_ATS_Ac_mapped.csv");

%%
distributions = {Par_wo.ranges, Par_Ac.ranges, ATS_wo.ranges, ATS_Ac.ranges};
legends       = {'Parental strain w/o acetic acid', 'Parental strain with acetic acid','ATS I strain w/o acetic acid', 'ATS I with acetic acid'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);
