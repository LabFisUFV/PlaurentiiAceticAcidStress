%%pFBA and Random_sampling

%Biomass_specific models with data from Almeida et al. (2023)
%First, it was necessary to adjust the dilution rate (ub to 0.055) and NGAM (lb to 0.05) to allow
%growth in such restriced uptake rates. This limitation is related to the
%model minimum growth allowance during reconstruction. 



%% Cleaning the workspace and the command window
clear;clc

%% ATS_Ac

%Load model 

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_ATS_Ac.mat');

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.3881;

% Fraction of enzymes in the model [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;



%Collect exchange reactions 

exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);

model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);

model.lb(strcmp(model.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model = setParam(model, 'eq', {'r_1714'}, 0); % close glucose
model = setParam(model, 'lb', {'r_4046'}, 0.05); %change lower bound NGAM

model   = setParam(model,'lb','r_2111',0.062);
%model   = setParam(model,'obj','r_2111',1);
model   = setParam(model,'obj','r_2111',0);
%model   = setParam(model,'obj','r_4046',1);


model = setParam(model, 'lb', {'r_1718'}, (-0.634-0.002)); % open xylose
model = setParam(model, 'ub', {'r_1718'}, (-0.634+0.002)); % open xylose

model = setParam(model, 'lb', {'r_1634'}, (-0.861-0.017)); % open acetate
model = setParam(model, 'ub', {'r_1634'}, (-0.861+0.017)); % open acetate

% sol  = solveLP(model,1);
% printFluxes(model, sol.x);

%% Running pFBA 

% Adaptaded from pFBA.m (COBRA Toolbox 3)

%% Configurations for MATLAB
changeCobraSolver('gurobi', 'LP');


%% Construct the LP problem to perform pFBA

%model = convertToIrreversible(model);


[~,nRxns] = size(model.S);
gprAssociated = find(sum(model.rxnGeneMat,2)>0);


LPproblem = buildLPproblemFromModel(model);

LPproblem.c = zeros(nRxns,1);
LPproblem.c(gprAssociated) = 1;

LPproblem.osense = 1;

%% Verify if the problem is a valid LP problem
fprintf('\n');
statusOK = verifyCobraProblem(LPproblem);
if statusOK == -1
    disp('Invalid LP problem')
end

%% Solve the problem
MinimizedFlux = solveCobraLP(LPproblem);

%% Get the solution(s)
MinimizedFlux.x = MinimizedFlux.full;

fprintf('\n');
fprintf('LP SOLUTION\n');

printFluxes(model, MinimizedFlux.x, true);
%printFluxes(model, MinimizedFlux.x, false);


%% Prepare for random sampling
%Adapted from https://dx.plos.org/10.1371/journal.pcbi.1011009

%Carbon uptake from Almeida et al. 2023 and oxygen, CO2 and prot_pool_exchange from
%the pFBA solution and grwoth adjusted from the FBA solution 

getpref('RAVEN');
setRavenSolver('gurobi');


%Load model 

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_ATS_Ac.mat');

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.3881;

% Fraction of enzymes in the model [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;



%Collect exchange reactions 

exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);

model = setParam(model, 'var', {'r_1992'}, -2.5135, 10);    % O2
model = setParam(model, 'var', {'r_1672'}, 2.828, 10);    % CO2

model.lb(strcmp(model.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model = setParam(model, 'eq', {'r_1714'}, 0); % close glucose
model = setParam(model, 'lb', {'r_4046'}, 0.05); %change lower bound NGAM

model = setParam(model, 'var', 'prot_pool_exchange', -97.025, 10);

model   = setParam(model,'var','r_2111',0.055,5);
%model   = setParam(model,'obj','r_2111',1);
%model   = setParam(model,'obj','r_2111',0);

model = setParam(model, 'lb', {'r_1718'}, (-0.634-0.002)); % open xylose
model = setParam(model, 'ub', {'r_1718'}, (-0.634+0.002)); % open xylose

model = setParam(model, 'lb', {'r_1634'}, (-0.861-0.017)); % open acetate
model = setParam(model, 'ub', {'r_1634'}, (-0.861+0.017)); % open acetate

%W/o good reactions - multiple carbon sources in the case of acetic acid

goodRxns = []; 

sol = rs(model,10000,true,true,true,goodRxns,true);


%% output: ec

out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model.rxns;

writetable(struct2table(out),'out_ATS_Ac.csv');

model_conv = readCbModel('papla-GEM.xml');
solMapped = mapRxnsToConv(model,model_conv,sol);

outMapped.rxns = model_conv.rxns;
outMapped.rxnNames = model_conv.rxnNames;
outMapped.median = full(median(solMapped,2));  %2 means for each row
outMapped.mean=full(mean(solMapped,2));
outMapped.std=full(std(solMapped,0,2));

writetable(struct2table(outMapped),'out_Mapped_ATS_Ac.csv');

%% calculate ATP, NADPH, and NADH balances from non-ec fluxes:

model_conv = ravenCobraWrapper(model_conv);

% adjust for each compound and repeat for each condition

for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.median,true);
    outATP_median.rxns    = model_conv.rxns(rxnIdx);
    outATP_median.rxnNames= model_conv.rxnNames(rxnIdx);
    outATP_median.rxnEqns = constructEquations(model_conv,rxnIdx);
    outATP_median.fluxes  = num2cell(fluxes);
    outATP_median = [outATP_median.rxns outATP_median.rxnNames outATP_median.rxnEqns outATP_median.fluxes];
end

writecell(outATP_median,'out_Mapped_ATS_Ac_ATP_median.csv');

for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.mean,true);
    outATP_mean.rxns    = model_conv.rxns(rxnIdx);
    outATP_mean.rxnNames= model_conv.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model_conv,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end

writecell(outATP_mean,'out_Mapped_ATS_Ac_ATP_mean.csv');


for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.median,true);
    outNADH_median.rxns    = model_conv.rxns(rxnIdx);
    outNADH_median.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end

writecell(outNADH_median,'out_Mapped_ATS_Ac_NADH_median.csv');

for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.mean,true);
    outNADH_mean.rxns    = model_conv.rxns(rxnIdx);
    outNADH_mean.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end

writecell(outNADH_mean,'out_Mapped_ATS_Ac_NADH_mean.csv');

for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.median,true);
    outNADPH_median.rxns    = model_conv.rxns(rxnIdx);
    outNADPH_median.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end

writecell(outNADPH_median,'out_Mapped_ATS_Ac_NADPH_median.csv');

for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.mean,true);
    outNADPH_mean.rxns    = model_conv.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end

writecell(outNADPH_mean,'out_Mapped_ATS_Ac_NADPH_mean.csv');



