%% PARROT 

% Modified from: 


%% Running pFBA 

% Adaptaded from pFBA.m (COBRA Toolbox 3)

%% Configurations for MATLAB
initCobraToolbox(false)
changeCobraSolver('gurobi', 'LP');


%% Prepare model
%Adapted from https://dx.plos.org/10.1371/journal.pcbi.1011009

%Carbon uptake from Almeida et al. 2023 and oxygen, CO2 and prot_pool_exchange from
%the pFBA solution and grwoth adjusted from the FBA solution 

getpref('RAVEN');
setRavenSolver('gurobi');


%Load model 

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Wo.mat');

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.4227;

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

model = setParam(model, 'var', {'r_1992'}, -1.5756, 10);    % O2
model = setParam(model, 'var', {'r_1672'}, 1.7964, 10);    % CO2

model.lb(strcmp(model.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
model = setParam(model, 'eq', {'r_1714'}, 0); % close glucose
model = setParam(model, 'lb', {'r_4046'}, 0.05) %change lower bound NGAM

model = setParam(model, 'var', 'prot_pool_exchange', -105.75, 10);

model   = setParam(model,'var','r_2111',0.043,5);
%model   = setParam(model,'obj','r_2111',1);
%model   = setParam(model,'obj','r_2111',0);

model = setParam(model, 'lb', {'r_1718'}, (-0.657-0.002)); % open xylose
model = setParam(model, 'ub', {'r_1718'}, (-0.657+0.002)); % open xylose


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
% printFluxes(model, MinimizedFlux.x, false);



