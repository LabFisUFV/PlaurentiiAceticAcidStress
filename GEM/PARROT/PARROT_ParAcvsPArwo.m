%% Parrot 

% initCobraToolbox(false)

% Par Wo vs. Par Ac
clear;clc
%% Par Wo

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Wo.mat');

modelREF = model;

% Total protein content in the cell [g protein/gDw]
Ptot  = 0.4227;

% Fraction of enzymes in the modelREFREF [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;

etotREF = Ptot*f*sigma;

%Collect exchange reactions 

exchangeRxns = modelREF.rxns(endsWith(modelREF.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


modelREF = setParam(modelREF, 'lb', exchangeRxns, 0);
modelREF = setParam(modelREF, 'lb', requiredRxns, -1000);

modelREF = setParam(modelREF, 'var', {'r_1992'}, -1.5756, 10);    % O2
modelREF = setParam(modelREF, 'var', {'r_1672'}, 1.7964, 10);    % CO2

modelREF.lb(strcmp(modelREF.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
modelREF = setParam(modelREF, 'eq', {'r_1714'}, 0); % close glucose
modelREF = setParam(modelREF, 'lb', {'r_4046'}, 0.05); %change lower bound NGAM

modelREF = setParam(modelREF, 'var', 'prot_pool_exchange', -105.75, 10);

modelREF   = setParam(modelREF,'var','r_2111',0.043,5);
%modelREF   = setParam(modelREF,'obj','r_2111',1);
%modelREF   = setParam(modelREF,'obj','r_2111',0);

modelREF = setParam(modelREF, 'lb', {'r_1718'}, (-0.657-0.002)); % open xylose
modelREF = setParam(modelREF, 'ub', {'r_1718'}, (-0.657+0.002)); % open xylose

%% Construct the LP problem to perform pFBA

%modelREF = convertToIrreversible(modelREF);


[~,nRxns] = size(modelREF.S);
gprAssociated = find(sum(modelREF.rxnGeneMat,2)>0);


LPproblem = buildLPproblemFromModel(modelREF);

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

printFluxes(modelREF, MinimizedFlux.x, true);
% printFluxes(modelREF, MinimizedFlux.x, false);

%% Par Ac

clear model

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Ac.mat');

modelALT = model; 


% Total protein content in the cell [g protein/gDw]
Ptot  = 0.3608;

% Fraction of enzymes in the modelALT [g enzyme/g protein]
f     = 0.5;

% Average enzyme saturation factor
sigma = 0.5;

etotALT = Ptot*f*sigma;

%Collect exchange reactions 

exchangeRxns = modelALT.rxns(endsWith(modelALT.rxnNames,'exchange'));


%Determine essential reaction to similuate minimum media

%Adjust parameters specific to the cultivation condition
% block all uptake and allow only required metabolites

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange


modelALT = setParam(modelALT, 'var', {'r_1992'}, -2.2179, 10);    % O2
modelALT = setParam(modelALT, 'var', {'r_1672'}, 2.4992, 10);    % CO2

modelALT.lb(strcmp(modelALT.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
modelALT = setParam(modelALT, 'eq', {'r_1714'}, 0); % close glucose
modelALT = setParam(modelALT, 'lb', {'r_4046'}, 0.05); %change lower bound NGAM

modelALT = setParam(modelALT, 'var', 'prot_pool_exchange', -90.2, 10);

modelALT   = setParam(modelALT,'var','r_2111',0.055,5);
%modelALT   = setParam(modelALT,'obj','r_2111',1);
%modelALT   = setParam(modelALT,'obj','r_2111',0);

modelALT = setParam(modelALT,'eq','r_1761',0);
modelALT = setParam(modelALT,'eq','r_1631',0);

modelALT = setParam(modelALT, 'lb', {'r_1718'}, (-0.554-0.009)); % open xylose
modelALT = setParam(modelALT, 'ub', {'r_1718'}, (-0.554+0.009)); % open xylose

modelALT = setParam(modelALT, 'lb', {'r_1634'}, (-0.757-0.007)); % open acetate
modelALT = setParam(modelALT, 'ub', {'r_1634'}, (-0.757+0.007)); % open acetate

%% Run PARROT

enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
enzymeIds(end,:) = [];

modelREF_enzlb = MinimizedFlux.x(enzymeIds);


for i = 1:length(modelREF_enzlb)
    if modelREF_enzlb(i) == 0
       modelREF_enzlb(i) = -1000;
    end   
end

modelREF.lb(enzymeIds) = modelREF_enzlb;


MinStr = 'Euclidean';
lambda = 0;

[solutionALT, solStatus] = PARROT(modelREF, modelALT, MinStr, lambda, etotREF, etotALT, MinimizedFlux);


printFluxes(modelALT, solutionALT, true);

modelALT_enz = solutionALT(enzymeIds);
modelREF_enz = MinimizedFlux.x(enzymeIds);

% adjust to copies/pgDW

N_A = 6.02214076e23;

modelREF_enz = -modelREF_enz*N_A/1e18;
modelALT_enz = -modelALT_enz*N_A/1e18;

enzymeRxnName = modelREF.rxnNames(enzymeIds);
geneAssociation = modelREF.grRules(enzymeIds);
Result = table();

Result.Name = enzymeRxnName;
Result.GrRule = geneAssociation;
Result.REF = modelREF_enz;
Result.ALT = modelALT_enz;

writetable(Result, "ParWOvsParAc.csv");

