function [solutionALT, solStatus] = PARROT(modelREF, modelALT, MinStr, lambda, etotREF, etotALT, MinimizedFlux)
% Calculates the enzyme usage distribution of an alternative
% growth condition given a known reference in a reference condition.
%
% USAGE:
%
%    [solutionALT, solStatus] = PARROT(modelREF, modelALT, MinStr, lambda, etotREF, etotALT)
%
% INPUTS:
%    modelREF:          Protein-constrained model with integrated protein
%                       measurements to be used as a reference condition.
%    modelALT:          Protein-constrained model, batch version.
%
% OPTIONAL INPUTS:
%    MinStr:            Minimization strategy, 'Euclidean' or 'Manhattan'.
%                       (Default = 'Manhattan')
%    lambda:            Weighting factor for metabolic fluxes. Must be a 
%                       rational and positive value. (Default = 0)
%    etotREF:           Total enzyme usage for the reference condition
%                       (calculated from ptot * sigma * f) (Default = 0.5)
%    etotALT:           Total enzyme usage for the alternative condition
%                       (calculated from ptot * sigma * f) (Default = 0.5)
%
% OUTPUTS:
%    solutionALT:       Enzyme usage distribution for alternative condition
%    solStatus:         Solution status
%
% NOTE: Code partially adapted from the COBRA Toolbox function MOMA.m 
%    written by Markus Herrgard (2006).
%
% Author: - Mauricio Ferreira (2023)

if (nargin<3)
    MinStr = "Manhattan";
end

if (nargin<4)
    lambda = 0;
end

if (nargin<5)
    etotREF = 0.5;
end

if (nargin<6)
    etotALT = 0.5;
end

if lambda < 0
    error("Invalid value for lambda. It must be a rational and positive value.")
end

% if lambda > 0
%     % Generating V1 and E1 (REF)
%     [~,nRxns] = size(modelREF.S);
%     enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
%     enzymeIds(end) = [];
% 
%     infNum = sum(isinf(modelREF.ub(enzymeIds)));
%     infBounds = find(isinf(modelREF.ub));
%     infBoundsProt = intersect(infBounds, enzymeIds);
% 
%     etotShare = modelREF.ub(nRxns)/infNum;
%     modelREF.ub(infBoundsProt) = etotShare;
% 
%     NormMin = buildLPproblemFromModel(modelREF);
% 
%     proteins = modelREF.enzymes(:);
% 
%     pIdx = [];
%     rxnIdx = {};
%     kcat = [];
% 
%     for i=1:numel(proteins)
%         pIdx(i,1) = find(strcmpi(modelREF.metNames,join(['prot_' char(proteins(i))],"")));
%         rxnIdx{i,1} = find(modelREF.S(pIdx(i),:) < 0);
%         kcat(i,1) = min(-1./modelREF.S(pIdx(i),rxnIdx{i}));
%     end
% 
%     NormMin.c = zeros(size(NormMin.A,2),1);
%     NormMin.c(enzymeIds) = 1;
% 
%     NormMin.osense = 1;
% 
%     MinimizedFlux = solveCobraLP(NormMin);
% end

if MinStr == "Euclidean"
    % Construct the QP problem to find Es
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsALT] = size(modelALT.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_'))); 
    [enzymeNum,~] = size(enzymeIds);
    rxnIds = find(~cellfun('isempty',modelREF.rxns(1:min(enzymeIds)-1)));
    
    infNum = sum(isinf(modelREF.ub(enzymeIds)));
    infBounds = find(isinf(modelREF.ub));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = modelREF.ub(nRxnsREF)/infNum;
    modelREF.ub(infBoundsProt) = etotShare;
    
    QPproblem = buildLPproblemFromModel(modelALT);
    
    if lambda > 0
        QPproblem.c(rxnIds) = MinimizedFlux.full(rxnIds);
        QPproblem.c(rxnIds) = -2*lambda*QPproblem.c(rxnIds);
        QPproblem.c(enzymeIds) = MinimizedFlux.full(enzymeIds)/etotREF;
        QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
        
        QPproblem.F = sparse(size(QPproblem.A,2));
        QPproblem.F(1:nRxnsALT,1:nRxnsALT) = 2*lambda*speye(nRxnsALT);
        QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
        
        QPproblem.osense = 1;
    else
        QPproblem.c = zeros(nRxnsALT,1);
        QPproblem.c(enzymeIds) = modelREF.ub(enzymeIds);
        QPproblem.c(enzymeIds) = QPproblem.c(enzymeIds)/etotREF;
        QPproblem.c(enzymeIds) = -2*QPproblem.c(enzymeIds);
        
        QPproblem.F = sparse(size(QPproblem.A,2));
        QPproblem.F(1:nRxnsALT,1:nRxnsALT) = speye(nRxnsALT);
        QPproblem.F(1:nRxnsALT,1:nRxnsALT) = 0;
        QPproblem.F(min(enzymeIds):max(enzymeIds),min(enzymeIds):max(enzymeIds)) = 2*eye(enzymeNum);
        
        QPproblem.osense = 1;
    end
    
    % Solve the problem (COBRA)
    QPsolution = solveCobraQP(QPproblem);
    
    % Get the solution(s)
    solutionALT = QPsolution.full;
    solutionALT(enzymeIds) = solutionALT(enzymeIds) * etotALT;
    
    solStatus = QPsolution.stat;

elseif MinStr == "Manhattan"
    % Construct the LP problem to find ES2
    [~,nRxnsREF] = size(modelREF.S);
    [~,nRxnsALT] = size(modelALT.S);
    enzymeIds = find(~cellfun('isempty',strfind(modelREF.rxnNames,'prot_')));
    [enzymeNum,~] = size(enzymeIds);
    
    infNum = sum(isinf(modelREF.lb(enzymeIds)));
    infBounds = find(isinf(modelREF.lb));
    infBoundsProt = intersect(infBounds, enzymeIds);
    
    etotShare = etotALT/infNum;
    modelREF.lb(infBoundsProt) = etotShare;
    
    modelDelta = struct();
    
    % Set the lower and upper bounds for the fluxes
    modelDelta.lb = [modelREF.lb; modelALT.lb];
    modelDelta.lb(modelDelta.lb==-Inf) = -1000;
    
    modelDelta.ub = [modelREF.ub; modelALT.ub];
    modelDelta.ub(modelDelta.ub==Inf) = 1000;
    
    % Concatenate the stoichiometric matrices and S-matrices of the two models
    modelDelta.S = [modelREF.S zeros(size(modelREF.S,1), size(modelALT.S,2)); zeros(size(modelALT.S,1), size(modelREF.S,2)) modelALT.S];
    
    % Set the RHS vector
    modelDelta.b = zeros(size(modelDelta.S,1),1);
    
    % Objective function
    [~,nVars] = size(modelDelta.S);
    % enzymeIds(modelREF.ub(enzymeIds)==Inf) = [];
    enzymeIds_cindex = enzymeIds + nRxnsREF;
    rxns_index = nRxnsREF - enzymeNum;

    if lambda > 0
        modelDelta.c = zeros(nVars,1);
        modelDelta.c(1:rxns_index) = lambda*MinimizedFlux.full(1:rxns_index);
        modelDelta.c(enzymeIds) = ones(length(modelREF.rxns(enzymeIds)),1);
        modelDelta.c(enzymeIds) = (MinimizedFlux.full(enzymeIds)/etotREF);
        modelDelta.c(enzymeIds_cindex) = -ones(length(modelALT.rxns(enzymeIds)),1);
    else
        modelDelta.c = zeros(nVars,1);
        modelDelta.c(enzymeIds) = ones(length(modelREF.rxns(enzymeIds)),1);
        modelDelta.c(enzymeIds) = modelDelta.c(enzymeIds)/etotREF;
        modelDelta.c(enzymeIds_cindex) = -ones(length(modelALT.rxns(enzymeIds)),1);
    end
    
    % Solve the problem
    LPsolution = solveCobraLP(modelDelta);
    
    % Get the solution(s)
    solutionALT = LPsolution.full(length(modelREF.rxns)+1:end);
    solutionALT(enzymeIds) = solutionALT(enzymeIds) * etotALT;

    solStatus = LPsolution.stat;
    
    else
    error("Invalid optimization strategy, use 'Manhattan' or 'Euclidean'");
    
end
end
