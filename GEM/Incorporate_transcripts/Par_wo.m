%% Par wo
%% Load model

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\GEM-biomass_Adj\ecpaplaGEM_stage3_BRENDA_wostd_ALE_Par_Wo.mat');

%% Load FVA data

path = "C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\FVA\FVA_fluxes_Par_wo.txt";

FVA_data = readtable(path, 'Delimiter', '\t');

%% Set boundaries 

for i=1:length(FVA_data.Rxns)
    if isequal(FVA_data.Rxns{i}, model.rxns{i})
        model = setParam(model, 'lb', model.rxns(i), FVA_data.minFlux(i));
        model = setParam(model, 'ub', model.rxns(i), FVA_data.maxFlux(i));
    end
end

%% Set Condition 

%% Incorporate transcripts

clear path i j

path = "C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\Incorporate_transcripts\AND_reactions_Norm_RAS_wo_zeros.csv";

Norm_Ras_data = readtable(path, 'Delimiter', '\t');

% AND reactions

model_1 = model;

Norm_Ras_data.Norm_Par_wo_1 = Norm_Ras_data.Norm_Par_wo_1;


for i=1:length(Norm_Ras_data.Rxn)
    for j=1:length(model.rxns)
        if isequal(Norm_Ras_data.Rxn{i}, model.rxns{j})
        model_1 = setParam(model_1, 'lb', model.rxns(j), model_1.lb(j)*Norm_Ras_data.Norm_Par_wo_1(i));
        model_1 = setParam(model_1, 'ub', model.rxns(j), model_1.ub(j)*Norm_Ras_data.Norm_Par_wo_1(i));
        end
    end
end


% Other reactions

clear path i j Norm_Ras_data 

path = "C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\p_laurentii_acetic_acid_stress\Papiliotrema_laurentii_acetic_acid_stress\GEM\Incorporate_transcripts\Expression_singlwo_zero.csv";

Norm_Ras_data = readtable(path, 'Delimiter', '\t');

for i=1:length(Norm_Ras_data.Gene_gattii)
    for j=1:length(model.rxns)
        if isequal(Norm_Ras_data.Gene_gattii{i}, model.grRules{j})
        model_1 = setParam(model_1, 'lb', model.rxns(j), model_1.lb(j)*Norm_Ras_data.Norm_Par_wo_1(i));
        model_1 = setParam(model_1, 'ub', model.rxns(j), model_1.ub(j)*Norm_Ras_data.Norm_Par_wo_1(i));
        end
    end
end


%%

goodRxns = []; 

sol = rs(model_1,10000,false,true,true,goodRxns,false);


%% output: ec

out.median=full(median(sol,2));  %2 means for each row
out.mean=full(mean(sol,2));
out.std=full(std(sol,0,2));
out.rxns = model_1.rxns;

writetable(struct2table(out),'out_Par_wo.csv');

model_conv = readCbModel('papla-GEM.xml');
solMapped = mapRxnsToConv(model_1,model_conv,sol);

outMapped.rxns = model_conv.rxns;
outMapped.rxnNames = model_conv.rxnNames;
outMapped.median = full(median(solMapped,2));  %2 means for each row
outMapped.mean=full(mean(solMapped,2));
outMapped.std=full(std(solMapped,0,2));

writetable(struct2table(outMapped),'out_Mapped_Par_wo.csv');

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

writecell(outATP_median,'out_Mapped_Par_wo_ATP_median.csv');

for i={'ATP'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.mean,true);
    outATP_mean.rxns    = model_conv.rxns(rxnIdx);
    outATP_mean.rxnNames= model_conv.rxnNames(rxnIdx);
    outATP_mean.rxnEqns = constructEquations(model_conv,rxnIdx);
    outATP_mean.fluxes  = num2cell(fluxes);
    outATP_mean = [outATP_mean.rxns outATP_mean.rxnNames outATP_mean.rxnEqns outATP_mean.fluxes];
end

writecell(outATP_mean,'out_Mapped_Par_wo_ATP_mean.csv');


for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.median,true);
    outNADH_median.rxns    = model_conv.rxns(rxnIdx);
    outNADH_median.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADH_median.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADH_median.fluxes  = num2cell(fluxes);
    outNADH_median = [outNADH_median.rxns outNADH_median.rxnNames outNADH_median.rxnEqns outNADH_median.fluxes];
end

writecell(outNADH_median,'out_Mapped_Par_wo_NADH_median.csv');

for i={'NADH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.mean,true);
    outNADH_mean.rxns    = model_conv.rxns(rxnIdx);
    outNADH_mean.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADH_mean.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADH_mean.fluxes  = num2cell(fluxes);
    outNADH_mean = [outNADH_mean.rxns outNADH_mean.rxnNames outNADH_mean.rxnEqns outNADH_mean.fluxes];
end

writecell(outNADH_mean,'out_Mapped_Par_wo_NADH_mean.csv');

for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.median,true);
    outNADPH_median.rxns    = model_conv.rxns(rxnIdx);
    outNADPH_median.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADPH_median.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADPH_median.fluxes  = num2cell(fluxes);
    outNADPH_median = [outNADPH_median.rxns outNADPH_median.rxnNames outNADPH_median.rxnEqns outNADPH_median.fluxes];
end

writecell(outNADPH_median,'out_Mapped_Par_wo_NADPH_median.csv');

for i={'NADPH'}
    [fluxes, rxnIdx] = getMetProduction(model_conv,i,outMapped.mean,true);
    outNADPH_mean.rxns    = model_conv.rxns(rxnIdx);
    outNADPH_mean.rxnNames= model_conv.rxnNames(rxnIdx);
    outNADPH_mean.rxnEqns = constructEquations(model_conv,rxnIdx);
    outNADPH_mean.fluxes  = num2cell(fluxes);
    outNADPH_mean = [outNADPH_mean.rxns outNADPH_mean.rxnNames outNADPH_mean.rxnEqns outNADPH_mean.fluxes];
end

writecell(outNADPH_mean,'out_Mapped_Par_wo_NADPH_mean.csv');

