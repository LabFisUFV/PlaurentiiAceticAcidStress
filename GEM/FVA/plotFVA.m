%
% Plot enzyme usage FVA
%
%%
clc;clear

%%
noBio_noUnd = readtable("FVA_fluxes_eciML1515u_CORAL_DLKcat_noBio_noUnd.csv");
noBio_Und = readtable("FVA_fluxes_eciML1515u_CORAL_DLKcat_noBio_Und.csv");
Bio_noUnd = readtable("FVA_fluxes_eciML1515u_CORAL_DLKcat_Bio_noUnd.csv");
Bio_Und = readtable("FVA_fluxes_eciML1515u_CORAL_DLKcat_Bio_Und.csv");

%%
distributions = {noBio_noUnd.ranges, noBio_Und.ranges, Bio_noUnd.ranges, Bio_Und.ranges};
legends       = {'No underground reactions, no fixed biomass', 'Underground reactions, no fixed biomass','No underground reactions, fixed biomass', 'Underground reactions, fixed biomass'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);
