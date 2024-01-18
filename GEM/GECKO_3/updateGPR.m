% Load your genome-scale metabolic model
load('papla-GEM.mat'); % Replace 'your_model.mat' with the actual filename of your model

% Define the two lists of genes
originalGenes = OrtPaplaCGattii(:,1); % Replace with your actual gene names from the model
updatedGenes = OrtPaplaCGattii(:,2); % Replace with the corresponding updated gene names

% Create a mapping between the original genes and updated genes
geneMap = containers.Map(originalGenes, updatedGenes);

% Update the gene names in model.genes and model.geneShortNames
for i = 1:numel(model.genes)
    if isKey(geneMap, model.genes{i})
        model.genes{i} = geneMap(model.genes{i});
    end
end

for i = 1:numel(model.geneShortNames)
    if isKey(geneMap, model.geneShortNames{i})
        model.geneShortNames{i} = geneMap(model.geneShortNames{i});
    end
end


% Iterate over the reactions in the model and update the GPR rules
for i = 1:numel(model.rxns)
    % Get the original GPR rule for the current reaction
    originalGPR = model.grRules{i};
    
    % Split the original GPR rule into individual gene associations
    geneAssociations = strsplit(originalGPR, ' or ');
    
    % Iterate over the gene associations and update the gene names
    updatedAssociations = cell(size(geneAssociations));
    for j = 1:numel(geneAssociations)
        % Split the gene association into individual genes
        genes = strsplit(geneAssociations{j}, ' and ');
        
        % Update the gene names using the mapping
        updatedGenes = cell(size(genes));
        for k = 1:numel(genes)
            % Check if the gene name touches "(" or ")"
            if startsWith(genes{k}, '(') && endsWith(genes{k}, ')')
                updatedGenes{k} = genes{k};
            else
                % Remove brackets and update the gene name using the mapping
                gene = regexprep(genes{k}, '[\(\)]', '');
                if isKey(geneMap, gene)
                    updatedGenes{k} = regexprep(genes{k}, gene, geneMap(gene));
                else
                    updatedGenes{k} = genes{k};
                end
            end
        end
        
        % Join the updated gene names back into a single string
        updatedAssociations{j} = strjoin(updatedGenes, ' and ');
    end
    
    % Join the updated gene associations back into a single GPR rule
    updatedGPR = strjoin(updatedAssociations, ' or ');
    
    % Update the GPR rule in the model
    model.grRules{i} = updatedGPR;
end


unique_genes = unique(model.genes);

model.genes=unique_genes;

model.geneShortNames=unique_genes;

% Save the updated model
save('updated_model.mat', 'model'); % Replace 'updated_model.mat' with the desired filename for the updated model

