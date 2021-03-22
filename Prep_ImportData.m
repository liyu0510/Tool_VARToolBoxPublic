%% Prep_ImportData Function
    % Goal:
    	% Read in csv files with the first rowing being variable names.
        % Save data,label, and map in a struct DATASET_VAR.
        % Save data in a matrix.
    % Input:
        % Input 1: Dataset name.
    
        
%% 


function Prep_ImportData(Dataset)
DATASET_VAR.TSERIES = readtable(Dataset);                               % Reads data, with names
DATASET_VAR.LABEL   = DATASET_VAR.TSERIES.Properties.VariableNames;     % Extract variable names
DATASET_VAR.VALUE   = 1:length(DATASET_VAR.LABEL);

DATASET_VAR.MAP = containers.Map(DATASET_VAR.LABEL,DATASET_VAR.VALUE);
    % Creates a mapping between labels and values

DATASET_VAR.TSERIES = table2array(DATASET_VAR.TSERIES);           % Reads data, with names
    
save('..\SVARToolBox\DATASET.mat','DATASET_VAR');                     % Save data into a matrix

end
% 
% DATASET_VAR.TSERIES = readtable('BQ1989_Data.xlsx');                               % Reads data, with names
% DATASET_VAR.LABEL   = DATASET_VAR.TSERIES.Properties.VariableNames;         % Extract variable names
% DATASET_VAR.VALUE   = 1:length(DATASET_VAR.LABEL);
% DATASET_VAR.MAP = containers.Map(DATASET_VAR.LABEL,DATASET_VAR.VALUE);
% DATASET_VAR.TSERIES = table2array(DATASET_VAR.TSERIES);                               % Reads data, with names
%     
% save('..\ProjVAR\DATASET.mat','DATASET_VAR');
%     % Save data into a matrix   