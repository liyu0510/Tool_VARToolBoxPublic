%% Mod_Setup Function
    % * Goal of this function
    % Take specified sample period and variable choice.
    % Choose the corresponding data.
    % Form lag matrix with the chosen variables and sample periods.
%%

function [VAR] = Mod_Setup(DATASET_VAR, Para, Start, End, VAR_cell, VAR_label_cell, VAR_label_short_cell)

%% 1.Locate time in the dataset

VAR.Start           =   Start;
VAR.End             =   End;

smpl_min_VAR_vec    =   VAR.Start;                  % [Year month;Year month]    
smpl_max_VAR_vec    =   VAR.End;                    % End-of-sample (change for robustness exercises)   

VAR.smpl_min_VAR = find(and(DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'Year'})))==smpl_min_VAR_vec(1,1), ...
    DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'Month'})))==smpl_min_VAR_vec(1,2)));
VAR.smpl_max_VAR = find(and(DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'Year'})))==smpl_max_VAR_vec(1,1), ...
    DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'Month'})))==smpl_max_VAR_vec(1,2)));
%%

%% 2.Extract chosen variables from dataset
VAR.select_vars               = VAR_cell;
VAR.select_vars_label         = VAR_label_cell;
VAR.select_vars_label_short   = VAR_label_short_cell;
%%

%% 3.Specify parameter 
VAR.p          =   Para.p;                                    % VAR lag length
VAR.const      =   Para.const;
VAR.irhor      =   Para.irhor;                                % Impulse Response Horizon
%%

%% 4.Setup model
VAR.vars             = DATASET_VAR.TSERIES(VAR.smpl_min_VAR:VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,VAR.select_vars)));
VAR.nobs             = size(VAR.vars,1);        % total number of observations
VAR.T                = length(VAR.vars)-abs(VAR.p);  % real number of observations that make it to the regression
VAR.n                = size(VAR.vars,2);        % dimension of multivariate (excluding constant)
VAR.year             = DATASET_VAR.TSERIES(VAR.smpl_min_VAR:VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,{'Year'})));
VAR.month            = DATASET_VAR.TSERIES(VAR.smpl_min_VAR:VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,{'Month'})));









