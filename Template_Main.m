%% Main Template
    % * Goal
        % 1.Estimate VAR;
        % 2.Identify structure shocks;
        % 3.Plot results:
            % 3.1.impulse response (default);
            % 3.2.variance decomposition (upon request);
            % 3.3.historical decomposition (upon request);
            % 3.4.conditional forecast (upon request).
    % * Instruction:
        % 1.Preparation:
            % 1.1.data to be used is labeled as "data.csv" and in the folder;
                % Data must be contain two columns of "Year" and "Month" 
            % 1.2.parameters in Prep_MainPara.m;
        % 2.Model:
        % 3.Plot:
            % 3.1.plot IRF;
            % 3.2.plot FEVD;
    
 
%% 0.Preparation.
% 0.1.clear space.
clear all;
close all;
clc;

%%
% 0.2.set work directory.
cd F:\GitHub\SVARToolBox;

%%
% 0.3.load data.
Prep_ImportData('Data_Template.xlsx');
load ..\SVARToolBox\DATASET DATASET_VAR;


%% 1.Model.
%
% 1.1.set parameter.
Template_Input_Para;

% 1.2.choose variables.
Template_Input_Estim;
VAR = Mod_Setup(DATASET_VAR, Para, Start, End, VAR_cell, VAR_label_cell, VAR_label_short_cell);


% 1.3.estimation.
% 1.3.1.estimation using OLS.
VAR_Est = Mod_EstOLS(VAR, [], [], 'non');
    VAR_Est_bs = Mod_EstOLS_bs(VAR_Est, Para);

% 1.4.identification.
% 1.4.1.identification with Mod_EstOLS.
Template_Input_Ident_basic;

% identification by short run.
VAR_Ident_short = Mod_EstOLS(VAR, short_order, short_shock, 'short');
    VAR_Ident_short_bs = Mod_EstOLS_bs(VAR_Ident_short, Para);
% Comment: When you specify shock position, you are specifying the ex post
% shock position. That is, the position after rearrangement.

% identification by long run.
VAR_Ident_long = Mod_EstOLS(VAR, long_order, long_shock, 'long');
    VAR_Ident_long_bs = Mod_EstOLS_bs(VAR_Ident_long, Para);

% 1.4.2.identification with Mod_EstOLS (require Mod_EstOLS output).
Template_Input_Ident_other;
    
% identification by proxy.
VAR_Ident_proxy = Mod_Ident(ident_input_proxy);
    VAR_Ident_proxy_bs = Ident_Proxy_bs(VAR_Ident_proxy, Para);

% identification by sign restriction.
VAR_Ident_sign = Mod_Ident(ident_input_sign);
    
% identification by max share.
VAR_Ident_Maxshare = Mod_Ident(ident_input_maxshare);
    VAR_Ident_Maxshare_bs = Ident_Maxshare_bs(VAR_Ident_Maxshare, Para);
% Comment: When you specify the variable to explain, you are specifying the
% ex post position, after rearrangment.
    


%% 2.Plot
% standardized plot.
Plot_IRF(VAR_Ident_short,VAR_Ident_short_bs,ParaPlot);
Plot_IRF_MultipleCI(VAR_Ident_short,VAR_Ident_short_bs,ParaPlot);
Plot_FEVD(VAR_Ident_short)

Plot_IRF(VAR_Ident_long,VAR_Ident_long_bs,ParaPlot);
Plot_IRF_MultipleCI(VAR_Ident_long,VAR_Ident_long_bs,ParaPlot);
Plot_FEVD(VAR_Ident_long)

Plot_IRF(VAR_Ident_proxy,VAR_Ident_proxy_bs,ParaPlot);
Plot_IRF_MultipleCI(VAR_Ident_proxy,VAR_Ident_proxy_bs,ParaPlot);
Plot_FEVD(VAR_Ident_proxy)

Plot_IRF(VAR_Ident_sign,VAR_Ident_sign,ParaPlot);
Plot_IRF_MultipleCI(VAR_Ident_sign,VAR_Ident_sign,ParaPlot);
Plot_FEVD(VAR_Ident_sign)

Plot_IRF(VAR_Ident_Maxshare,VAR_Ident_Maxshare_bs,ParaPlot);
Plot_IRF_MultipleCI(VAR_Ident_Maxshare,VAR_Ident_Maxshare_bs,ParaPlot);
Plot_FEVD(VAR_Ident_Maxshare)
