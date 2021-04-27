%% BQ 1989
        
%% 

%% 0.Preparation
% 0.1.clear space
clear all;
close all;
clc;


% 0.2.set work directory
cd F:\GitHub\SVARToolBox;


% 0.3.load data
Prep_ImportData('Data_BQ1989.xlsx')    
load ..\SVARToolBox\DATASET DATASET_VAR;


% 0.4.set parameter
% VAR Parameters
Para.p      = 8;                  % VAR lag length
Para.const  = 1;  
Para.irhor  = 40;                 % Impulse response horizon

% Simulation Parameters
Para.nboot  = 1000;               % Number of bootstrap samples (equals 10000 in the paper)
Para.clevel = 95;   % 68;         % Bootstrap percentile shown
Para.CI     = [5 16 84 95];       % Bootstrap multiple percentile shown

% Plot Paramters
ParaPlot.nCol = 2;
ParaPlot.nRow = 3;
ParaPlot.M    = 1;
ParaPlot.fontsize  =   14;            % Fontsize in figures

%%

%% 1.Model

% choose variables
Start = [1948 10];
End   = [1987 10];
VAR_cell        =    {'GDPGrowth', 'Unemployment'};  
VAR_label_cell  =    {'FD Log GDP', 'Unemployment'}; 
VAR_label_short_cell  =    {'GDP', 'Unemployment'};

VAR = Mod_Setup(DATASET_VAR, Para, Start, End, VAR_cell, VAR_label_cell, VAR_label_short_cell);

% input for long run identification
long_order = [1, 2];
long_shock  = 2;


% long run identification
VAR_Ident_long = Mod_EstOLS(VAR, long_order, long_shock, 'long');
    VAR_Ident_long_bs = Mod_EstOLS_bs(VAR_Ident_long, Para);
        
VAR_Ident_long_total = Mod_EstOLS(VAR, long_order, [], 'long');
    VAR_Ident_long_total_bs = Mod_EstOLS_bs(VAR_Ident_long_total, Para);

%%

%% 2.Plot 
% standardized plot
Plot_IRF_MultipleCI(VAR_Ident_long,VAR_Ident_long_bs,ParaPlot);
Plot_IRF_Cumulative(VAR_Ident_long,VAR_Ident_long_bs,ParaPlot);

% exported plot
subplot(1,2,1)
plot(cumsum(VAR_Ident_long_total.IRF_total(:,1,1)),'LineWidth',2) % only cumsum here because it is FD log
hold on
plot(VAR_Ident_long_total.IRF_total(:,2,1),'--r','LineWidth',2)
hold on
plot(zeros(VAR_Ident_long_total.irhor),'-k')
title('Supply shock')
legend({'GDP Level';'Unemployment'})

subplot(1,2,2)
plot(cumsum(-VAR_Ident_long_total.IRF_total(:,1,2)),'LineWidth',2)
hold on
plot(-VAR_Ident_long_total.IRF_total(:,2,2),'--r','LineWidth',2)
hold on
plot(zeros(VAR_Ident_long_total.irhor),'-k')
title('Demand shock')
legend({'GDP Level';'Unemployment'})


%%


