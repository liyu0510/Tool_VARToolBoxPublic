%% Uhlig 2005
        
%% 

%% 0.Preparation
% 0.1.clear space
clear all;
close all;
clc;


% 0.2.set work directory
cd F:\GitHub\SVARToolBox;


% 0.3.load data
Prep_ImportData('Data_Uhlig2005.xlsx')    
load ..\SVARToolBox\DATASET DATASET_VAR;


% 0.4.set parameter
% VAR Parameters
Para.p      = 12;                  % VAR lag length
Para.const  = 0;  
Para.irhor  = 60;                 % Impulse response horizon
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
Start = [1965 1];
End   = [2003 12];
					

VAR_cell              =    {'RealGDP', 'CPI','CommodityPrice','FedFund','NonBorrReserves','TotalReserves'};  
VAR_label_cell        =    {'RealGDP', 'CPI','CommodityPrice','FedFund','NonBorrReserves','TotalReserves'};   
VAR_label_short_cell  =    {'RGDP', 'CPI','CP','EFFR','NBReser','TReser'};  

VAR = Mod_Setup(DATASET_VAR, Para, Start, End, VAR_cell, VAR_label_cell, VAR_label_short_cell);

% short run identification 
VAR_Ident_short = Mod_EstOLS(VAR, [1,2,3,4,5,6], 4, 'short');


% sign restriction identification
ndraw = 1000;
R(:,:,1) = [ 0     0     0    % Real GDP
             1     6    -1    % Inflation
             1     6    -1    % Commodity Price
             1     6     1    % Fed Fund
             1     6    -1    % NonBorr. Reserves
             0     0     0];  % Total Reserves

ident_input_sign     = {'sign', Para, VAR_Ident_short, R};

VAR_Ident_sign = Mod_Ident(ident_input_sign);
  
%%

%% 2.Plot 
% standardized plot
Plot_IRF_MultipleCI(VAR_Ident_sign,VAR_Ident_sign,ParaPlot);

%%


