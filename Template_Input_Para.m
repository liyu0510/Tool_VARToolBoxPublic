%% Prep_MainPara Module
    % * VAR parameter contains parameter for model setup and results.
    % * Bootstrap parameter contains parameter regarding bootstrapping for CI.
    % * Plot parameter contains parameters regarding single or multiple plots.
%%

%% 1.VAR Parameters
Para.p      = 4;                  % VAR lag length
Para.const  = 0;                  % constant term
Para.irhor  = 40;                 % Impulse response horizon


%% 2.Simulation Parameters
Para.nboot  = 1000;               % Number of bootstrap samples (equals 10000 in the paper)
Para.clevel = 95;   % 68;         % Bootstrap percentile shown
Para.CI     = [5 16 84 95];       % Bootstrap multiple percentile shown
% !If change length of CI, need to adjust Plot_IRF_MultipleCI accordingly.




%% 3.Plot Paramters
ParaPlot.nCol = 2;
ParaPlot.nRow = 3;
ParaPlot.fontsize  =   14;            % Fontsize in figures
