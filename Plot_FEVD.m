function Plot_FEVD(VAR)
%% 1.Import parameter

FEVD = VAR.FEVD;

% Initialize FEVD matrix
[nsteps]  = size(FEVD,1);
[nshocks] = size(FEVD,2);   % Because proxy VAR only identified one shock
[nvars]   = size(FEVD,3);

% Define the nRows and nColumns for the subplots
nRow = round(sqrt(nvars));
nCol = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);


%% Plot 1: by shock
%================================================
for jj=1:nshocks
    for ii=1:nvars
% FEVD(t,j,k): matrix with 't' steps, the FEVD due to 'j' shock for 'k' variable
% 哪个是变量，哪个是 reaction？？
        subplot(nRow,nCol,ii);
        plot(steps,FEVD(:,jj,ii),'LineStyle','-','color',[0.01 0.09 0.44],'LineWidth',2);
        hold on
        plot(x_axis,'k','LineWidth',0.5)
        if isfield(VAR,'FEVDH') && isfield(VAR,'FEVDL')
            plot(steps,VAR.FEVDH(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
            hold on
            plot(steps,VAR.FEVDL(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
        end
        
        if isfield(VAR,'pos_shock')
            title([VAR.select_vars_label_short_order{ii} ' to ' VAR.select_vars_label_short_order{VAR.pos_shock}], 'FontWeight','bold','FontSize',7); 
        else
            title([VAR.select_vars_label_short_order{ii}], 'FontWeight','bold','FontSize',7); 
        end
        
        
        xlim([1 nsteps]); ylim([0 1]);
       
    end
end

save2pdf(['newfigures/FEVD_hor' num2str(VAR.irhor)])
% close
end

