function Plot_FEVD(VAR)

% FEVD(t,j,k): 
% Matrix with 't' steps, the FEVD due to 'j' shock for 'k' variable.
% Careful about who's shock variable and who's responding variable.

%% 1.Parameter

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


%% 2. Plot
for jj=1:nshocks
    for ii=1:nvars
        
        subplot(nRow,nCol,ii);
        plot(steps,FEVD(:,jj,ii),'LineStyle','-','color',[0.01 0.09 0.44],'LineWidth',2);
        
        hold on
        plot(x_axis,'k','LineWidth',0.5)
        
        if isfield(VAR,'FEVDH') && isfield(VAR,'FEVDL')
            plot(steps,VAR.FEVDH(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
            hold on
            plot(steps,VAR.FEVDL(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
        end
        
        title([VAR.select_vars_label_order{ii}], 'FontWeight','bold','FontSize',7,'interpreter','latex','FontSize', 13);
%         if isfield(VAR,'pos_shock')
%             title([VAR.select_vars_label_short_order{ii} ' to ' VAR.select_vars_label_short_order{VAR.pos_shock}],...
%                 'FontWeight','bold','FontSize',7,'interpreter','latex','FontSize', 13); 
%         else
%             title([VAR.select_vars_label_short_order{ii}], 'FontWeight','bold','FontSize',7,...
%                 'interpreter','latex','FontSize', 13);
%         end
        
        xlabel ('Periods after shock','interpreter','latex','FontSize', 13);
        
%         ylabel ('Forecast Error Variance Explained','interpreter','latex','FontSize', 13);
        if isfield(VAR,'pos_shock')
            ylabel (['Forecast Error Variance Explained by ' VAR.select_vars_label_short_order{VAR.pos_shock}],'interpreter','latex','FontSize', 13);
        else
            ylabel ('Forecast Error Variance Explained','interpreter','latex','FontSize', 13);

        end       
        
 
        set(gca,'YGrid','on','XGrid','on');
        set(gcf, 'units', 'inches', 'position', [.1 .1 10 8])
        set(gcf, 'Color', 'w');
        xlim([1 nsteps]); ylim([0 1]);
       
    end
end

save2pdf(['newfigures/FEVD_hor' num2str(VAR.irhor)])
% close
end

