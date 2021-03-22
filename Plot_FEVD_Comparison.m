function Plot_FEVD_Comparison(VAR, VAR2)
%% 1.Import parameter

FEVD  = VAR.FEVD;
FEVD2 = VAR2.FEVD;

% Initialize FEVD matrix
[nsteps]  = size(FEVD,1);
[nshocks] = size(FEVD,2);   % Because proxy VAR only identified one shock
[nvars]   = size(FEVD,3);

% Define the rows and columns for the subplots
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);


%% Plot 1: by shock
%================================================
for jj=1:nshocks 
    for ii=1:nvars
        subplot(row,col,ii);
        
        plot(steps,FEVD(:,jj,ii),'LineStyle','-','Color',[0.01 0.09 0.44],'LineWidth',2);
        hold on
        
        if isfield(VAR,'FEVDH') && isfield(VAR,'FEVDL')
            plot(steps,VAR.FEVDH(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
            hold on
            plot(steps,VAR.FEVDL(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
        end
        
        plot(steps,FEVD2(:,jj,ii),'LineStyle','-','Color',[0.39 0.58 0.93],'LineWidth',2);
        hold on
        
        plot(x_axis,'k','LineWidth',0.5)
        if isfield(VAR2,'FEVDH') && isfield(VAR2,'FEVDL')
            plot(steps,VAR2.FEVDH(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
            hold on
            plot(steps,VAR2.FEVDL(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
        end
        
        if isfield(VAR,'pos_shock')
            title([VAR.select_vars_label_order{ii} ' to ' VAR.select_vars_label_order{VAR.pos_shock}], 'FontWeight','bold','FontSize',10); 
        else
            title([VAR.select_vars_label_order{ii}], 'FontWeight','bold','FontSize',10); 
        end
        
        xlim([1 nsteps]); ylim([0 1]);
       
    end
end
end
