function Plot_IRF_Comparison(VAR,VAR2,VARbs,VARbs2,para)

nCol     = round(sqrt(VAR.n));
nRow     = ceil(sqrt(VAR.n));
fontsize = para.fontsize;            % Fontsize in figures

set(gcf,'DefaultAxesFontSize',fontsize);
set(gcf,'DefaultTextFontSize',fontsize);

for nvar=1:VAR.n
    subplot(nCol,nRow,nvar); plot(1:VAR.irhor,VAR.IRF(:,nvar),'LineWidth',2,'Color',[0.01 0.09 0.44]); 
    title(VAR.select_vars_label_order(1,nvar)); set(gcf, 'Color', 'r');
    hold on
    	
    if isfield(VARbs,'IRFH') && isfield(VARbs,'IRFL')
    subplot(nCol,nRow,nvar); plot(1:VAR.irhor,VARbs.IRFH(:,nvar),'LineWidth',1,'Color',[0.5098 0.4863 0.4863],'LineStyle','--');
    hold on
    subplot(nCol,nRow,nvar); plot(VARbs.IRFL(:,nvar),'LineWidth',1,'Color',[0.5098 0.4863 0.4863],'LineStyle','--');
    hold on
    end
    
    subplot(nCol,nRow,nvar); plot(1:VAR2.irhor,VAR2.IRF(:,nvar),'LineWidth',2,'Color',[0.39 0.58 0.93]); 
    title(VAR2.select_vars_label_order(1,nvar)); set(gcf, 'Color', 'w');
    hold on
    
    if isfield(VARbs2,'IRFH') && isfield(VARbs2,'IRFL')
    subplot(nCol,nRow,nvar); plot(1:VAR2.irhor,VARbs2.IRFH(:,nvar),'LineWidth',1,'Color',[0.0745 0.6235 1.0000],'LineStyle','--');
    hold on
    subplot(nCol,nRow,nvar); plot(VARbs2.IRFL(:,nvar),'LineWidth',1,'Color',[0.0745 0.6235 1.0000],'LineStyle','--');
    hold on
    end

    set(gca,'YGrid','off','XGrid','on');    
    set(gcf, 'units', 'inches', 'position', [2 2 15 8])
    set(gcf, 'Color', 'w');
    grid
end

save2pdf(['newfigures/IRF_Comp_hor' num2str(VAR.irhor)])


