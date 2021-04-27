function Plot_IRF(VAR, VARCI, para)

%% 1.Import parameter
nCol     = round(sqrt(VAR.n));
nRow     = ceil(sqrt(VAR.n));
fontsize = para.fontsize;            % Fontsize in figures


set(gcf,'DefaultAxesFontSize',fontsize);
set(gcf,'DefaultTextFontSize',fontsize);

for nvar=1:VAR.n
    
    subplot(nCol,nRow,nvar); 
    plot(1:VAR.irhor,VAR.IRF(:,nvar),'LineWidth',2,'Color',[0.01 0.09 0.44]); 
    title(VAR.select_vars_label_order(1,nvar)); set(gcf, 'Color', 'w');
    hold on
    
	if isfield(VARCI,'IRFH') && isfield(VARCI,'IRFL')
        subplot(nCol,nRow,nvar); plot(1:VAR.irhor,VARCI.IRFH(:,nvar),'LineWidth',1,'Color',[0.5098 0.4863 0.4863],'LineStyle','--','HandleVisibility','off');
        hold on
        subplot(nCol,nRow,nvar); plot(1:VAR.irhor,VARCI.IRFL(:,nvar),'LineWidth',1,'Color',[0.5098 0.4863 0.4863],'LineStyle','--');
        hold on
    end
    set(gca,'YGrid','off','XGrid','on');
    grid

end

save2pdf(['newfigures/IRF_hor' num2str(VAR.irhor)])
close

