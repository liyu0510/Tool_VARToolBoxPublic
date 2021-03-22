function Plot_IRF_Comparison_Shade(VAR, VAR2, VARbs, VARbs2, para)

nCol     = round(sqrt(VAR.n));
nRow     = ceil(sqrt(VAR.n));
fontsize = para.fontsize;            % Fontsize in figures

set(gcf,'DefaultAxesFontSize',fontsize);
set(gcf,'DefaultTextFontSize',fontsize);

for nvar=1:VAR.n
    
	subplot(nCol,nRow,nvar); 
    plot(1:VAR.irhor,VAR.IRF(:,nvar),'LineWidth',2,'Color',[0, 0.4470, 0.7410]); 
    title(VAR.select_vars_label_order(1,nvar)); set(gcf, 'Color', 'w');
    hold on
    
    if isfield(VARbs,'A') 
        
        A = VARbs.A;
        subplot(nCol,nRow,nvar);
        plot(A(:,nvar,1),'LineStyle','--','Color',[0.3010, 0.7450, 0.9330],'LineWidth',1);
        plot(A(:,nvar,2),'LineStyle','-','Color',[0.3010, 0.7450, 0.9330],'LineWidth',1);
        plot(A(:,nvar,3),'LineStyle','-','Color',[0.3010, 0.7450, 0.9330],'LineWidth',1);
        plot(A(:,nvar,4),'LineStyle','--','Color',[0.3010, 0.7450, 0.9330],'LineWidth',1);
        xlim([1 VAR.irhor]);
        plot(zeros(1,VAR.irhor),'k','LineWidth',0.5)
    
    elseif isfield(VARbs,'IRFH') && isfield(VARbs,'IRFL')        
        ciplot(VARbs.IRFL(:,nvar),VARbs.IRFH(:,nvar),1:VAR.irhor,[0.5 0.5 0.5],0.2);
        hold on    

        subplot(nCol,nRow,nvar); plot(1:VAR.irhor,VARbs.IRFH(:,nvar),'LineWidth',1,'Color',[0.3010, 0.7450, 0.9330],'LineStyle','--');
        hold on
        subplot(nCol,nRow,nvar); plot(VARbs.IRFL(:,nvar),'LineWidth',1,'Color',[0.3010, 0.7450, 0.9330],'LineStyle','--');
        hold on
    
    end
    
    subplot(nCol,nRow,nvar);

    plot(1:VAR2.irhor,VAR2.IRF(:,nvar),'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
    title(VAR2.select_vars_label_order(1,nvar)); set(gcf, 'Color', 'w');
    hold on
    
    if isfield(VARbs2,'A') 
        
        A = VARbs2.A;
        
        plot(A(:,nvar,1),'LineStyle','--','Color',[0.9290, 0.6940, 0.1250]	,'LineWidth',1);
        plot(A(:,nvar,2),'LineStyle','-','Color',[0.9290, 0.6940, 0.1250]	,'LineWidth',1);
        plot(A(:,nvar,3),'LineStyle','-','Color',[0.9290, 0.6940, 0.1250]	,'LineWidth',1);
        plot(A(:,nvar,4),'LineStyle','--','Color',[0.9290, 0.6940, 0.1250]	,'LineWidth',1);
        xlim([1 VAR.irhor]);
        plot(zeros(1,VAR.irhor),'k','LineWidth',0.5)
    
    elseif isfield(VARbs2,'IRFH') && isfield(VARbs2,'IRFL')
            
        subplot(nCol,nRow,nvar); plot(1:VAR2.irhor,VARbs2.IRFH(:,nvar),'LineWidth',1,'Color',[0.9290, 0.6940, 0.1250]	,'LineStyle','--');
        hold on
        subplot(nCol,nRow,nvar); plot(VARbs2.IRFL(:,nvar),'LineWidth',1,'Color',[0.9290, 0.6940, 0.1250]	,'LineStyle','--');
        hold on
    
    end 

    set(gca,'YGrid','off','XGrid','on');    
    set(gcf, 'units', 'inches', 'position', [2 2 15 8])
    set(gcf, 'Color', 'w');
    grid
end

save2pdf(['newfigures/IRF_Comp_Shade_hor' num2str(VAR.irhor)])


