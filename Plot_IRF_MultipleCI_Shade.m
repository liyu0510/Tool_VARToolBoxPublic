function Plot_IRF_MultipleCI_Shade(VAR, VARCI, para)

%% 1. Parameter
nCol     = round(sqrt(VAR.n));
nRow     = ceil(sqrt(VAR.n));
fontsize = para.fontsize;            % Fontsize in figures

set(gcf,'DefaultAxesFontSize',fontsize);
set(gcf,'DefaultTextFontSize',fontsize);

%% 2. Plot
for nvar=1:VAR.n
       
    subplot(nCol,nRow,nvar); 
    plot(1:VAR.irhor,VAR.IRF(:,nvar),'LineWidth',2,'Color',[0.01 0.09 0.44]); 
%     title(VAR.select_vars_label_order(1,nvar),'interpreter','latex','FontSize', 13);
%     title(['Resp. of ' char(VAR.select_vars_label_order(1,nvar))],'interpreter','latex','FontSize', 13); 
    title([char(VAR.select_vars_label_order(1,nvar))],'interpreter','latex','FontSize', 13); 

    hold on

    if isfield(VARCI,'A') 
        
        A = VARCI.A;
        
        ciplot(A(:,nvar,1),A(:,nvar,4),1:VAR.irhor,[0.5 0.5 0.5],0.2);
        hold on
        
        plot(A(:,nvar,1),'LineStyle','--','Color',[0.5098 0.4863 0.4863],'LineWidth',1);
        plot(A(:,nvar,2),'LineStyle','-','Color',[0.5098 0.4863 0.4863],'LineWidth',1);
        plot(A(:,nvar,3),'LineStyle','-','Color',[0.5098 0.4863 0.4863],'LineWidth',1);
        plot(A(:,nvar,4),'LineStyle','--','Color',[0.5098 0.4863 0.4863],'LineWidth',1);
        xlim([1 VAR.irhor]);
        plot(zeros(1,VAR.irhor),'k','LineWidth',0.5)
    end
    
	xlabel ('Periods after shock','interpreter','latex','FontSize', 13);
    ylabel ('Impulse responese ','interpreter','latex','FontSize', 13);
%     ylabel (['Resp. of '  char(VAR.select_vars_label_order(1,nvar))],'interpreter','latex','FontSize', 13);
    
    set(gca,'YGrid','on','XGrid','on');
    set(gcf, 'units', 'inches', 'position', [.1 .1 10 8])
    set(gcf, 'Color', 'w');


end

save2pdf(['newfigures/IRF_Multi_Shade_hor' num2str(VAR.irhor)])
% close
