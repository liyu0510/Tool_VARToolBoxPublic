%% Res_IRF Function
    % Goal:
        % Given identified result, compute IRF
    % Input:
        % Input 1: VAR estimation results.
        % Input 2: Identified columns (one or multiple) of B matrix.
    % Output:
        % Output 1: IRF(j,t,i): response of variable j to shock i at period t.
        % Output 2: VAR struct with IRF field.
        

%% 
function [IRF, VAR] = Res_IRF(VAR, IdentifiedB)

    ncols = size(IdentifiedB, 2);                           % Number of shocks identified.
    IRF   = zeros(VAR.irhor, length(VAR.select_vars), ncols); % Initial IRF.
    impulse = zeros(ncols,1);               % Initial shock.

    for c = 1:ncols
        
	% shock size
%         % one stdev shock
%             impulse(c,1) = 1;     
        % unitary shock (not applicable for proxy?)
        if isfield(VAR,'pos_shock')  && isfield(VAR,'Ident_column')==0     % short ident, long ident
            impulse(c,1) = 1/IdentifiedB(VAR.pos_shock); 
        elseif isfield(VAR,'pos_shock')&& isfield(VAR,'Ident_column')  % maxshare ident
            impulse(c,1) = 1/IdentifiedB(VAR.N); 
        elseif isfield(VAR,'b1')            % proxy ident
            impulse(c,1) = 1/IdentifiedB(1); 
        elseif isfield(VAR,'Bnew')          % sign ident
            impulse(c,1) = 1/IdentifiedB(c);       
        end
                    
        irs = zeros((VAR.irhor+VAR.p), length(VAR.select_vars));    
        irs(VAR.p+1,:) = IdentifiedB(:,c)*impulse(c,1); 
        
        for jj=2:VAR.irhor    % VAR.p+1 has been filled
            lvars = (irs(VAR.p+jj-1:-1:jj,:))'; % for a given period, choose previous p period values¡£
            irs(VAR.p+jj,:) = lvars(:)'*VAR.bet(1+VAR.const:VAR.n*VAR.p+VAR.const,:);  
        end
        
        IRF(:,:,c) = irs(VAR.p+1:end,:); 
    end
   
%         no_fig = plot_figure_onlyChol(VAR,VARChol,VARCholbs,4,2,M,switch_extern);
%         no_fig2 = plot_figure_onlyChol_nonConfi(VAR,VARChol,4,2,M,switch_extern);

    VAR.IRF = IRF;


end