%% Mod_EstOLS Function
    % Goal:
        % Use OLS to estimate the VAR model and provide VAR results.
        % If specified, run short run restriction and provide SVAR results.
    % Input:
        % Input 1: Model setup from Mod_setup.
        % Input 2: Order of variables. Write [] if doesn't want to identify.
        % Input 3: Position of shock. Write [] if want all shocks.
        % Input 4: Identification scheme: 'short','long', 'non'.
    % Output:
        % Output 1: identified B matrix (if specify short run restriction)(in VAR.B).
        % Output 2: IRF (in VAR.irs if specify one shock; in VAR.irs_total if not).
        % Output 3: FEVD (in VAR.FEVD if specify one shock; in VAR.FEVD_total if not).
%%

function [VAR] = Mod_EstOLS_restricted(VAR, var_order, restr, chi2, restrF, pos_shock, ident)

%% 1.Organize variables.
    if isempty(var_order)
        VAR.var_order       = cumsum(ones(1,size(VAR.vars, 2)));
                
        VAR.vars_order      = VAR.vars;
        VAR.select_vars_label_order = VAR.select_vars_label;
        VAR.select_vars_label_short_order = VAR.select_vars_label_short;
        VAR.select_vars_label_short_order = VAR.select_vars_label_short;
       
    else
        VAR.var_order       = var_order;
        VAR.vars_order      = VAR.vars(:,VAR.var_order); 
        VAR.select_vars_label_order = VAR.select_vars_label(:,VAR.var_order);
        VAR.select_vars_label_short_order = VAR.select_vars_label_short(:,VAR.var_order);
    end
    
    VAR.vars_order(isnan(VAR.vars_order)) = 0;
%%    

%% 2.Restricted estimation.

% Set data
    T = VAR.T; % length of time series
    k = VAR.n; % number of variables
    p = VAR.p; % number of lags
    cons = VAR.const;
    Yobs = VAR.vars_order';
    
    VAR.restr  = restr;
    VAR.chi2   = chi2;
    VAR.restrF = restrF;
    
    Y = zeros(k,T);
    if cons == 1      % if constant is included
        Z = zeros(k*p+1,T);
    else              % otherwise
        Z = zeros(k*p,T);
    end
    
    for j = 1:T
        Y(:,j) = Yobs(:,j+p); % set the regressand matrix for each period j
        dummy  = Yobs(:,(j+p-1):(-1):(j)); % set the regressor matrix for the corresponding period j
        dummy  = dummy(:); % line it up
        if cons == 1  % if constant is included
            Z(:,j) = [1;dummy];
        else          % otherwise
            Z(:,j) = dummy;
        end
    end
        
% Set estimated coefficients
    VAR.Bet = zeros(k,k,p); % estimated coefficient (withVAR constant)
    if cons == 1
        VAR.BCONS = zeros(k,1); % estimated coefficient on a constant
    end

% Do regression line by line, incorporate restrictions
    if cons==1
        dummy=zeros(k,k*p+1); 
    else
        dummy=zeros(k,k*p); 
    end
    
    dummyalt=dummy;
    restr_big=0*dummyalt;
    
%     if ~exist('restrF','var')
%         restrF=0*restr; 
%     end
%     keyboard
    
    for rownum=1:k
        off  = find(restr(rownum,:) == 0);
        on   = find(restr(rownum,:) == 1);
        offF = find(restrF(rownum,:) == 0);
        onF  = find(restrF(rownum,:) == 1);
        onC  = find(restr(rownum,:) == 1 & restrF(rownum,:) == 0);
        % we want to do a kronecker sum, so have to do the product and take
        % logs. 
%         select=log2(kron(2.^([0:p-1]*k),2.^off));
        if cons==1
            Z1=[Z(1,:); Z(1+off,:)];
            dummy(rownum,[1 off+1]) = Y(rownum,:)*Z1'/(Z1*Z1');
        else
            Z1=[Z(off,:)];
            dummy(rownum,off) = Y(rownum,:)*Z1'/(Z1*Z1');
        end
        
        if chi2==1
            
            if cons==1
                Z1=[Z(1,:); Z(1+offF,:)];
                dummyalt(rownum,[1 offF+1]) = Y(rownum,:)*Z1'/(Z1*Z1');
                restr_big(rownum,onC+1)=1;
            else
                Z1=[Z(offF,:)];
                dummyalt(rownum,[offF]) = Y(rownum,:)*Z1'/(Z1*Z1');
                restr_big(rownum,onC)=1;
            end
        end
    end
    
    if cons ~= 1
        for kk=1:p
            VAR.Bet(:,:,kk) = dummy(:,((kk-1)*k+1):(kk*k));
        end
    end

    if cons == 1
        for kk=1:p
            VAR.Bet(:,:,kk) = dummy(:,(1+(kk-1)*k+1):(1+kk*k));
        end
        VAR.BCONS = dummy(:,1);
    end
    
% Get estimated residuals and covariance matrix.
    VAR.res   = (Y-dummy*Z)'; % residuals
%     VAR.Sigma = (1/(T-k*p-1))*(VAR.res*VAR.res'); % covariance matrix
        % BDG code version (wrong?):     
        VAR.Sigma = (1/(T+p-k*p-1))*(VAR.res'*VAR.res); % covariance matrix
    VAR.bet   = dummy';
    VAR.betT  = dummy;


    % get companion matrix.
    VAR.betcomp = [VAR.betT(:,1+VAR.const:VAR.n*VAR.p+VAR.const); eye(VAR.n*(VAR.p-1)) zeros(VAR.n*(VAR.p-1),VAR.n)];  

    % get moving average coefficients.
        % OD:   order of the MA representation to be computed
    if isfield(VAR,'OD')
        OD = VAR.OD;
    else
        OD = 1e2;
    end

    VAR.BB = zeros(VAR.n,VAR.n,OD);
    VAR.BB(:,:,1) = eye(VAR.n,VAR.n);
    coef = VAR.betT(:,1+VAR.const:VAR.n*VAR.p+VAR.const);
    COEF = zeros(VAR.n,VAR.n, VAR.p);

for kk=1:VAR.p
        COEF(:,:,kk) = coef(:,((kk-1)*VAR.n+1):(kk*VAR.n));
end

for ii = 1:OD
    dummy = zeros(VAR.n,VAR.n);
    for jj=1:min(ii,VAR.p) 
        % previous p block. If ii < p, previous ii block.
        % min to make sure ii-jj+1 >= 1.
        % written as sum of jj matrix
        dummy = dummy + VAR.BB(:,:,ii-jj+1)*COEF(:,:,jj);
    end
    VAR.BB(:,:,ii+1) = dummy;
end

% Recursive structure: jj = 1 ~ p (p=8)
% ii - 1 + 1
% ii - 1
% ii - 2
% ii - 3
% ii - 4
% ii - 5
% ii - 6
% ii - 7



%%

%% 3.Run identification.
    if isempty(var_order)
        % Do not impose short run restriction. Run VAR.
        VAR.B   =   eye(size(VAR.Sigma));   
        VAR.ident = 'non';

    else
        if (strcmp(ident, 'short'))
            VAR.ident = 'short';
            VAR.B     =   chol(VAR.Sigma,'lower');
            
        elseif (strcmp(ident, 'long'))            
            betinf_big = inv(eye(length(VAR.betcomp))-VAR.betcomp); % From the companion.
            Finf = betinf_big(1:VAR.n,1:VAR.n);
            b  = chol(Finf*VAR.Sigma*Finf')'; % Identification: u2 has no effect on y1 in the long run.
            B  = Finf\b;   
            VAR.ident = 'long';
            VAR.B     = B;
            
        elseif (strcmp(ident, 'non')) 
            % disp('No identification imposed. Run VAR.');       
            VAR.ident = 'non';
            VAR.B   =   eye(size(VAR.Sigma));
            
        else
            
            disp('Not applicable. Use Mod_Ident function for other identification schemes.');
        end
    end
%%

%% 4.Obtain result: IRF and FEVD.
    if ~isempty(pos_shock) && pos_shock > size(VAR.Sigma, 2)
        disp('Shock position exceeds number of variables.');
   
    elseif ~isempty(pos_shock) && pos_shock <= size(VAR.Sigma, 2) &&  pos_shock > 0 
        VAR.pos_shock = pos_shock;
        % initial shock: eps(1,1)=1
        
        IRF      = Res_IRF(VAR, VAR.B(:,VAR.pos_shock));
        VAR.IRF  = IRF;
        FEVD     = Res_FEVD(VAR, VAR.B(:,VAR.pos_shock));
        VAR.FEVD = FEVD;
            
        % Calculate irs directly.     
%             irs(VAR.p+1,:) = VAR.B(:,VAR.pos_shock); 
%                 % Automatically create a matrix called 'irs', includes p+1 rows, 
%                 % with the p+1 row being the pos_shock row of VAR.B.
%                 % Therefore, ncol of irs equals nrows of VAR.B: irs(nsteps, nvar).
% 
%                 for jj=2:VAR.irhor    % Row p+1 has already been filled£¡
%                 lvars = (irs(VAR.p+jj-1:-1:jj,:))'; 
%                     % For period jj, the P periods before jjth period are all useful.
%                     % Take the last p rows before jjth row.
%                         % Transpose it in order to be multiplied by the beta matrix,
%                         % because beta matrix is ordered backward.
% 
%                 irs(VAR.p+jj,:) = lvars(:)'*VAR.bet(1+VAR.const:VAR.n*VAR.p+VAR.const,:);  
%                     % Take VAR lvars matrix and extent it as a vector so as to multiplied.
%                     % it by the beta matrix. Each column of the beta matrix correspond
%                     % to one variable of y.
%                     % VAR.bet(1:VAR.p*VAR.n,:) is a np x n vector, omitting the constant coefficients (why?)
%                         % Why doesn't impulse response consider constants and trends?
%                         % Because constant and trends doesn't change?
%                 end
%                 % After the loop, irs is a matrix with 'p+VAR.irhor' rows, 
%                 % containing the IRF of 'j' variable to 'pos_shock'th shock. 
%             VAR.irs = irs(VAR.p+1:end,:); 
%             % Only keep the rows that are non zero.

    elseif isempty(pos_shock)
    % Shock not specified. Want all shocks.
    
        % 1) Calculate using function.
            IRF_total      = Res_IRF(VAR, VAR.B);
            VAR.IRF_total  = IRF_total;
            FEVD_total     = Res_FEVD(VAR, VAR.B); 
            VAR.FEVD_total = FEVD_total;
        
        % 2) Calculate it directly.     
%             % Build three dimensional cube to store reactions to all shocks.
%             % (t,j,k): matrix with 't' steps, containing the IRF of 'j' variable to 'k' shock.
%                 % We can write two loops:
%                 irs_total   = zeros(VAR.irhor+VAR.p,length(VAR.B),length(VAR.B));  
%                 % Allocate space in advance.
% 
%                 for ii=1:length(VAR.B)
%                     % loop through shocks (# number of shocks equals # of variables)
%                     irs_total(VAR.p+1,:,ii) = VAR.B(:,ii); 
%                     for jj=2:VAR.irhor
%                        lvars_total = (irs_total(jj+VAR.p-1:-1:jj,:,ii))';
%                        irs_total(VAR.p+jj,:,ii) = lvars_total(:)'*VAR.bet(1+VAR.const:VAR.n*VAR.p+VAR.const,:);   
%                     end
%                 end
%                 VAR.irs_total = irs_total(VAR.p+1:end,:,:);

    end
