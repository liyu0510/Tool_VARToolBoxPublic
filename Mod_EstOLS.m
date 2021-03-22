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

function [VAR] = Mod_EstOLS(VAR, var_order, pos_shock, ident)

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
    
%% 2.Make lag matrix.
X      = lagmatrix(VAR.vars_order,1:VAR.p);    % Note that here is 1:VAR.p
X      = X(VAR.p+1:end,:);               % The first VAR.p period is not complete, drop them
Y      = VAR.vars_order(VAR.p+1:end,:);
%%

%% 3.Add constant.
%   const : 0, no constant, no trend
%           1, constant, no trend
%           2, constant, trend
%           3, constant, trend, trend^2

    if VAR.const == 0
            X = X;    
    elseif VAR.const == 1                       % [constant, X]
            X = [ones(size(X,1),1) X];
    elseif VAR.const == 2                       % [time trend and constant, X]
            trend = 1:size(X,1);
            X = [ones(size(X,1),1) trend' X];
    elseif VAR.const == 3                       % [squared time trend, linear time trend, constant, X]
            trend = 1:size(X,1);
            X = [ones(size(X,1),1) trend' (trend').^2 X];
    end

VAR.X  = X;
VAR.Y  = Y;

VAR.sizeX = size(VAR.X);
VAR.sizeY = size(VAR.Y);
%%

%% 4.Run OLS estimation.
VAR.bet  = X\Y; 
% VAR.bet = (Y'*X/(X'*X))'; % estimated coefficient

VAR.betT    = VAR.bet';
VAR.res     = Y-X*VAR.bet;
VAR.Sigma   = (VAR.res'*VAR.res)/(VAR.T-(VAR.n*VAR.p + VAR.const));

% get different representation
    if VAR.const  ~= 1
        for kk=1:VAR.p
            VAR.Bet(:,:,kk) = VAR.betT(:,((kk-1)*VAR.n+1):(kk*VAR.n));
        end
    end

    if VAR.const  == 1
        for kk=1:VAR.p
            VAR.Bet(:,:,kk) = VAR.betT(:,(1+(kk-1)*VAR.n+1):(1+kk*VAR.n));
        end
        VAR.BCONS = VAR.betT(:,1);
    end

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

%% 5.Run identification.
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

%% 6.Obtain result: IRF and FEVD.
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
%                     % Take out lvars matrix and extent it as a vector so as to multiplied.
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
