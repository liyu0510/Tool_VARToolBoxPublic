%% Ident_Maxshare Function
    % Goal:
        % Use max share identification method to identify the VAR model
    % Input(varagin):
        % input(must): 
            % Input 1: VAR with short run identification results.
            % Input 2: Position of variables to explain.
            % Input 3: The horizon to explain (end).
        % input(optional):
            % Input 4: The horizon to explain (beginning) (default to be one).
            % Input 5: Starting vector (default to be ones).
            % Input 6: Maximum number of iteration (default to be 10 000).
    % Output:
            % Output 1: identification: one column of identified B matrix.
            % Output 2: IRF.
            % Output 3: FEVD.
%% 

function [VAR] = Ident_Maxshare(varargin)
	
    %% 1.Set optimization parameters.
        
        % Extract struct from input.
        VAR         = varargin{1};      % VAR after estimation and short run identification.
        VAR.N       = varargin{2};      % The position of variable whose variation you want to explain.
        VAR.k       = varargin{3};      % The number of horizons to explain.    

        if nargin < 3
            disp('Not enough inputs.');
            
            elseif nargin == 3
            VAR.k0        = 1;                         % The beginning number of horizons to explain is default to be 1
            VAR.q0        = ones(VAR.n,1)/sqrt(VAR.n); % The beginning vector.
            VAR.m         = 1500;                      % Maximum number of iteration.
            VAR.NofInputs = 3;
            
            elseif nargin == 4
            VAR.k0        = varargin{4};  
            VAR.q0        = ones(VAR.n,1)/sqrt(VAR.n);
            VAR.m         = 1500;  
            VAR.NofInputs = 4;
            
            elseif nargin == 5
            VAR.k0        = varargin{4};  
            VAR.q0        = varargin{5};
            VAR.m         = 1500;  
            VAR.NofInputs = 5;
        
            elseif nargin == 6
            VAR.k0        = varargin{4};  
            VAR.q0        = varargin{5};
            VAR.m         = varargin{6};      
            VAR.NofInputs = 6;
            

            
            else
            disp ('Too many inputs.')  
        end
 
        
    %% 2.Run maximization problem.
    
	options =  optimoptions(@fmincon,'Algorithm','sqp' ...
        , 'MaxFunctionEvaluations', 1e7 ...
        , 'MaxIterations', VAR.m ...
        , 'StepTolerance', 1.0000e-08 ...
        , 'ConstraintTolerance', 1.0000e-06);      

       %%
        % 2.1. First way: using external function Res_FEVD.
        
%         function [FEVD_total] = Res_FEVD_sum(VAR, q)
%             Ident_Column_q = VAR.B*q; 
%             [FEVD_q,~] = Res_FEVD(VAR, Ident_Column_q);
% 
%             % Compute cumulative sum.
%             FEVD_interest = FEVD_q(:,:,VAR.N);   
%                 % We are interested in the FEVD of the [Nth] variable.
%             FEVD_total =  sum(FEVD_interest((1:VAR.k),(1:size(Ident_Column_q, 2))),'all');     
%                 % We care about the total FEVD of first [k] periods.
%         end
        
        VAR.maxshare = fmincon(@(q) -Res_FEVD_sum (VAR, q),VAR.q0,[],[],[],[],[],[], @confuneq_general,options);

       %%
        % 2.2. Second way: define function directly.
        
%         function [FEVD_total] = FEVDshare(VAR, q)
% 
% 
%             % 2.2.1.Extract variables from struct.
%                 nsteps = VAR.irhor ;          % Number of steps for IRFs and FEVDs.
%                 nvar   = VAR.n;               % Number of variables.
%                 sigma  = VAR.Sigma;           % Adjusted var-covar of Cholesky.
%                     % sigma  = VAR.Sigma_m;       % Adjusted var-covar of Proxy,with limited sample.
%                 N      = VAR.N;               % The position of variable whose variation you want to explain.
%                 k      = VAR.k;               % The number of horizons to explain.
% 
%             % 2.2.2.Define the matrix to be filled.
%                 MSE   = zeros(nvar,nvar,nsteps);     % Mean sqaured error that stores VAR-COVAR matrix.
%                 SE    = zeros(nsteps,nvar);          % Standard error, only the diagonals pf MSE.
%                 MSE_j = zeros(nvar,nvar,nsteps);     % Mean squared error caused by variable j.
%                 PSI   = zeros(nvar,nvar,nsteps);     % Multiplication matrix of shock.
%                 FEVD  = zeros(nsteps,nvar,nvar);     % At period t, the MSE of variable i, explained by variable j.
% 
%             % 2.2.3.Calculate MSE.
% 
%             % The calculation of MSE does not require identification.
%             % Compute MSE multiplier
%                 IRFjunk = VAR.irs_total_unitshock;   % VAR impulse response.
%                 for mm = 1:nvar
%                     PSI(:,mm,:) = reshape(IRFjunk(:,:,mm)',1,nvar,nsteps);  % Note that IRFjunk is transposed.
%                         % redefine the sidze：
%                             % IRFjunk(:,:,mm) is nstep x nvar matrix;
%                             % IRFjunk(:,:,mm)'is nvar x nstep matrix;
%                         % after reshape, PSI(:,mm,:)is also nvar x nstep matrix;    
%                             % ready to be filled in.
%                         % What does 1 stand for？
%                      % Recall
%                         % PSI   = zeros(nvar,nvar,nsteps); 
%                         % IRFjunk  = zeros(nsteps,nvar,nvar); 
%                         % IRFjunk(t,j,k): matrix with 't' steps, containing the IRF of 'j' variable to 'k' shock;
%                         % PSI(j, k, t): matrix with t steps, containing the multiplier of j variable to 'k' shock;
%                         % Save every impulse response to shock k.
%                 end
% 
%             % Calculate total MSE.
%                 MSE(:,:,1) = sigma; 
%                 for kk = 2:nsteps
%                    MSE(:,:,kk) = MSE(:,:,kk-1) + PSI(:,:,kk)*sigma*PSI(:,:,kk)';
%                 % PSI(:, :, kk): the multiplier of all variables react to all shocks at period kk.
%                 end
%                 % Get the column of invA corresponding to the mm_th shock.
% 
%              % 2.2.4.Compute MSE due to particular shocks (i.e, FEVD).
% 
%                 % Number of shocks you want to identify:
%                 mm = 1; 
%                     % At this stage, we only identify one shock. So mm is set to 1.
%                     % If we want to identify more than one shocks:
%                         % q should be a multiple-column vector, part of a orthonormal matrix;
%                         % loop through mm to fill FEVD.
%                 invA = VAR.B;     % A matrix is the inverse of B matrix
%                 column = invA*q; 
%                     % invA is nxn matrix，q is nx1 vector.
%                     % This correspond to the column in matrix B that amplifies jth shock
%                     % The actual position of j doesn't matter.
%                     % Why? because this column can be on any position.
% 
%                 MSE_j(:,:,1) = column*column'; 
%                     % Why is it a matrix? Because it is the VAR-coVAR induced by variable j's changes
%                 for kk = 2:nsteps
%                     MSE_j(:,:,kk) = MSE_j(:,:,kk-1) + PSI(:,:,kk)*(column*column')*PSI(:,:,kk)';   
%                 end
% 
%                 % Compute the Forecast Error Covariance Decomposition
%                 FECD = MSE_j./MSE; % total covariance below, covariance due to j above
%                 % 点除：两个立方体一一对应的input挨个相除？？
%  
%             % 2.2.5.Obtain cumulative FEVD share of desired variable.
% 
%             % Select only the variance terms.               
%                 for nn = 1:nsteps
%                     for ii = 1:nvar
%                         FEVD(nn,mm,ii) = FECD(ii,ii,nn);  % 在第 n 期，第 i 个变量被第 m 个变量影响
%                         SE(nn,:) = sqrt(diag(MSE(:,:,nn))');  % 为什么必须写在循环里面？
%                     end
%                 end
%                 VAR.FEVD = FEVD;
%             % Compute cumulative sum.
%                 FEVD_interest = FEVD(:,:,N);                   % We are interested in the FEVD of the [Nth] variable
%                 FEVD_total =  sum(FEVD_interest((1:k),(1:size(column, 2))),'all');     % We care about the total FEVD of first [k] periods
%                 % Note that here it's just 1st column, because we only
%                 % have one shock.
%                 % However, if we have multiple shocks, then is it to
%                 % maximize the sum of all these shocks' reaction, or
%                 % separately adding up two columns?
%         end
%     
%         VAR.maxshare = fmincon(@(q) -FEVDshare (VAR, q),VAR.q0,[],[],[],[],[],[], @confuneq_general,options);
        
    %% 3.Obtain result: IRF and FEVD.
        VAR.Ident_column = VAR.B*VAR.maxshare;   
            % Derived the identified column: the multiplication of the 
            % identified orthonormal columns and the Cholesky matrix.
            
        % initial shock: eps(1,1)=1       
            IRF      = Res_IRF(VAR, VAR.Ident_column);
            VAR.IRF  = IRF;
            FEVD     = Res_FEVD(VAR, VAR.Ident_column);
            VAR.FEVD = FEVD;
	
end