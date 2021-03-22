%% Res_FEVD Function
    % Goal:
        % Given identified result, compute FEVD.
    % Input:
        % Input 1: VAR estimation results.
        % Input 2: Identified columns (one or multiple) of B matrix.
    % Output:
        % Output 1: FEVD(t,j,i): matrix with 't' steps, the FEVD due to 'j' shock for 'i' variable.
        % Output 2: VAR struct with FEVD.

%% 

function [FEVD, VAR] = Res_FEVD(VAR, Ident_Columns)

    ncols = size(Ident_Columns,2);
    % Number of shocks identified.
        
                 
       %%
        % 2.1.Extract variables from struct.
            nsteps = VAR.irhor ;          % Number of steps for IRFs and FEVDs.
            nvar   = VAR.n;               % Number of variables.
            sigma  = VAR.Sigma;           % Adjusted var-covar of Cholesky.
            % sigma  = VAR.Sigma_m;       % Adjusted var-covar of Proxy,with limited sample.
   
       %%
        % 2.2.Define the matrix to be filled.
            MSE   = zeros(nvar,nvar,nsteps);     % MSE(i,j,t) eaction of i to j.
            SE    = zeros(nsteps,nvar);          % Standard error, only the diagonals pf MSE.
            MSE_i = zeros(nvar,nvar,nsteps);     % Mean squared error caused by variable j.
            PSI   = zeros(nvar,nvar,nsteps);     % Multiplication matrix of shock.
        
       %%
        % 2.3.Calculate MSE.
 
            % Compute MSE multiplier
                % The calculation of MSE does not require identification.
                % We need VAR impulse response.
                % Why? 研究 shock 是如何被放大的。
                % 为什么研究shock的放大需要PSI？不是应该需要 companion matrix？
                % PSI 本质上就是 companion matrix？
                % PSI 是从哪个 code 里来的？Bianchi..
            
            
                irs_total_unitshock   = zeros(VAR.irhor+VAR.p,length(VAR.B),length(VAR.B));  
                % irs_total_unitshock(t,j,i):t period, j impulse response from i.
                % Allocate space in advance.
                
                pseudoB = eye(size(VAR.B));
                % We take pseudoB to be identity matrix because we want pure multiplicative effect from the coefficients.
                
                for i = 1:length(pseudoB)
                    % loop through shocks (# number of shocks equals # of variables)
                    % j is the jth shock
                    irs_total_unitshock(VAR.p+1,:,i) = pseudoB(:,i); 
                        % jj is the jjth period reaction to that shock
                        for jj=2:VAR.irhor
                           lvars_total = (irs_total_unitshock(jj+VAR.p-1:-1:jj,:,i))';
                           irs_total_unitshock(VAR.p+jj,:,i) = lvars_total(:)'*VAR.bet(1+VAR.const:VAR.n*VAR.p+VAR.const,:);   
                        end
                end

                IRFjunk = irs_total_unitshock(VAR.p+1:end,:,:);   % VAR impulse response multiplier.
                    % IRFjunk(t,j,i):t period,j impulse response from i.
                    % What's the difference btw IRFjunk and VAR.betcomp?
                    % IRFjunk is computed IRF.
                    % But isn't IRF just multiplication of VAR.betcomp?
                    
                for c = 1:size(VAR.B,1)    
                  % IRFjunk(t,j,i):t period,j impulse response from i.
                  % PSI(j,i,t): j impulse response from i, t period.
                    PSI(:,c,:) = reshape(IRFjunk(:,:,c)',1,nvar,nsteps);  % Note that IRFjunk is transposed.
                        % redefine the sidze：
                            % IRFjunk(:,:,mm) is nstep x nvar matrix;
                            % IRFjunk(:,:,mm)'is nvar x nstep matrix;
                        % after reshape, PSI(:,mm,:)is also nvar x nstep matrix;    
                            % ready to be filled in.
                        % What does 1 stand for？
                     % Recall
                        % PSI   = zeros(nvar,nvar,nsteps); 
                        % IRFjunk  = zeros(nsteps,nvar,nvar); 
                        % IRFjunk(t,j,k): matrix with 't' steps, containing the IRF of 'j' variable to 'k' shock;
                        % PSI(j, k, t): matrix with t steps, containing the multiplier of j variable to 'k' shock;
                        % Save every impulse response to shock k.
                end
            
            % Calculate total MSE.
                MSE(:,:,1) = sigma; 
                % MSE(j,i,t): j impulse response to i, period t.
                    for kk = 2:nsteps
                       MSE(:,:,kk) = MSE(:,:,kk-1) + PSI(:,:,kk)*sigma*PSI(:,:,kk)';
                    % PSI(:, :, kk): the multiplier of all variables react to all shocks at period kk.
                    end

            
       %%
        % 2.4.Compute MSE due to particular shocks (i.e, FEVD).   
        FEVD  = zeros(nsteps,ncols,nvar);     
        % FEVD(t,i,j):period t, the MSE of variable i, explained by variable j.
        % FEVD
            % The number of rows is the number of forecast periods.
            % The number of columns identified should be equal to the columns of FEVD.
            % The number of pages is the number of reactions
            
        for i = 1:ncols  % loop for the shocks
            column = Ident_Columns(:,i); 
                
            MSE_i(:,:,1) = column*column'; 
            % MSE_i(j,i,t): j impulse response to i, period t.
                % Why is it a matrix? Because it is the VAR-coVAR induced by variable j's changes
            for kk = 2:nsteps
                MSE_i(:,:,kk) = MSE_i(:,:,kk-1) + PSI(:,:,kk)*(column*column')*PSI(:,:,kk)';   
            end
            
            % Compute the Forecast Error Covariance Decomposition
            FECD = MSE_i./MSE; 
            % FECD(j,i,t): j to i, period t.
            % total var-covar below, var-covar due to j above
            % 点除：两个立方体一一对应的input挨个相除？？
    
       %%     
        % 2.5.Obtain cumulative FEVD share of desired variable.        
        % Select only the variance terms.               
            for t = 1:nsteps
                for j = 1:nvar
                    FEVD(t,i,j) = FECD(j,j,t);  
                    % FECD(j,i,t): j to i, period t.
                    % FEVD(t,i,j): j to i, period t
                    % 在第 t 期，第 j 个变量被第 i 个变量影响
                    SE(t,:) = sqrt(diag(MSE(:,:,t))');  % 为什么必须写在循环里面？
                end
            end
            
        end
        
    VAR.FEVD = FEVD;
    
end

