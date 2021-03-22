%% Mod_EstOLS_bs Function
     % Goal:
        % Use wild bootstrap to obtain confidence interval for Cholesky estimation.
     % Input:
        % Input 1: Estimated VAR from Mod_EstOLS.
        % Input 2: Number of draws.
        % Input 3: Confidence level.
     % Output:
        % Output 1: VAR.bs
     % Comment:
        % Why in the Bootstrap results, the shock variable's IRF does not
        % fall into the CI?

function VARbs = Mod_EstOLS_bs(VAR, Para)
    
	%% 1. Specify shock position  
        if ~isfield(VAR,'pos_shock')
                disp ("Must specify a shock. Default to be the first.");
                pos_shock_bs = 1;
        else
                pos_shock_bs = VAR.pos_shock;
        end
        
        nboot  = Para.nboot;
        clevel = Para.clevel;
        CI     = Para.CI;

    %% 2. Perform wild bootstrap

       %%
        % 2.1. Removes a constant.
            res = VAR.res;
%             res = detrend(VAR.res,'constant');
                % Removes the mean value from each column of the matrix.
                % Why do you need to do this?


                
       %%
        % 2.2. Generate simulated residual with (-1, 1) sampling.
            jj=1;
            
            while jj < nboot+1
                rr = 1-2*(rand(VAR.T,1)>0.5);
                % T x 1 -1 or 1 random numbers with prob 0.5-0.5
                    % But this will make CI super large?
                    % Why VAR.T? Because VAR.T = length(VAR.vars)-VAR.p;
                    % Therefore, we need to generate this number of shock
                    % sampling for the whole sample.

                % T x n randomly choose the sign of the time T shocks (all)
                resb = (res.*(rr*ones(1,VAR.n)))';
                    % .* means element by element multiplication.
                    % By this operation we thereby obtain a "sampled"
                    % residual, with which we can build new dataset.
                  
                varsb = zeros(VAR.p+VAR.T,VAR.n);
                
                

       %%
        % 2.3. Generate simualted data using these residuals.
                % Initial values: first p periods are true value.
                    varsb(1:VAR.p,:) = VAR.vars_order(1:VAR.p,:);


                % Loop to fill in simulated value.
                for j = VAR.p+1:VAR.p+VAR.T
                    lvars = (varsb(j-1:-1:j-VAR.p,:))';
                    % Put in constant terms.
                        %   const : 0, no constant, no trend
                        %           1, constant, no trend
                        %           2, constant, trend
                        %           3, constant, trend, trend^2
                    trend = 1:size(VAR.vars,1);
                    if VAR.const == 0
                        linevars = lvars(:)';
                        elseif VAR.const == 1
                            linevars = [1 lvars(:)'];
                        elseif VAR.const == 2
                            linevars = [1 trend(j) lvars(:)'];
                        elseif VAR.const == 3
                            linevars = [1 trend(j) trend(j).^2 lvars(:)'];
                    end
                    
                    varsb(j,:) = linevars*VAR.bet + resb(:,j-VAR.p)';
         
                end 
                
                VARBS = VAR;
                VARBS.vars = varsb;    
                VARBS.var_order = cumsum(ones(1,size(VARBS.vars_order, 2)));
                VARBS.pos_shock = pos_shock_bs;

       %%
        % 2.4. Run Cholesky with simulated data and obtain IRF.
                if isfield(VAR,'restr')
                    VARBS = Mod_EstOLS_restricted(VARBS, VARBS.var_order, ...
                        VARBS.restr, VARBS.chi2, VARBS.restrF, ...
                        VARBS.pos_shock, VARBS.ident);
                else
                    VARBS = Mod_EstOLS(VARBS, VARBS.var_order, VARBS.pos_shock, VARBS.ident);
                end
        
                IRF = VARBS.IRF(:,:);

       %%
        % 2.5. Save IRF on jjth and loop on.
                VARbs.IRF(:,jj) = IRF(:);   % Fill the jjth time of wild bootstrap
                % dimension
                    % IRF [n x irhor ]
                    % VARbs.IRF [n x irhor x nboot]
                jj=jj+1;   
            end  
            
	%% 3. Save estimation at different quantile.
        for j=1    
             VARbs.IRFL(:,:,j) = reshape(quantile(VARbs.IRF(:,:,j)',(1-clevel/100)/2),VAR.irhor,size(VARBS.IRF,2));
             VARbs.IRFH(:,:,j) = reshape(quantile(VARbs.IRF(:,:,j)',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.IRF,2));
             
             jj = 1;
             while jj < size(CI,2)+1 
             VARbs.A(:,:,jj,j) = reshape(prctile((VARbs.IRF(:,:,j))',CI(jj)),VAR.irhor,size(VARBS.IRF,2));
             
             jj = jj+1;
             end
        end
 
end

