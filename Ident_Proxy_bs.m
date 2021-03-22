%% Ident_Proxy_bs Function
    % Use wild bootstrap to obtain confidence interval for Cholesky estimation.

function VARbs = Ident_Proxy_bs(VAR, Para)

     nboot  = Para.nboot;
     clevel = Para.clevel;
     CI     = Para.CI;

	%% 1. Perform wild bootstrap
    
       %%
        % 1.1. Removes a constant.
            % What if there is no constant to begin with?
            % Need to ask.


            res = detrend(VAR.res,'constant');
            % Up until now, we only consider the case of no constant and
            % with constant. Not the case with time trend or squared time
            % trend.

                
       %%
        % 1.2. Generate simulated residual with (-1, 1) sampling.
            jj=1;
            
            while jj < nboot + 1
                rr = 1-2*(rand(VAR.T,1)>0.5);
                %T x 1 -1 or 1 random numbers with prob 0.5-0.5
                       %T x n randomly choose the sign of the time T shocks (all)
                       
                resb = (res.*(rr*ones(1,VAR.n)))';
                varsb = zeros(VAR.p+VAR.T,VAR.n);

       %%
        % 1.3. Generate simualted data using these residuals.    
                % Initial values
                varsb(1:VAR.p,:) = VAR.vars_order(1:VAR.p,:);
               
                for j=VAR.p+1:VAR.p+VAR.T
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
                
                VARBS = Mod_EstOLS(VARBS, [], [], 'non');
                % Bootstrap ”Î identification scheme Œﬁπÿ£ø£ø£ø
                

       %%
        % 1.4. Generate proxy data.
                % Proxy variable also needs to be rearranged in the same manner.
                VARBS.proxies = [VAR.m.*(rr(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1)*ones(1,size(VAR.proxies,2)))];
       
       %%
        % 1.5. Assemble proxy dataset.
                DATASET_FACTOR_bs.TSERIES(VARBS.smpl_min_FACTORS:VARBS.smpl_max_FACTORS,:) = VARBS.proxies;        % Reads data, with names
                DATASET_FACTOR_bs.TSERIES(VARBS.smpl_min_FACTORS:VARBS.smpl_max_FACTORS,(size(DATASET_FACTOR_bs.TSERIES,2)+1)) = ...
                    VARBS.DATASET_FACTOR.TSERIES(VARBS.smpl_min_FACTORS:VARBS.smpl_max_FACTORS,1);
                DATASET_FACTOR_bs.TSERIES(VARBS.smpl_min_FACTORS:VARBS.smpl_max_FACTORS,(size(DATASET_FACTOR_bs.TSERIES,2)+1)) = ...
                    VARBS.DATASET_FACTOR.TSERIES(VARBS.smpl_min_FACTORS:VARBS.smpl_max_FACTORS,2);

                DATASET_FACTOR_bs.LABEL   = [VARBS.select_factors,'Year','Month'];     % Extract variable names
                DATASET_FACTOR_bs.VALUE   = 1:length(DATASET_FACTOR_bs.LABEL);

                DATASET_FACTOR_bs.MAP = containers.Map(DATASET_FACTOR_bs.LABEL,DATASET_FACTOR_bs.VALUE);
%                  DATASET_FACTOR_bs.TSERIES = table2array(DATASET_FACTOR_bs.TSERIES);           % Reads data, with names
                    
       %%
        % 1.6. Run proxy estimation with simulated data and obtain IRF.
                
                % Define proxy input.
                if VARBS.NofInputs == 4
                    ident_input_proxy_bs = {VARBS, DATASET_FACTOR_bs, VARBS.select_factors, VARBS.proxyposit};
                elseif VARBS.NofInputs == 6
                    ident_input_proxy_bs = {VARBS, DATASET_FACTOR_bs, VARBS.select_factors, VARBS.proxyposit ...
                       , VARBS.FactorStart, VARBS.FactorEnd};
                elseif VARBS.NofInputs == 8
                    ident_input_proxy_bs = {VARBS, DATASET_FACTOR_bs, VARBS.select_factors, VARBS.proxyposit ...
                        , VARBS.FactorStart, VARBS.FactorEnd ...
                        , VARBS.Factor_label_cell, VARBS.Factor_label_short_cell};
                end  

                % Run proxy identification
                VARBS = Ident_Proxy(ident_input_proxy_bs{1:length(ident_input_proxy_bs)});
                
                % What is VAR.k and why it's in here but not in OLS_bs?
                    % VAR.k is the number of variables being proxied.
                    % Always 1 at this stage.
                    for j=1:VAR.k
                        IRF = VARBS.IRF(:,:,j);
                        VARbs.IRF(:,jj,j) = IRF(:);
                    end
              clear DATASET_FACTOR_bs
              jj=jj+1;   

            end  

        %% 2. Save lower and upper estimation.
         for j=1:VAR.k  
             VARbs.IRFL(:,:,j)=reshape(quantile(VARbs.IRF(:,:,j)',(1-clevel/100)/2),VAR.irhor,size(VARBS.IRF,2));
             VARbs.IRFH(:,:,j)=reshape(quantile(VARbs.IRF(:,:,j)',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.IRF,2));
             
             jj = 1;
             while jj < size(CI,2)+1 
             VARbs.A(:,:,jj,j) = reshape(prctile((VARbs.IRF(:,:,j))',CI(jj)),VAR.irhor,size(VARBS.IRF,2));
             
             jj = jj+1;
             end
         end