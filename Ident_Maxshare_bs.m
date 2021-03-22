%% Ident_maxshare_bs Function
    % Use wild bootstrap to obtain confidence interval for Cholesky estimation.

function VARbs = Ident_Maxshare_bs(VAR, Para)

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
                VARBS.var_order = cumsum(ones(1,size(VARBS.vars_order, 2)));
                
                VARBS = Mod_EstOLS(VARBS, VARBS.var_order, VARBS.pos_shock, 'short');
                % Bootstrap ”Î identification scheme Œﬁπÿ£ø£ø£ø
                 
       %%
        % 1.4. Run maxshare estimation with simulated data and obtain IRF.
                
                % Define maxshare input.
                if VARBS.NofInputs == 3
                    ident_input_maxshare_bs = {VARBS, VARBS.N, VARBS.k};
                elseif VARBS.NofInputs == 4
                    ident_input_maxshare_bs = {VARBS, VARBS.N, VARBS.k, VARBS.k0};
                elseif VARBS.NofInputs == 5
                    ident_input_maxshare_bs = {VARBS, VARBS.N, VARBS.k, VARBS.k0, VARBS.q0};
                elseif VARBS.NofInputs == 6
                    ident_input_maxshare_bs = {VARBS, VARBS.N, VARBS.k, VARBS.k0, VARBS.q0, VARBS.m};
                end  
                

                % Run maxshare identification
                VARBS = Ident_Maxshare(ident_input_maxshare_bs{1:length(ident_input_maxshare_bs)});
                
                % What is VAR.k and why it's in here but not in OLS_bs?
                    % VAR.k is the number of variables being proxied.
                    % Always 1 at this stage.
                    for j=1:size(VARBS.Ident_column,2)

                        IRF = VARBS.IRF(:,:,j);
                        VARbs.IRF(:,jj,j) = IRF(:);
                    end
  
              jj=jj+1;   

            end  

        %% 2. Save lower and upper estimation.
         for j=1:size(VARBS.Ident_column,2)
             VARbs.IRFL(:,:,j)=reshape(quantile(VARbs.IRF(:,:,j)',(1-clevel/100)/2),VAR.irhor,size(VARBS.IRF,2));
             VARbs.IRFH(:,:,j)=reshape(quantile(VARbs.IRF(:,:,j)',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.IRF,2));
             
             jj = 1;
             while jj < size(CI,2)+1 
             VARbs.A(:,:,jj,j) = reshape(prctile((VARbs.IRF(:,:,j))',CI(jj)),VAR.irhor,size(VARBS.IRF,2));
             
             jj = jj+1;
             end
             
             
         end