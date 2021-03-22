%% Ident_Proxy Function
    % Goal:
        % Use proxy identification method to identify the VAR model.
        % Proxy method as in Mertens and Ravn (2013).
        % Take instrument(factor) as give. 
        % Adjust if instrument is not of the same length as other variables.
    % Input(varagin):
        % input(must):
            % Input 1: Estimated VAR. 
            % Input 2: Factor dataset. 
            % Input 3: Factor name.
            % Input 4: Position of proxied.
        % input(optional):
            % Input 5: Factor label.
            % Input 6: Factor short label.
            % Input 7: Factor start time.
            % Input 8: Factor end time.
    % Output:
        % Output 1: Identification: one column of identified B matrix.
        % Output 2: IRF.
        % Output 3: FEVD.
%%

function [VAR] = Ident_Proxy(varargin)

    %% 1.Determine the proxy inputs.
        
        if nargin == 4
            % Case 1: Only VAR and Factor name.            
                VAR                          = varargin{1};
                VAR.DATASET_FACTOR           = varargin{2};
                VAR.select_factors           = varargin{3};
                VAR.proxyposit               = varargin{4};
                VAR.Factor_label_cell        = VAR.select_factors;
                VAR.Factor_label_short_cell  = VAR.select_factors;
                
                VAR.NofInputs = 4;

            elseif nargin == 6
            % Case 2: VAR, Factor name, start & end time.
                VAR                          = varargin{1};
                VAR.DATASET_FACTOR           = varargin{2};
                VAR.select_factors           = varargin{3};
                VAR.proxyposit               = varargin{4};
                VAR.FactorStart              = varargin{5};
                VAR.FactorEnd                = varargin{6};
                VAR.Factor_label_cell        = VAR.select_factors;
                VAR.Factor_label_short_cell  = VAR.select_factors;
                
                VAR.NofInputs = 6;

            elseif nargin == 8
            % Case 3: VAR, Factor name, start & end time, label & short label.
                VAR                           = varargin{1};
                VAR.DATASET_FACTOR            = varargin{2};
                VAR.select_factors            = varargin{3};
                VAR.proxyposit                = varargin{4};
                VAR.FactorStart               = varargin{5};
                VAR.FactorEnd                 = varargin{6};
                VAR.Factor_label_cell         = varargin{7};
                VAR.Factor_label_short_cell   = varargin{8};
                
                VAR.NofInputs = 8;
        else
            disp ('Wrong number of inputs.');
        end

        % Adjust to put the proxied variable to the first.
            proxy_order     = [cumsum(ones(1,size(VAR.vars, 2)))];
            proxy_order(proxy_order == VAR.proxyposit) = [];
            VAR.proxy_order = [VAR.proxyposit, proxy_order];
            VAR.vars_proxy  =   VAR.vars(:,VAR.proxy_order); 
        
    %% 2.Extract chosen factor variables from dataset.
    
       %%
        % 2.1.Factors sample starts minimum p periods after the VAR sample.
            % Now, if the data is in quarterly frequency, p periods means 4p months.
            %	VAR.p/4
            %If the data is in monthly frequency, p periods means p months.
            %	VAR.P/12

            if isfield(VAR,{'FactorStart'})
                smpl_min_factors_vec  = VAR.FactorStart;
                    if  smpl_min_factors_vec(1,1) >= VAR.Start(1,1) + VAR.p/4
                        smpl_min_FACTORS =  smpl_min_factors_vec;
                    else
                        smpl_min_FACTORS(1,1) = VAR.Start(1,1) + floor(VAR.p/4);
                        smpl_min_FACTORS(1,2) = VAR.Start(1,2) + 3*(VAR.p-4*floor(VAR.p/4));
                    end

            else
                smpl_min_FACTORS(1,1) = VAR.Start(1,1) + floor(VAR.p/12);
                smpl_min_FACTORS(1,2) = VAR.Start(1,2) + (VAR.p-12*floor(VAR.p/12));
            end

       %%    
        % 2.2.Factors sample ends maximum minimum before the VAR sample.
            if isfield(VAR,{'FactorEnd'})
                smpl_max_factors_vec  = VAR.FactorEnd;
                    if (smpl_max_factors_vec(1,1) > VAR.End(1,1))
                        smpl_max_FACTORS    =   VAR.End;
                    else
                        smpl_max_FACTORS = smpl_max_factors_vec; %smpl_max_VAR; %[2012 6]; %this has crisis in it!
                    end
            else
                smpl_max_FACTORS = VAR.End;
            end

       %% 
        % 2.3.Find the dates in the dataset.    
            VAR.smpl_min_FACTORS = find(and(VAR.DATASET_FACTOR.TSERIES(:,cell2mat(values(VAR.DATASET_FACTOR.MAP,{'Year'})))==smpl_min_FACTORS(1,1), ...
                VAR.DATASET_FACTOR.TSERIES(:,cell2mat(values(VAR.DATASET_FACTOR.MAP,{'Month'})))==smpl_min_FACTORS(1,2)));
            VAR.smpl_max_FACTORS = find(and(VAR.DATASET_FACTOR.TSERIES(:,cell2mat(values(VAR.DATASET_FACTOR.MAP,{'Year'})))==smpl_max_FACTORS(1,1), ...
                VAR.DATASET_FACTOR.TSERIES(:,cell2mat(values(VAR.DATASET_FACTOR.MAP,{'Month'})))==smpl_max_FACTORS(1,2)));

       %% 
        % 2.4.Find variables in the dataset.
    
        VAR.proxies = VAR.DATASET_FACTOR.TSERIES(VAR.smpl_min_FACTORS:VAR.smpl_max_FACTORS,cell2mat(values(VAR.DATASET_FACTOR.MAP,VAR.select_factors))); 

        % number of proxies (at this stage default to be one).
            % Idea: Multiple high frequency time series as proxy?
        VAR.k  = 1;
        
        % Assuming proxies start at least p periods later
        VAR.m  = VAR.proxies(1:end,:);
            % Why do we need this?

        [VAR.T_m,VAR.n_m] = size(VAR.proxies);
        VAR.T_m_end       = VAR.smpl_max_VAR-VAR.smpl_max_FACTORS;      
            %The difference between where the VAR and the factors end in the VAR sample

       %% 
        % 2.5.Locate the position of the variable to be instrumented.

    %% 3.Proxy identification.
    
        % Only the restricted sample is used for identification.
        VAR.Sigma_m = (VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:)'*VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:))/(VAR.T_m-VAR.n*VAR.p-1);
        Phib        = [ones(VAR.T_m,1) VAR.m]\VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:);
            % What is this Phib? It is the coefficient of regress all residuals on exogenous proxy.
 
        Res_m     =  VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:)-[ones(VAR.T_m,1) VAR.m]*Phib;
        Res_const =  VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1)-ones(VAR.T_m,1)*(ones(VAR.T_m,1)\VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1));
        XX_m      =  [ones(VAR.T_m,1) VAR.m];

        % Calculate robust standard errors
        SS_m  =   zeros(VAR.n_m+1,VAR.n_m+1);
        for ii=1:VAR.T_m
            SS_m  =   SS_m+1/VAR.T_m*XX_m(ii,:)'*XX_m(ii,:)*Res_m(ii,1)^2;
        end

        Avarb_m  =   inv(1/VAR.T_m*XX_m'*XX_m)*SS_m*inv(1/VAR.T_m*XX_m'*XX_m);
        RR_m     =   [zeros(VAR.n_m,1) eye(VAR.n_m)];
        WW_m     =   VAR.T_m*(RR_m*Phib(:,1))'*inv(RR_m*Avarb_m*RR_m')*(RR_m*Phib(:,1));
        VAR.F_m_rob     =   WW_m/VAR.n_m;

        SST_m       = Res_const'*Res_const;
        SSE_m       = Res_m(:,1)'*Res_m(:,1);
        VAR.F_m     = ((SST_m-SSE_m)/VAR.n_m)/(SSE_m/(length(VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1))-(VAR.n_m+1)));
        VAR.R2_m    = (1-SSE_m/(SST_m));
        VAR.R2adj_m = VAR.R2_m-(VAR.n_m/((length(VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1))-(VAR.n_m+1))))*(1-VAR.R2_m); 
        
        % Mertens and Ravn (2013) method of proxy identification:
            % 1) Two stage regression:
                % FIRFt stage: obtain a forecast of u1.
                uhat1           =   [ones(VAR.T_m,1) VAR.m]*Phib(:,1);
                    % Here we only take the first column, which we believe
                    % to be u^p_t.
                    
                    % However, what if the proxied variable is not ranked
                    % at first?
                    % In this case, we might have to change Phib(:,1) to
                    % Phib(:,proxied), where proxied is the location of the
                    % variable that we want to instrument.
                    
                    % ...Or do we?
                    % An easier way might be to reordere the variables, as
                    % we do in Cholesky. When doing proxy estimation, we
                    % always reorder the variables to rank the proxied
                    % variable to be the first variable in line.
                    % So this is what we will do.
                
                % Second stage: obtain consistent estimator of sq/sp.
                b21ib11_TSLS    =   [ones(VAR.T_m,1) uhat1]\VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,2:end);  
                    % This step is to regress all the variables
                b21ib11_TSLS    =   b21ib11_TSLS(2:end,:)';
                b21ib11         =   b21ib11_TSLS;

            % 2) Matrix calculation: Identification of b11 and b12 from the covariance matrix of the VAR
                % Define matrix block.
                Sig11   = VAR.Sigma_m(1:VAR.k,1:VAR.k);
                Sig21   = VAR.Sigma_m(VAR.k+1:VAR.n,1:VAR.k);
                Sig22   = VAR.Sigma_m(VAR.k+1:VAR.n,VAR.k+1:VAR.n);
                ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;  % Q matrix in GK
                
                % Calculate true sp.
                b12b12p  = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11)); 
                b11b11p  = Sig11-b12b12p;
                b11      = sqrt(b11b11p);
                VAR.b1   = [b11; b21ib11*b11];  % VAR.b1 as a vector
                VAR.Phib = Phib;

    %% 4.Obtain result: IRF and FEVD.
        % initial shock: eps(1,1)=1
       
            IRF      = Res_IRF(VAR, VAR.b1);
            VAR.IRF  = IRF;
            FEVD     = Res_FEVD(VAR, VAR.b1);
            VAR.FEVD = FEVD;

%         % 2) Calculate it directly.  
%             IRF(VAR.p+1,:) = VAR.b1(:,1);   
%                
% 
%             for jj=2:VAR.irhor
%                 lvars = (IRF(VAR.p+jj-1:-1:jj,:))';
%                 IRF(VAR.p+jj,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
%             end
% 
%             VAR.IRF   = IRF(VAR.p+1:VAR.p+VAR.irhor,:);

end