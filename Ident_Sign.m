%% Ident_Sign Function
    % Goal:
        % Use sign restriction identification method to identify the VAR model.
        % Sign restriction method as in Uhlig (2005) rejection method.
    % Input(varagin):
        % input(must):
            % Input 1: Estimated VAR.
            % Input 2: Sign restriction matrix.
        % input(optional):
            % Input 5: .
    % Output:
        % Output 1: Identification B matrix.
        % Output 2: IRF.
        % Output 3: FEVD.
    % Comment:
        % 111: Why? Shouldn't it be n - nzar?

%%

function [VAR] = Ident_Sign(varargin)
% =======================================================================
% Compute IRFs, FEVDs, and HDs for a VAR model estimated with VARmodel and 
% identified with sign restrictions
% =======================================================================
% VAR = SR(VAR,R,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - VAR: structure, result of VARmodel function
%   - R: 3-D matrix containing the sign restrictions (n,3,nshocks)
%       described below
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% -----------------------------------------------------------------------
% OUTPUT
%   - VAR
%       * IRFall  : 4-D matrix of IRFs  (irhor,n,nshocks,ndraws)
%       * FEVDall : 4-D matrix of FEVDs (irhor,n,nshocks,ndraws)
%       * HDall   : 4-D matrix of HDs (irhor,n,nshocks,ndraws)
%       * Ball : 4-D matrix of Bs (n,n,nshocks,ndraws)
%       * IRFmed  : median of IRFall
%       * FEVDmed : median of FEVDall
%       * Bmed : median of Ball
%       * HDmed   : median of HDall
%       * IRFinf  : 16th percentile of IRFall
%       * FEVDinf : 16th percentile of FEVDall
%       * IRFinf  : 16th percentile of HDall
%       * IRFsup  : 84th percentile of IRFall
%       * FEVDsup : 84th percentile of FEVDall
%       * HDsup   : 84th percentile of HDall
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% Note. This code follows the notation as in the lecture notes available at
% https://sites.google.com/site/ambropo/MatlabCodes

% The R matrix is a 3-D matrix with dimension (n,3,nshocks). Assume you
% have 3 variables and you want to identify just one shock. Then:
% 
%              from        to          sign
%  R(:,:,1) = [ 1           4           1          % VAR1
%               1           4          -1          % VAR2
%               0           0           0];        % VAR3
% 
%   - The first column defines the first period from which the restriction is imposed 
%   - The second column defines the last period to which the restriction is imposed 
%   - The last column defines the sign of the restriction: positive (1) or negative (-1)
%   - To leave unrestricted set all three columns to zero
% 
% In the above example we set the following restrictions. VAR1 must respond
% with positive sign from period 1 to period 4; VAR2 must respond with 
% negative sign from period 1 to period 4; VAR3 is left unrestricted
%
% An additional shock could be defined as follows:
%              from        to          sign
%  R(:,:,2) = [ 1           4           1          % VAR1
%               1           4           1          % VAR2
%               1           4          -1];        % VAR3

    

    %% 1.Determine the sign inputs.
    
        if nargin == 3
                Para                         = varargin{1};
                VAR                          = varargin{2};
                VAR.R                        = varargin{3};
                VAR.ndraws                   = 1000;
        elseif nargin == 4
                Para                         = varargin{1};
                VAR                          = varargin{2};
                VAR.R                        = varargin{3};
                VAR.ndraws                   = varargin{4}; 
        end

       % retrieve parameters
            % VAR parameters
                n     = VAR.n;
                irhor = VAR.irhor;
            % sign restriction parameters
                R       = VAR.R;
                ndraws  = VAR.ndraws;
                nshocks = size(R,3);
            % confidence interval parameters
                clevel = Para.clevel;
                CI     = Para.CI;
            % number of zero restrictions (nzr) % what is zero restriction anyway
                nzr = 0;      
                for ss=1:nshocks
                    if sum(sum(R(:,:,ss))) == 0 % if the sum is zero, there is no sign restriction
                        nzr = nzr +1;
                    end    
                end
                
            % number of sign restrictions (nsr)  
%                  nsr = nshocks - nzr;
                % why?? 不是应该 n - nzar吗？

                   nsr = n - nzr;

        % initialize empty matrix for the IRF draws
            Bstore    = nan(n,nshocks,ndraws); 
            IRFstore  = nan(irhor,n,nshocks,ndraws); % total draws? not accepted draws?
            FEVDstore = nan(irhor,n,nshocks,ndraws); 
%             HDstore   = nan(T+p,n,n,ndraws); 

        % initialize the vectors for the candidate IRF
            irhor_check = max(max(R(:,2,:)));            % maximal length of IRF function to be checked
            irhor_R     = max(max(R(:,2,:)-R(:,1,:)))+1; % number of periods with restriction



    %% 2.Sign restriction with Uhlig (2005) rejection method.

            jj = 0; % accepted draws
            kk = 0; % total draws
            ww = 1; % index for printing on screen
            
            while jj < ndraws
                kk = kk+1;
                if jj == 0 && kk>4000
                    error('Max number of iterations reached. Try different sign restrictions') 
                end

            % 2.1.draw a random orthonormal shifting matrix
                % flat prior for a
                    % randn gives you random variables normally distributed
                if nzr == 0 % no zero restrictions
                    S = OrthNorm(n);
                else % with zero restrictions
                    % so basically, you want to maintain the shifting to
                    % lower right box of the original B matrix
                    % this means that you are only shifting the n-nsr+1
                    % onward columns, and not the first n-nsr columns
                    S = eye(n);
                    auxS = OrthNorm(nsr); 
                    % But if nsr=1 then it's always 1? No shifting??????????
                    % Now, this shifting matrix approach is a bit weird.
                    % Not just a normalized vector, but a matrix?
                    S(n-nsr+1:n,n-nsr+1:n) = auxS;
                    % Rotate the matrix but only alter the last nsr columns
                    clear auxS
                end

            % 2.2.draw a F and sigma from the posterior
                % normal-wishart prior for Σ and beta
                [sigma_draw, bet_draw] = Mod_EstBayesian_drawpost(VAR);


            % 2.3.set up VAR_draw
                VAR_draw       = VAR;
                VAR_draw.Ft    = bet_draw;
                VAR_draw.sigma = sigma_draw;
                VAR_draw.S     = S;
                VAR_draw.irhor = irhor_check;                
                VAR_draw.Bnew  = VAR.B*S;

            % 2.4.compute IRFs only for the restricted periods
                % IRF_draw has dimension [irhor_check x n]
                [IRF_draw, VAR_draw] = Res_IRF(VAR_draw,VAR_draw.Bnew);


            % 2.5.check whether sign restrictions are satisfied
                checkall = ones(n,nshocks);
                checkall_flip = ones(n,nshocks); 
                % why do you need a second one? In case of a negative shock?
                
                for ss=1:nshocks 
                    % how many shock? how many columns
                    for ii = 1:n
                        if R(ii,1,ss) ~= 0
                            
                            if R(ii,3,ss) == 1     % if the sign restriction is positive
                                check           = IRF_draw((R(ii,1,ss)):R(ii,2,ss),ii,ss) > 0;
                                checkall(ii,ss) = min(check);
                                % Check flipped signs
                                check_flip           = IRF_draw((R(ii,1,ss)):R(ii,2,ss),ii,ss) < 0;
                                checkall_flip(ii,ss) = min(check_flip);
                                
                            elseif R(ii,3,ss) == -1 % if the sign restriction is negative
                                check           = IRF_draw((R(ii,1,ss)):R(ii,2,ss),ii,ss) < 0;
                                checkall(ii,ss) = min(check);
                                % Check flipped signs
                                check_flip           = IRF_draw((R(ii,1,ss)):R(ii,2,ss),ii,ss) > 0;
                                checkall_flip(ii,ss) = min(check_flip);
                            end
                        end
                    end
                    clear aux aux_flip
                end
    
                % If restrictions are satisfied, compute IRFs for the desired periods (nstep)
                VAR_draw.irhor = irhor;
                                
                    if min(min(checkall))==1 
                        jj = jj+1;
                        % IRF
                        [aux_irf, VAR_draw] = Res_IRF(VAR_draw, VAR_draw.Bnew(:,1:nshocks));
                        IRFstore(:,:,:,jj)  = aux_irf;
                        % FEVD
                        aux_fevd = Res_FEVD(VAR_draw, VAR_draw.Bnew(:,1:nshocks));
                        FEVDstore(:,:,:,jj) = aux_fevd;
                        % B
                        aux_B = VAR_draw.Bnew(:,1:nshocks);
                        Bstore(:,:,jj) = aux_B;
                        
                        % Display number of loops
                        if jj==10*ww
                            disp(['Loop: ' num2str(jj) ' / ' num2str(kk) ' draws']);
                            ww=ww+1;
                        end
                        
                    elseif min(min(checkall_flip)) == 1 
                        jj = jj+1 ;
                        % IRF (change the sign!)
                        [aux_irf, VAR_draw] = Res_IRF(VAR_draw, VAR_draw.Bnew(:,1:nshocks));
                        IRFstore(:,:,:,jj)  = -aux_irf;
                        % FEVD
                        aux_fevd = Res_FEVD(VAR_draw, VAR_draw.Bnew(:,1:nshocks));
                        FEVDstore(:,:,:,jj)  = aux_fevd;
                        % B
                        aux_B = VAR_draw.Bnew(:,1:nshocks);
                        Bstore(:,:,jj) = aux_B;
                        
                        % Display number of loops
                            if jj==10*ww
                                disp(['Loop: ' num2str(jj) ' / ' num2str(kk) ' draws']);
                                ww=ww+1;
                            end
                    end
            end


    %% 3.Store results
        VAR = rmfield(VAR,'IRF');  % remove previous results.

        % 3.1.Store all accepted IRFs and FEVDs
            VAR.IRFall  = IRFstore;
            VAR.FEVDall = FEVDstore;
            VAR.Ball    = Bstore;

        % 3.2.Compute and save median impulse response
            VAR.IRF(:,:,:)     = median(IRFstore,4);
            VAR.FEVDmed(:,:,:) = median(FEVDstore,4);
            VAR.Bmed(:,:,:)    = median(Bstore,3);

        % 3.3.Compute lower and upper bounds
            aux      = prctile(IRFstore,[(100-clevel) clevel],4);
            VAR.IRFL = aux(:,:,:,1);
            VAR.IRFH = aux(:,:,:,2);
            
            jj = 1;
            while jj < size(CI,2)+1 
            VAR.A(:,:,jj) = reshape(prctile(IRFstore,CI(jj),4),VAR.irhor,size(VAR.IRF,2));
             
            jj = jj+1;
            end

            aux         = prctile(FEVDstore,[(100-clevel) clevel],4);
            VAR.FEVDinf = aux(:,:,:,1);
            VAR.FEVDsup = aux(:,:,:,2);

end
