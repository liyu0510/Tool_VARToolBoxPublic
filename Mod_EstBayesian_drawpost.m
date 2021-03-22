%% Mod_EstBayesian_drawpost Function
    % Goal:
        % Draw from the posterior distribution of a VAR model.
    % Input:
        % Input 1: estimated VAR.
    % Output:
        % Output 1: posterior sigma matrix.
        % Output 2: posterior coefficient matrix.
    % Comment:
        % 27: inv_sigma_draw = wishrnd(inv_sigma_hat/T,T);   % T inside? Why? 


function [sigma_draw, bet_draw] = Mod_EstBayesian_drawpost(VAR)


%% 1.Get relevant information from VAR structure.
%     k    = [];   % what is k?

    T             = VAR.T;   
    X             = VAR.X;
    bet_hat       = VAR.bet;
    sigma_hat     = VAR.Sigma;
    inv_sigma_hat = inv(sigma_hat);
    

%% 2.Draw the VCV matrix (sigma).
    inv_sigma_draw = wishrnd(inv_sigma_hat/T,T);     % No S0? Because we assume S0=0
%     inv_sigma_draw = wishrnd(inv_sigma_hat,T)/T;     % Also work?
%     inv_sigma_draw = wishrnd(inv_sigma_hat,T);       % Not correct
        
    sigma_draw     = inv(inv_sigma_draw);


%% 3.Draw the coeffiecient matrix (bet).
    aux1       = inv(X'*X);
    aux2       = kron(sigma_draw,aux1);
    bethat_vec = bet_hat(:);
    betdraw    = mvnrnd(bethat_vec,aux2);             % No B0? Beacuse we assume N0=0
    bet_draw   = reshape(betdraw,size(bet_hat));
    
end










