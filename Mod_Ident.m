%% Mod_Ident Function
    % Goal:
        % Perform further identification with estimated VAR model.
    % Input:
        % Input 1: ident{'name of identification', parameters}
    % Output:
        % Output 1: Identification: B matrix (whole or partial).
        % Output 2: IRF.
        % Output 3: FEVD.
%%

function [VAR_Out] = Mod_Ident(ident)

    %% 1.Define proxy parameter.
    if (strcmp(ident{1}, 'proxy')) 
        disp('Proxy identification');
        VAR_Out = Ident_Proxy(ident{2:length(ident)});
    %%
    
    %% 2.Define sign restriction parameter.
    elseif (strcmp(ident{1}, 'sign')) 
        disp('Sign restriction identification.');
        VAR_Out = Ident_Sign(ident{2:length(ident)});
    %%
    
    %% 3.Define max share parameter.
    elseif (strcmp(ident{1}, 'maxshare')) 
        disp('Max share identification.');
        VAR_Out = Ident_Maxshare(ident{2:length(ident)});
    %%
    
    %% 4.Identification scheme not specified.
    else
        disp('Identification not correctly specified.');
    %%
    end
end



