%% Res_FEVD_sum Function
    % Goal:
        % Use max share identification method to identify the VAR model
    % Input:
        % Input 1: VAR with short run identification results.
        % Input 2: Horizon to calculate.
    % Output:
        % Output 1: Sum of FEVD of VAR.Nth variable over VAR.k periods.

%%
function [FEVD_total] = Res_FEVD_sum(VAR, q)
            Ident_Column_q = VAR.B*q;   
            [FEVD_q,~]     = Res_FEVD(VAR, Ident_Column_q);

            % Compute cumulative sum.
            % FEVD(t,i,j):period t, the MSE of variable i, explained by variable j.

            FEVD_interest = FEVD_q(:,:,VAR.N);   
                % We are interested in the FEVD of the [Nth] variable.
            FEVD_total =  sum(FEVD_interest((VAR.k0:VAR.k),(1:size(Ident_Column_q,2))),'all');     
                % We care about the total FEVD of [k0:k] periods.
        end