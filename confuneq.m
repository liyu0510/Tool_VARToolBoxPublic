        function [c, ceq] = confuneq(q)
        % Nonlinear inequality constraints
        % Nonlinear equality constraints
        c   = [];
        ceq = norm(q)-1;
        