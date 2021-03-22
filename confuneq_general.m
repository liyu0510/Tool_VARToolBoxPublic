        function [c, ceq] = confuneq_general(q)
        % Nonlinear inequality constraints
        % Nonlinear equality constraints
        c   = [];
        ceq = q'*q-eye(size(q,2)); 
            % Orthonormal matrix
            % 减这个零矩阵是干啥？？
            % 以及这个 case 难道不是更 general？包含了 confuneq？？
      
        