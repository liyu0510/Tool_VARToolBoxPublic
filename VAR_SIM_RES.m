function out = VAR_SIM_RES(A,B,SIGMA,T,seed,Y0,res)

% This program simulates T realizations from a stochastic process of the form:
% Y_t = A_1*Y_{t-1} + ... + A_P*Y_{t-p} + B + u_t
% where:
% u_t is normal and iid over time
% E(u_t) = 0
% E(u_t*u_t') = SIGMA
%
% Initial conditions used are Y0, if Y0 is specified, or 0 otherwise.
%
% Simulation is done using the seed "seed" for the
%
% The program returns a structure, "out", containing:
% out.Y : the stochastic process Y_t, t=1,...,T (k+1 by T)
% out.u : sequence of shocks (k by T)
% out.EV : eigenvalues of the companion matrix associated to the polynomial
% A
% out.Y0 : the initial conditions
%
% June 2011

k          = size(A,1);
p          = size(A,3);
COMPMATRIX = zeros(k*p,k*p);
for j=1:p
    COMPMATRIX(1:k,(1+(j-1)*k):(j*k)) = A(:,:,j);
    if j < p
    COMPMATRIX((j*k+1):((j+1)*k),(1+(j-1)*k):(j*k)) = eye(k,k);
    end;
end;
    
out.EV     = eig(COMPMATRIX);

%if max(abs(out.EV)) > 1
%    fprintf('\n ... stability condition not met ...\n');
%else
%    fprintf('\n ... stability condition met ...\n');
%end;

if nargin <= 5 
    out.Y0 = zeros(k,p);
else
    out.Y0 = Y0;
end;

%randn('state',seed);

%if k > 1
%out.u = (mvnrnd(zeros(1,k),SIGMA,T))';
%else
%out.u = sqrt(SIGMA).*randn(1,T);
%end;

out.u = res;

out.Y = zeros(k,T);

dummy = zeros(k,T+p);
dummy(:,1:p) = out.Y0(:,1:p);

for t=1:T
    dummyY       = dummy(:,(t-1+p):(-1):(t-p+p));
    dummyY       = dummyY(:);
    dummy(:,t+p) = B + COMPMATRIX(1:k,:)*dummyY + out.u(:,t);
end;

out.Y = dummy(:,(p+1):(p+T));




        



