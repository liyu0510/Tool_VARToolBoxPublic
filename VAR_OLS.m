function out = VAR_OLS(Yobs,p,cons,OD,restr,chi2,restrF)
% out = VAR_OLS(Yobs,p,cons,OD)
% 
% Program to estimate a VAR(p) by OLS
%--------------- Input ----------------%
% Yobs: k x T matrix of data 
% p:    VAR(p)
% cons: 1 if constant should be included in the regression
% OD:   order of the MA representation to be computed
%--------------------------------------%
%
% Uses the general matrix notation:
% Y = B*Z + U
% Y is k x (T-p)
% B is (k) x (k*p+1) if no constant, (k) x (k*(p+1)) if constant
% Z is (k*p) x (T-p+1) if no constant, (k*(p+1))*(T-p+1) if constant
%
%--------------- Output ----------------%
% out.*
% - B: (k x k x p) estimated parameter of the VAR specification without constant
% - BCONS: (k x 1) estimated constant parameter (if cons == 1)
% - u: (k x (T-p)) estimated residuals
% - SIGMA: estimates of the variance-covariance matrix of the residuals (with
% small-sample degree-of-freedom adjustment)
% - BB: The associated moving-average representation
% - mu: The associated long-run mean
%--------------- Restrictions ----------%
% The matris restr imposes restrictions on the VAR coefficients. It is
% equal to 0 for the unrestricted coefficients and 1 elsewhere.
% If chi2=1, the function will display the chi2 statistic for the
% restriction. 

if nargin<6; chi2=0; end

warning('OFF','stats:pvaluedw:ExactUnavailable');

T = size(Yobs,2); % length of time series
k = size(Yobs,1); % number of variables

% Set data
Y = zeros(k,T-p);
if cons == 1      % if constant is included
    Z = zeros(k*p+1,T-p);
else              % otherwise
    Z = zeros(k*p,T-p);
end

for j = 1:(T-p)
    Y(:,j) = Yobs(:,j+p); % set the regressand matrix for each period j
    dummy  = Yobs(:,(j+p-1):(-1):(j)); % set the regressor matrix for the corresponding period j
    dummy  = dummy(:); % line it up
    if cons == 1  % if constant is included
        Z(:,j) = [1;dummy];
    else          % otherwise
        Z(:,j) = dummy;
    end
end
    
% Set estimated coefficients
out.B = zeros(k,k,p); % estimated coefficient (without constant)
if cons == 1
    out.BCONS = zeros(k,1); % estimated coefficient on a constant
end
    
if nargin<5

    dummy = Y*Z'/(Z*Z'); % estimated coefficient
    
else
    % case with restrictions. Need to get zeros into dummy

    % do regression line by line, incorporate restrictions
    if cons==1; 
        dummy=zeros(k,k*p+1); 
    else
        dummy=zeros(k,k*p); 
    end
    
    dummyalt=dummy;
    restr_big=0*dummyalt;
    
    if ~exist('restrF','var'); restrF=0*restr; end
%     keyboard
    
    for rownum=1:k
        off=find(restr(rownum,:)==0);
        on=find(restr(rownum,:)==1);
        offF=find(restrF(rownum,:)==0);
        onF=find(restrF(rownum,:)==1);
        onC = find(restr(rownum,:)==1 & restrF(rownum,:)==0);
        % we want to do a kronecker sum, so have to do the product and take
        % logs. 
%         select=log2(kron(2.^([0:p-1]*k),2.^off));
        if cons==1
            Z1=[Z(1,:); Z(1+off,:)];
            dummy(rownum,[1 off+1]) = Y(rownum,:)*Z1'/(Z1*Z1');
        else
            Z1=[Z(off,:)];
            dummy(rownum,off) = Y(rownum,:)*Z1'/(Z1*Z1');
        end
        
        if chi2==1
            
            if cons==1
                Z1=[Z(1,:); Z(1+offF,:)];
                dummyalt(rownum,[1 offF+1]) = Y(rownum,:)*Z1'/(Z1*Z1');
                restr_big(rownum,onC+1)=1;
            else
                Z1=[Z(offF,:)];
                dummyalt(rownum,[offF]) = Y(rownum,:)*Z1'/(Z1*Z1');
                restr_big(rownum,onC)=1;
            end
        end
    end
end

if cons ~= 1
    for kk=1:p
        out.B(:,:,kk) = dummy(:,((kk-1)*k+1):(kk*k));
    end
end

if cons == 1
    for kk=1:p
        out.B(:,:,kk) = dummy(:,(1+(kk-1)*k+1):(1+kk*k));
    end
    out.BCONS = dummy(:,1);
end

% keyboard

% Get estimated residuals and covariance matrix
out.u = Y-dummy*Z; % residuals
out.SIGMA = (1/(T-k*p-1))*(out.u*out.u'); % covariance matrix

% Test restrictions
if chi2==1
    u_alt=Y-dummyalt*Z;
    SIGMA_alt=(1/(T-k*p-1))*(u_alt*u_alt');
    cv=kron(SIGMA_alt,inv(cov(Z(2:end,:)')))/size(Y,2);
    restr_big=restr_big(:,2:end)';
    dummyalt=dummyalt(:,2:end)';

    R=diag(restr_big(:));
    R=R(sum(R,2)==1,:);
    chi2_=(R*dummyalt(:))'*inv(R*cv*R')*(R*dummyalt(:));
    nrest=size(R,1);

    disp('chi2_stat      DF      p-value');
    disp([chi2_ nrest  1-chi2cdf(chi2_,nrest)]);
end

% Get sums of coefficients and p-value
bcv=kron(inv(cov(Z(2:end,:)')),out.SIGMA)/size(Y,2);
dd=dummy(:,2:end);dd=dd(:);
sum_restr=repmat(eye(k^2),1,p);
out.sumcoeff=sum_restr*dd;
out.sumcoeff_cov=sum_restr*bcv*sum_restr';
    
% Get moving average coefficients
out.BB = zeros(k,k,OD);
out.BB(:,:,1) = eye(k,k);

for ii = 1:OD
    dummy = zeros(k,k);
    for jj=1:min(ii,p)
        dummy = dummy + out.BB(:,:,ii-jj+1)*out.B(:,:,jj);
    end
    out.BB(:,:,ii+1) = dummy;
end

% Get long-run average coefficients
if cons ~=1
    out.mu = zeros(k,1);
end

if cons == 1
    dummy = eye(k,k);
    for s=1:p
        dummy = dummy - out.B(:,:,s);
    end
    out.mu = dummy^(-1)*out.BCONS;
end
    

out.X=Z;
out.Y=Y;
