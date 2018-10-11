function [beta,  rawBeta, R2, actInd, betaAll] = get_FWA_l1LineSearch(Y,X, lb,ub, maxJ, standardizeFlag, nVar)
% estimates the frank wolfe algorithm with l1 constraint: fully vectorized. For best
% performance, use XX with diagonal equal to 1.




if nargin <3
    
    lb = -1;
    
end


if nargin <4
    
    ub = 1;
    
end

if nargin <5
    
   maxJ = 500;
    
end

if nargin < 6
   
    standardizeFlag =0;
    
end


XX                  = X'*X;
XY                  = X'*Y;
YY                  = Y'*Y;
mY                  = sum(Y);
n                   = size(X,1);
if nargin <7
    
    
    nVar = size(XX,1);
    
end

K                 = size(XX,1);

if standardizeFlag==1

      indZero        = (diag(XX)==0);
      rootXX         = sqrt(diag(XX)/n);
      rootXX(indZero) = 1;
      D              = diag(1./rootXX)*(XX/n)*diag(1./rootXX);
      XY             = bsxfun(@rdivide, XY/n, rootXX);      
      X              = X*diag(1./rootXX);
else
    
      rootXX         = ones(K,1);
      D              = XX/n;
      XY             = XY/n;      
    
    
    
end

beta             = zeros(K,1);
betaAll          = cell(maxJ,1);
obj              = 0;
j                = 1;
nActive          = 0;
actInd           = zeros(floor(maxJ)+2,1);
stop             = false;
while (j<=maxJ) && (nActive<= nVar) && ~stop

A                = (XY-D*beta);
[ma sj]          = max(abs(A));

a                = zeros(K,1);
a(sj)            = ub*(A(sj)>0)+lb*(A(sj)<0);



% beta             = (1-2/(1+j))*beta+a*2/(1+j);
% 

func             = @(gamma) objective(X,Y,beta, a, gamma);
gamma            = fminbnd(func,0,1);
beta             = (1-gamma)*beta+gamma*a;

betaAll{j}       = beta./rootXX;

actInd(j)        = sj;

j                = j+1;
nActive          = sum(beta~=0);

stop            = abs(func(gamma)-obj)<1e-6*obj ;%& j >.25*maxJ;
obj             = func(gamma);



end
  
if (nActive> nVar)

    beta        = beta/(1-gamma);
    beta(sj)    = 0;
    
end    
   
rawBeta         = beta;

beta            = beta./rootXX;      
R2              = 1-(YY-2*beta'*(n*rootXX.*XY)+beta'*XX*beta)/(YY-(mY^2/n));
end




function out = objective(X,Y,B,bNew,gamma)

BNew         = (1-gamma)*B+gamma*bNew;
out          = (Y-X*BNew)'*(Y-X*BNew);


end












