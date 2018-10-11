function [beta  rawBeta R2 actInd betaAll] = get_RGA(YY,XY, XX, mY, n,lb,ub,...
                                                     maxJ, standardizeFlag, nVar, decayRate)
% estimates the RGA with constraints: fully vectorized. For best
% performance, use XX with diagonal equal to 1.
%copyright Alessio Sancetta

if nargin <6
    
   lb = -Inf;
    
end


if nargin <7
    
    ub = Inf;
    
end

if nargin <8
    
   maxJ = 500;
    
end

if nargin < 9
   
    standardizeFlag =0;
    
end

if nargin <10
    
    
    nVar = size(XX,1);
    
end

if nargin <11
    
    
    decayRate = 1;
    
end



K                 = size(XX,1);

if standardizeFlag==1

      indZero        = (diag(XX)==0);
      rootXX         = sqrt(diag(XX)/n);
      rootXX(indZero) = 1;
      D              = diag(1./rootXX)*(XX/n)*diag(1./rootXX);
      XY             = bsxfun(@rdivide, XY/n, rootXX);      
    
    
else
    
      rootXX         = ones(K,1);
      D              = XX/n;
      XY             = XY/n;      
    
    
    
end

w                = 1./(1:round(maxJ+1)).^decayRate;

beta             = zeros(K,1);
j                = 1;
nActive          = 0;
actInd           = zeros(floor(maxJ)+2,1);
while (j<=maxJ) && (nActive<= nVar)

A                = XY-(1-w(j))*D*beta;
[ma sj]          = max(abs(A));

a                = zeros(K,1);
a(sj)            = max(min(A(sj)/D(sj,sj)/w(j),ub),lb);
beta             = (1-w(j))*beta+a*w(j);
betaAll{j}       = beta./rootXX;

actInd(j)        = sj;

j                = j+1;
nActive          = sum(beta~=0);
end
  
if (nActive> nVar) && j>1

    beta        = beta/(1-w(j-1));
    beta(sj)    = 0;
    
end    
    
rawBeta         = beta;
beta            = beta./rootXX;
%R2              = 1-(YY-2*beta'*(n*rootXX.*XY)+beta'*XX*beta)/(YY-(mY^2/n));
R2              = (beta'*XX*beta)/(YY-(mY^2/n));
      
end




