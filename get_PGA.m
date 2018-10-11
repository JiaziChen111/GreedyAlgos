function [beta  rawBeta R2 actInd betaAll] = get_PGA(YY,XY, XX, mY, n, shrink, maxJ, standardizeFlag, nVar)
% estimates the RGA with constraints: fully vectorized. For best
% performance, use XX with diagonal equal to 1.


if nargin <6
    
   shrink = .1;
    
end

if nargin <7
    
   maxJ = 500;
    
end

if nargin < 8
   
    standardizeFlag =0;
    
end

if nargin <9
    
    
    nVar = size(XX,1);
    
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

beta             = zeros(K,1);
j                = 1;
nActive          = 0;
actInd           = zeros(floor(maxJ)+2,1);

while (j<=maxJ) && (nActive<= nVar)

A                = XY-D*beta;
[ma sj]          = max(abs(A));

a                = zeros(K,1);
a(sj)            = A(sj)/D(sj,sj);
beta             = beta+shrink*a;

betaAll{j}       = beta./rootXX;


actInd(j)        = sj;

j                = j+1;
nActive          = sum(beta~=0);
end
  
if (nActive> nVar)

    beta(sj)    = 0;
    
end    
    
rawBeta         = beta;
beta            = beta./rootXX;
%R2              = 1-(YY-2*beta'*(n*rootXX.*XY)+beta'*XX*beta)/(YY-(mY^2/n));
R2              = (beta'*XX*beta)/(YY-(mY^2/n));
       
end















