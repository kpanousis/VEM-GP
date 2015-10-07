function Lambda = buildPriorPrecisionMatrixFromLDS_Custom(params,T)
%
% Lambda = buildPrecisionMatrixFromLDS(params,T)
%
% construct the precision matrix of the prior across all time points and
% dimensions, as described in Paninski et al, A new look at state-space
% models for neural data, 2009
%
% c/o L Buesing and J Macke, 01/2014


xDim   = size(params.A,1);
invQ   = pinv(params.Q);
invQ0  = pinv(params.Q0);
AinvQ  = params.A'*invQ;
AinvQA = AinvQ*params.A;


Lambda = sparse(T*xDim,T*xDim);
Lambda(1:xDim,1:xDim) = invQ0;

for t=1:T-1
  xidx = ((t-1)*xDim+1):(t*xDim);
%   size(Lambda(xidx,xidx))
%   size(AinvQA)
  Lambda(xidx,xidx) = Lambda(xidx,xidx)+AinvQA(xidx);
  Lambda(xidx,xidx+xDim) = -AinvQ(xidx);
  Lambda(xidx+xDim,xidx) = -AinvQ(xidx)';
  Lambda(xidx+xDim,xidx+xDim) = Lambda(xidx+xDim,xidx+xDim)+invQ(xidx);
end
Lambda = sparse((Lambda+Lambda')/2);
