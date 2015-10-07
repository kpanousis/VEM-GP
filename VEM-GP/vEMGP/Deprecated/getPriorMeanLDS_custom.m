function Mu = getPriorMeanLDS_custom(params,T,varargin)
%
%
%

seq  = [];
A    = params.A;
x0   = params.x0;
xAdd = [];

assignopts(who,varargin);

xDim = size(params.A,1);

Mu = zeros(xDim,T);
Mu(:,1) = x0;

if params.notes.useB
  if ~isempty(seq)
      disp('useB');
%       size(seq.u)
%       size(params.gp.B)
    Mu(:,1) = Mu(:,1)+params.B*seq.u(:,1);
  else
    error('params.model.notes.useB == true   but no seq given!')
  end
end

if ~isempty(xAdd)
  Mu(:,1) = Mu(:,1)+xAdd(:,1);
end


for t=2:T
  Mu(:,t) = A*Mu(:,t-1);

  if ~isempty(seq) && params.notes.useB
    Mu(:,t) = Mu(:,t)+params.B*seq.u(:,t);
  end
  if ~isempty(xAdd)
    Mu(:,t) = Mu(:,t)+xAdd(:,t);
  end
  
end