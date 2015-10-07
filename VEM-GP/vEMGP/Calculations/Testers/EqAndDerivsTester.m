function [f,varargout]=EqAndDerivsTester(chypers,times,clickTimes,clickSigns,mu,Vsm,VVsm,entryToReturnRow,varargin)
%Wrapper validating EqAndDerivs
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%5 June 2015

%Input:
%   chypers: the hyperparameters of the model (can put 1:2,1:3..., but the
%   rest must be in varargin)
%   times: the times to consider
%   clickTimes: the times of the clicks
%   clickSigns: the signs of the clicks
%   mu: the posterior mean
%   Vsm: the posterior variance
%   VVsm: the posterior covariance
%   entryToReturnRow: the row to return
%   varargin: the rest of the hypers
%Output:
%   f: The value of Eq at a specified point in time
%   derivs: the gradient with respect to the hyperparameters

%% Concatenate varargin to hypers
if (length(varargin)>=1)
    hypers=[chypers; varargin{1}];
else
    hypers=chypers;
end

%% Input Checking
assert(numel(hypers)==6,'Something Wrong With the Number of Hyperparameters');

%% Input handling:
t_prime=[times(1:end-1)];
times=times(2:end);

%% Calculate Eq and its derivatives
[f,derivs]=EqAndDerivatives(times,t_prime,clickTimes,clickSigns,mu',Vsm,VVsm,hypers);

%% For derivative Check
% Return the row specified by the entryToReturnRow variable
f = f(entryToReturnRow,1);
assert(numel(f) == 1, 'Somehow messing up number of times to examine!')
assert( all(size(derivs) == [6,1]), 'Somehow messing up the size of derivs!');

derivsOut = zeros(size(derivs));
for hyperIdx = 1:size(derivs,1)
    if ~isempty(derivs{hyperIdx,1})
        derivsOut(hyperIdx) = derivs{hyperIdx,1}(entryToReturnRow);
    end
end

%% Assign varargout
varargout{1}=derivsOut(1:numel(chypers));

end