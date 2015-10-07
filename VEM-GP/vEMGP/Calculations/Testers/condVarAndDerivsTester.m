function [ f,derivs] = condVarAndDerivsTester( chypers,times,clickTimes,entry, varargin )
%CONDVARANDDERIVSTESTER Wrapper to validate the derivatives of the
%conditional variance
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%18 June 2015
%Inputs:  chypers: hyperparameters to check the derivatives (number can
%         vary but the rest must be in varargin)
%         times: the times to consider
%         clickTimes: the times of the clicks
%         entry: the entry to return
%         varargin: if numel(chypers)<6, the rest of the hyperparameters
%Outputs: f: the conditional Variance at time entry
%         derivs: the derivs of the conditional variance

%% Check if we have extra inputs and concatenate
if (length(varargin)>=1)
    hypers=[chypers; varargin{1}];
else
    hypers=chypers;
end

%% Input Checking
assert(numel(hypers)==6,'Something Wrong With the Number of Hyperparameters');

%% Variable Initilization
times_prime=times(1:end-1);
times=times(2:end);
condVar=zeros(numel(times),1);
derivs=zeros(numel(times),numel(chypers));

%% Position of the derivatives
loc_sigmai=1;
loc_sigmaa=2;
loc_sigmas=3;
loc_lambdaVal=4;
loc_phi=5;
loc_logTauPhi=6;

%% Calculate the conditional Varianace and its derivatives for each time step
for t=1:numel(times)
    
    condVar(t)=condVar_calc(hypers,times(t),times_prime(t),clickTimes);
    derivs(t,loc_sigmai)=condVar_deriv_s_i(hypers,times(t),times_prime(t));
    derivs(t,loc_sigmaa)=condVar_deriv_s_a(hypers,times(t),times_prime(t));
    derivs(t,loc_sigmas)=condVar_deriv_s_s(hypers,times(t),times_prime(t),clickTimes);
    derivs(t,loc_lambdaVal)=condVar_deriv_lambda(hypers,times(t),times_prime(t),clickTimes);
    derivs(t,loc_phi)=condVar_deriv_phi(hypers,times(t),times_prime(t),clickTimes);
    derivs(t,loc_logTauPhi)=condVar_deriv_logTauPhi(hypers,times(t),times_prime(t),clickTimes);
    
    
end

%% Assign Output
f=condVar;

%% Uncomment for derivative check
% Return the specified by entry variable, row of the vectors
% and the derivs only for the variables in chypers
entryToReturnRow = entry;
f = condVar(entryToReturnRow,1);
assert(numel(f) == 1, 'Somehow messing up number of times to examine!')
derivs=derivs(entryToReturnRow,1:numel(chypers))';

end

