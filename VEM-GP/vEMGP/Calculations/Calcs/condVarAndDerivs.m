function [ condVar,derivs ] = condVarAndDerivs( hypers,times,times_prime,ciVals,ciDerivs,clickTimes)
%CONDVARANDDERIVS Calculate and return the conditional variance 
%Konstantinos Panagiotis Panousis
%Gatbsby Computational Neuroscience Unit
%University College London
%16 June 2015
%Inputs: hypers: the vector of hyperparameters (6x1)
%        times: the times to consider
%        clickTimes: the times of the clicks
%Outputs: condVar: the conditional variance calculated at each timestep so
%         size(condVar)=numel(t),1
%        Derivs: the derivatives of the conditional variance with respect
%        to the 5 hyperparameters (numel(t)xnumel(hypers))

%% Input Checking
assert(numel(hypers)==6,'Something Wrong With the Number of Hyperparameters');

%% Initialize needed variables
condVar=zeros(numel(times),1);
derivs=zeros(numel(times),numel(hypers));

%% Positions of the hyperparameters for the derivatives
loc_sigmai=1;
loc_sigmaa=2;
loc_sigmas=3;
loc_lambdaVal=4;
loc_phi=5;
loc_logTauPhi=6;

%% Calculate conditional variance and its derivatives for each time step
for t=1:numel(times)
    
    condVar(t)=condVar_calc(hypers,times(t),times_prime(t),ciVals,clickTimes);
    derivs(t,loc_sigmai)=condVar_deriv_s_i(hypers,times(t),times_prime(t));
    derivs(t,loc_sigmaa)=condVar_deriv_s_a(hypers,times(t),times_prime(t));
    derivs(t,loc_sigmas)=condVar_deriv_s_s(hypers,times(t),times_prime(t),ciVals,clickTimes);
    derivs(t,loc_lambdaVal)=condVar_deriv_lambda(hypers,times(t),times_prime(t),ciVals,clickTimes);
    derivs(t,loc_phi)=condVar_deriv_phi(hypers,times(t),times_prime(t),ciVals,ciDerivs,clickTimes);
    derivs(t,loc_logTauPhi)=condVar_deriv_logTauPhi(hypers,times(t),times_prime(t),ciVals,ciDerivs,clickTimes);
    
    
end

% End function
end

