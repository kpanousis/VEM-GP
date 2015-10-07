function [ deriv ] = condVar_deriv_phi( hypers,t,t_prime,ciVals,ciDerivs,clickTimes)
%CONDVAR_DERIV_PHI Calculates the derivative of the conditional
%variance with respect to parameters phi
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%20 June 2015
%Inputs:  hypers: the hyperparameter vector
%         t: the time t to consider
%         t_prime: the previous time step
%         clickTimes: the times of the clicks
%Outputs: deriv: the deriv of the conditional variance with respect to log
%         tauPhi at time step t

%% Input Checking
assert(numel(hypers)==6,'Expected 6 hyperparameters');

%% Get hypers and location of derivative of phi in ciDerivs
sigma_s_sq=exp(hypers(3));
lambda=hypers(4);

loc_phi=2;

%% Sum Term Calculation
value=0;
for i=1:numel(clickTimes)
    if (clickTimes(i)>=t_prime && clickTimes(i)<t)
        value=value+ciVals(i)*exp(2*lambda*(t-clickTimes(i)))*ciDerivs(i,loc_phi);
    end
end

%% Final Calculation
deriv=2*sigma_s_sq*value;

end

