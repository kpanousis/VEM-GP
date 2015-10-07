function [ deriv_sigma_s ] = condVar_deriv_s_s( hypers,t, t_prime ,ciVals,clickTimes)
%CONDVAR_DERIV_S_S Derivative of conditional variance with respect to log
%sigma_s
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%15 June 2015
%This function is used by F_GP_deriv and calculates the derivative for only
%one time step.
%Inputs: hypers: the (6,1) hyperaparameter vector
%        t: the current time t_i, a scalar
%        t_prime: the t_{i-1}
%Output: deriv_sigma_s: the derivative of the conditional variance
%evaluated at time t

%% Input Checking
if (numel(hypers)~=6)
    error('Something wrong with the number of hyperparameters');
end

%% Hypers
sigma_s_sq=exp(hypers(3));
lambda=hypers(4);

%% Calculation of the sum term
value=0;
for i=1:numel(clickTimes)
    
    if (clickTimes(i)>=t_prime && clickTimes(i)<t)
        
        value=value+ciVals(i)^2*exp(2*lambda*(t-clickTimes(i)));
    end
end

%% Final Calculation
deriv_sigma_s=sigma_s_sq*value;

end

