function [ deriv_sigma_i ] = condVar_deriv_s_i( hypers,t,t_prime )
%CONDVAR_DERIV_S_I Derivatives of the conditional variance with respect
%to \log sigma_i^2
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%15 June 2015
%This function is used by F_GP_deriv and calculates the derivative for only
%one time step.
%Inputs: hypers: the (6,1) hyperaparameter vector
%        t: the current time t_i, a scalar
%        t_prime: the t_{i-1}
%Output: deriv_sigma_i: the derivative of the conditional variance
%evaluated at time t

%% Input Checking
if (numel(hypers)~=6)
    error('Something wrong with the number of hyperparameters');
end

%% Hypers
sigma_i_sq=exp(hypers(1));
lambda=hypers(4);

%% Final Calculation
deriv_sigma_i=sigma_i_sq*exp(2*lambda*t)*(t_prime==0);

end

