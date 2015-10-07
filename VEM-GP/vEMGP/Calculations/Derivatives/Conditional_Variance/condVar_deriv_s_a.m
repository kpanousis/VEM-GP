function [ deriv_sigma_a ] = condVar_deriv_s_a(hypers,t,t_prime)
%CONDVAR_DERIV_S_A Derivative of conditional variance with respect to log
%sigma alpha^2
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience
%University College London
%15 June 2015
%This function is used by F_GP_deriv and calculates the derivative for only
%one time step.
%Inputs: hypers: the (6,1) hyperaparameter vector
%        t: the current time t_i, a scalar
%        t_prime: the t_{i-1}
%Output: deriv_sigma_a: the derivative of the conditional variance
%evaluated at time t

%% Input Checking
if (numel(hypers)~=6)
    error('Something wrong with the number of hyperparameters');
end

%% Hypers
sigma_a_sq=exp(hypers(2));
lambda=hypers(4);

%% Final Calculation
if (lambda==0)
    deriv_sigma_a=sigma_a_sq*(t-t_prime);
else
    deriv_sigma_a=sigma_a_sq*expm1(2*lambda*(t-t_prime))/(2*lambda);
end

end

