function [ deriv_lambda ] = condVar_deriv_lambda( hypers,t,t_prime,ciVals,clickTimes)
%CONDVAR_DERIV_LAMBDA Derivative of the conditional variance with respect
%to the parameter lambda
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%15 June 2015
%This function is used by F_GP_deriv and calculates the derivative for only
%one time step.
%Inputs: hypers: the (6,1) hyperaparameter vector
%        t: the current time t_i, a scalar
%        t_prime: the t_{i-1}
%        ClickTimes: the times of the clicks
%Output: deriv_lambda: the derivative of the conditional variance
%evaluated at time t

%% Input Checking
if (numel(hypers)~=6)
    error('Something wrong with the number of hyperparameters');
end

%% Get Hypers
sigma_i_sq=exp(hypers(1));
sigma_a_sq=exp(hypers(2));
sigma_s_sq=exp(hypers(3));
lambda=hypers(4);


%% Individual Derivatives

%% Sigma_i_sq
% derivative with respect to log sigma_i_sq
deriv_sigma_i_term=2*t*sigma_i_sq*exp(2*lambda*t)*(t_prime==0);

%% Sigma_a term
% Derivative with repsect to log sigma_a_sq
if (lambda==0)
    diff=t-t_prime;
    deriv_sigma_a_term=sigma_a_sq*diff^2;
else
    deriv_sigma_a_term=sigma_a_sq*(exp(2*lambda*(t-t_prime))*(2*lambda*(t-t_prime)-1)+1)/(2*lambda^2);
end

%% Sigma_s_sq term
% derivative with respect to log sigma_s_sq
% Will need a small modification for clickTime in [t',t) and not
% clickTime approximately equal to t'
value=0;
for i=1:numel(clickTimes)
    
    if (clickTimes(i)>=t_prime && clickTimes(i)<t)
        value=value+ciVals(i)^2*exp(2*lambda*(t-clickTimes(i)))*(t-clickTimes(i));
    end
end
deriv_sigma_s_term=2*sigma_s_sq*value;


%% Final Calculation
deriv_lambda=deriv_sigma_i_term+deriv_sigma_s_term+deriv_sigma_a_term;

end

