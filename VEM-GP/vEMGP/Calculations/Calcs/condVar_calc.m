function [ condVar] = condVar_calc( hypers,t,t_prime,ciVals,clickTimes)
%CONDVAR_CALC Calculate the conditional variance for t and t'
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%15 June 2015
%This function is used by F_GP_deriv and F_GP_calc and calculates the conditional
%variance for only one time step.
%Inputs: hypers: the (6,1) hyperaparameter vector
%        t: the current time t_i, a scalar
%        t_prime: the t_{i-1}
%Outputs: condVar: the conditional variance, calculated at time t


%% Get Hypers
sigma_i_sq=exp(hypers(1));
sigma_a_sq=exp(hypers(2));
sigma_s_sq=exp(hypers(3));
lambda=hypers(4);



%% Calculations/ Individual Terms
%% Sigma _i term
sigma_i_term=sigma_i_sq*exp(2*lambda*t)*(t_prime==0);

%% Sigma_a term
if (lambda==0)
    sigma_a_term=sigma_a_sq*(t-t_prime);
else
    sigma_a_term=sigma_a_sq*(expm1(2*lambda*(t-t_prime)))/(2*lambda);
end

%% Sigma_s term
value=0;
for i=1:numel(clickTimes)
    if (clickTimes(i)>=t_prime && clickTimes(i)<t)
            value=value+ciVals(i)^2*exp(2*lambda*(t-clickTimes(i)));
    end
end
sigma_s_term=sigma_s_sq*value;


%% Final Calculation
condVar=sigma_i_term+sigma_a_term+sigma_s_term;

end

