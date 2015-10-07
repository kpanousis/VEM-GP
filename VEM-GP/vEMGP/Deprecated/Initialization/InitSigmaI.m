function [sigma_i]=InitSigmaI(hypers,Eq,times,clickTimes)
%INITSIGMAI Initialize sigma_i^2
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%9 June 2015
%Inputs: hypers: The hyperparameters of the GP. Should have size 6.
%        seq: the struct containing the posterior mu,sigma as well as the
%             click times, click signs and considered times
%        trial: the currrent trial we are looking at
%Output: sigma_s: the calculated initial value of sigma_i^2

%% Get Hypers

if (numel(hypers)==6)
    % log_sigma_i=hypers(1);
    log_sigma_a=hypers(2);
    log_sigma_s=hypers(3);
    lambda=hypers(4);
    phi=hypers(5);
    log_tauPhi=hypers(6);
else
    disp('Something wrong with the number of hyperparameters');
end


%% Numerator Individual Terms

%the second term of the numerator containing \sigma_a^2 term
if (lambda~=0)
    s_a_term=exp(log_sigma_a)*(1-exp(2*lambda*times))./(2*lambda);
else
    s_a_term=exp(log_sigma_a)*(times-[0 times(1:numel(times)-1)]);
end

%the final term with the \sigma_s^2 term
sum_ci=0;
cZero=0;
times_prime=[0 times(1:numel(times)-1)];
for j=1:numel(clickTimes)
    if (clickTimes(j)>times(2))
        break;
    end
    if (clickTimes(j)>=times_prime(j) && clickTimes(j)<times(j))
        ciVals=cAndDerivs(clickTimes(j),cZero,exp(log_tauPhi),phi);
        sum_ci=sum_ci+ ciVals.^2*exp(2*lambda*(times-clickTimes(j)));
    end
end
s_s_term=exp(log_sigma_s)*sum_ci;


%% Numerator Final Calculation
numerator=Eq-s_s_term-s_a_term;

%% Denominator Final Calculation
denom=exp(2*lambda*times);

%% Final Calculation of \sigma_i^2

sigma_i=numerator/denom;

end

