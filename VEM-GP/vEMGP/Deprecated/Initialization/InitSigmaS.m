function [sigma_s]=InitSigmaS( hypers,times,times_prime,clickTimes,expect )
%INITSIGMAS Initialize sigma_s^2
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%9 June 2015
%Inputs: hypers: The hyperparameters of the GP. Should have size 6.
%        seq: the struct containing the posterior mu,sigma as well as the
%             click times, click signs and considered times
%        trial: the currrent trial we are looking at
%Output: sigma_s: the calculated initial value of sigma_s^2


%% Hyperparameter Count Check

if (numel(hypers)~=6)
    disp('Something wrong with the number of hyperparameters');
end

%% Calculate sigma_s
sigma_s=(1/numel(clickTimes))*EqPerTimeStepCalc(times,times_prime,clickTimes,hypers,expect,true);

end

