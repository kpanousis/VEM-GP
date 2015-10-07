function [ params ] = ParamsGPInit(curr_params,seq)
%PARAMSGPINIT Initialize the parameters of the through the calculations
%in "Variational Gaussian Process EM for Neural Accumulators"
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%9 June 2015
%Inputs: curr_params: the current set of params, from which we need the gp
%        hypers struct
%        seq: Containing the clickTimes, clickSigns, etc
%Output params: the curr_params with added fields for the initial values of
%       the log sigmas

%% Input checking
if (~isfield(curr_params,'gp'))
    error('The curr_params do not have the GP field!');
end

%% Get the hypers of the GP
hypers=curr_params.model.hypers;

if (numel(hypers)~=6)
    error('Wrong number of hyperparameters!');
end

%% Assign output

params=curr_params;

%% Initialize log sigma_i^2,log sigma_s^2 and log sigma_a^2 and add fields
% sum_sigma_i=0;
sum_sigma_s=0;
sum_sigma_a=0;

for trial=1:numel(seq);
    %% Get the needed variables per trial
    t=seq(trial).times;
    t_prime=t(1:end-1);
    t=t(2:end);
    clickTimes=seq(trial).clickTimes;
    clickSigns=seq(trial).clickSigns;
    mu=seq(trial).posterior.xsm;
    Vsm=seq(trial).posterior.Vsm;
    VVsm=seq(trial).posterior.VVsm;
    %% Calculate the Eq term
    [ciVals,ciDerivs]=cAndDerivs(clickTimes,0,exp(params.model.hypers(6)),params.model.hypers(5));
    expect=EqAndDerivatives(t,t_prime,ciVals,ciDerivs,clickTimes,clickSigns,mu',Vsm,VVsm,hypers);
    
    %% Initialize
%     sum_sigma_i=sum_sigma_i+log(InitSigmaI(hypers,expect',t,clickTimes));
    sum_sigma_s=sum_sigma_s+InitSigmaS(hypers,t,t_prime,clickTimes,expect');
    sum_sigma_a=sum_sigma_a+InitSigmaA(hypers,t,t_prime,clickTimes,expect');
end
N=numel(seq);
params.gp.hypersPrior(2:3)=[log(sum_sigma_a/N);log(sum_sigma_s/N)];
params.gp.hypersPrior(2:3)

