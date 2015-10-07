function [ VarInfparams] = setParameters( params,seq )
%SETPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

%% Add the times and the gp hypers to calculate the conditional variacne in VAriationalInferenceDualCost
% Konstantinos Panagiotis Panousis
%3 July 2015
%this thing right here seems a bit slow, cause for each trial we must
%calculate the conditional variance,
%at least calculate the conditional variance and use it also for the M
%Step
disp('Setting Parameters');
VarInfparams=params;
for tr=1:numel(seq)
    VarInfparams.times=seq(tr).times;
    VarInfparams.clickTimes=seq(tr).clickTimes;
    VarInfparams.clickSigns=seq(tr).clickSigns;
    VarInfparams.hypers=params.gp.hypersPrior;
    VarInfparams.Q0=exp(VarInfparams.hypers(1));
    VarInfparams.x0=1e-6;
    VarInfparams.B=1;
    VarInfparams.A=exp(VarInfparams.hypers(4)*(1/params.gp.samplingRate));
    VarInfparams.model.notes.useB=true;
    VarInfparams.Q=condVarAndDerivs(VarInfparams.hypers,VarInfparams.times(2:end),VarInfparams.times(1:end-1),VarInfparams.clickTimes);
end

end

