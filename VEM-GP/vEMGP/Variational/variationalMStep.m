function [ NEWparams,sessionDataOut] = variationalMStep(params,sessionDataIn,nParams)
%VARIATIONALMSTEP Main Function for Implementing the Variational MStep
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%17 June 2015

%Inputs:  params: the struct containing the parameters of the GLM and the
%          GP
%         sessionDataIn: Struct that contains the necessary info for each
%         session in order to perform the E Step. This info includes C and
%         d of the GLM and seq, where
%           seq: struct containing the y,times, clickTimes etc

%Outputs: NEWparams: the updated struct
%         sesionDataOut: the updated sessions struct

%% Some Messages
fprintf('-------------------------------------------\n')
disp('Variational M Step');
fprintf('-------------------------------------------\n')

%% Input Checking
assert(isfield(params,'gp'),'Current Parameters Do Not Contain GP field!');
assert(numel(params.model.hypers)==6,'Wrong Number of Hyperparameters');
assert(isfield(sessionDataIn(1).seq,'posterior'),' No Posterior Field on the Seq struct. Please Perform E-Step First!');

%% Assign Variables
sessionDataOut=sessionDataIn;
curr_hypers=params.model.hypers;

%% Update for the GLM Variables C and d
params.model.A=exp(curr_hypers(4)*(1/params.model.samplingRate));
NEWparams=params;

%For each session, get the parameters C and d, optimize and update the
%sessionDataOut Struct
%The PLDSMStepObservationCost is a function of the popSpike_dyn package
for i=1:numel(sessionDataOut)
    NEWparams.model.C=sessionDataOut(i).C;
    NEWparams.model.d=sessionDataOut(i).d;
    NEWparams = PLDSMStepObservation(NEWparams,sessionDataOut(i).seq);
    NEWparams.model.C
    sessionDataOut(i).C=NEWparams.model.C;
    sessionDataOut(i).d=NEWparams.model.d;
end

%% Options for minFunc
%Change some of the default values for the options for minFunc to get some
%messages for the optimization
options = params.opts.algorithmic.MStepObservation.minFuncOptions;
options.display='full';


%% GP Parameters
%%% Maximize the F_GP using minFunc and get the new Hypers
if (nParams<6)
    newHypers=minFunc(@F_GP_AndDerivsTester,curr_hypers(1:nParams),options,sessionDataOut,curr_hypers(nParams+1:6))
    NEWparams.model.hypers=[newHypers; curr_hypers(nParams+1:6)];
else
    newHypers=minFunc(@F_GP_AndDerivsTester,curr_hypers(1:nParams),options,sessionDataOut)
    NEWparams.model.hypers=newHypers;
end

% In case we optimize more than 3 parameters (lambda, phi, log_tau_phi) we
% need to recalculate the input.
if (nParams>3)
    fprintf('--------------------------------------------------------------------------\n')
    disp('Recalculating the input U, due to changes in lambda and/or other params...');
    fprintf('--------------------------------------------------------------------------\n');
    for sessionIdx=1:numel(sessionDataOut)
        for trialIdx=1:numel(sessionDataOut(sessionIdx).seq)
            sessionDataOut(sessionIdx).seq(trialIdx).u=calculate_Input(sessionDataOut(sessionIdx).seq(trialIdx),NEWparams.model.hypers);
        end
    end
end

end