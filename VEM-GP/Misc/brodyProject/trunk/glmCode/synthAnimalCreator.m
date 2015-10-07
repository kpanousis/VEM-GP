%synthAnimalCreator.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%31 July, 2014
%Script which creates a synthetic animal and experimental data from it for
%use with the glmBrody fitting routines.

%% Clear the decks
% clear
% close all
% clc

%% Set paths to code:
% addpath('~/SVNrep/brodyProject/trunk')
% addpath('~/SVNrep/brodyProject/trunk/glmCode')  %Where this file should be
% addpath('~/SVNrep/brodyProject/trunk/utils')
% addpath('~/SVNrep/brodyProject/trunk/kernel')
% addpath('~/SVNrep/brodyProject/trunk/fittingTheta')
% addpath('~/SVNrep/brodyProject/trunk/dataLoading')  %Need this because it contains brody_PPCdataDirectory.mat
% addpath('~/SVNrep/brodyProject/trunk/testingTools')
addpath(genpath('~/Documents/MATLAB/minFunc_2012/'))

%% Reset the RNG to the default (start MATLAB) settings:
rng('default')

%% Define animal designator:
animalDesignator = 'synth1';
synthesizingModelNum = 1;
synthesizingModelType = 'glmBrodyLaplaceSelfOnly';
%groundTruthGPHypers = [log(0.0001),log(0.01),log(0.00001),-0.05,0.8,log(0.05)]'; %These are: [log(sigmaI^2), log(sigmaA^2), log(sigmaS^2), lambda, phi, log(tauPhi)]
groundTruthGPHypers = [log(0.0001),log(0.01),log(0.25),-0.05,0.8,log(0.05)]'; %These are: [log(sigmaI^2), log(sigmaA^2), log(sigmaS^2), lambda, phi, log(tauPhi)]
% groundTruthGPHypers(3)=0;
%For line search:
maxNumLineSearchEvals = 40;

%% Will we use momentum in the optimization?
useMomentum = false;
% useMomentum = true;
if (exist('trialsInput','var')==1)
    trials=trialsInput;
else
    disp('Default Trials 100');
    trials=100;
end
if (exist('neuronsInput','var')==1)
    neurons=neuronsInput;
else
    disp('Default Neurons: 1');
    neurons=1;
end
%% Trial creation parameters
trialsPerSessionVec  = [trials]; %Column vector of per-session numbers of trials to synthesize
neuronsPerSessionVec = [neurons];  %Column vector of per-session numbers of neurons

synthesisStruct.lapseRate = 0;

synthesisStruct.baseClickRate = 40; %Hz
synthesisStruct.startTime = 0; %Seconds, global time
synthesisStruct.assumed_state_0_exit_ToCPokeStart = 2.5; %Seconds, ~exponential distribution.
%synthesisStruct.assumedEndStimToEndFix = 0.1; %Seconds
synthesisStruct.assumedCPokeEndToCPokeOut = 0.1; %Seconds
%synthesisStruct.assumedEndFixToSidePoke = 1; %Seconds  %Typ 0.5 to 1.5
%synthesisStruct.assumedSidePokeToStartNextTrial = 0.5; %Seconds
synthesisStruct.assumedCPokeOutToStartNextTrial = 8 ; %seconds


synthesisStruct.totalCPokeStartToCPokeEndTime = 1; %Seconds: time between initiation and end of fixation.  This seems to vary by session, but be fixed for an individual session.  We can perhaps play with this later.
synthesisStruct.startDrawingClicksTimeAfterCPokeStart = 0.5; %seconds  %Time the animal must center-fixate before any clicks are played.
% % synthesisStruct.trialTimeMin = 0.1; %Seconds
% % synthesisStruct.trialTimeMax = 1.2; %Seconds
% % synthesisStruct.lapseRate = 0.2;  %Proportion of time the animal chooses by
% % %flipping a coin, rather than using the accumulator value to solve the task

%Above: Hanks et al.: "For neural recording sessions, the duration of
%stimulus presentation was varied randomly from trial to trial,
%with durations ranging from 0.1 to 1.2 s."  However, not sure how to use
%this.

%Notes: 
%At least in the session which produced file phys_data_1874.mat, the
%interval cpoke_start to cpoke_end was ALWAYS 1 second.  Surveying more
%files, it appears that this parameter varies for different animals, but is
%the same across all sessions with that animal.
%
% So what it looks like is that cpoke_end - cpoke_start is fixed, the
% clicks are drawn over this interval, and stim_start is set equal to the
% time of the first click.  Should be relatively easy to generate.

synthesisStruct.basalFiringRate = 10; %Hz: spikes per second: sets the generating d weight values.
%For generating the appropriate self-weights:
% candidateFunction = 0.05 * (0.01 * normpdf(-5:0.5:5) + 0.01 * normpdf(-1:0.5:9) + 0.01 * normpdf(-3:0.5:7)  - 2.5 * poisspdf(0:20,0.5)) ;
% candidateFunction = fliplr(candidateFunction(1:end-1));
% candidateFunction = candidateFunction';

%Moved definitions of these parameters up here:
if (exist('samplingInput','var')==1)
    samplingRate=samplingInput;
else
    disp('Default Sampling Rate 100Hz');
    samplingRate = 100; %Hz
end
relevantSelfHistory_glmBrodyLaplace = 0.020; %Seconds

candidateFunction = candidateFunctionCreator('positiveBump',relevantSelfHistory_glmBrodyLaplace,samplingRate);


synthesisStruct.selfWeightTemplate = candidateFunction;
clear candidateFunction;

synthesisStruct.wRegularizationWeight = 0.1;




%% Set some key model parameters
bupLength = 0.003; %Seconds
%samplingRate = 1000; %Hz
bleedTime = 0.300; %Time to construct representations before and after first and last bup of trial.  May wish to expand at end to try to understand time of decision, with reference to the event time of decision?
maxLagConsidered = 0.000; %Seconds: 400 ms (what I had before 24 April) is probably too long: 300 ms would be better; we're trying to avoid movement-tuned activity at the end of the trial.


%relevantSelfHistory_GLM = 0.100; %Seconds
%relevantSelfHistory_glmBrodyLaplace = 0.020; %Seconds
relevantStimHistory = 0.300; %Seconds

%Basis specification:
numPoissonCompsSelfBasis = 10;

%In lassoGLM
alpha = 0.01;  %Weights L2 regularlization more heavily than L1.

%Assume a particular basal firing rate in the prior:
assumedBasalFiringRate_prior = 10; %Hz

%% Create the appropriate fittingParameters structures
%For speed, set fewer than full number of trials
%maxNumTrialsForFittingGLMs = 20;
trialNumbersFor_glmBrodyModelFitting = inf;%[20,inf]; %For glmBrodyLaplace model.
%numberOfWorkers_glmBrody = 4;

fittingParameters.alpha = alpha;
fittingParameters.saveLocalRepresentations = false;
fittingParameters.maxCycles = 10;
fittingParameters.maxLatentNewtonIters = 40;
%fittingParameters.numberOfParallelWorkers = numberOfWorkers_glmBrody;
fittingParameters.maxTrialNumber = 30; 
%fittingParameters.gamma_aStep = 0.5; %Goes into the naive newton step
fittingParameters.gamma_thetaStep = 1;
% Momentum
if useMomentum
    fittingParameters.massMomentum = 0.1;
end

%For line search:
fittingParameters.maxNumEvalsLineSearch = maxNumLineSearchEvals;
fittingParameters.alphaNewtonLS = 1e-4;
fittingParameters.objectiveFunctionIdx = 2;
fittingParameters.minimumAllowedSlope = 0.001;

%% Create the modelGlobalParameters structures
%For GLM Brody model with Laplace approximation:
modelGlobalParameters_glmBrodyLaplace.samplingRate = samplingRate;
modelGlobalParameters_glmBrodyLaplace.bleedTime = bleedTime;  %used only for the fitting of W and d using lassoglm.
modelGlobalParameters_glmBrodyLaplace.relevantSelfHistory = relevantSelfHistory_glmBrodyLaplace;
%modelGlobalParameters_glmBrodyLaplace.gpHypers = [log(1e-4),log(1e-4),log(0.36),0,0.2,log(0.1)]';
%modelGlobalParameters_glmBrodyLaplace.gpHypers = [log(1e-2),log(1e-2),log(0.36),0,0.2,log(0.1)]';
fprintf('Using hyperparameter values which should give approximately equal weight to initial variance, \n drift, and click incorporation in each trial, under the assumption that it lasts approximately 1/2 second.\n')
modelGlobalParameters_glmBrodyLaplace.gpHypers = [log(1),log(4),log(0.1),0,0.2,log(0.1)]';
%Above: gpHypers entries are log(sigmaI^2), log(sigmaA^2), log(sigmaS^2), lambda, phi, log(tauPhi).
modelGlobalParameters_glmBrodyLaplace.gpHypersPriorMean = modelGlobalParameters_glmBrodyLaplace.gpHypers;
modelGlobalParameters_glmBrodyLaplace.gpHypersPriorVariance = eye(numel(modelGlobalParameters_glmBrodyLaplace.gpHypers));  %Diagonal, unit covariance; should change the variances in the future.
modelGlobalParameters_glmBrodyLaplace.glmHypersPriorMean = [0, log(assumedBasalFiringRate_prior/samplingRate), 0]';  %Sets the mean for C_i, d_i, W_{i,j}
modelGlobalParameters_glmBrodyLaplace.glmHypersPriorVariance = [0.25;0.25;0.25];
modelGlobalParameters_glmBrodyLaplace.maxLagConsidered = maxLagConsidered;

%% Create the basis creation structures
selfBasisStruct.nPoissonCompsBasis = numPoissonCompsSelfBasis;

%% Create the glmBrody class object

synthAnimalStruct = glmBrody(animalDesignator);

%% Add the desired model

synthAnimalStruct.addModel('glmBrodyLaplaceSelfOnly')


%% Set the model global parameters
synthAnimalStruct.modelGlobalParameters{1} = modelGlobalParameters_glmBrodyLaplace;
synthAnimalStruct.modelGlobalParameters{1}.gpHypers = groundTruthGPHypers; 
%% Create the required bases for the model
synthAnimalStruct.createBasis(1,'hybridPoissonGaussian','self',selfBasisStruct)

%% Synthesize the data
synthAnimalStruct.synthesizeData(synthesizingModelNum,trialsPerSessionVec,neuronsPerSessionVec,synthesisStruct);

%% Plot the initial self-weights and save other components

startingValues.gpHypers = synthAnimalStruct.modelGlobalParameters{1}.gpHypers;
startingValues.W{1} = synthAnimalStruct.modelSessionParameters{1,1}.W;
startingValues.d{1} = synthAnimalStruct.modelSessionParameters{1,1}.d;
startingValues.C{1} = synthAnimalStruct.modelSessionParameters{1,1}.C;

%startingValues.W{2} = synthAnimalStruct.modelSessionParameters{2,1}.W;
%startingValues.d{2} = synthAnimalStruct.modelSessionParameters{2,1}.d;
%startingValues.C{2} = synthAnimalStruct.modelSessionParameters{2,1}.C;

% for sessionIdx = 1:numel(trialsPerSessionVec)
%     figureNum = sessionIdx;
%     synthAnimalStruct.plotSessionWeights(1,sessionIdx,figureNum);
% end


%% Do the fitting
%synthAnimalStruct.fitModel(1,fittingParameters);

%% Plot the final self-weights
%for sessionIdx = 1:numel(trialsPerSessionVec)
%    figureNum = sessionIdx + numel(trialsPerSessionVec);
%    synthAnimalStruct.plotSessionWeights(1,sessionIdx,figureNum);
%end


%% Save out a final state
%save('finalSynthAnimalState')
