%loadRealData.m
%Konstantinos Panagiotis Panousis 
%Gatsby Computational Neuroscience Unit
%University College London
%28 AUgust 2015
%Script that loads the brody Data
%Not yet for full use, maybe find target experimental data that correspond
%to some restrictions
%Many parts from code from Dr. Thomas Desautels
%-No permission to publish or redistribute-
clear;
clc;

%The data are not publicly available either way.
addpath '~/Desktop/Matlab_Gatsby/Matlab/brodyData/brodyData'
addpath '~/Desktop/Matlab_Gatsby/Matlab/brodyData/brodyData/phys_data'

animalIdx=5;

%% Set some key model parameters
bupLength = 0.003; %Seconds
samplingRate = 1000; %Hz
bleedTime = 0.300; %Time to construct representations before and after first and last bup of trial.  May wish to expand at end to try to understand time of decision, with reference to the event time of decision?
maxLagConsidered = 0.000; %Seconds: 400 ms (what I had before 24 April) is probably too long: 300 ms would be better; we're trying to avoid movement-tuned activity at the end of the trial.


relevantSelfHistory_GLM = 0.100; %Seconds
relevantSelfHistory_glmBrodyLaplace = 0.020; %Seconds
relevantStimHistory = 0.300; %Seconds

%Assume a particular basal firing rate in the prior:
assumedBasalFiringRate = 10; %Hz

%Basis specification:
numPoissonCompsSelfBasis = 10;
numPoissonCompsStimBasis = 20;

%% Create the basis creation structures
selfBasisStruct.nPoissonCompsBasis = numPoissonCompsSelfBasis;
stimBasisStruct.nPoissonCompsBasis = numPoissonCompsStimBasis;


%For GLM Brody model with Laplace approximation:
modelGlobalParameters_glmBrodyLaplace.samplingRate = samplingRate;
modelGlobalParameters_glmBrodyLaplace.bleedTime = bleedTime;  %used only for the fitting of W and d using lassoglm.
modelGlobalParameters_glmBrodyLaplace.relevantSelfHistory = relevantSelfHistory_glmBrodyLaplace;
%modelGlobalParameters_glmBrodyLaplace.gpHypers = [log(1e-4),log(1e-4),log(0.36),0,0.2,log(0.1)]';  
modelGlobalParameters_glmBrodyLaplace.gpHypers = [log(1e-2),log(1e-2),log(0.36),0,0.8,log(0.05)]';  
%Above: gpHypers entries are log(sigmaI^2), log(sigmaA^2), log(sigmaS^2), lambda, phi, log(tauPhi).
modelGlobalParameters_glmBrodyLaplace.gpHypersPriorMean = modelGlobalParameters_glmBrodyLaplace.gpHypers;
modelGlobalParameters_glmBrodyLaplace.gpHypersPriorVariance = eye(numel(modelGlobalParameters_glmBrodyLaplace.gpHypers));  %Diagonal, unit covariance; should change the variances in the future.
modelGlobalParameters_glmBrodyLaplace.glmHypersPriorMean = [0, log(assumedBasalFiringRate/samplingRate), 0]';  %Sets the mean for C_i, d_i, W_{i,j}
modelGlobalParameters_glmBrodyLaplace.glmHypersPriorVariance = [1;1;1];
modelGlobalParameters_glmBrodyLaplace.maxLagConsidered = maxLagConsidered;  


%% Instantiate the structure
loadingStructure    = load('brody_PPCdataDirectory','cellsCellArray','animalsList','sessionNumbersCellArray','nCellsLoaded','animalNumCells');
%sessionsToLoad      = loadingStructure.sessionNumbersCellArray{animalIdx};
%cellsToLoad         = loadingStructure.cellsCellArray{animalIdx};
sessionsToLoad      = loadingStructure.sessionNumbersCellArray{animalIdx}(7);
cellsToLoad         = loadingStructure.cellsCellArray{animalIdx}(7);
animalStruct        = glmBrody(loadingStructure.animalsList{animalIdx},  sessionsToLoad,   cellsToLoad,  'phys_data_');


animalStruct.addModel('glmBrodyLaplaceSelfOnly');
animalStruct.modelGlobalParameters{1} = modelGlobalParameters_glmBrodyLaplace;  
animalStruct.createBasis(1, 'hybridPoissonGaussian', 'self', selfBasisStruct)

nModel = 1;

nParamsToOptimize=4;

%% For the diary/other transcripts
% fprintf('Now fitting run %d.\n\n',run)

%% Data Transform
%Now fix the data in a way that is compatible with MBS code
%so use the data_transform.m function
fprintf('\n-------------------------------------------\n')
disp('Data Transform');
fprintf('-------------------------------------------\n')

%hyperparameters random initialization
%all negative except phi parameter
log_sigma_i_sq=10*rand-10;
log_sigma_a_sq=5*rand-5;
log_sigma_s_sq=5*rand-5;
lambda=4*rand-4;
phi=rand(1,1);
log_tauPhi=6.*rand-6;
hyp=[log_sigma_i_sq;log_sigma_a_sq;log_sigma_s_sq;lambda;phi];
hyp=hyp(1:nParamsToOptimize);

hypersPrior=[hyp;animalStruct.modelGlobalParameters{1,1}.gpHypers(nParamsToOptimize+1:6)]
c_glm=zeros(numel(cellsToLoad{1}),1);
useWeights=false;
[sessionData,params]=data_transform(animalStruct,nModel,hypersPrior,c_glm,useWeights);
%clearvars -except outputFileNamePrefix terminalStateDirectory runnum hypersPrior sessionData params nParamsToOptimize synthAnimalStruct;
% groundTruthHypers=synthAnimalStruct.modelGlobalParameters{1,1}.gpHypers;

%print out some statistics of the artificial data
fprintf('\n-------------------------------------------\n')
disp('Statistics Session 1');
fprintf('-------------------------------------------\n')
fprintf('Max spike count:    %i \n', max(vec([sessionData(1).seq.y])))
fprintf('Mean spike count:   %d \n', mean(vec([sessionData(1).seq.y])))
fprintf('Freq non-zero bin:  %d \n', mean(vec([sessionData(1).seq.y])>0.5))
pause(1);

%% Variational EM
fprintf('\n-------------------------------------------\n')
disp('Variational Expectation Maximization');
fprintf('-------------------------------------------\n')

maxIter=100;           
% setting max no of EM iterations
params.opts.algorithmic.EMIterations.maxIter     = maxIter;	
% setting max CPU time for EM to inf
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;		

tic; [params_vEM, sessionData_vEM, varBound_vEM, EStepTimes, MStepTimes]=vEM(params,sessionData,nParamsToOptimize); toc;

maxIter=200;
params=[];
params_popSpike = PLDSInitialize(sessionData(1).seq,1,'params',params);
params_popSpike.model.C=c_glm;
params_popSpike.model.B=1;
params_popSpike.model.inferenceHandle = @PLDSVariationalInference;
params_popSpike.opts.algorithmic.EMIterations.maxIter     = maxIter;	
params_popSpike.model.notes.useS=useWeights;
for i=1:numel(sessionData)
    sessionData(i).C=params_popSpike.model.C;
    sessionData(i).d=params_popSpike.model.d;
end
tic; [params_popSpike, sessionData_PopSpikeEM,varBound_PopSpikeEM,EStepTimes,MStepTimes]=PopSpikeEM_SessionsMod(params_popSpike,sessionData); toc;

params=[];
params_popSpikeLaPlace = PLDSInitialize(sessionData(1).seq,1,'params',params);
params_popSpikeLaPlace.model.C=c_glm;
params_popSpikeLaPlace.model.B=1;
params_popSpikeLaPlace.opts.algorithmic.EMIterations.maxIter     = maxIter;	
params_popSpikeLaPlace.model.inferenceHandle = @PLDSLaplaceInference;
params_popSpikeLaPlace.model.notes.useS=useWeights;


for i=1:numel(sessionData)
    sessionData(i).C=params_popSpikeLaPlace.model.C;
    sessionData(i).d=params_popSpikeLaPlace.model.d;
end
tic; [params_popSpikeLaPlace, sessionData_PopSpikeEM_LaPlace,varBound_PopSpikeEM_LaPlace,EStepTimes,MStepTimes]=PopSpikeEM_SessionsMod(params_popSpikeLaPlace,sessionData); toc;