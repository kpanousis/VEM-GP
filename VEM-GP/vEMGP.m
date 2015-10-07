%vEMGP.m
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%Variational Gaussian Process EM implementation
%8 June 2015

%% Clear & Close Everything
clear;
close all;
clc;

%% Set Path
run('set_path');

%% USER INPUT
%Request user input to read data if file exists or create data if it does not
[data_create,sessionInput,trialsInput,neuronsInput,samplingInput,nParamsToOptimize]=user_input();

%% Load Case
if (~data_create)
        %Create filename var and try to load data
        filename=['vEMGP/Data/Sessions',num2str(sessionInput),'_Trials',num2str(trialsInput,'_%u'),'_',num2str(1/samplingInput*1e3),'ms_Neurons'...
            num2str(neuronsInput,'_%u'),'.mat'];
        fprintf('\n-------------------------------------------\n')
        disp('Data Load');
        fprintf('-------------------------------------------\n')
        %Check existence of the file and load
        if exist(filename,'file')==2
            fprintf('\nLoading %s\n',filename);
            load(filename);
            pause(1);
            %if it does not exist create new file
        else
            data_create=true;
            fprintf('404 File Not Found.\nCreating new Data for %d Sessions,%d Trials, %d sampling Rate(Hz) and %d Neurons'...
                ,sessionInput,trialsInput,samplingInput,neuronsInput);
        end
end

filenameResults=['Sessions',num2str(sessionInput),'_Trials',num2str(trialsInput,'_%u'),'_',num2str(1/samplingInput*1e3),'ms_Neurons'...
        num2str(neuronsInput,'_%u'),'Results.mat'];

%% Create Data
if (data_create)
    %Synthetic Data Creation
    %first create the data using Thomas Desautels 'synthAnimalCreator.m'
    %function
    fprintf('\n-------------------------------------------\n')
    disp('Data Creation');
    fprintf('-------------------------------------------\n')
    
    %SynthAnimalCreator script created by Dr. Thomas Desautels 
    %-no permission for publication or redistribution-
    run('synthAnimalCreator.m');
    
    filename=['vEMGP/Data/Sessions',num2str(sessionInput),'_Trials',num2str(trialsInput,'_%u'),'_',num2str(1/samplingInput*1e3),'ms_Neurons'...
        num2str(neuronsInput,'_%u'),'.mat'];
    
    save(filename,'synthAnimalStruct','trialsPerSessionVec','synthesizingModelNum');
    clearvars -except filenameResults synthAnimalStruct synthesizingModelNum trialsPerSessionVec nParamsToOptimize;
    
end
%%%%%%%%%%%%%%%%%%%%%Data Creation/Load End%%%%%%%%%%%%%%%%%%%%

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

%create the hyper vector
hyp=[log_sigma_i_sq;log_sigma_a_sq;log_sigma_s_sq;lambda;phi];

%only the parameters we want to optimize
hyp=hyp(1:nParamsToOptimize);

%for the rest use the true generating parameters
hypersPrior=[hyp;synthAnimalStruct.modelGlobalParameters{1,1}.gpHypers(nParamsToOptimize+1:6)]

%initialize the c_glm to zero 
c_glm=zeros(numel(synthAnimalStruct.modelSessionParameters{1,1}.C),1);

%call data_transform function
useWeights=true;
[sessionData,params]=data_transform(synthAnimalStruct,synthesizingModelNum,hypersPrior,c_glm,useWeights);

%clear all the variables except the ones we need for vEM
clearvars -except filenameResults hypersPrior sessionData params nParamsToOptimize synthAnimalStruct c_glm useWeights;

%assign the ground truth hypers to workspace for check
groundTruthHypers=synthAnimalStruct.modelGlobalParameters{1,1}.gpHypers;

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

%The Max Number of Iterations
maxIter=100;           
% setting max no of EM iterations
params.opts.algorithmic.EMIterations.maxIter     = maxIter;	
% setting max CPU time for EM to inf
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;		

%Call the main vEM_GP function that implements the VEM_GP algorithm
tic; [params_vEM, sessionData_vEM, varBound_vEM, EStepTimes, MStepTimes]=vEM_GP(params,sessionData,nParamsToOptimize); toc;

%Re-initialize the parameters for calling vEM of Macke et al.
params=[];
params_popSpike = PLDSInitialize(sessionData(1).seq,1,'params',params);
params_popSpike.model.C=c_glm;
params_popSpike.model.B=1;
params_popSpike.model.notes.useS=useWeights;
params_popSpike.model.inferenceHandle = @PLDSVariationalInference;
params_popSpike.opts.algorithmic.EMIterations.maxIter     = maxIter;	

%Set the initial C and d for each session
for i=1:numel(sessionData)
    sessionData(i).C=params_popSpike.model.C;
    sessionData(i).d=params_popSpike.model.d;
end

%Main call to popSpikeEM to perform the vEM algorithm
tic; [params_popSpike, sessionData_PopSpikeEM,varBound_PopSpikeEM,EStepTimes,MStepTimes]=PopSpikeEM_SessionsMod(params_popSpike,sessionData); toc;

%Re-initialize the parameters for calling the lEM algorithm 
params=[];
params_popSpikeLaPlace = PLDSInitialize(sessionData(1).seq,1,'params',params);
params_popSpikeLaPlace.model.C=c_glm;
params_popSpikeLaPlace.model.B=1;
params_popSpike.model.notes.useS=useWeights;
params_popSpikeLaPlace.opts.algorithmic.EMIterations.maxIter     = maxIter;	
params_popSpikeLaPlace.model.inferenceHandle = @PLDSLaplaceInference;

%Set the initial C and d for each session
for i=1:numel(sessionData)
    sessionData(i).C=params_popSpikeLaPlace.model.C;
    sessionData(i).d=params_popSpikeLaPlace.model.d;
end

%Main call to popSpikeEM to perform the lEM algorithm
tic; [params_popSpikeLaPlace, sessionData_PopSpikeEM_LaPlace,varBound_PopSpikeEM_LaPlace,EStepTimes,MStepTimes]=PopSpikeEM_SessionsMod(params_popSpikeLaPlace,sessionData); toc;

%Save the results
save(filenameResults);
