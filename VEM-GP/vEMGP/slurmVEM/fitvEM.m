function fitvEM(synthesizingModelNum,fileLoadSaveStructure,runnum,gpHypers)
%FITVEM function called by slurmVEM in order to submit jobs to the slurm
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%14 August 2015
%Function which runs the core of vEMGP intended to be used with SLURM.

%% Establish the locations and names of various important files and objects:
modelRootLocation           = fileLoadSaveStructure.modelRootLocation;
% initialFileName             = fileLoadSaveStructure.initialFileName;
% nameOfClassObjInThatFile    = fileLoadSaveStructure.nameOfInitialClassObj;
terminalStateDirectory      = fileLoadSaveStructure.terminalStateDirectory;
%intermediateStateDirectory  = fileLoadSaveStructure.intermediateStateDirectory;
outputFileNamePrefix        = fileLoadSaveStructure.outputFileNamePrefix;
sourceName=fileLoadSaveStructure.sourceName;

run('set_path_slurm');
 
% %% Turn on diary, with a unique name via the current time code.
diaryName = ['diaryvEM_run',num2str(runnum),'_',datestr(now,'HHMMSSddmmyyyy'),'.log'];
diary(fullfile(terminalStateDirectory,diaryName))

%Set parameters we want to optimize and load the data
nParamsToOptimize=4;
synthAnimalStruct=load(sourceName,'synthAnimalStruct');

%% For the diary/other transcripts
% fprintf('Now fitting run %d.\n\n',run)

%% Data Transform
%Now fix the data in a way that is compatible with MBS code
%so use the data_transform.m function
fprintf('\n-------------------------------------------\n')
disp('Data Transform');
fprintf('-------------------------------------------\n')

% Fix the hyperparameters vector
hyp=gpHypers(1:nParamsToOptimize);
hyp=-abs(hyp);
if (nParamsToOptimize>5)
    hyp(5)=abs(hyp(5));
end

%for the rest of the parameters use the true generating parameters
hypersPrior=[hyp;synthAnimalStruct.modelGlobalParameters{1,1}.gpHypers(nParamsToOptimize+1:6)]

%transform the data
[sessionData,params]=data_transform(synthAnimalStruct,synthesizingModelNum,hypersPrior);

clearvars -except outputFileNamePrefix terminalStateDirectory runnum hypersPrior sessionData params nParamsToOptimize synthAnimalStruct;
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
%The Max Number of Iterations
maxIter=100;           
% setting max no of EM iterations
params.opts.algorithmic.EMIterations.maxIter     = maxIter;	
% setting max CPU time for EM to inf
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;		

%main call to vEM_Gp that implements the algorithm
tic; [params_vEM, sessionData_vEM, varBound_vEM, EStepTimes, MStepTimes]=vEM_GP(params,sessionData,nParamsToOptimize); toc;

%Re-initialize params to call vEM
params=[];
params_popSpike = PLDSInitialize(sessionData(1).seq,1,'params',params);
params_popSpike.model.C=c_glm;
params_popSpike.model.B=1;
params_popSpike.model.inferenceHandle = @PLDSVariationalInference;
params_popSpike.opts.algorithmic.EMIterations.maxIter     = maxIter;	

for i=1:numel(sessionData)
    sessionData(i).C=params_popSpike.model.C;
    sessionData(i).d=params_popSpike.model.d;
end

%Main call to popSpike vEM implementation
tic; [params_popSpike, sessionData_PopSpikeEM,varBound_PopSpikeEM,EStepTimes,MStepTimes]=PopSpikeEM_SessionsMod(params_popSpike,sessionData); toc;

%Re-Initialize params for the lEM
params=[];
params_popSpikeLaPlace = PLDSInitialize(sessionData(1).seq,1,'params',params);
params_popSpikeLaPlace.model.C=c_glm;
params_popSpikeLaPlace.model.B=1;
params_popSpikeLaPlace.opts.algorithmic.EMIterations.maxIter     = maxIter;	
params_popSpikeLaPlace.model.inferenceHandle = @PLDSLaplaceInference;

for i=1:numel(sessionData)
    sessionData(i).C=params_popSpikeLaPlace.model.C;
    sessionData(i).d=params_popSpikeLaPlace.model.d;
end

%Main call to popSpike lEM implementation
tic; [params_popSpikeLaPlace, sessionData_PopSpikeEM_LaPlace,varBound_PopSpikeEM_LaPlace,EStepTimes,MStepTimes]=PopSpikeEM_SessionsMod(params_popSpikeLaPlace,sessionData); toc;

%% Save the model state
modelSaveFileName = [outputFileNamePrefix,'_fittingRun',num2str(runnum),'_',datestr(now,'HHMMSSddmmyyyy')];
save(fullfile(terminalStateDirectory,modelSaveFileName))
fprintf(['Saved model state to file ',modelSaveFileName,'.\n\n'])


%% Turn off the diary
diary off
end

