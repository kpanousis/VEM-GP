%VEMGP_SLURM run vEM on slurm
%Konstantinos Panagiotis Panousis
%14 August 2015
%Gatsby Computational Neuroscience Unit
%University College London
%
%Script to create runs for SLURM

%% Clear everything
clear;
close all;
clc;

%% set paths
run('set_path');
pathTo_fit='~/Desktop/Matlab/vEMGP/slurmVEM/';


%% set constants 
modelRootLocation='~/Desktop/Matlab/';

totalRuns=30;

sourceDataName='~/Desktop/Matlab/vEMGP/Data/Sessions1_Trials_2000_10ms_Neurons_3.mat';
runFolder='vEMGP/Runs/Run_s1_2000tr_10ms_Neurons3_cg_S_true_fixed/';

%load some stuff, 
%this loads synthAnimalStruct, synthesizingModel num and trials per session
%vector
% load(sourceDataName);
synthesizingModelNum=1;
%% Are we making live SLURM runs, or testing?
liveSlurmRunsLogical = true;

timeEstimateString = '0-23:59:59';  %Format is days-hours:minutes:seconds
qosString = 'normal';

matlabStartFlagsString = '-nodisplay -nosplash -singleCompThread';  %This
matlabLocationString = '/opt/matlab-R2014b/bin/matlab';  %Setting to non-default value

%% Start the diary:
diaryName = ['masterDiary_vEMGP_',datestr(now,'HHMMddmmyyyy'),'.log'];
diaryDir=fullfile(modelRootLocation,'/diary');
if (~exist(diaryDir,'dir'))
    mkdir(diaryDir);
end
diary(fullfile(diaryDir,diaryName))

%% Create storage cell arrays
startTimeString         = datestr(now,'ddmmyyyy_HHMMSS');


lowLevelScriptNames     = cell(1,totalRuns);
lowLevelScriptPaths     = cell(1,totalRuns);
midLevelScriptNames     = cell(1,totalRuns);
midLevelScriptPaths     = cell(1,totalRuns);
midLevelFullPathsCell   = cell(1,totalRuns);
evalStrings              = cell(1,totalRuns);
endMatName='endModelRun';

sessionDir=fullfile(modelRootLocation,runFolder);
if (~exist(sessionDir,'dir'))
    mkdir(sessionDir);
end
for i=1:totalRuns
    
    gpHypersString=['gpHypers=[',num2str(randn(1,1)),';',num2str(randn(1,1)),';',num2str(randn(1,1)),';',num2str(randn(1,1)),'];'];
    mainDirectoryName=['vEMGP_Run_',num2str(i)];
    dir=fullfile(sessionDir,mainDirectoryName);
    if (~exist(dir,'dir'))
        mkdir(dir);
    end
    mkdir(fullfile(dir,'intermediateStates'));
    
    thisRunDirectory=dir;
    midLevelScriptNames{i} = ['vEMGP_', startTimeString, '_run_',num2str(i), '_mid'];
    midLevelScriptPaths{i} = thisRunDirectory;
    lowLevelScriptNames{i} = ['vEMGP_', startTimeString, '_run_',num2str(i), '_low'];
    lowLevelScriptPaths{i} = thisRunDirectory;
        
    midLevelFullPathsCell{i}  = fullfile(midLevelScriptPaths{i},midLevelScriptNames{i});
    
    outputFileNamePrefixString = 'vFitState';
    
    createFLSSString = ['fLSS.modelRootLocation = ''',modelRootLocation,...
        '''; fLSS.initialFileName = ''',endMatName,''';fLSS.sourceName = ''',sourceDataName,'''; fLSS.nameOfInitialClassObj = ''synthAnimalStruct''; fLSS.terminalStateDirectory = ''',fullfile(midLevelScriptPaths{i}),...
        '''; fLSS.intermediateStateDirectory = ''',fullfile(midLevelScriptPaths{i},'intermediateStates'),'''; fLSS.outputFileNamePrefix = ''',outputFileNamePrefixString,'''; fLSS.saveWithinCycleIntermediates = false; '];

    addpathString = ['addpath(''',pathTo_fit,'''); '];
    
    %This line will be the one actually used to run the fitting, using the
    %structure we're creating above:    
    runFittingString = ['fitvEM(',num2str(synthesizingModelNum),', fLSS, ',num2str(i),',gpHypers',');'];
    
    %Concatenate the strings to be executed by MATLAB: 
    executionString = horzcat(gpHypersString,createFLSSString,addpathString,runFittingString);
    
    startupString = '';
    shutdownString = 'exit;';
    evalStrings{i}             = [startupString,' ',executionString,' ',shutdownString];  %This is the complete command MATLAB will run for a given run.
end

%% Create the low- and mid-level scripts for each seed
for i = 1:totalRuns
    makeMidLevelSlurmScript(midLevelScriptNames{i}, midLevelScriptPaths{i}, fullfile(lowLevelScriptPaths{i}, lowLevelScriptNames{i}), timeEstimateString, qosString)
    %makeLowLevelSlurmScript(evalString{seedIdx}, lowLevelScriptPaths{seedIdx}, lowLevelScriptNames{seedIdx})
    makeLowLevelSlurmScript(evalStrings{i}, lowLevelScriptPaths{i}, lowLevelScriptNames{i}, matlabStartFlagsString, matlabLocationString)
end

%% Create the batch script
batchScriptPath = modelRootLocation;
batchScriptName = ['fullvEM_MASTER_',startTimeString];
makeTopLevelSlurmScript(batchScriptName, batchScriptPath, midLevelFullPathsCell)


%% For testing purposes, kick to user
%keyboard

%% Run the batch script, starting the jobs with SLURM
system(['sh ',fullfile(batchScriptPath,batchScriptName)]);


%% Save this workspace for future use in reassembling the runs
dir=(fullfile(modelRootLocation,'workspaces_fitting'))'
if (~exist(dir,'dir'))
    mkdir(dir);
end
save(fullfile(dir,['fittingvEMmultipleRuns_SourceWorkspace_',startTimeString,'_Source_',sourceDataName]))

pause(10)

system('ls > /dev/null')  %Trying to make the slurm job submissions print to MATLAB; this may or may not succeed.

%% Print the final time for the diary:
fprintf(['Finished creating SLURM jobs at time ',datestr(now,'HH:MM:SS, dd/mm/yyyy'),'.\n']);
%% Turn off the diary
diary off

