%set_path.m
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%8 June 2015
%Script which sets the path for executing the variational EM procedure
%using PopSpikeDyn,brodyProject and vEMGP packages

%% Setting Path for vEMGP
fprintf(['\n-----------------------------------------------\n'])
fprintf(['Setting up paths for code-package vEMGP!\n'])
addpath('vEMGP/');
addpath('vEMGP/Data_Transform/');
addpath('vEMGP/Variational/');
addpath('vEMGP/Calculations/');
addpath('vEMGP/Calculations/Derivatives/');
addpath('vEMGP/Calculations/Derivatives/Conditional_Variance/');
addpath('vEMGP/Calculations/Derivatives/F_GP/');
addpath('vEMGP/Calculations/Testers/');
addpath('vEMGP/Calculations/Calcs/');
addpath('vEMGP/Data/');
addpath('vEMGP/slurmVEM/');



%% Setting Path for brodyProject
fprintf(['Setting up paths for code-package brodyProject!\n'])
addpath 'Misc/brodyProject/trunk'
addpath('Misc/brodyProject')
addpath('Misc/brodyProject/trunk/glmCode/')  %Where this file should be
addpath('Misc/brodyProject/trunk/utils/')
addpath('Misc/brodyProject/trunk/kernel/')
addpath('Misc/brodyProject/trunk/fittingTheta/')
addpath('Misc/brodyProject/trunk/dataLoading/')  %Need this because it contains Misc/brody_PPCdataDirectory.mat
addpath('Misc/brodyProject/trunk/testingTools/')
addpath('Misc/brodyProject/trunk/testingTools/minFunc/')
addpath('Misc/brodyProject/trunk/testingTools/minFunc/minFunc/');
addpath('Misc/brodyProject/trunk/testingTools/minFunc/logisticExample/')
addpath('Misc/brodyProject/trunk/testingTools/minFunc/autoDif/');


%% Setting path for Pop Spike Dyn
fprintf(['Setting up paths for code-package PopSpikeDyn!'])
addpath('Misc/MBS_code/');
addpath('Misc/MBS_code/pop_spike_dyn/');
addpath Misc/MBS_code/pop_spike_dyn/core/VariationalInferenceDualLDS
addpath Misc/MBS_code/pop_spike_dyn/core/PLDS
addpath Misc/MBS_code/pop_spike_dyn/core/LDS
addpath Misc/MBS_code/pop_spike_dyn/core
addpath Misc/MBS_code/pop_spike_dyn/utils
addpath Misc/MBS_code/pop_spike_dyn/examples

%% Path Complete
fprintf(['\nCompleted!'])
fprintf(['\n-----------------------------------------------\n'])

%% MinFunc Use
%use_our_minfunc=input(['\nIf you do not have a working version of minFunc by Mark Schmidt installed, \n' ...
   % 'and want to use the one packaged with this code (version 2012), please type "Y"...\n'],'s');
   use_our_minfunc='y';
use_our_minfunc=strcmpi(use_our_minfunc,'y');

if use_our_minfunc
    addpath Misc/MBS_code/pop_spike_dyn/utils/minFunc/minFunc
    addpath Misc/MBS_code/pop_spike_dyn/utils/minFunc/minFunc/compiled
    try
        mcholC(eye(3));
        disp('Execution of mex-file mcholC worked, so it seems as if functions are compiled for your system.')
    catch
        disp('Execution of mex-file mcholC failed, seems as if the functions are not compiled for your system.')
        disp('Go to ./utils/minFunc/mex and compile the mex-files therein');
    end
end





