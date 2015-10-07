%% Setting Path for vEMGP
fprintf(['\n-----------------------------------------------\n'])
fprintf(['Setting up paths for code-package vEMGP!\n'])
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Data_Transform/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Variational/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Initialization/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Calculations/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Calculations/Derivatives/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Calculations/Derivatives/Conditional_Variance/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Calculations/Derivatives/F_GP/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Calculations/Testers/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Calculations/Calcs/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Data/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/slurmVEM/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/vEMGP/Runs/');



%% Setting Path for brodyProject
fprintf(['Setting up paths for code-package brodyProject!\n'])
addpath '~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk'
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/glmCode/')  %Where this file should be
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/utils/')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/kernel/')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/fittingTheta/')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/dataLoading/')  %Need this because it contains Misc/brody_PPCdataDirectory.mat
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/testingTools/')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/testingTools/minFunc/')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/testingTools/minFunc/minFunc/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/testingTools/minFunc/logisticExample/')
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/brodyProject/trunk/testingTools/minFunc/autoDif/');


%% Setting path for Pop Spike Dyn
fprintf(['Setting up paths for code-package PopSpikeDyn!'])
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/');
addpath('~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/');
addpath ~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/core/VariationalInferenceDualLDS
addpath ~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/core/PLDS
addpath ~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/core/LDS
addpath ~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/core
addpath ~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/utils
addpath ~/Desktop/Matlab_Gatsby/Matlab/Misc/MBS_code/pop_spike_dyn/examples

%% Path Complete
fprintf(['\nCompleted!'])
fprintf(['\n-----------------------------------------------\n'])

%% MinFunc Use
%use_our_minfunc=input(['\nIf you do not have a working version of minFunc by Mark Schmidt installed, \n' ...
   % 'and want to use the one packaged with this code (version 2012), please type "Y"...\n'],'s');
   use_our_minfunc='y';
use_our_minfunc=strcmpi(use_our_minfunc,'y');

