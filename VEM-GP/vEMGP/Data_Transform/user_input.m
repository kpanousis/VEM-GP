function [data_create,sessionNum,trialsNum,neuronsNum,samplingRate,nParamsToOptimize]=user_input()
%USER_INPUT function to read some parameters for the user in order to
%create or load data
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%August 2015

%Outputs: data_create: flag to create or load data
%         sessionNum: the number of sessions
%         trialsNum: the number of trials
%         neuronsNum: the number of neurons
%         samplingRate: the samplingRate
%         nParamsToOptimize: the number of params to optimize

%Create or load data?
read_data=input('1)Load Data from file\n2)Create new Data\nInput: ');
data_create=false;

if (read_data==1)
    message='load';
    %If directory is empty request input for data creation
    if (isempty(dir('vEMGP/Data/*.m*')))
        disp('Empty Data Directory.');
        message='create';
        data_create=true;
        %else display file and request input for loading the appropriate data
    else
        dir vEMGP/Data/*.m*
    end
else
    message='create';
    data_create=true;
end

%% User Input
fprintf('-----------------\nSession Num\n-----------------\n');
sessionNum=input(['Please specify the number of sessions to ',message,' data.\nInput: ']);
trialsNum=zeros(1,sessionNum);
neuronsNum=zeros(1,sessionNum);

%For each session ask the user for the specifics
for ses=1:sessionNum
    fprintf('-----------------\nSession %d Inputs\n-----------------\n',ses);
    message=['Please give number of trials for session ',num2str(ses),'.\nInput: '];
    trialsNum(ses)=input(message);
    message=['Please give number of neurons for session ',num2str(ses),'.\nInput: '];
    neuronsNum(ses)=input(message);
end

%sampling Rate input
fprintf('--------------------\nSampling Rate Input\n--------------------\n');
samplingRate=input('Please specify the sampling rate (Hz).\nInput: ');

%User Input of the number of params to optimize
fprintf('---------------------------\nParameter Optimization Num\n---------------------------\n')
nParamsToOptimize=input('Please specify the number of parameters to optimize (1-6)\nInput: ');



end