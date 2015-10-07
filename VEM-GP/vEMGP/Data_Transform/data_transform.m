function [sessionData, params] = data_transform(synthAnimalStruct,synthesizingModelNum,hypers,c_glm,useWeights)
%DATA_TRANSFORM Function to transform the brody project data to a
%compatible form for the pop_Spike_dyn package
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University Collge London
%8 June 2015
%Inputs: synthAnimaStruct object containing all the needed data
%        synthesizingModelNum: the number of the model
%        hypers: the hyperparameters vector
%        c_glm: the initial c for the algorithm (usually zero)
%        useWeights: use the self history flag
%Outputs: SessionData that contain the necessary info for each trial. This
%         includes C and d for the GLM and seq, 
%             seq: Containing fields
%                 seq.y: NxT for a spike raster from N neurons
%                 seq.T: the 2nd dimension of seq.Y
%                 and fields for the GP (clickTimes,Signs etc)
%         params: initialized by the modified pop_spike_dyn code


%% Bleed Times
%Get the bleed times from the data
bleedTimes=[synthAnimalStruct.modelGlobalParameters{synthesizingModelNum}.relevantSelfHistory,...
    synthAnimalStruct.modelGlobalParameters{synthesizingModelNum}.maxLagConsidered];

%% Firing Rate
%Get the firing from the trial Parser function of brody project package
%Checked trial Parser, everything seems fine
samplingRate=synthAnimalStruct.modelGlobalParameters{synthesizingModelNum}.samplingRate;


%% Firing and Times
%get number of sessions
nSessionsIndices=1:size(synthAnimalStruct.trialData,1);

%Initialize the seq struct
sessionData(nSessionsIndices(end)).seq=[];
for i=nSessionsIndices
    sessionData(i).seq(numel(synthAnimalStruct.trialData{i,1})).y=rand(1,1);
end


%For each session and each trial get the needed data
for sessionIdx=nSessionsIndices
    for trialIdx=1:numel(synthAnimalStruct.trialData{sessionIdx,1})
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate the stimulus drive for each trial
        %%% Code from Dr. T. Desautels         %%%%%
        %%% No permission to be made public or %%%%%
        %%% redistribute                       %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [firing, startFrameNumber,clickTimes,clickSigns]=...
            synthAnimalStruct.trialParser(sessionIdx,trialIdx,bleedTimes,...
            1/samplingRate);
        
        if (useWeights)
            [latentTrajectory] = meanAndDerivs((1/samplingRate*(0:1:(size(firing,1)-startFrameNumber)))',clickTimes,clickSigns,...
                synthAnimalStruct.modelGlobalParameters{synthesizingModelNum}.gpHypers);
            nFrames = size(latentTrajectory,1);
            nNeurons = size(firing,2);
            [nHistoryFrames, nCompsPseudoBasis] = size(synthAnimalStruct.modelGlobalParameters{synthesizingModelNum}.pseudoBasisSelf);
            nHistoryFeatures = nNeurons * nCompsPseudoBasis;  %The number of features available at each time point.
            firstEndFrame   = startFrameNumber - 1;
            firstStartFrame = firstEndFrame - nHistoryFrames + 1;
            featureValues = zeros(nHistoryFeatures, nFrames); %TOTALNUMBEROFHISTORYFEATURES x nFramesOfInterest
            for frameIdx = 1:nFrames
                chunkOfFiringWeNeed = firing(firstStartFrame + (frameIdx - 1) : firstEndFrame + (frameIdx - 1),:);  %nHistoryFrames x nNeurons
                %obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf
                transformedIntoNewBasis = chunkOfFiringWeNeed' * synthAnimalStruct.modelGlobalParameters{synthesizingModelNum}.pseudoBasisSelf; %Now nNeurons (source) x nFeatures
                reshapedInNewBasis = reshape(transformedIntoNewBasis',[],1);
                %Above: This row translates to reading across the rows of the above matrix (as in normal English text) and filling in a column vector with the entries in the order encountered.
                %The result is a column vector with nNeurons (number of source neurons) blocks of nFeatures.
                %Store this in the appropriate column, where featureValues will
                %eventually be TOTALNUMBEROFHISTORYFEATURES x nFramesOfInterest
                featureValues(:,frameIdx) = reshapedInNewBasis;
            end
            W=synthAnimalStruct.modelSessionParameters{sessionIdx,1}.W;
            if ~all(size(W) == [nNeurons , nNeurons * nCompsPseudoBasis])
                if all(size(W) == [nNeurons , nCompsPseudoBasis]) && any(strcmp(synthAnimalStruct.modelTypes{synthesizingModelNum},{'glmBrodyLaplaceSelfOnly','glmBrodyFixedSelfOnly'}))  %If we're working with a self-weights only model and it's using the sparse representation:
                    
                    %Reshape W such that the rows of W make the "blocks" of a new,
                    %block-diagonal W.
                    newW = zeros(nNeurons, nNeurons * nCompsPseudoBasis);
                    for targetNeuron = 1:nNeurons
                        newW(targetNeuron, (targetNeuron - 1) *  nCompsPseudoBasis + 1 : targetNeuron *  nCompsPseudoBasis) = W(targetNeuron, :);
                    end
                    
                    W = newW; %reshaped version
                else
                    error('Unrecognized shape of W!')
                end
            end
            sessionData(sessionIdx).seq(trialIdx).s=W*featureValues;
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        %calculate all the rest
        %%%%%%%%%%%%%%%%%%%%%%%
        firing=firing(startFrameNumber:end,:)';
        sessionData(sessionIdx).seq(trialIdx).y=firing;
        sessionData(sessionIdx).seq(trialIdx).T=size(firing,2);
        %Added this on the 28/6 for getting Vsm and VVsm augmented (and mu)
        sessionData(sessionIdx).seq(trialIdx).times=1/samplingRate*(0:1:(size(firing,2)-1));
        sessionData(sessionIdx).seq(trialIdx).clickTimes=clickTimes;
        sessionData(sessionIdx).seq(trialIdx).startFrameNumber=startFrameNumber;
        sessionData(sessionIdx).seq(trialIdx).clickSigns=clickSigns;
        sessionData(sessionIdx).seq(trialIdx).sessionNum=sessionIdx;
        %put input u (click Times) for the E Step
        sessionData(sessionIdx).seq(trialIdx).u=calculate_Input(sessionData(sessionIdx).seq(trialIdx),hypers);
        
    end
end

%% Initialize Params
params=[];
xDim=1;
params.gp.useGP=true;

%% Initial and ground truth hypers
params.model.hypers=hypers;
params.model.initialHypers=hypers;
params.model.groundTruthGP=synthAnimalStruct.modelGlobalParameters{1,1}.gpHypers;

%% Markovian representation
params.model.Q0=exp(params.model.hypers(1));
params.model.x0=0;
params.model.A=exp(params.model.hypers(4)*(1/samplingRate));
params.model.B=1;
params.model.samplingRate=samplingRate;

%% Params for the GLM
%The PLDSInitialize is a function of popSPike_dyn package
%It uses PLDSsetDefaultParameters, again of the same package that was
%modified to initialize the GP params a little different
params=PLDSInitialize(sessionData(1).seq,xDim,'params',params);
if (~useWeights)
    params.model.notes.useS=false;
end
params.model.C=c_glm;
for i=1:numel(sessionData)
    sessionData(i).C=c_glm;
    sessionData(i).d=params.model.d;
    if (useWeights)
        sessionData(i).W=synthAnimalStruct.modelSessionParameters{i,1}.W;
    end
end

end




