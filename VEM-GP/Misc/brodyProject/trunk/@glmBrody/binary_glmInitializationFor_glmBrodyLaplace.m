function binary_glmInitializationFor_glmBrodyLaplace(obj,modelNumber, sessionIdx, fittingParameters)
%GLMINITIALIZATIONFOR_GLMBRODYLAPLACE Method which initializes W, d for a
%glmBrodyLaplace model.
%glmInitializationFor_glmBrodyLaplace.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%2 May, 2014
%
%SETS THE FIELDS obj.modelSessionParameters{sessionIdx, modelNumber}.d
%AND obj.modelSessionParameters{sessionIdx, modelNumber}.W
%via calculations done here.

fprintf('Running binary r180 version of glmInitializationFor_glmBrodyLaplace.\n')
compressedWSelfWeightRepresentation = true;


assert(any(strcmp(obj.modelTypes{modelNumber},{'glmBrodyLaplace','glmBrodyLaplaceSelfOnly'})),...
    'Model is not of type glmBrodyLaplace or glmBrodyLaplaceSelfOnly!')

if strcmp(obj.modelTypes{modelNumber}, 'glmBrodyLaplaceSelfOnly')
    selfOnlyFlag = true;
    fprintf('Fitting a glmBrodyLaplaceSelfOnly model with only self-weights: no cross!\n\n')
elseif strcmp(obj.modelTypes{modelNumber}, 'glmBrodyLaplace')
    selfOnlyFlag = false;
    fprintf('Fitting a glmBrodyLaplace model with self- and cross-weights.\n\n')
else
    error('Unrecognized model type!')  %We shouldn't get here, since this should fail the assertion immediately above.
end

%% Set some constants we'll need
bleedTime = obj.modelGlobalParameters{modelNumber}.bleedTime;
nCompsPseudoBasis = size(obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf,2);
samplingRate = obj.modelGlobalParameters{modelNumber}.samplingRate;
relevantSelfHistory = obj.modelGlobalParameters{modelNumber}.relevantSelfHistory;
relevantSelfSamples = relevantSelfHistory * samplingRate;

nNeurons = numel(obj.cellsIDCell{sessionIdx,1});


%% Begin the fitting procedure
fprintf(['GLM parameters d and W are absent for session number ',num2str(sessionIdx),'; using traditional GLM \n initialization.\n\n']);

[relevantFrames, spikesInFrames] = obj.getSpikesAndBupsInFrames(modelNumber, sessionIdx, fittingParameters, bleedTime,'trialBased');

allRelevantInputs = zeros(sum(relevantFrames(:,2) - relevantFrames(:,1) + 1), nNeurons * nCompsPseudoBasis);
allRelevantSpikeLogicals = false(sum(relevantFrames(:,2) - relevantFrames(:,1) + 1), nNeurons);
previousFrames = 0;
for windowIdx = 1:size(relevantFrames,1)
    
    for relativeFrameIdx = 1:(relevantFrames(windowIdx,2) - relevantFrames(windowIdx,1) + 1)  %index within the window, as described by frame number
        globalFrameNumber = relevantFrames(windowIdx,1) + (relativeFrameIdx - 1);  %Match this to the global frame number
        globalSelfFrameNumbers = globalFrameNumber - relevantSelfSamples : globalFrameNumber - 1;
        
        %Some of the below could perhaps use previous timestep's computations
        %Compute the recent stimulus and spiking history as vectors
        if relativeFrameIdx == 1
            spikeBinHistory = false(relevantSelfSamples,nNeurons);
            
            %Create the spiking history
            for localSelfFrameNumber = 1:relevantSelfSamples
                globalSelfFrame = globalSelfFrameNumbers(localSelfFrameNumber);
                for targetNeuronIdx = 1:nNeurons
                    if any(spikesInFrames{targetNeuronIdx} == globalSelfFrame)
                        spikeBinHistory(localSelfFrameNumber,targetNeuronIdx) = true;
                        %This puts the oldest at the top, and the newest at the bottom.
                    end
                end
            end
            
        else  %We've just computed these for the previous frame: just update with the newest information
            %Use the (now former) spike state to update
            %spikeBinHistory
            spikeBinHistory = [spikeBinHistory(2:end,:);spikeStateNow];
        end
        
        
        %Create the self-history
        spikeStateNow = false(1,nNeurons);
        for targetNeuronIdx = 1:nNeurons
            if any(spikesInFrames{targetNeuronIdx} == globalFrameNumber)
                spikeStateNow(targetNeuronIdx) = true;
            end
        end
        
        %Compute the feature-transformed spiking history.
        %obj.modelGlobalParameters(modelNumber).pseudoBasisSelf has
        %been constructed such that it is an nFeatures x nSamples
        %matrix, with the most recent components in the bottom rows and
        %the oldest at the top. Further, the basis functions are
        %ordered such that the most recent (i.e., that which puts the
        %most weight on the most recent time-steps) is in the left-most
        %column and the oldest (weight farther back in time) is on the
        %right. pseudoBasis' thus indexes from most recent features to
        %oldest features (by row) and by frame (column).
        transformedSpikingHistory = obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf' * spikeBinHistory;  %nCompsPseudoBasis x nNeurons
        
        
        %Write these values out.
        allRelevantInputs(previousFrames + relativeFrameIdx,:) = reshape(transformedSpikingHistory,1,[]);
        %Above reordering reads down the columns of transformedSpikingHistory, and when it hits the end of a column, moves to the next column, continuing to write horizontally.
        %This means that the row vector created is organized in large
        %blocks by source neuron and within the blocks by feature of the
        %source neuron's recent spiking.
        allRelevantSpikeLogicals(previousFrames + relativeFrameIdx,:) = spikeStateNow;
        
        
    end
    previousFrames = previousFrames + relativeFrameIdx;
end

%Print status.
fprintf('Finished assembling local feature representations.\n');

%% Initialize the weight matrices by using lassoglm
%to perform fitting, getting back the matrix of regularlized weights
%and the stats structure, which includes the regularization parameter.

%Preallocate the weight matrices
d = zeros(nNeurons,1);
if ~selfOnlyFlag || ~compressedWSelfWeightRepresentation
    W = zeros(nNeurons,nNeurons * nCompsPseudoBasis);  %Rows: target neuron; Columns: large blocks are source neuron, then at finer scale, features of source neuron.
elseif compressedWSelfWeightRepresentation
    W = zeros(nNeurons, nCompsPseudoBasis);  %Rows: target neuron; Columns: features of same neuron's history.
else
    error('Unrecognized condition!')
end

successfullyFitANeuronFlag = false;

for targetNeuronIdx = 1:nNeurons  
    %Note that targetNeuronIdx indexes the TARGET of the glm, the firing we're trying to explain given a particular collection of features from the history.
    
    %Define source and target columns:
    if selfOnlyFlag && ~compressedWSelfWeightRepresentation
        relevantColumnsInFeatures = (targetNeuronIdx - 1) * nCompsPseudoBasis + 1 : targetNeuronIdx * nCompsPseudoBasis;
        relevantColumnsInW = relevantColumnsInFeatures;
    elseif selfOnlyFlag
        relevantColumnsInW = 1:nCompsPseudoBasis;
        relevantColumnsInFeatures = relevantColumnsInW + (targetNeuronIdx - 1) * nCompsPseudoBasis;
    else
        relevantColumnsInFeatures = 1:size(allRelevantInputs,2);
        relevantColumnsInW = relevantColumnsInFeatures;
    end
    
    
    if ~successfullyFitANeuronFlag
        
        if any(allRelevantSpikeLogicals(:,targetNeuronIdx)) && ~all(allRelevantSpikeLogicals(:,targetNeuronIdx))  %That is, the neuron fired a spike during the period of interest, and did not fire a spike in every window, such that we have something to fit, potentially
            %Create the rows of the corresponding matrices.
            [parameterizedWeights,FitInfo] = lassoglm(allRelevantInputs(:,relevantColumnsInFeatures),   allRelevantSpikeLogicals(:,targetNeuronIdx), 'poisson','Alpha',fittingParameters.alpha,'NumLambda',17);
            
            %For now, choosing a value in the middle.
            colNum = ceil(size(parameterizedWeights,2)/2);
            lambdaValue = FitInfo.Lambda(colNum);
            
            %Pull out the relevant components into C, d, W.
            d(targetNeuronIdx,1) = FitInfo.Intercept(colNum);
            W(targetNeuronIdx,relevantColumnsInW) = parameterizedWeights(:, colNum)';
            successfullyFitANeuronFlag = true;
        elseif ~any(allRelevantSpikeLogicals(:,targetNeuronIdx))
            %The cell didn't fire, so set all of its weights to zero: this
            %makes the most sense.
            d(targetNeuronIdx,1) = 0;
            W(targetNeuronIdx,relevantColumnsInW) = zeros(1,numel(relevantColumnsInFeatures));
        else %The cell fired in every window; very, very unlikely.
            d(targetNeuronIdx,1) = 1;  %Bad choice, but could be OK, since it gives us something to work with, anyway.
            W(targetNeuronIdx,relevantColumnsInW) = zeros(1,numel(relevantColumnsInFeatures));
        end
    else
        %Use the same Lambda value as for the first run
        [parameterizedWeights,FitInfo] = lassoglm(allRelevantInputs(:,relevantColumnsInFeatures),   allRelevantSpikeLogicals(:,targetNeuronIdx), 'poisson','Alpha',fittingParameters.alpha,'Lambda',lambdaValue);
        d(targetNeuronIdx,1) = FitInfo.Intercept(1);
        W(targetNeuronIdx,relevantColumnsInW) = parameterizedWeights(:, 1)';
    end
    %So if this isn't the first neuron fit, we're choosing a consistent value.
    %Choose and set the appropriate weights from the lassoglm fitting
    %Key question: how should we determine which set of weights in the
    %regularlization schedule is the correct one?  Lassoglm implements
    %cross-validation, but I don't know if we should bother.  Also,
    %can't really keep that consistent across multiple independent fits
    %of the individual neurons.
    
    %At the moment, just choosing a moderate lambda value from the
    %first neuron's run, and then using the same lambda value for all
    %of the rest of the neurons.
end

%Choose and set the appropriate weights from the lassoglm fitting
%Key question: how should we determine which set of weights in the
%regularlization schedule is the correct one?

%Print a status message indicating this step is complete.
fprintf(['Completed initializing weights for session number ',num2str(sessionIdx),'.\n\n']);


%THE BELOW IS A PROBLEM FOR FULL PARALLELIZATION:
%One option would be to NOT set the fields in
%obj.modelSessionParameters until the end, instead just using d and
%W here and below.
obj.modelSessionParameters{sessionIdx, modelNumber}.d = d;
obj.modelSessionParameters{sessionIdx, modelNumber}.W = W;
% 
% clear d
% clear W
