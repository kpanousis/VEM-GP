function generatedSpikesInFrames = drawSyntheticSpikeFramesBinary(obj,synthesizingModelNum,sessionIdx,aOfTInGlobalFrameNumbers)
%DRAWSYNTHETICSPIKEFRAMESBINARY For a session, draw the frame numbers of spikes.
%
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%1 August, 2014
%
%Method which uses the GLM model to run forward in time and draws spikes
%for all cells in a particular session.  The resulting frames with spikes
%are given in the 1 x nCells cell array spikeFrames.

%% Input checking
assert(strcmp(obj.modelTypes{synthesizingModelNum},'glmBrodyLaplaceSelfOnly'),'Invalid model for spike generation.')
assert(numel(obj.trialData) >= sessionIdx,'Session number invalid')

%% Set constants
nNeurons = numel(obj.cellsIDCell{sessionIdx,1});
relevantSelfSamples = size(obj.modelGlobalParameters{1,synthesizingModelNum}.pseudoBasisSelf,1);  %The number of frames into the past the GLM will use to calculate spike rates
nCompsPseudoBasis   = size(obj.modelGlobalParameters{1,synthesizingModelNum}.pseudoBasisSelf,2);  %The number of basis components in the self-history representation.
samplingRate = obj.modelGlobalParameters{1,synthesizingModelNum}.samplingRate;

%% Pull down the weights

W = obj.modelSessionParameters{sessionIdx, synthesizingModelNum}.W;

% From calculateFiringRates_BrodyGLM.m: Sparse W mode switching:
if ~all(size(W) == [nNeurons , nNeurons * nCompsPseudoBasis])
    if all(size(W) == [nNeurons , nCompsPseudoBasis]) && strcmp(obj.modelTypes{synthesizingModelNum},'glmBrodyLaplaceSelfOnly')  %If we're working with a self-weights only model and it's using the sparse representation:

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



%% Generate some past history
spikeBinHistory = false(relevantSelfSamples,nNeurons);

%Initialize by independent poisson randoms; this should be fine for now.
for neuronIdx = 1:nNeurons
   %spikeBinHistory(:,neuronIdx) = rand(relevantSelfSamples,1) > poisspdf(0,obj.modelSessionParameters{sessionIdx,synthesizingModelNum}.d(neuronIdx)/samplingRate);
   spikeBinHistory(:,neuronIdx) = rand(relevantSelfSamples,1) > poisspdf(0,exp(obj.modelSessionParameters{sessionIdx,synthesizingModelNum}.d(neuronIdx)));
end

poissonValsGreaterThanOne = 0;

%% Draw the remainder of the history, recording spike times
%Reworked guts from drawSimulatedNeurons.m, from 7 February, 2014.
generatedSpikesInFrames = cell(1, nNeurons);
for globalFrameNum = 1:size(aOfTInGlobalFrameNumbers)
   %Use spikeBinHistory and a[t] to calculate the values of log(lambda_i)
   
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
   transformedSpikingHistory = obj.modelGlobalParameters{synthesizingModelNum}.pseudoBasisSelf' * spikeBinHistory;  %now nBasisCompsSelf x nNeurons
   
   logLambda = obj.modelSessionParameters{sessionIdx, synthesizingModelNum}.d + obj.modelSessionParameters{sessionIdx, synthesizingModelNum}.C * aOfTInGlobalFrameNumbers(globalFrameNum,1) + W * reshape(transformedSpikingHistory,[],1);
   
   
   %draw spikes from a poisson random with lambdas = lambda_i
   initialDraw = poissrnd(exp(logLambda'));  %since we need a row vector.
   
   %Add to poissonValsGreaterThanOne number of values > 1;
   poissonValsGreaterThanOne = poissonValsGreaterThanOne + sum(initialDraw > 1);
   spikeStateNow = logical(initialDraw);
   
   
   for neuronIdx = 1:nNeurons
       if spikeStateNow(neuronIdx)
           generatedSpikesInFrames{neuronIdx} = vertcat(generatedSpikesInFrames{neuronIdx},globalFrameNum);
       end
   end
           
    
    
   spikeBinHistory = [spikeBinHistory(2:end,:);spikeStateNow]; 
   %Note: The above implies that new values of spikeBinHistory are added
   %onto the bottom, and old values get bumped off the top; this must be
   %consistent with the way the basis is chosen.
end



%% Alert the user if our rates appear to be too high.
if poissonValsGreaterThanOne > 0
    fprintf('During drawing process, had %d frames in which a cell fired more than one spike.\n\n',poissonValsGreaterThanOne)
end