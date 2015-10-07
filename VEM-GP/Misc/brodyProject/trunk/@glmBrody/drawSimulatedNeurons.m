function simNeuronFiring = drawSimulatedNeurons(obj,modelNumber,sessionIdx,varargin)
%Draw the activity of a simulated neuron, using the existing weights for a
%particular model.
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%7 February, 2014
%
%Note: since we only run the GLM model inside certain windows, our
%artificial neurons will be completely silent outside of the relevantFrames
%windows; SO LONG AS THESE ARE CONSISTENTLY DEFINED HERE AND IN FITTING
%ROUTINES this will not be a problem for using the simulated cells for
%testing.

%% Set key constant:
doubleStart = true;

%% Input handling
nExtraArgsIn = max(nargin-3,0);
%Determine which neurons we would like to have for our simulated
%recording:  Note that we're going to have to draw all of them anyway
%(since we're assuming cross-linking weights) and so this really boils down
%to which ones we're going to keep at the end.
if nExtraArgsIn == 1
    fittingParameters = varargin{1};
    sourceNeurons = 1:numel(obj.cellsIDCell{sessionIdx});  %Draw every neuron
elseif nExtraArgsIn == 2
    fittingParameters = varargin{1};  %Chiefly interested in the .maxTrialNumber field here, which we're going to pass down
    sourceNeurons = varargin{2};  %Draw for the chosen neurons (designated numerically)
    assert(all(sourceNeurons >0 & sourceNeurons <= numel(obj.cellsIDCell{sessionIdx})),'Source neurons out of range for this session!')
elseif nExtraArgsIn == 0
    fittingParameters.maxTrialNumber = numel(obj.trialData{sessionIdx});
    sourceNeurons = 1:numel(obj.cellsIDCell{sessionIdx});  %Draw every neuron
else
    error('Neurons improperly specified!')
end

if isfield(fittingParameters,'maxTrialNumber')
    fprintf(['Simulating firing for at most the first ',num2str(fittingParameters.maxTrialNumber),' trials of session ',num2str(sessionIdx), '.\n\n']);
end

%% Set some constants
%Find the existing number of neurons
nNeurons = numel(obj.cellsIDCell{sessionIdx});
poissonValsGreaterThanOne = 0;

stimulusWidth = 2;  %Left and right bups.
%Shouldn't be necessary, but just for sanity:
assert(stimulusWidth == 2,'Stimulus width is not 2!')

samplingRate = obj.modelGlobalParameters{modelNumber}.samplingRate;
relevantStimHistory = obj.modelGlobalParameters{modelNumber}.relevantStimHistory;
relevantSelfHistory = obj.modelGlobalParameters{modelNumber}.relevantSelfHistory;
relevantStimSamples = relevantStimHistory * samplingRate;
relevantSelfSamples = relevantSelfHistory * samplingRate;
bupFrames = ceil(samplingRate * obj.modelGlobalParameters{modelNumber}.bupLength);



%% Begin key computations
%Determine the relevant windows for the stimulus.

[relevantFrames, realSpikesInFrames, leftBupsInFrames, rightBupsInFrames, refTime] = obj.getSpikesAndBupsInFrames(modelNumber, sessionIdx, fittingParameters, obj.modelGlobalParameters{modelNumber}.bleedTime);

%For each window, we will look at (the previous spiking history and) the stimulus
%history, and use this to draw spikes for the population; these spikes will
%be organized by frames.

generatedSpikesInFrames = cell(1,nNeurons);
for windowIdx = 1:size(relevantFrames,1)
    %Create an initial history state for the neurons, using the
    %activity of the real neurons in the preceeding frames as
    %initialization.
    
    %For each time step in the window (Guts taken from the GLM routines)
    %Question: do we need to go back 2x the relevantSelfSamples, such that
    %we generate freshly up until relevantSelfSamples before the start of
    %the window for future fitting?
    
    if doubleStart
        startValue = -relevantSelfSamples + 1;
    else
        startValue = 1;
    end
    
    
    for relativeFrameIdx = startValue:(relevantFrames(windowIdx,2) - relevantFrames(windowIdx,1) + 1)  %index within the window, as described by frame number
        globalFrameNumber = relevantFrames(windowIdx,1) + (relativeFrameIdx - 1);  %Match this to the global frame number
        globalSelfFrameNumbers = globalFrameNumber - relevantSelfSamples : globalFrameNumber - 1;
        globalStimFrameNumbers = globalFrameNumber - relevantStimSamples + 1 : globalFrameNumber;
        
        %Some of the below could perhaps use previous timestep's computations
        %Compute the recent stimulus and spiking history as vectors
        if relativeFrameIdx == startValue
            error('Have not yet implemented non-binary spiking in this routine!')
            spikeBinHistory = false(relevantSelfSamples,nNeurons);
            stimBinHistory = false(relevantStimSamples,stimulusWidth);
            
            %Initialize the spiking history with the real neurons' spiking
            %history in the previous period
            for localSelfFrameNumber = 1:relevantSelfSamples
                globalSelfFrame = globalSelfFrameNumbers(localSelfFrameNumber);
                for neuronIdx = 1:nNeurons
                    if any(realSpikesInFrames{neuronIdx} == globalSelfFrame)
                        error('This code needs to change for non-binary spiking!')
                        spikeBinHistory(localSelfFrameNumber,neuronIdx) = true;
                        %This puts the oldest at the top, and the newest at the bottom.
                    end
                end
            end
            
            %Create the bup history
            for localStimFrameNumber = 1:relevantStimSamples
                globalStimFrame = globalStimFrameNumbers(localStimFrameNumber);
                %for stimIdx = 1:stimulusWidth
                %Below: very inefficient implementation; not using
                %the fact that the lists are ordered at all.
                if any(leftBupsInFrames + bupFrames-1 >= globalStimFrame & leftBupsInFrames <= globalStimFrame)
                    stimBinHistory(localStimFrameNumber,1) = true;
                end
                if any(rightBupsInFrames + bupFrames-1 >= globalStimFrame & rightBupsInFrames <= globalStimFrame)
                    stimBinHistory(localStimFrameNumber,2) = true;
                end
                %end
            end
            
        else  %We've just computed these for the previous frame: just update with the newest information
            %Use the (now former) spike state to update
            %spikeBinHistory
            spikeBinHistory = [spikeBinHistory(2:end,:);spikeStateNow];
            
            %And similarly, compute only the latest bups to update the
            %bup history.
            globalStimFrame = globalStimFrameNumbers(end);
            lastBupLine = false(1,stimulusWidth);
            if any(leftBupsInFrames + bupFrames-1 >= globalStimFrame & leftBupsInFrames <= globalStimFrame)
                lastBupLine(1,1)  = true;
            end
            
            if any(rightBupsInFrames + bupFrames-1 >= globalStimFrame & rightBupsInFrames <= globalStimFrame)
                lastBupLine(1,2)  = true;
            end
            stimBinHistory = [stimBinHistory(2:end,:); lastBupLine];
            
        end
        
        %Now, we're going to draw this, not read it out of the data
        % %Create the self-history
        % spikeStateNow = false(1,nNeurons);
        % for neuronIdx = 1:nNeurons
        %     if any(spikesInFrames{neuronIdx} == globalFrameNumber)
        %         spikeStateNow(neuronIdx) = true;
        %     end
        % end
        %
        % %Create the stimulus history: this is a binary matrix with a height equal to the
        % %number of stimulus frames to be examined, and a width equal to the number of stimuli.
        
        
        
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
        transformedSpikingHistory = obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf' * spikeBinHistory;  %now nBasisCompsSelf x nNeurons
        
        %Same procedure with the stimulus history.
        transformedStimHistory = obj.modelGlobalParameters{modelNumber}.pseudoBasisStim' * stimBinHistory;  %now nBasisCompsStim x nStimuli
    
        %Weight matrices are C (which is (nNeurons x (nBasisCompsStimuli x
        %nStimuli)), d (nNeurons x 1), and W (nNeurons x (nBasisCompsSelf x
        %nNeurons)).  Reshape transformedSpikingHistory,
        %transformedStimHistory appropriately to get the right outputs to
        %logLambda
        
        %Apply the weights to compute the log rates for this time step
        
        logLambda = obj.modelSessionParameters{sessionIdx, modelNumber}.d + obj.modelSessionParameters{sessionIdx, modelNumber}.C *reshape(transformedStimHistory,[],1)   + obj.modelSessionParameters{sessionIdx, modelNumber}.W * reshape(transformedSpikingHistory,[],1);
        
        %Draw from those rates to obtain the spiking for this bin
        initialDraw = poissrnd(exp(logLambda'));  %since we need a row vector.
        
        %Add to poissonValsGreaterThanOne number of values > 1;
        poissonValsGreaterThanOne = poissonValsGreaterThanOne + sum(initialDraw > 1);
        error('This code (above and below) needs to change for non-binary spiking!')
        spikeStateNow = logical(initialDraw);
        
        
        %Record this spiking activity.
        for neuronIdx = 1:nNeurons
            if spikeStateNow(neuronIdx) > 0
                generatedSpikesInFrames{neuronIdx} = vertcat(generatedSpikesInFrames{neuronIdx},globalFrameNumber);
                error('This code needs to change for non-binary spiking!')
            end
        end
    end
end

%Use refTime from above to convert the frames back into global times

%Using obj.trialData{sessionIdx}(trialIdx).state_0_exits, we may obtain the
%time at which any trial started: thus, we can segment the spike times we
%generate above back into the appropriate trials and put them in aan
%appropriate output cell array.

keepingSpikesInFrames = generatedSpikesInFrames(sourceNeurons);


keepingSpikesInTimes = cell(size(keepingSpikesInFrames));
for neuronIdx = 1:numel(sourceNeurons)
    keepingSpikesInTimes{neuronIdx} = (keepingSpikesInFrames{neuronIdx} - 1) / samplingRate + refTime;
end

%The above are now back in global times at the start of the respective
%frames.  %NOTE (7 FEB 2014): NOT SURE IF FRAME NUMBERS ARE STABLE UNDER REPEATED
%APPLICATIONS OF THIS PROCESS. TRY TO MAKE SURE THIS IS REVERSIBLE.



%Create simNeuronFiring, a cell array, where each entry in the cell array
%is the firing time for a cell in the same format as in the original data,
%i.e., just as if this were an entry in obj.trialData{sessionIdx}
simNeuronFiring = cell(1,numel(obj.trialData{sessionIdx}));
for trialIdx = 1:size(simNeuronFiring,2)
    simNeuronFiring{1,trialIdx} = cell(1,numel(sourceNeurons));
    %Check if any spikes were fired by each virtual cell during this trial
    
    
    %Look at global times marking the start and end of the trial
    tStart = obj.trialData{sessionIdx}(trialIdx).state_0_exits;
    if trialIdx ~= numel(obj.trialData{sessionIdx})
        tEnd = obj.trialData{sessionIdx}(trialIdx + 1).state_0_exits;
    else
        tEnd = inf;
    end
    
    %For each neuron, store any spikes in this window
    for neuronIdx = 1:numel(sourceNeurons)
        simNeuronFiring{1,trialIdx}{neuronIdx} = keepingSpikesInTimes{neuronIdx}(keepingSpikesInTimes{neuronIdx} >= tStart & keepingSpikesInTimes{neuronIdx} < tEnd) - obj.trialData{sessionIdx}(trialIdx).state_0_exits;
    end
end    



if poissonValsGreaterThanOne > 0
    fprintf(['Set ',num2str(poissonValsGreaterThanOne),' instances of > 1 spike fired in a bin to 1.\n\n'])
end