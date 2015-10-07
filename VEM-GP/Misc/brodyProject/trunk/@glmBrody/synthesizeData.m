function [varargout] = synthesizeData(obj,synthesizingModelNum,trialsPerSessionVec,neuronsPerSessionVec,synthesisStruct)
%SYNTHESIZEDATA Synthesize data for a glmBrody model class object.
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%31 July, 2014
%
%Method which can be used to synthesize data for a glmBrody class object.
%Creates fields trialData, sessionIDCell, cellsIDCell,
%modelSessionParameters and modelTrialParameters for the synthetic
%sessions.
%Assumes that the model and basis have been created.  Most key
%constants are passed in within synthesisStruct.  trialsPerSessionVec and
%neuronsPerSessionVec are column vectors with the numbers of trials and
%neurons per session respectively.

%% Assertions about class object and inputs
%First, assert that there are no non-synthetic sessions and trials.
assert(numel(obj.trialData) == 0,'Trial data is non-empty!')

%Assert that the class of the model from which we are synthesizing the data
%is a glmBrodyLaplaceSelfOnly model
assert(strcmp(obj.modelTypes{synthesizingModelNum},'glmBrodyLaplaceSelfOnly'),'Synthesizing model must currently be a glmBrodyLaplaceSelfOnly.')

%Assert that the number of elements of trialsPerSessionVec and
%neuronsPerSessionVec are equal and every one is >= 0
assert(numel(trialsPerSessionVec) == numel(neuronsPerSessionVec) && all(trialsPerSessionVec >= 1)...
    && all(neuronsPerSessionVec >= 1),...
    'The number of trials per session and the number of neurons per session must be >= 1.\n Also, these vectors must have the same number of elements.')

%If there is no lapse rate in the generating process, indicate this to the
%user:
if ~isfield(synthesisStruct,'lapseRate')
    fprintf('No lapse rate (field lapseRate) detected in synthesisStruct!\n')
else
    assert(synthesisStruct.lapseRate <= 1 && synthesisStruct.lapseRate >= 0, 'Invalid (not within [0,1]) lapse rate!')
    fprintf('Lapse rate used for synthesis is %1.3f.\n',synthesisStruct.lapseRate)
end

%% Pull down the bleed times and sampling rate
bleedTimes = [obj.modelGlobalParameters{synthesizingModelNum}.relevantSelfHistory,obj.modelGlobalParameters{synthesizingModelNum}.maxLagConsidered];
samplingRate = obj.modelGlobalParameters{1,synthesizingModelNum}.samplingRate;
%% Create the data
obj.trialData = cell(numel(trialsPerSessionVec),1);
latentVarStateAllSessions = cell(numel(trialsPerSessionVec),1);
latentTimeAllSessions = cell(numel(trialsPerSessionVec),1);
for sessionIdx = 1:numel(trialsPerSessionVec)
    obj.sessionIDCell{sessionIdx,1} = sprintf('synthSess%d',sessionIdx);
    obj.cellsIDCell{sessionIdx,1} = cell(1,neuronsPerSessionVec(sessionIdx));
    for cellIdx = 1:neuronsPerSessionVec(sessionIdx)
        obj.cellsIDCell{sessionIdx,1}{1,cellIdx} = sprintf('sS%d_Cell%d',cellIdx);
    end
    
    
    previous_cpoke_out_InGlobalTime = 0;
    for trialIdx = 1:trialsPerSessionVec(sessionIdx)
        %Generate the stimuli
        %Fields state_0_exits (global time), gamma, right_bups, left_bups, bup_diff,
        %stim_start, cpoke_end, cpoke_out.
        %Set rR
        %rR = rand * synthesisStruct.baseClickRate;
        if rand > 0.2
            rR = synthesisStruct.baseClickRate * min(1,max(0,1/2 + sign(rand-1/2) .* (0.1 + 0.1 * randn)));
        else
            rR = synthesisStruct.baseClickRate * rand;
        end
        
        obj.trialData{sessionIdx,1}(1,trialIdx).gamma = log(rR/(synthesisStruct.baseClickRate - rR));  %Where rR is the generating poisson rate on the right side, rL is the left, and rR + rL = baseClickRate;
        
        
        %Determine when the animal pokes its nose into the center port
        obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_start = floor(samplingRate * random('exp',synthesisStruct.assumed_state_0_exit_ToCPokeStart))/samplingRate;
        
        %Draw the click times
        cumulativeTimeSinceCPokeIn = floor(samplingRate * synthesisStruct.startDrawingClicksTimeAfterCPokeStart)/samplingRate;  %Start the stimulus only after a specified time post-cpoke_start
        clickTimes = [];
        while cumulativeTimeSinceCPokeIn < synthesisStruct.totalCPokeStartToCPokeEndTime;
            
            clickTimes = horzcat(clickTimes, cumulativeTimeSinceCPokeIn + floor(samplingRate * random('exp',1/synthesisStruct.baseClickRate))/samplingRate);
            cumulativeTimeSinceCPokeIn = clickTimes(end);
            if ~any(clickTimes < synthesisStruct.totalCPokeStartToCPokeEndTime) %That is, no clicks occurred
                %Mulligan: draw a new click time, such that every trial has
                %at least the first (double-sided) click.
                clickTimes = [];
                cumulativeTimeSinceCPokeIn = synthesisStruct.startDrawingClicksTime;
            end
        end
        
        %Remove any over-time elements of the stimulus.
        clickTimes = clickTimes(clickTimes < synthesisStruct.totalCPokeStartToCPokeEndTime);
        
        %Now, assign the clicks to either right or left
        left_bups_build = clickTimes(1);
        right_bups_build = clickTimes(1);
        for clickIdx = 2:numel(clickTimes)
            if rand < rR / synthesisStruct.baseClickRate
                right_bups_build = [right_bups_build;clickTimes(clickIdx)];
            else
                left_bups_build = [left_bups_build;clickTimes(clickIdx)];
            end
        end
        %right_bups_build = unique(right_bups_build);  %Remove any clicks which occurred at exactly the same time
        %left_bups_build  = unique(left_bups_build);  %Remove any clicks which occurred at exactly the same time
        obj.trialData{sessionIdx,1}(1,trialIdx).left_bups   = unique(floor(samplingRate * (obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_start + left_bups_build))  /samplingRate);
        obj.trialData{sessionIdx,1}(1,trialIdx).right_bups  = unique(floor(samplingRate * (obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_start + right_bups_build)) /samplingRate);
        obj.trialData{sessionIdx,1}(1,trialIdx).stim_start  = floor(samplingRate * (obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_start + clickTimes(1)))    /samplingRate;
        
        clear('clickTimes','left_bups_build','right_bups_build')
        
        %Now, generate the other time stamps
        
        obj.trialData{sessionIdx,1}(1,trialIdx).bup_diff = numel(obj.trialData{sessionIdx,1}(1,trialIdx).right_bups) - numel(obj.trialData{sessionIdx,1}(1,trialIdx).left_bups);
        
        obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_end = floor(samplingRate * (obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_start + synthesisStruct.totalCPokeStartToCPokeEndTime))/samplingRate;
        obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_out = floor(samplingRate * (obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_end + random('exp',synthesisStruct.assumedCPokeEndToCPokeOut)))/samplingRate;
        
        %Set below via value of a[t]
        %obj.trialData{sessionIdx,1}(1,trialIdx).pokedR = ;
        
        obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits = previous_cpoke_out_InGlobalTime + synthesisStruct.assumedCPokeOutToStartNextTrial * (1 + max(-0.4,min( 0.3 * randn,0.4)));  %Selects an additional interval from a random distribution centered on synthesisStruct.assumedCPokeOutToStartNextTrial (8 seconds), and ranging no more than 0.6 to 1.4 of this value.
        obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits = floor(samplingRate * obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits)/samplingRate;
        previous_cpoke_out_InGlobalTime = obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits  + obj.trialData{sessionIdx,1}(1,trialIdx).cpoke_out;
        
        obj.trialData{sessionIdx,1}(1,trialIdx).spikes = cell(1,neuronsPerSessionVec(sessionIdx));
        for neuronIdx = 1:neuronsPerSessionVec(sessionIdx)
            obj.trialData{sessionIdx,1}(1,trialIdx).spikes{1,neuronIdx} = [];
        end
    end
    
    
    
    
    
    %Draw the corresponding values for the latent variable state
    latentVarState = cell(1,trialsPerSessionVec(sessionIdx));
    for trialIdx = 1:trialsPerSessionVec(sessionIdx)
        %Draw a[t] using the stimulus
        %clickTimes = ;
        %clickSigns = ;
        %In the fitting methods:
        [firing,startFrameNumber,clickTimes,clickSigns] = obj.trialParser(sessionIdx,trialIdx,bleedTimes,1/obj.modelGlobalParameters{synthesizingModelNum}.samplingRate);
        
        
        %times = ;
        %in the fitting methods:
        times = (1/obj.modelGlobalParameters{synthesizingModelNum}.samplingRate * (0:1:(size(firing,1) - startFrameNumber)))';
        clear firing
        clear startFrameNumber
        
        priorMeanTrajectory = meanAndDerivs(times,clickTimes,clickSigns,obj.modelGlobalParameters{synthesizingModelNum}.gpHypers);  %Column vector
        kernelMatrix = kernelAndDerivs(times,clickTimes,clickSigns,obj.modelGlobalParameters{synthesizingModelNum}.gpHypers);
        
        latentVarState{1,trialIdx} = mvnrnd(priorMeanTrajectory',kernelMatrix)';
        latentTimes{1,trialIdx} = times;
        
        %Draw the animal's action based on the terminal value of
        %latentVarState
        if isfield(synthesisStruct,'lapseRate')
            coinFlipToUseAccumulator = rand;
            if coinFlipToUseAccumulator >= synthesisStruct.lapseRate
                obj.trialData{sessionIdx,1}(1,trialIdx).pokedR = latentVarState{1,trialIdx}(end) > 0;
            else %Flip a coin
                obj.trialData{sessionIdx,1}(1,trialIdx).pokedR = rand > 0.5;
            end
        else
            obj.trialData{sessionIdx,1}(1,trialIdx).pokedR = latentVarState{1,trialIdx}(end) > 0;
        end
        
    end
    
    
    %     %Stitch the latentVarStates in the cell array into a continuous
    %     %transcript, where the values in between the sessions first
    %     %exponentially decay to zero and then linearly ramp up to the starting
    %     %value of the next session (during at most the last second of the
    %     %change, and no more than half).
    
    globalFrameNumber1IsAtTime = 0;
    assert(floor(globalFrameNumber1IsAtTime * samplingRate) / samplingRate == globalFrameNumber1IsAtTime,'The globalFrameNumber1IsAtTime value must be a integer number of frames from 0 (typically zero).')
    
    %Preallocate
    startFrameNumbers           = zeros(1,trialsPerSessionVec(sessionIdx));
    endFrameNumbers             = zeros(1,trialsPerSessionVec(sessionIdx));
    state_0_exitsFrameNumbers   = zeros(1,trialsPerSessionVec(sessionIdx));
    %Find values
    for trialIdx = 1:trialsPerSessionVec(sessionIdx)
        state_0_exitsFrameNumbers(trialIdx) = round(samplingRate * obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits);
        assert(obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits == state_0_exitsFrameNumbers(trialIdx) / samplingRate,'A value of state_0_exits is not precisely on a global frame!')
        assert(obj.trialData{sessionIdx,1}(1,trialIdx).stim_start == round(samplingRate * obj.trialData{sessionIdx,1}(1,trialIdx).stim_start) / samplingRate,'A value of stim_start is not precisely on a global frame!')
        
        startFrameNumbers(trialIdx) = round(samplingRate * (obj.trialData{sessionIdx,1}(1,trialIdx).state_0_exits - globalFrameNumber1IsAtTime + obj.trialData{sessionIdx,1}(1,trialIdx).stim_start));  %This can have numerical issues due to scale.  Not happy about the rounding, but it should work.
        endFrameNumbers(trialIdx)   = startFrameNumbers(trialIdx) + numel(latentVarState{1,trialIdx}) - 1;
    end
    
    %Above looks complicated.  Instead, we'll just treat a[t] as sparse.
    aOfTInGlobalFrameNumbers = sparse(startFrameNumbers(1):endFrameNumbers(1),1,latentVarState{1,1},endFrameNumbers(end),1);
    %Note: Above may need to be extended a bit after this in time, such
    %that we have some bleed after the end of the last stim.
    
    for trialIdx = 2:trialsPerSessionVec(sessionIdx)
        aOfTInGlobalFrameNumbers(startFrameNumbers(trialIdx):endFrameNumbers(trialIdx),1) = latentVarState{1,trialIdx};
    end
    
    %Using this discrete-time transcript in a[t], generate spike times
    
    
    % Create the generating weights
    %Note: we already have
    %obj.modelGlobalParameters{1,modelIdx}.pseudoBasisSelf, and we should
    %use this basis to approximate a specified shape, e.g., via least
    %squares.
    obj.createGeneratingWeights(synthesizingModelNum,sessionIdx,synthesisStruct)
    
    spikeFramesCell = obj.drawSyntheticSpikeFrames(synthesizingModelNum,sessionIdx,aOfTInGlobalFrameNumbers);
    
    %Assertion that the generated spikes do not lie in the same cell:
    %that would translate to non-binary spiking.
    
    %for neuronIdx = 1:numel(spikeFramesCell)
    %    assert(numel(spikeFramesCell{neuronIdx}) == numel(unique(spikeFramesCell{neuronIdx})), 'More than one spike occurred during a single frame!')
    %end
            
    %Break this back down into trials, finding the spike times and saving
    %the individual trials' data into the
    %trialData{sessionIdx,1}(1,trialIdx).spikes (row) cell arrays
    for trialIdx = 1:trialsPerSessionVec(sessionIdx)
        %Determine start and end frames
        startGlobalFrame = state_0_exitsFrameNumbers(trialIdx);
        if trialIdx < trialsPerSessionVec(sessionIdx)
            endGlobalFrame = state_0_exitsFrameNumbers(trialIdx+1) - 1;
        else
            endGlobalFrame = inf;
        end
              
        
        for cellIdx = 1:neuronsPerSessionVec(sessionIdx)
            framesWithSpikes = spikeFramesCell{1,cellIdx}(spikeFramesCell{1,cellIdx} >= startGlobalFrame & spikeFramesCell{1,cellIdx} <= endGlobalFrame) - startGlobalFrame + 1;
            
            %Convert these to times at which spikes were fired
            timesAtWhichSpikesWereFired = (framesWithSpikes-1) / obj.modelGlobalParameters{synthesizingModelNum}.samplingRate;  %Note: I think it's probably best to just leave the times identical, such that they always definitely end up in the same bin.
            %store this in
            obj.trialData{sessionIdx,1}(1,trialIdx).spikes{cellIdx} = timesAtWhichSpikesWereFired;
        end
    end
    
    latentVarStateAllSessions{sessionIdx,1} = latentVarState;
    latentTimeAllSessions{sessionIdx,1} = latentTimes;
    
end

%% Output handling
if nargout > 0
    if nargout > 2 || nargout == 1
        error('Unsupported number of outputs!')
    else
        varargout{1} = latentTimeAllSessions;
        varargout{2} = latentVarStateAllSessions;
    end
end


%% Make a log entry
obj.logEntry(sprintf('Using model number %d, made synthetic data for %d sessions and corresponding cells (trial numbers and cell numbers included).',synthesizingModelNum,numel(trialsPerSessionVec)),{trialsPerSessionVec, neuronsPerSessionVec})

