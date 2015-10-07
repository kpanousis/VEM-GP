function [varargout] = accumulatorVsRate(obj,modelNumber, sessionIdx, cellIdx, varargin)
%ACCUMULATORVSRATE Plots a series of smoothed tuning curves for a cell
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%24 March, 2014
%
%For a large set of lags, calculates the correlation of the smoothed firing
%rate with the evidence accumulator value, plots the tuning curve, and (if
%requested) returns the lags and the correlation coefficients, such that
%the value can be maximized.

%% Output checking
nExtraArgsOut = max(nargout,0);
assert(nExtraArgsOut == 0 || nExtraArgsOut == 2, 'Improper number of outputs of accumulatorVsRate.')

%% Input checking
nExtraArgsIn = max(nargin-4,0);
assert(nExtraArgsIn <= 1);
if nExtraArgsIn == 1
    priorOrPosterior = varargin{1};
else
    priorOrPosterior = 'priorUnderCurrentGPHypers';
    fprintf('Assuming that you intended to use the current GP prior hyperparameters to compute the accumulator trajectories.\n\n')
end

%% Set key values
lagMin = 0;  %Seconds
%lagMax = 0.400;
lagMax = obj.modelGlobalParameters{modelNumber}.maxLagConsidered;  %In seconds
dLag = 0.010;

nCells = numel(obj.cellsIDCell{sessionIdx,1});
startTime = 0;

assert((lagMax - lagMin)/dLag == round((lagMax - lagMin)/dLag), 'Non-integer number of potential lags selected.')


dt = 1/obj.modelGlobalParameters{modelNumber}.samplingRate;

%filteringScheme = 'stepwiseOld';
%filteringScheme = 'stepwiseNew';
filteringScheme = 'halfGaussian';
filterParams = [0,0.2,3];  %How should the third entry in filterParams be
%set?


%% Create and store the per-trial representations
accumulatorValuesStor = cell(numel(obj.trialData{sessionIdx}),1);
smoothedFiringStor = cell(nCells,numel(obj.trialData{sessionIdx}));

%% At global level, get the times of all spikes in the whole session:
[allSpikes, ~,~,perTrialTimeOffsets] = obj.listAllBandS_Method(sessionIdx,numel(obj.trialData{sessionIdx}),startTime);  %Returns a 1 x nCells cell array, where the ith entry is the set of spike times for all trials.
%Note: perTrialTimeOffsets = obj.trialData{sessionIdx}(trialIdx).state_0_exits - startTime;
%allSpikes{cellIdx} are column vectors.
for cellIdx = 1:nCells
    assert(issorted(allSpikes{cellIdx}),'A cell''s spiking is not sorted!\n')
end

switch priorOrPosterior
    case {'priorUnderCurrentGPHypers','priorUnderHyperpriorMean'}
        maxTrialNum = numel(obj.trialData{sessionIdx});
    case 'posterior'
        maxTrialNum = numel(obj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories);
        if maxTrialNum ~= numel(obj.trialData{sessionIdx})
            fprintf('Using %s, the number of trials which have been fitted, as the max trial number.\n\n',num2str(maxTrialNum))
        end
end

for trialIdx = 1:maxTrialNum
    
    %Parse the trial:
    [firing, startFrameNumber, clickTimes, clickSigns] = obj.trialParser(sessionIdx,trialIdx,[0,lagMax],dt);
    
    switch filteringScheme
        case 'halfGaussian'
            if ~islogical(firing)
                fprintf('Assuming halfGaussian filtering scheme is OK for non-binary firing variable.\n')
            end
        otherwise
            assert(islogical(firing),'Not debugged for non-binary firing (i.e., spike logical) object.  This looks as if it will cause substantial problems below!')
    end
    %clickTimes and clickSigns are relative to
    %obj.trialData{sessionIdx}(trialIdx).stim_start, where
    %obj.trialData{sId}(trialIdx).stim_start is the time at the beginning
    %of frame startFrameNumber (now set to 0).
    switch priorOrPosterior
        case {'priorUnderCurrentGPHypers','priorUnderHyperpriorMean'}
            %Calculate times:
            %Want last entry to be at least lagMax before the last entry in firing.
            maxFrameNumber = size(firing,1) - startFrameNumber - ceil(lagMax/dt);
            %The above follows because the condition that there is a frame which
            %starts at time maxFrameNumber * dt + dt * ceil(lagmax/dt) (so we have a corresponding
            %rate to do a correlation with) implies that the value of
            %maxFrameNumber must allow this to be so; we've already fixed the
            %length of firing and lagMax.  Thus, size(firing,1) - (startFrameNumber
            %- 1) (which is the frame number of the last element of firing) must be
            %at least (maxFrameNumber + ceil(lagMax/dt) + 1), since the value of
            %the accumulator at maxFrameNumber * dt will be compared with the
            %firing rate in this frame.  Solving the inequality gives the equation
            %above.
            
            times = dt * (0:1:maxFrameNumber)';
            %Above: instants at which we're going to calculate the prior mean;
            %correspond to the start of frame startFrameNumber and integer numbers
            %of frames thereafter.
            
            %Calculate the prior mean accumulator trajectory (stimulus, no neural
            %data, NO RESPONSE) (The last is a key difference from the analysis of
            %T. Hanks)
            switch priorOrPosterior
                case 'priorUnderCurrentGPHypers'
                    priorOrPosteriorLatentTrajectory = meanAndDerivs(times,clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers);
                case 'priorUnderHyperpriorMean'
                    priorOrPosteriorLatentTrajectory = meanAndDerivs(times,clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypersPriorMean);
            end
            %Calculate the smoothed empirical firing rate for this cell, with the
            %same starting frame  (this will have to go at least lagmax longer than
            %the end of the priorMeanTrajectory vector.
            
        case 'posterior'
            times = obj.modelAltDataRepresentations{sessionIdx,modelNumber}.timesForLatentsRelativeToStimStart{trialIdx};  %Has all of the entries: in reference to stimStart.
            priorOrPosteriorLatentTrajectory = obj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories{trialIdx};
    end
    
    for cellIdx = 1:nCells
        thisCellAllSpikesInLocalTime = allSpikes{cellIdx} - perTrialTimeOffsets(trialIdx) - startTime - obj.trialData{sessionIdx}(trialIdx).stim_start;
        %Above: allSpikes{cellIdx} -(obj.trialData{sessionIdx}(trialIdx).state_0_exits - startTime) - startTime - obj.trialData{sessionIdx}(trialIdx).stim_start
        %     = allSpikes{cellIdx} - obj.trialData{sessionIdx}(trialIdx).state_0_exits - obj.trialData{sessionIdx}(trialIdx).stim_start
        %     = allSpikesNowInSameReferenceAsTheClicks
        %
        %Thus time zero in this representation corresponds to the same time
        %as the beginning of the (startFrameNumber)th frame.
        %
        %state_0_exits is the t = 0 reference for all of the array_data
        %fields.
        
        
        %Figure out frame start/stop times for firing
        frameNumbers = (1 - (startFrameNumber - 1)) : 1 : (size(firing, 1) - (startFrameNumber -1) );
        startFrameTimes = (frameNumbers - 1) * dt;
        endFrameTimes = (frameNumbers) * dt;
        %Above: Ascending order
        
        switch filteringScheme
            case 'stepwiseOld'
                %% Old rate calculation:
                %Should we be picking up some bleed to the beginning and end here, such
                %that the filtering will have a bit of extra length?
                %Alternatively, we can cut down the ammount of history we're using to
                %do the fitting.
                
                %Let's use 1/gap as the estimate of the firing rate.  This is pretty
                %easy and just requires knowing when the last spike was.
                
                %As a hack, we'll use 2x the gap on the ends.  Not the best, but seems
                %reasonable.
                
                %For the period up until the next spike, the assumed rate will be the
                %1/(next-this): so all indices this:next-1 will have this value
                
                rates = zeros(size(firing,1),1);
                spikeIndicesInFiring = find(firing(:,cellIdx));
                %reconstitutedSpikeTimes = spikeIndicesInFiring * dt + obj.trialData{sessionIdx}(trialIdx).stim_start;
                %spikeIndicesInFiringAugmented = [-spikeIndicesInFiring(1)+1,spikeIndicesInFiring,size(firing,1) + (size(firing,1) - spikeIndicesInFiring(end)) -1
                startIndices = vertcat(1,spikeIndicesInFiring);
                endIndices = vertcat(spikeIndicesInFiring-1,size(firing,1));
                for rateIdx = 1:numel(startIndices)
                    if endIndices(rateIdx) == startIndices(rateIdx)
                        rates(startIndices(rateIdx),1) = 1;
                        
                    elseif rateIdx == 1 || rateIdx == numel(startIndices)
                        rates(startIndices(rateIdx):endIndices(rateIdx),1) = 1/(2*(endIndices(rateIdx) - startIndices(rateIdx)));  %The factor of two here is to account for the uniform expectation of where the break falls in an interval (correct?)
                    else
                        rates(startIndices(rateIdx):endIndices(rateIdx),1) = 1/(  (endIndices(rateIdx) - startIndices(rateIdx)));
                    end
                end
                
                
            case 'stepwiseNew'
                %% New Firing rate calculation:
                %USES THE GLOBAL TIMESTAMPS OF THE SPIKE TIMES TO CALCULATE THE FIRING
                %RATE OF THE NEURON
                
                
                
                rates = nan(size(firing,1),1);
                
                %Localize spikes:
                relevantSpikes = thisCellAllSpikesInLocalTime(max(1,find(thisCellAllSpikesInLocalTime > startFrameTimes(1),1,'first') - 1) : find(thisCellAllSpikesInLocalTime < endFrameTimes(end),1,'last') + 1);
                minFlag = false;
                if find(thisCellAllSpikesInLocalTime > startFrameTimes(1),1,'first') == 1
                    minFlag = true;
                end
                spikesAndIntervalsLogical = startFrameTimes(ones(numel(relevantSpikes),1),:)' <= relevantSpikes(:,ones(size(firing,1),1))' & endFrameTimes(ones(numel(relevantSpikes),1),:)' > relevantSpikes(:,ones(size(firing,1),1))';
                %i,j is true if spike j occurred within time interval i [start,end).
                anySpikesInThisInterval = any(spikesAndIntervalsLogical,2);
                
                %         lastSpike = nan(size(firing,1),1);
                %         containsSpikes = cell(size(firing,1),1);
                %         nextSpike = nan(size(firing,1),1);
                
                %% 1/interval rate calculation:
                
                if ~any(anySpikesInThisInterval)
                    assert(numel(relevantSpikes) == 2)
                    rates(:) = 1/(relevantSpikes(end) - relevantSpikes(1));
                else
                    
                    %First pass:
                    framesWithSpikes = find(anySpikesInThisInterval);
                    for idxInFWS = 1: numel(framesWithSpikes)-1
                        if idxInFWS == 1
                            if ~minFlag
                                %Fill in the start
                                previousStartingSpikeOccurredAtTime             = relevantSpikes(find(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS),:),1,'first')        - 1);
                                previousIntervalEndingSpikeOccurredAtTime       = relevantSpikes(find(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS),:),1,'first')           );
                                rates(1:framesWithSpikes(idxInFWS)-1,1) = 1/ (previousIntervalEndingSpikeOccurredAtTime - previousStartingSpikeOccurredAtTime);
                            else %minFlag
                                
                            end
                        end
                        startingSpikeOccurredAtTime = relevantSpikes(find(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS),:),1,'last'));
                        endingSpikeOccurredAtTime = relevantSpikes(find(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS + 1),:),1,'first'));
                        rates(framesWithSpikes(idxInFWS) + 1 : framesWithSpikes(idxInFWS + 1) - 1,  1) = 1 / (endingSpikeOccurredAtTime - startingSpikeOccurredAtTime )  ;
                        
                    end
                    %Fill in the end
                    lastSpikeIdx = find(spikesAndIntervalsLogical(framesWithSpikes(end),:),1,'last');
                    lastSpikeOccurredAtTime = relevantSpikes(lastSpikeIdx);
                    nextSpikeOccurredAtTime = relevantSpikes(lastSpikeIdx + 1);
                    rates(framesWithSpikes(end) + 1 : end,  1) = 1/(nextSpikeOccurredAtTime -  lastSpikeOccurredAtTime);
                    
                    %Second Pass:
                    for idxInFWS = 1:numel(framesWithSpikes)
                        if minFlag && idxInFWS == 1
                            %figure out how to calculate the rate otherwise
                            fprintf('NEED TO FIGURE OUT HOW TO CALCULATE RATE IN A SPECIAL CASE!\n\n')
                        else
                            intervalBreaksWithinFrame = [startFrameTimes(framesWithSpikes(idxInFWS)), relevantSpikes(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS),:)), endFrameTimes(framesWithSpikes(idxInFWS))];
                            proportions = diff(intervalBreaksWithinFrame) ./ dt;
                            
                            previousSpikeTime = relevantSpikes(find(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS),:),1,'first')        - 1);
                            nextSpikeTime = relevantSpikes(find(spikesAndIntervalsLogical(framesWithSpikes(end),:),1,'last') + 1);
                            intervalBreaksGlobal = [previousSpikeTime, relevantSpikes(spikesAndIntervalsLogical(framesWithSpikes(idxInFWS),:)), nextSpikeTime];
                            ratesInIntervals = 1./diff(intervalBreaksGlobal);
                            
                            
                            rates(framesWithSpikes(idxInFWS),1) = ratesInIntervals * proportions';
                        end
                    end
                end
            case 'halfGaussian'
                %% Filtering calculation:
                %Here, looking at rates at starts of frames: not exactly
                %equivalent with what I was doing above (which took account
                %of spikes WITHIN the frame) but should be OK.
                rates = brodyFilterConv(startFrameTimes,thisCellAllSpikesInLocalTime,filterParams,'fastHalfGauss');
            otherwise
                error('Unrecognized filtering scheme!')
        end
        
        
        
        
        
        
        %Store
        if any(isnan(rates))
            fprintf('Some rates are NaN! Fix this!\n\n')
        end
        smoothedFiringStor{cellIdx,trialIdx} = rates;
        
    end
    
    %Store the representations
    accumulatorValuesStor{trialIdx} = priorOrPosteriorLatentTrajectory;
    
end


%% Do the calculations and make the requested plot
switch  priorOrPosterior
    case {'priorUnderCurrentGPHypers','priorUnderHyperpriorMean'}
        lags = lagMin:dLag:lagMax;
    case 'posterior'
        lags = 0;
        assert(lagMin == 0,'lagMin is non-zero!');
end
corrCoeffsOut = zeros(size(lags));
inputVector = -10:0.25:10;  %values for the accumulator in the regression


accumulatorValuesConcatenated = [];

for trialIdx = 1:numel(obj.trialData{sessionIdx})
    accumulatorValuesConcatenated = vertcat(accumulatorValuesConcatenated,accumulatorValuesStor{trialIdx});
end
fxnToFit = @(betaVals,x) betaVals(1) + betaVals(2) ./ (1 + exp(-betaVals(3) .* (x-betaVals(4))));


for lagIdx = 1:numel(lags)
    tic
    smoothedFiringRatesConcatenated = zeros(size(accumulatorValuesConcatenated));
    assemblyIdx = 0;
    numFramesLag = (dLag / dt) * (lagIdx - 1);
    
    %Get all of the inputs (accumulator states) and outputs (smoothed firing rates) and compute the correlation.
    
    %The different part here is that we're going to be doing this for
    %different lags: so the entries in smoothedFiringRatesConcatenated are
    %going to be different.
    
    for trialIdx = 1:numel(obj.trialData{sessionIdx})
        smoothedFiringRatesConcatenated(assemblyIdx + 1: assemblyIdx + numel(accumulatorValuesStor{trialIdx})) = smoothedFiringStor{trialIdx}(1+numFramesLag:numel(accumulatorValuesStor{trialIdx}) + numFramesLag);
        assemblyIdx = assemblyIdx + numel(accumulatorValuesStor{trialIdx});
    end
    
    
    bothNotNan = ~isnan(accumulatorValuesConcatenated) & ~isnan(smoothedFiringRatesConcatenated);
    
    rTemp = corrcoef(accumulatorValuesConcatenated(bothNotNan),smoothedFiringRatesConcatenated(bothNotNan));
    corrCoeffsOut(lagIdx) = rTemp(2,1);
    
    %Fit a sigmoid curve to this data
    %%Old Method:
    %fittedCurveObject = fit(accumulatorValuesConcatenated,smoothedFiringRatesConcatenated,'poly1');
    maxSmoothCat = max( smoothedFiringRatesConcatenated);
    beta0 = [1/4 * maxSmoothCat,rand * 1/4 * maxSmoothCat,randn * 10,randn * 0.01];
    %New Method:
    betaOut = nlinfit(accumulatorValuesConcatenated,smoothedFiringRatesConcatenated,fxnToFit,beta0,statset('Display','iter'));
    
    
    %Plot the scatterplot and the curve
    figure(1)
    clf
    hold all
    scatter(accumulatorValuesConcatenated,smoothedFiringRatesConcatenated)
    
    %%Old plotting:
    %plot(fittedCurveObject)
    
    %New Plotting
    inputVector = min(accumulatorValuesConcatenated):(1/100 *  (max(accumulatorValuesConcatenated) - min(accumulatorValuesConcatenated)) ):max(accumulatorValuesConcatenated);
    plot(inputVector,fxnToFit(betaOut,inputVector));
    switch priorOrPosterior
        case 'priorUnderCurrentGPHypers'
            xlabelString = 'Current GP hyperparameters Prior Accumulator Value';
        case 'priorUnderHyperpriorMean'
            xlabelString = 'Hyperprior mean GP hyperparameters Prior Accumulator Value';
        case 'posterior'
            xlabelString = 'MAP Accumulator Value';
    end
    xlabel(xlabelString)
    ylabel('Approximate firing rate')
    title(strcat('Accumulator value for cell ', obj.cellsIDCell{cellIdx},' at lag ',num2str(lags(lagIdx)*1000),' ms. r = ',num2str(corrCoeffsOut(lagIdx)),'.'))
    
    figure(2)
    clf
    mSFRC = max(smoothedFiringRatesConcatenated) + 0.001;
    binCenters = {-7:0.25:7,mSFRC/50 : mSFRC/25 : mSFRC * 49/50};
    hist3([accumulatorValuesConcatenated,smoothedFiringRatesConcatenated],binCenters)
    heatMapValues = hist3([accumulatorValuesConcatenated,smoothedFiringRatesConcatenated],binCenters);
    xlabel(xlabelString)
    ylabel('Approximate firing rate')
    
    figure(3)
    clf
    [X,Y] = meshgrid(binCenters{1},binCenters{2});
    contourf(X',Y',heatMapValues)
    xlabel(xlabelString)
    ylabel('Approximate firing rate')
    
    disp(num2str(lags(lagIdx)/dt))
    disp(betaOut)
    pause(5-toc);
    
end



%% Return the possible outputs
if nExtraArgsOut == 2
    varargout{1} = lags;
    varargout{2} = corrCoeffsOut;
end