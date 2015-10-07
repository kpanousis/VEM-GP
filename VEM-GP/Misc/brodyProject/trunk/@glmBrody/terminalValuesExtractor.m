function [terminalTime, terminalGPPriorMeanAndVariance, terminalCurrentPriorMeanAndVariance, terminalCurrentMAPValue, terminalFiringRate, pokedR] = terminalValuesExtractor(obj,modelNumber,sessionIdx,maxTrialNumber, firingRateConvolutionFilterType,filterParams)
%TERMINALVALUESEXTRACTOR Returns terminal values for regression.
%terminalValuesExtractor
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%23 April, 2014
%
%Method of glmBrody class which, for a given model number and session,
%extracts the final timestamp (relative to stimulus start) terminal values
%of the GP prior under the mean hyperparameters of the original hyperprior,
%the terminal value of the prior under the current GP hyperparameters, the
%terminal value of the MAP latent trajectory, and the final firing rate
%(smoothed) for each trial.
%
%Inputs: 
%   modelNumber: The model number of glmBrodyLaplace model examined in the
%       the glmBrody class object.
%   sessionIdx: Which session is being examined.  The method could be used
%       once for each of several sessions.
%   maxTrialNumber: The largest number of trials to consider; may use fewer
%       trials if fewer have been used in iterative training.
%   firingRateConvolutionFilterType: the type of convolution function being
%       used in filtering the impulse train.
%   filterParams: Parameters of this filter to be passed to
%       brodyFilterConv.
%   
%Outputs:
%   terminalTime: nTrials x 1 array of terminal timestamps (with respect to
%       obj.trialData{sessionIdx}(trialNumber).stim_start for each trial;
%       for the following, let this value be T.
%   terminalGPPriorMeanAndVariance: nTrials x 2 array: The mean and 
%       variance of the prior on a(T), DEPENDENT ONLY ON STIMULUS
%       INFORMATION, where the Theta_GP = the mean values of the
%       hyperprior.
%   terminalCurrentPriorMeanAndVariance: nTrials x 2 array: Mean and 
%       variance of the prior on a(T), using Theta_GP = the CURRENT values
%       for the GP hyperparameters; this estimate is thus dependent on the
%       neural firing through them (because they have been fitted
%       iteratively, such that they are dependent on firing activity).
%   terminalCurrentMAPValue: nTrials x 1 array: MAP latent value a(T), 
%       where the prior is defined according to the current GP
%       hyperparameters.
%   terminalFiringRate: nTrials x nCells array (where cells are ordered
%       according to their appearance in obj.cellsIdCell{sessionIdx}): an
%       estimated firing rate of the cells at T.
%   pokedR: nTrials x 1 vector of choices (0 = left, +1 = right).

plottingWholeTrajectories = false;

nTrials = min(maxTrialNumber,numel(obj.trialData{sessionIdx}));
fprintf('Extracting data from the first %s trials.',num2str(nTrials))

assert(strcmp(obj.modelTypes{1,modelNumber},'glmBrodyLaplace'),'Model selected is not of type glmBrodyLaplace!')


%Preallocate storage arrays:
terminalTime = zeros(nTrials,1);
terminalGPPriorMeanAndVariance = zeros(nTrials,1);
terminalCurrentPriorMeanAndVariance = zeros(nTrials,1);
terminalCurrentMAPValue = zeros(nTrials,1);
terminalFiringRate = zeros(nTrials,1);
pokedR = false(nTrials,1);


%% Find all spike times for (the non-pruned data in) this session
[allSpikes, ~,~,perTrialTimeOffsets] = obj.listAllBandS_Method(sessionIdx,numel(obj.trialData{sessionIdx}),0);  %Returns a 1 x nCells cell array, where the ith entry is the set of spike times for all trials.


%% Trial-specific calculations
for trialIdx = 1:nTrials
    %Calculate the times of interest
    times = obj.modelAltDataRepresentations{sessionIdx,modelNumber}.timesForLatentsRelativeToStimStart{trialIdx};  %Has all of the entries: in reference to stim_start.
    terminalTime(trialIdx) = times(end);
    
    pokedR(trialIdx,1) = obj.trialData{sessionIdx}(trialIdx).pokedR;

    
    if plottingWholeTrajectories
              aMAP  = obj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories{trialIdx};
        
        [~, ~, clickTimes, clickSigns] = obj.trialParser(sessionIdx,trialIdx,[0,0],1/obj.modelGlobalParameters{modelNumber}.samplingRate);
        %clickTimes and clickSigns are relative to
        %obj.trialData{sessionIdx}(trialIdx).stim_start, where
        %obj.trialData{sId}(trialIdx).stim_start is the time at the beginning
        %of frame startFrameNumber (now set to 0).
        
        
        %Compute the current values of the GP mean and variance, along with the
        %value of the MAP estimate.
        GPPriorMean             = meanAndDerivs(         times,   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypersPriorMean);
        GPPriorVariance         = diag(kernelAndDerivs(       times,   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypersPriorMean));
        
        CurrentPriorMean        = meanAndDerivs(    times,   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers);
        CurrentPriorVariance    = diag(kernelAndDerivs(  times,   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers));
        for cellIdx = 1:numel(allSpikes)
            thisCellAllSpikesInLocalTime = allSpikes{cellIdx} - perTrialTimeOffsets(trialIdx) - obj.trialData{sessionIdx}(trialIdx).stim_start;  %Puts the spike times in the same time reference (i.e., stim_start is zero)
            
            rates = brodyFilterConv(times,thisCellAllSpikesInLocalTime,filterParams,firingRateConvolutionFilterType);
        end
        
        assert(numel(allSpikes) == 1,'Hacky code: assumes only one cell. Fix!')
        
        figure(1)
        clf
        subplot(2,1,1)
        hold all
        plot(times,rates,'r-')
        axis([0,max(times),0,20])

        if pokedR(trialIdx,1)
            scatter(times(end),5,'rx')
        else
            scatter(times(end),10,'bo')
        end
        
        subplot(2,1,2)
        hold all
        plot(times,GPPriorMean,'g-')
        plot(times,GPPriorMean + sqrt(GPPriorVariance),'g--')
        plot(times,GPPriorMean - sqrt(GPPriorVariance),'g--')
        plot(times,aMAP,'b-')
        axis([0,max(times),-10,10])

        pause(3)
        
        
    else
        terminalCurrentMAPValue(trialIdx,1)  = obj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories{trialIdx}(end);
        
        [~, ~, clickTimes, clickSigns] = obj.trialParser(sessionIdx,trialIdx,[0,0],1/obj.modelGlobalParameters{modelNumber}.samplingRate);
        %clickTimes and clickSigns are relative to
        %obj.trialData{sessionIdx}(trialIdx).stim_start, where
        %obj.trialData{sId}(trialIdx).stim_start is the time at the beginning
        %of frame startFrameNumber (now set to 0).
        
        
        %Compute the current values of the GP mean and variance, along with the
        %value of the MAP estimate.
        terminalGPPriorMeanAndVariance(trialIdx,1) = meanAndDerivs(         terminalTime(trialIdx),   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypersPriorMean);
        terminalGPPriorMeanAndVariance(trialIdx,2) = kernelAndDerivs(       terminalTime(trialIdx),   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypersPriorMean);
        
        terminalCurrentPriorMeanAndVariance(trialIdx,1) = meanAndDerivs(    terminalTime(trialIdx),   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers);
        terminalCurrentPriorMeanAndVariance(trialIdx,2) = kernelAndDerivs(  terminalTime(trialIdx),   clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers);
        for cellIdx = 1:numel(allSpikes)
            thisCellAllSpikesInLocalTime = allSpikes{cellIdx} - perTrialTimeOffsets(trialIdx) - obj.trialData{sessionIdx}(trialIdx).stim_start;  %Puts the spike times in the same time reference (i.e., stim_start is zero)
            
            rates = brodyFilterConv(terminalTime(trialIdx),thisCellAllSpikesInLocalTime,filterParams,firingRateConvolutionFilterType);
            terminalFiringRate(trialIdx,cellIdx)  = rates(end);
        end
    
    end
    
    
    
end