function [klDivsVec,varargout] = calculateKLTrialsPrior(obj,modelNumber,referenceGPParameters,sessionsIndices, trialIndices)
%CALCULATEKLTRIALSPRIOR calculates the KL divergences between glmBrodyLaplace models 
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%28 November, 2014
%
%Inputs:
%   modelNumber: indexing into glmBrodyLaplace object which model will
%       be examined.
%   referenceGPParameters: Define the reference distribution from which we
%       will calculate the KL divergence
%   sessionsIndices: vector of sessions to use in calculating the KL divergences; []
%       --> all.
%   trialIndices:  If sessionsIndices is [], this must also be [] and all 
%       trials will be used, else error.  If sessionsIndices is non-empty, 
%       cell array of same size as sessionsIndices.   
%   
%Outputs: 
%   klDivsVec: a vector of KL divergences between the two induced
%   distributions, one entry per trial used.
%   varargout: if two entries, returns column vectors of corresponding sessionIdx,
%       trialIdx to the entries in klDivsVec

				    error('The usage of kAISInverse2.m below is out of date; it no longer calculates the logdet.')

%% Model checking
assert(any(strcmp(obj.modelTypes{modelNumber},{'glmBrodyLaplace','glmBrodyLaplaceSelfOnly'})),'Unsupported model type!')
assert(isfield(obj.modelGlobalParameters{modelNumber},'gpHypers'),'The chosen glmBrodyLaplace(SelfOnly) model does not have a gpHypers field!')

%% Output checking
nExtraArgsOut = max(0,nargout - 1);
if ~any(nExtraArgsOut == [0,2])
    error('Improper number of outputs from calculateKLTrialsPrior!')
end

%% Input handling
assert(isnumeric(referenceGPParameters) && all(size(referenceGPParameters) == [6,1]),'Input referenceGPParameters is not correctly formatted! Must be 6x1 vector!')

if isempty(sessionsIndices)
    sessionsIndices = (1:numel(obj.trialData))';
    trialIndices = cell(size(sessionsIndices));
    for sessionIdx = 1:numel(sessionsIndices)
        trialIndices{sessionIdx} = 1:numel(obj.trialData{sessionIdx});
    end
else
    assert(all(size(sessionsIndices) == size(trialIndices)),'Non-matching session and trial indices array sizes!')
    assert(iscell(trialIndices),'trialIndices must be a cell array!')
    for sessionIdx = 1:numel(sessionsIndices)
        assert(~any(trialIndices{sessionIdx} > numel(obj.trialData{sessionIdx})),sprintf('An entry in session %s of trialIndices is larger than the number of trials in the session!',sessionIdx))
        assert(size(trialIndices{sessionIdx},1) == 1, 'Each cell of the trialIndices cell array must contain a ROW vector of indices!')
    end
end

%% Preallocate storage
nTrialsToCompute = 0;
for sessionIdx = sessionsIndices'
    nTrialsToCompute = nTrialsToCompute + numel(trialIndices{sessionIdx});
end
    
klDivsVec       = zeros(nTrialsToCompute,1);
sessionsOutVec  = zeros(nTrialsToCompute,1);
trialsOutVec    = zeros(nTrialsToCompute,1);
%% Compute the KL divergences

idxInKlDivsVec = 0;
for sessionIdx = sessionsIndices'
    for trialIdx = trialIndices{sessionIdx}
        idxInKlDivsVec = idxInKlDivsVec + 1;
        %Overall:
        %Compute the per-trial KL divergence between the reference
        %distribution and model fitted distribution at discretized times
        
        %Specific:
        %Compute the click times and signs
        bleedTimes = [obj.modelGlobalParameters{modelNumber}.relevantSelfHistory,obj.modelGlobalParameters{modelNumber}.maxLagConsidered];
        [firing,startFrameNumber,clickTimes,clickSigns] = obj.trialParser(sessionIdx,trialIdx,bleedTimes,1/obj.modelGlobalParameters{modelNumber}.samplingRate);
        
        %Compute the times of relevance for input into the mean and kernel
        %functions
        timesForLatentCalculation = (1/obj.modelGlobalParameters{modelNumber}.samplingRate * (0:1:(size(firing,1) - startFrameNumber)))';
        
        %Compute mu, sigma, and sigmaInv for the reference parameters
        mu = meanAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,referenceGPParameters);
        SigmaMat = kernelAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,referenceGPParameters); %weird that clickSigns is an input to K... shouldn't be.
        %SigmaInv = kAISInverse2(timesForLatentCalculation,clickTimes,referenceGPParameters);
        [SigmaInv,~,logDetSigmaInv] = kAISInverse2(timesForLatentCalculation,clickTimes,referenceGPParameters);
        
        %Compute m, Kinv for the chosen parameters
        mVec = meanAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers);
        %KMat = kernelAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers); %weird that clickSigns is an input to K... shouldn't be.
        [KInv,~,logDetKInv] = kAISInverse2(timesForLatentCalculation,clickTimes,obj.modelGlobalParameters{modelNumber}.gpHypers);
        
        %The KL divergence is then given by
        %1/2 * [ log|sigmaInv| - log|Kinv|] + 1/2 tr(Kinv*sigma) - T/2 +
        %1/2 (m - mu)' * Kinv * (m-mu);
        %klDivsVec(idxInKlDivsVec) =   1/2 * sum(log(eig(SigmaInv)))...
        %                            - 1/2 * sum(log(eig(Kinv)))...
        %                            + 1/2 * trace(Kinv * SigmaMat)...
        %                            - numel(timesForLatentCalculation)/2 ...
        %                            + 1/2 * (mVec - mu)' * Kinv * (mVec - mu);
        klDivsVec(idxInKlDivsVec) =   1/2 * logDetSigmaInv...
                                    - 1/2 * logDetKInv...
                                    + 1/2 * trace(KInv * SigmaMat)...
                                    - numel(timesForLatentCalculation)/2 ...
                                    + 1/2 * (mVec - mu)' * KInv * (mVec - mu);
                                
        sessionsOutVec(idxInKlDivsVec) = sessionIdx;
        trialsOutVec(idxInKlDivsVec) = trialIdx;
    end
end


if nExtraArgsOut == 2
    varargout{1} = sessionsOutVec;
    varargout{2} = trialsOutVec;
end

end
