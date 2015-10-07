function [klDivsVec,varargout] = calculateKLTrialsPosterior(obj,modelNumber,referenceObj,sessionsIndices, trialIndices)
%CALCULATEKLTRIALSPOSTERIOR calculates the KL divergences between glmBrodyLaplace models 
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%28 November, 2014
%
%Inputs:
%   modelNumber: indexing into glmBrodyLaplace object which model will
%       be examined.
%   referenceObj: Another glmBrody object, which will serve as the
%       reference distribution in the KL divergence.
%   sessionsIndices: vector of sessions to use in calculating the KL divergences; []
%       --> all.
%   trialIndices:  If sessionsIndices is [], this must also be [] and all 
%       trials will be used, else error.  If sessionsIndices is non-empty, 
%       cell array of same size as sessionsIndices.   
%   
%Outputs: 
%   klDivsVec: a vector of KL divergences between the two induced
%       distributions, one entry per trial used.
%   varargout: if two entries, returns column vectors of corresponding sessionIdx,
%       trialIdx to the entries in klDivsVec
				    error('The usage of kAISInverse2.m below is out of date; it no longer calculates the logdet.')

%% Model checking
assert(isa(referenceObj,'glmBrody'),'referenceObj input is not a glmBrody model!')

modelTypeThisStruct = obj.modelTypes{modelNumber};
modelTypeReferenceStruct = referenceObj.modelTypes{modelNumber};
assert(any(strcmp(modelTypeThisStruct,{'glmBrodyLaplace','glmBrodyLaplaceSelfOnly'})),'Unsupported model type!')
assert(strcmp(modelTypeReferenceStruct,modelTypeThisStruct),'The two models must have the same model type!')
assert(isfield(obj.modelGlobalParameters{modelNumber},'gpHypers'),'The chosen glmBrodyLaplace(SelfOnly) model does not have a gpHypers field!')
assert(isfield(referenceObj.modelGlobalParameters{modelNumber},'gpHypers'),'The reference glmBrodyLaplace(SelfOnly) model does not have a gpHypers field!')

referenceGPParameters = referenceObj.modelGlobalParameters{modelNumber}.gpHypers;

%% Output checking
nExtraArgsOut = max(0,nargout - 1);
if ~any(nExtraArgsOut == [0,2])
    error('Improper number of outputs from calculateKLTrialsPrior!')
end
%% Input handling
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

%Check that the sizes of the modelAltDataRepresentations fields
%agree.

%Check # sessions is the same (and # models, actually)
assert(all(size(obj.modelAltDataRepresentations) == size(referenceObj.modelAltDataRepresentations)),'Non-matching sizes of the modelAltDataRepresentations fields of the current and reference objects.')
%for sessionIdx = sessionsIndices
%    assert(all(size(obj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories)...
%        == size(referenceObj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories)),...
%        'Current and reference objects do not have matching sized modelAltDataRepresentations.posteriorLatentTrajectories fields!')
%end
clear sessionIdx


%% Preallocate storage
nTrialsToCompute = 0;
for sessionIdx = sessionsIndices'
    nTrialsToCompute = nTrialsToCompute + numel(trialIndices{sessionIdx});
end
    
klDivsVec       = zeros(nTrialsToCompute,1);
sessionsOutVec  = zeros(nTrialsToCompute,1);
trialsOutVec    = zeros(nTrialsToCompute,1);

%%Compute the KL divergences

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
        
        
        %%%% Compute the prior precisions and covariances
        %mu = meanAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,referenceGPParameters);
        SigmaMat = kernelAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,referenceGPParameters); %weird that clickSigns is an input to K... shouldn't be.
        %SigmaInv = kAISInverse2(timesForLatentCalculation,clickTimes,referenceGPParameters);
        [SigmaInv,~,logDetSigmaInv] = kAISInverse2(timesForLatentCalculation,clickTimes,referenceGPParameters);
        
        %Compute m, Kinv for the chosen parameters
        %mVec = meanAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers);
        KMat = kernelAndDerivs(timesForLatentCalculation,clickTimes,clickSigns,obj.modelGlobalParameters{modelNumber}.gpHypers); %weird that clickSigns is an input to K... shouldn't be.
        [KInv,~,logDetKInv] = kAISInverse2(timesForLatentCalculation,clickTimes,obj.modelGlobalParameters{modelNumber}.gpHypers);
        
        %Pull the posterior means from their computed values at the current
        %fitting state.
        mu_aMAP = referenceObj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories{trialIdx};
        m_aMAP  =          obj.modelAltDataRepresentations{sessionIdx,modelNumber}.posteriorLatentTrajectories{trialIdx};
        
        %Calculate the latent firing rate values for each of the trials for
        %each model.
        %These should be of dimension nNeurons
        [latentFiringRate_reference ,~,~] = referenceObj.calculateFiringRates_BrodyGLM(firing, mu_aMAP, startFrameNumber, modelNumber, sessionIdx);
        [latentFiringRate_thisObj   ,~,~] =          obj.calculateFiringRates_BrodyGLM(firing, m_aMAP , startFrameNumber, modelNumber, sessionIdx);
        
        %Calculate the induced precision from the firing rate calculateion
        LambdaMat_reference = diag(sum(   (referenceObj.modelSessionParameters{sessionIdx,modelNumber}.C(:,ones(size(KMat,1),1))') .^ 2 .* latentFiringRate_reference, 2));
        LambdaMat_thisObj   = diag(sum(   (         obj.modelSessionParameters{sessionIdx,modelNumber}.C(:,ones(size(KMat,1),1))') .^ 2 .* latentFiringRate_thisObj,   2));
        
        %The KL divergence is then given by
        %1/2 * [ log|sigmaInv| - log|Kinv|] + 1/2 tr(Kinv*sigma) - T/2 +
        %1/2 (m - mu)' * Kinv * (m-mu);
%         %klDivsVec(idxInKlDivsVec) =   1/2 * sum(log(eig(SigmaInv)))...
%         %                            - 1/2 * sum(log(eig(Kinv)))...
%         %                            + 1/2 * trace(Kinv * SigmaMat)...
%         %                            - numel(timesForLatentCalculation)/2 ...
%         %                            + 1/2 * (mVec - mu)' * Kinv * (mVec - mu);
%         klDivsVec(idxInKlDivsVec) =   1/2 * logDetSigmaInv...
%                                     - 1/2 * logDetKInv...
%                                     + 1/2 * trace(KInv * SigmaMat)...
%                                     - numel(timesForLatentCalculation)/2 ...
%                                     + 1/2 * (mVec - mu)' * KInv * (mVec - mu);
        %Unfortunately, we now have to compute this with the posterior
        %precision, not the prior precision, and further, this is going to
        %royally screw up the trace term.
        
        precisionReference = SigmaInv + LambdaMat_reference;  
        precisionThisObj   = KInv     + LambdaMat_thisObj;
        covariancePosteriorReference = inv(precisionReference);%Typically, this should be sufficiently regularized by the addition of the PSD lambda that we're not going to have too hard of a time inverting it.
        
        
        klDivsVec(idxInKlDivsVec) =   1/2 * sum(log(eig(precisionReference)))...
                                    - 1/2 * sum(log(eig(precisionThisObj)))...
                                    + 1/2 * trace(precisionThisObj * covariancePosteriorReference)...
                                    - numel(timesForLatentCalculation)/2 ...
                                    + 1/2 * (m_aMAP - mu_aMAP)' * precisionThisObj * (m_aMAP - mu_aMAP);

                                
        sessionsOutVec(idxInKlDivsVec) = sessionIdx;
        trialsOutVec(idxInKlDivsVec) = trialIdx;
    end
end



%% Extra outputs
if nExtraArgsOut == 2
    varargout{1} = sessionsOutVec;
    varargout{2} = trialsOutVec;
end
