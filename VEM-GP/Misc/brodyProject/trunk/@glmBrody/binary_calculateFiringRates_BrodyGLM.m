function [firingRates,varargout] = binary_calculateFiringRates_BrodyGLM(obj,firing,latentTrajectory,startFrameNumber, modelNumber, sessionIdx,varargin)
%BINARY_CALCULATEFIRINGRATES_BRODYGLM Returns the inferred firing rates (lambdas).
%binary_calculateFiringRates_BrodyGLM.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%12 March, 2014
%
%Method which calculates the firing rate of the GLM model, using the latent
%trajectory given.  The basic expression here is log(lambda) = C * a[t] + d
%+ W * h(y[t-k:t-1]);  We thus need the current value of a[t] (at the
%beginning of the frame), k steps of past history of firing, the
%self-basis, and the weights (C,d,W).
%Inputs: 
%   firing: (T_calculate + (startFrameNumber-1)) x nNeurons BINARY matrix of spiking
%       activity.
%   latentTrajectory: T_calculate x 1 latent GP model capturing evidence
%       accumulation.
%   startFrameNumber: Index into firing which lets us know where to start.
%       Somewhat redundant with the implicit value given by the sizes of
%       firing and latentTrajectory
%   modelNumber: Index into the class object's structures, such that we can
%       access stored parameters.
%Outputs:
%   firingRates: The T_calculate x nNeurons matrix of firing rates
%
% NOTE: From my previous work, it appears that lambda has been selected as
% the rate of firing PER BIN, and thus is dependent on the choice of bin
% length.  Clearly, we can just rescale everything by the bin length, but
% this will be necessary for a data interpretation step.
nExtraArgsOut = max(0,nargout - 1);
assert(nExtraArgsOut <= 2, 'Too many outputs!')

nExtraArgsIn = max(0,nargin-6);
assert(nExtraArgsIn == 0 || nExtraArgsIn == 3,'Unrecognized input format; wrong number of inputs.')



%% Input checking
assert(size(firing,1) - startFrameNumber + 1 == size(latentTrajectory,1),'Sizes of inputs don''t match! Check this!')
assert(size(latentTrajectory,2) == 1, 'Passing more than one latent state! Not currently equipped to deal with this case!')
assert(startFrameNumber > 0,'Start frame must be properly specified.')
assert(numel(obj.modelTypes) >= modelNumber  && any(strcmp(obj.modelTypes(modelNumber),{'glmBrodyLaplace','glmBrodyLaplaceSelfOnly'})), 'The model number given does not exist or is not a glmBrodyLaplace model.')

assert(islogical(firing),'Firing is non-logical!  This shouldn''t happen in this version of the code!')
%% Pull down to local variables for notational simplicity
%(e.g., the weights)
if nExtraArgsIn == 0
    C = obj.modelSessionParameters{sessionIdx, modelNumber}.C;  %nNeurons x 1
    d = obj.modelSessionParameters{sessionIdx, modelNumber}.d;  %nNeurons x 1
    W = obj.modelSessionParameters{sessionIdx, modelNumber}.W;  %nNeurons,nNeurons * nCompsPseudoBasis (i.e., nNeurons x nHistoryFeatures)
elseif nExtraArgsIn == 3
    C = varargin{1};
    d = varargin{2};
    W = varargin{3};
end


%% Create some important constants
nFrames = size(latentTrajectory,1);
nNeurons = size(firing,2);
[nHistoryFrames, nCompsPseudoBasis] = size(obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf);

nHistoryFeatures = nNeurons * nCompsPseudoBasis;  %The number of features available at each time point.

%% Sparse W mode switching:
if ~all(size(W) == [nNeurons , nNeurons * nCompsPseudoBasis])
    if all(size(W) == [nNeurons , nCompsPseudoBasis]) && strcmp(obj.modelTypes{modelNumber},'glmBrodyLaplaceSelfOnly')  %If we're working with a self-weights only model and it's using the sparse representation:

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


%% Calculate weighted basis representation of the spiking history

%In the glm fitters, we do the following:
% %Compute the feature-transformed spiking history.
% %obj.modelGlobalParameters(modelNumber).pseudoBasisSelf has
% %been constructed such that it is an nFeatures x nSamples
% %matrix, with the most recent components in the bottom rows and
% %the oldest at the top. Further, the basis functions are
% %ordered such that the most recent (i.e., that which puts the
% %most weight on the most recent time-steps) is in the left-most
% %column and the oldest (weight farther back in time) is on the
% %right. pseudoBasis' thus indexes from most recent features to
% %oldest features (by row) and by frame (column).
% transformedSpikingHistory = obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf' * spikeBinHistory;
%
% %Write these values out.
% allRelevantInputs(previousFrames + relativeFrameIdx,:) = reshape(transformedSpikingHistory,1,[]);
% %Above reordering reads down the columns of transformedSpikingHistory, and when it hits the end of a column, moves to the next column, continuing to write horizontally.
%
%
%The reason we want a vector representation of the inputs of the individual
%cells is because this will enable us to use a weight matrix which is
%nNeurons x TOTALNUMBEROFHISTORYFEATURES and do this in one line, rather
%than needing individual weight matrices for each cell.  Then we can
%present the history TOTALNUMBEROFHISTORYFEATURES x nFramesOfInterest
%matrix and do everything in one line.
%Note that we expect that the "diagonal" blocks of this weight matrix are
%going to be similar, since they represent the true self-weights of each
%neuron.

firstEndFrame   = startFrameNumber - 1;
firstStartFrame = firstEndFrame - nHistoryFrames + 1;
featureValues = zeros(nHistoryFeatures, nFrames); %TOTALNUMBEROFHISTORYFEATURES x nFramesOfInterest
for frameIdx = 1:nFrames
    chunkOfFiringWeNeed = firing(firstStartFrame + (frameIdx - 1) : firstEndFrame + (frameIdx - 1),:);  %nHistoryFrames x nNeurons
    %obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf
    transformedIntoNewBasis = chunkOfFiringWeNeed' * obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf; %Now nNeurons (source) x nFeatures
    reshapedInNewBasis = reshape(transformedIntoNewBasis',[],1);  
    %Above: This row translates to reading across the rows of the above matrix (as in normal English text) and filling in a column vector with the entries in the order encountered.
    %The result is a column vector with nNeurons (number of source neurons) blocks of nFeatures.
    %Store this in the appropriate column, where featureValues will
    %eventually be TOTALNUMBEROFHISTORYFEATURES x nFramesOfInterest
    featureValues(:,frameIdx) = reshapedInNewBasis;
end


weightContribution = W * featureValues;  
%Above:This should be a nNeurons (target neurons) x nFrames matrix,
%capturing the contribution of the self history term to log(lambda) (the
%target rate column vector) at each time step.


%% Calculate the firing rates using these values
firingRates = exp(C * latentTrajectory' + d(:,ones(1,size(latentTrajectory,1))) + weightContribution)';  
%Above: nNeurons x nFrames: the latent lambda values for each neuron at each time step under the current model

if nExtraArgsOut == 1
    varargout{1} = featureValues;
elseif nExtraArgsOut == 2
        varargout{1} = C;  %These two components are the ones we'd need to calculate the inferred rate without doing all of this trial parsing again.
        varargout{2} = (d(:,ones(1,size(latentTrajectory,1))) + weightContribution);  
end
