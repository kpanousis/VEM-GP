function createGeneratingWeights(obj,synthesizingModelNum,sessionIdx,synthesisStruct)
%CREATEGENERATINGWEIGHTS Creates glm weights 
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%1 August, 2014
%
%Need to set C, d, and W.

%% Input checking
assert(strcmp(obj.modelTypes{synthesizingModelNum},'glmBrodyLaplaceSelfOnly'),'Unsupported synthesizing model type!')
assert(isfield(synthesisStruct,'selfWeightTemplate') && size(synthesisStruct.selfWeightTemplate,2) == 1,'The field synthesisStruct.selfWeightTemplate must be present and a column vector.')
assert(numel(obj.modelGlobalParameters{1,synthesizingModelNum}.glmHypersPriorMean) == 3 && numel(obj.modelGlobalParameters{1,synthesizingModelNum}.glmHypersPriorVariance) == 3,'The glm weights hyperprior is misspecified.')
%% Set constants

nNeurons = numel(obj.cellsIDCell{sessionIdx,1});
pseudoBasisSelf = obj.modelGlobalParameters{1,synthesizingModelNum}.pseudoBasisSelf;
nCompsPseudoBasis = size(pseudoBasisSelf,2);  %The number of columns is the number of pseudobasis components.

%% Set d
%We're going to do this slightly wrong, but it should be OK
for neuronIdx = 1:nNeurons
    obj.modelSessionParameters{sessionIdx,synthesizingModelNum}.d(neuronIdx,1) = log(synthesisStruct.basalFiringRate / obj.modelGlobalParameters{1,synthesizingModelNum}.samplingRate);
end

%% Set C
%Draw these from NARROWED prior
for neuronIdx = 1:nNeurons
    obj.modelSessionParameters{sessionIdx,synthesizingModelNum}.C(neuronIdx,1) = obj.modelGlobalParameters{1,synthesizingModelNum}.glmHypersPriorMean(1) + 1/4 * sqrt(obj.modelGlobalParameters{1,synthesizingModelNum}.glmHypersPriorVariance(1)) * randn ;
end

%% Set W
if isfield(synthesisStruct,'wRegularizationWeight')
    leastSquaresSolnToWeightCurve = inv(pseudoBasisSelf' * pseudoBasisSelf + synthesisStruct.wRegularizationWeight * eye(size(pseudoBasisSelf,2)) ) * pseudoBasisSelf' * synthesisStruct.selfWeightTemplate;
else
    leastSquaresSolnToWeightCurve = inv(pseudoBasisSelf' * pseudoBasisSelf ) * pseudoBasisSelf' * synthesisStruct.selfWeightTemplate;
end
%Above: synthesisStruct.selfWeightTemplate should be lengthHistory x 1:
%result should be nCompsPseudoBasis x 1 

%Draw randomized weights around those of the least squares solution, with
%relatively small variance.
%obj.modelSessionParameters{sessionIdx,synthesizingModelNum}.W = ...
%    0.01 * obj.modelGlobalParameters{1,synthesizingModelNum}.glmHypersPriorVariance(3)^(1/2) * randn(nNeurons,nCompsPseudoBasis)...
%    + repmat(leastSquaresSolnToWeightCurve',nNeurons,1);

obj.modelSessionParameters{sessionIdx,synthesizingModelNum}.W = repmat(leastSquaresSolnToWeightCurve',nNeurons,1);