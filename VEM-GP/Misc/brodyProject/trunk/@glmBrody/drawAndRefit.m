function [newParams, oldParams] = drawAndRefit(obj,modelNumber,sessionIdx,varargin)
%DRAWANDREFIT.m draws from a model, creates a set of simulated neurons, adds these to the real neurons, and fits the model again.
%drawAndRefit.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%10 February, 2014
%
%Method of the glmBrody class which draws from a speficied model, takes the
%output simulatedSpikes structure and concatenates it back onto the
%original trial data, with appropriate 'simulated####' designations for the
%neurons.

%Handle the input arguments appropriately:
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
    error('Improper specification of inputs!\n\n')
end




%% Drawing using drawSimulatedNeurons

simulatedSpikes = drawSimulatedNeurons(obj,modelNumber,sessionIdx,fittingParameters,sourceNeurons);


%% Data handling to put these back into the original structure

nSimulatedNeurons = numel(sourceNeurons);
nOriginalNeurons = numel(obj.cellsIDCell{sessionIdx});
obj.cellsIDCell{sessionIdx} = horzcat(obj.cellsIDCell{sessionIdx},cell(1,nSimulatedNeurons));

%Create new cell designators and add the new data to the original
%structure.
for simNeuronIdx = 1:nSimulatedNeurons
    destIdx = nOriginalNeurons + simNeuronIdx;
    obj.cellsIDCell{sessionIdx}{destIdx} = strcat('sim',obj.cellsIDCell{sessionIdx}{sourceNeurons(simNeuronIdx)});
end

for trialIdx = 1:min(fittingParameters.maxTrialNumber, numel(obj.trialData{sessionIdx}))
    %Preallocate the storage cell array
    obj.trialData{sessionIdx}(trialIdx).spikes = horzcat(obj.trialData{sessionIdx}(trialIdx).spikes,cell(nSimulatedNeurons,1));
    for simNeuronIdx = 1:nSimulatedNeurons
        %Move the data into the trialData structure
        obj.trialData{sessionIdx}(trialIdx).spikes{nOriginalNeurons + simNeuronIdx} = simulatedSpikes{trialIdx}{simNeuronIdx, 1};
    end
end




%% Save the old parameters for comparison
oldParams.modelGlobalParameters = obj.modelGlobalParameters{modelNumber};
oldParams.modelSessionParameters = obj.modelSessionParameters{sessionIdx,modelNumber};
oldParams.modelTrialParameters = obj.modelTrialParameters{sessionIdx,modelNumber};


%% Now refit the model
%Note: this currently uses all sessions, not just the one for which we've
%simulated data.  Change?
obj.fitModel(modelNumber,fittingParameters)

newParams.modelGlobalParameters = obj.modelGlobalParameters{modelNumber};
newParams.modelSessionParameters = obj.modelSessionParameters{sessionIdx,modelNumber};
newParams.modelTrialParameters = obj.modelTrialParameters{sessionIdx,modelNumber};


%% Plot the new parameters for testing and examination
obj.plotSessionWeights(modelNumber,sessionIdx)
end