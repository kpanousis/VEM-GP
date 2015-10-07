function [varargout] = empiricalfiringRateCalculator(obj,sessionIdx,modelNum,figureNum)
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%21 October, 2014
%
%Method which tries a range of values of the evidence accumulator, produces
%a period of firing resulting from the accumulator value being held
%constant, and returns a resulting firing rate curve for each cell in the
%session.

%% Warning, since this is an initial implementation and not really well-refined.
fprintf('Using a version of empiricalFiringRateCalculator which does no post-processing on the drawn spiking; compare with drawSyntheticSpikeFrames.m.\n\n')


%% Input Checking
%Simple checks
assert(size(obj.trialData,1) >= sessionIdx,'No session for sessionIdx is present!')
assert(size(obj.modelTypes,2) >= modelNum,'No model for modelNum is present!')

%% Prepare the values of the accumulator:
accumulatorValuesRange = [-40,-30,-20:2:-10,-9:1:-3,-2.5:0.5:2.5,3:1:9,10:2:20,30,40];
maxTime = 10; %Seconds

dt = 1/obj.modelGlobalParameters{modelNum}.samplingRate;

nNeurons = numel(obj.cellsIDCell{sessionIdx});

%% Preallocate storage:
ratesPerCell = nan(nNeurons,numel(accumulatorValuesRange));
%% For each accumulator value, calculate the resulting rate
for aValIdx = 1:numel(accumulatorValuesRange)
    valOfA = accumulatorValuesRange(aValIdx);
    
    
    %and a transcript of the accumulator at the indexed value
    timesForAccumulator = (0:dt:maxTime)';
    accumulatorTranscript = valOfA * ones(size(timesForAccumulator));
    
    %Note: this uses initial period spiking according to d
    try
        spikeFramesCell = obj.drawSyntheticSpikeFrames(modelNum,sessionIdx,accumulatorTranscript);
        
        
        
        for cellIdx = 1:numel(spikeFramesCell)
            spikeFramesCell{cellIdx} = spikeFramesCell{cellIdx} - 1;  %to reset from frame idx to relative frame idx.
            assert(all(spikeFramesCell{cellIdx} <= maxTime / dt ) && all(spikeFramesCell{cellIdx} >= 0),'Spike times fall outside of the desired range! NOT GOOD!')
            ratesPerCell(cellIdx,aValIdx) = numel(spikeFramesCell{cellIdx})/maxTime;
        end
    catch errStructFromDraw
        fprintf('A Cell failed at rate %2.2f!\n',valOfA)
        %keyboard
    end
end

%% Create a plot
figure(figureNum)
clf
hold all
%For each cell, plot the resulting rate curve
for cellIdx = 1:size(ratesPerCell,1)
    plot(accumulatorValuesRange,ratesPerCell(cellIdx,:))
end
xlabel('Accumulator Value')
ylabel('Simulated net firing rate')

%% Output handling

if nargout > 0
    varargout{1} = ratesPerCell;
end