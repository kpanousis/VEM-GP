function [overall,inStimPeriod] = measureEmpiricalFiringRates(obj, sessionIdx)
%MEASUREEMPIRICALFIRINGRATES calculate firing rate in trialData
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%6 March 2015
%
%Method of glmBrody class which calculates the empirical firing rate of
%cells in a session.
%Inputs: 
%   sessionIdx: session to examine
%Outputs: 
%   overall: row vector of firing rates for all time in the session.
%   inStimPeriod: row vector containing firing rate only during stimulus period. 


overallSpikes = zeros(size(obj.cellsIDCell{sessionIdx,1}));
overallTime = 0;

inStimPeriodSpikes = zeros(size(obj.cellsIDCell{sessionIdx,1}));
inStimPeriodTime = 0;


for trialIdx = 1:numel(obj.trialData{sessionIdx,1})
    %Calculate the overall time
    if trialIdx < numel(obj.trialData{sessionIdx,1})
    allTimeForThisTrial = obj.trialData{sessionIdx,1}(trialIdx + 1).state_0_exits -  obj.trialData{sessionIdx,1}(trialIdx).state_0_exits;
    else %We're on the last trial; we're going to be a little off, but this should be OK.
        allTimeForThisTrial = max([max(obj.trialData{sessionIdx}(trialIdx).left_bups), max(obj.trialData{sessionIdx}(trialIdx).right_bups), obj.trialData{sessionIdx}(trialIdx).cpoke_end]);
        for cellIdx = 1:numel(obj.cellsIDCell{sessionIdx,1})
            if numel(obj.trialData{sessionIdx,1}(trialIdx).spikes{cellIdx}) > 0
                allTimeForThisTrial = max(allTimeForThisTrial,obj.trialData{sessionIdx,1}(trialIdx).spikes{cellIdx}(end));
            end
            %Note: since the above should essentially always be the
            %determining factor on the end time, I'm assuming that the
            %spikes fire at some rate so we have a hope of making a
            %reasonable (under-) estimate of the required time.
        end
    end
    
    %Calculate the inStimPeriod time
    inStimPeriodStart = obj.trialData{sessionIdx,1}(trialIdx).stim_start;
    inStimPeriodEnd = obj.trialData{sessionIdx,1}(trialIdx).cpoke_end;
    
    %Add these both to the appropriate running totals.
    overallTime = overallTime + allTimeForThisTrial;
    inStimPeriodTime = inStimPeriodTime + (inStimPeriodEnd - inStimPeriodStart);
    
    for cellIdx = 1:numel(obj.cellsIDCell{sessionIdx,1})
    %Count how many spikes occurred during the stimulus period and add to
    %the appropriate entry
    spikesInThisTrial = obj.trialData{sessionIdx,1}(trialIdx).spikes{cellIdx};
    inStimPeriodSpikes(cellIdx) = inStimPeriodSpikes(cellIdx) + sum(spikesInThisTrial>=inStimPeriodStart & spikesInThisTrial < inStimPeriodEnd);
    
    %Add the number of spikes which occurred in total to the appropriate
    %entry.
    overallSpikes(cellIdx) = overallSpikes(cellIdx) + numel(spikesInThisTrial);
    end
end


overall = overallSpikes / overallTime;
inStimPeriod = inStimPeriodSpikes / inStimPeriodTime;