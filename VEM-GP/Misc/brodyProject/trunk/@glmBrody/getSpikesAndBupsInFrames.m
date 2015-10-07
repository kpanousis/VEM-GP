function [relevantFrames, spikesInFrames, varargout] = getSpikesAndBupsInFrames(obj,modelNumber,sessionIdx,fittingParameters,bleedTime,varargin)
%GETSPIKESANDBUPSINFRAMES Returns frames to examine, & frames with spikes.
%
%Method which examines a full session (concatenating together trials) to
%come up with frame indices (in a global reference) for  contiguous periods
%in which bups are happening (relevantFrames) and the frame numbers in
%which all spikes occur. May also return bup frames for left and right,
%along with the chosen global reference timestamp.
%
%Notes: Importantly, ignores trial structure beyond the bups;
%relevantFrames is set solely with reference to the bups.

nExtraArgsIn = max(0,nargin - 5);
assert(nExtraArgsIn <=1,'Too many inputs to getSpikesAndBupsInFrames!')

if nExtraArgsIn == 1
    windowType = varargin{1};
else
    fprintf('Using click-based segmentation of relevant windows in history generation (getSpikesAndBupsInFrames).\n\n')
    windowType = 'clickBased';
end

%NEED:
switch windowType
    case 'trialBased'
        fprintf('In getSpikesAndBupsInFrames, using bleedTimeCenterFixation = obj.modelGlobalParameters{modelNumber}.maxLagConsidered !\n')
        bleedTimeCenterFixation =  obj.modelGlobalParameters{modelNumber}.maxLagConsidered;
end


nExtraArgsOut = max(nargout - 2,0);
assert(nExtraArgsOut ~= 1  && nExtraArgsOut <= 3, 'Incorrect output syntax!')


%First, compile the list of all time stamps at which events occurred.
%These times are all in a global reference.
if isfield(fittingParameters,'maxTrialNumber')
    maxTrials = min(fittingParameters.maxTrialNumber,numel(obj.trialData{sessionIdx}));
    [allSpikes, allLeftBups,allRightBups] = obj.listAllBandS_Method(sessionIdx,maxTrials);
else
    maxTrials = numel(obj.trialData{sessionIdx});
    [allSpikes, allLeftBups,allRightBups] = obj.listAllBandS_Method(sessionIdx);
end

nNeurons = numel(allSpikes);
samplingRate = obj.modelGlobalParameters{modelNumber}.samplingRate;


%For sanity, set the reference time to a value which is the nearest
%whole second before the first event (probably a spike).
minBupTime = min(min(allLeftBups),min(allRightBups));

minSpikeTime = inf;
for cellIdx = 1:nNeurons
    if ~(numel(allSpikes{cellIdx}) == 0)  %if this cell fired any spikes at all
        minSpikeTime = min(minSpikeTime,min(allSpikes{cellIdx}));
    end
end
if ~isinf(minSpikeTime)
    minRefTime = floor(min(minBupTime,minSpikeTime));  %largest whole second time value before any of the events.
else
    minRefTime = floor(minBupTime); %In the extremely unlikely case that the only (or all) cell(s) fired no spikes in the portion of the session being examined, set the reference time here to the minBupTime.
end

allLeftBups = allLeftBups - minRefTime;
allRightBups = allRightBups - minRefTime;
allBups = sort(vertcat(allLeftBups,allRightBups),'ascend');  %Used below for determining relevant windows.

for cellIdx = 1:numel(allSpikes)
    allSpikes{cellIdx} = allSpikes{cellIdx} - minRefTime;
end

switch windowType
    case 'clickBased'
        relevantWindows = [];  %This will contain in column 1 bupTimes - bleedTime and in column 2, bupTimes+bleedTime
        for bupIdx = 1:numel(allBups)
            if isempty(relevantWindows) || allBups(bupIdx) > relevantWindows(end,2) + bleedTime
                %That is, if the bup is outside of the last window:
                %Start a new window
                relevantWindows = vertcat(relevantWindows, [allBups(bupIdx) - bleedTime, allBups(bupIdx) + bleedTime]);
            else
                %If the bup is inside of the last window, extend the window
                %further.
                relevantWindows(end,2) = allBups(bupIdx) + bleedTime;
            end
        end
    case 'trialBased'
        relevantWindows = zeros(maxTrials,2);
        for trialIdx = 1:maxTrials
            %Create the relevant window
            tStart = obj.trialData{sessionIdx}(trialIdx).state_0_exits;  %As in listAllBandS_Method: lets us set the window to a global time.
            relevantWindows(trialIdx,1) = (tStart - minRefTime) + min(obj.trialData{sessionIdx}(trialIdx).left_bups(1), obj.trialData{sessionIdx}(trialIdx).right_bups(1)) - bleedTime;
            relevantWindows(trialIdx,2) = (tStart - minRefTime) + min(obj.trialData{sessionIdx}(trialIdx).cpoke_end + bleedTimeCenterFixation, obj.trialData{sessionIdx}(trialIdx).cpoke_out);
        end
        assert(all(relevantWindows(2:end,1) > relevantWindows(1:(end-1),2)),'In getSpikesAndBupsInFrames, trial-based windows overlap! Check this!')
        
end
%Now, use the bups matrices to set the appropriate time windows,




%Now, convert everything to frames:
%Take the first column to be the floor (by frames) and the second to the the ceil (by frames)

%7 Feb, 2014: setting it such that relevant windows remains un-rounded, but
%relevantFrames is rounded correctly: relevant windows isn't used
%elsewhere, so this should be OK.
%relevantWindows = [floor(relevantWindows(:,1)*samplingRate)/samplingRate, ceil(relevantWindows(:,2)*samplingRate)/samplingRate];
%relevantFrames = round(relevantWindows * samplingRate - [zeros(size(relevantWindows,1),1),ones(size(relevantWindows,1),1)]);
relevantFrames = round(...
    [floor(relevantWindows(:,1)*samplingRate)/samplingRate, ceil(relevantWindows(:,2)*samplingRate)/samplingRate]...
    * samplingRate...
    - [zeros(size(relevantWindows,1),1),ones(size(relevantWindows,1),1)]);



leftBupsInFrames = round(floor(allLeftBups * samplingRate) + 1);
rightBupsInFrames = round(floor(allRightBups * samplingRate) + 1);

%Note: these frame numbers are the numbers of the frame in which the event
%occurs, where the start of frame 1 is at global time minRefTime, and each
%frame is 1/samplingRate long.


%Put spikes in frames as well
spikesInFrames = cell(1,nNeurons);
for cellIdx = 1:numel(spikesInFrames)
    spikesInFrames{cellIdx} = round(floor(allSpikes{cellIdx} * samplingRate)+1);
end


if nExtraArgsOut == 2
    varargout{1} = leftBupsInFrames;
    varargout{2} = rightBupsInFrames;
elseif nExtraArgsOut == 3
    varargout{1} = leftBupsInFrames;
    varargout{2} = rightBupsInFrames;
    varargout{3} = minRefTime;
end


end