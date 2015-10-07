function [spikingHistoryInt,startFrameNumber,clickTimes,clickSigns] = trialParser(obj,sessionIdx,trialIdx,bleedTimes,frameLength)
%TRIALPARSER returns GLM-centric information on a designated trial
%trialParser.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%11 March, 2014
%Method which looks at a particular session and trial and returns the
%appropriate representation of the trial's events for the glmBrody model.
%Inputs: obj: glmBrody class object
%        sessionIdx, trialIdx: indices into obj.trialData: we are
%           interested in obj.trialData{sessionIdx}{trialIdx}
%        bleedTimes: 
%           contains [backwardBleedTime, forwardBleedTime]: 
%               backwardBleedTime: An amount of
%               time in seconds, typically on the order of 0.1 to 0.5.
%               Typically a value from
%               obj.modelGlobalParameters{modelIdx}.relevantSelfHistory or
%               .bleedTime.  Since we're doing this without reference to a
%               particular model, we allow this parameter to be free here.
%               forwardBleedTime: Again, an amount of time in seconds,
%               which corresponds to the interval after the last bup which
%               should be kept.
%        frameLength: The length of a single frame (i.e., time bin for
%           spiking) in seconds.  Typically 0.001 or 0.01 seconds.
%Outputs:
%       spikingHistory: NON-BINARY spiking representation for all of the
%           neurons recorded: TxnNeurons.  
%       startFrameNumber: integer designator of which row (frame) in spikingHistory
%           is the one which began with the start signal (i.e., obj.trialData{sessionIdx}(trialIdx).stim_start).
%       clickTimes: nClicks x 1: times (in seconds) of the clicks with respect to this
%           moment.
%       clickSigns: nClicks x 1: (+/-1) signs of these clicks.


%% Input processing
%assert(exist(obj.trialData{sessionIdx}(trialIdx),'var'),'Data requested does not exist.')
%Above: apparently you can't make this assertion about a structure.
%fprintf('Hello world!  I''m the new version of trialParser!\n')

if numel(bleedTimes) == 1
    backwardBleedTime = bleedTimes;
    forwardBleedTime = 0;
elseif numel(bleedTimes) == 2
    backwardBleedTime = bleedTimes(1);
    forwardBleedTime = bleedTimes(2);
else
    error('Wrong number of elements of bleedTimes: must be either 1 or 2 elements.')
end

%fprintf('Using bleedTimeCenterFixation = forwardBleedTime !\n')
bleedTimeCenterFixation = forwardBleedTime;


assert(backwardBleedTime >= 0, 'Negative bleed time given!  Should be positive!')
assert(forwardBleedTime >= 0, 'Negative bleed time given!  Should be positive!')
if frameLength >= 1; fprintf('The frame length should be a value which is a fraction of a second; check this!\n\n'); end
if frameLength < 5e-4; fprintf('The chosen frame length is infeasibly short; check this!\n\n'); end
if backwardBleedTime / frameLength ~= ceil(backwardBleedTime / frameLength); fprintf('Ideally, the backward bleed time should be an integer multiple of the frame length.\n\n'); end

%% Some basic checks on the data from this trial
%The stimulus start is greater than the bleed time
assert(obj.trialData{sessionIdx}(trialIdx).stim_start >= backwardBleedTime, 'The start time of the stimulus is less than the backward bleed time from the end of the previous trial; we will have to stitch together data from previous trials.')
%NOTE: there is no corresponding checking for the forward bleed times
%(which are typically set to be zero)
assert(sum(obj.trialData{sessionIdx}(trialIdx).left_bups == obj.trialData{sessionIdx}(trialIdx).stim_start) == 1, 'The number of left bups at the start time is not one.')
assert(sum(obj.trialData{sessionIdx}(trialIdx).right_bups == obj.trialData{sessionIdx}(trialIdx).stim_start) == 1, 'The number of right bups at the start time is not one.')

%% Create some constants
nNeurons =  numel(obj.trialData{sessionIdx}(trialIdx).spikes);

%% Reset the time references for the l/r bups, double-click as the start
%Grab the other bups
leftBupsKeep = obj.trialData{sessionIdx}(trialIdx).left_bups(obj.trialData{sessionIdx}(trialIdx).left_bups ~=  obj.trialData{sessionIdx}(trialIdx).stim_start);
rightBupsKeep = obj.trialData{sessionIdx}(trialIdx).right_bups(obj.trialData{sessionIdx}(trialIdx).right_bups ~=  obj.trialData{sessionIdx}(trialIdx).stim_start);

%Convert everything to the new time reference, where zero is obj.trialData{sessionIdx}(trialIdx).stim_start
leftBupsKeep  = sort(leftBupsKeep  - obj.trialData{sessionIdx}(trialIdx).stim_start);  %Just in case they're not sorted; they should be, but we're going to use this property below
rightBupsKeep = sort(rightBupsKeep - obj.trialData{sessionIdx}(trialIdx).stim_start);
spikeTimesKeep = obj.trialData{sessionIdx}(trialIdx).spikes;
for cellIdx = 1:nNeurons
    spikeTimesKeep{cellIdx} = spikeTimesKeep{cellIdx} - obj.trialData{sessionIdx}(trialIdx).stim_start;
end

%% Prepare the click times and click signs
%This is the easier bit
clickTimes = zeros(numel(leftBupsKeep) + numel(rightBupsKeep),1);
clickSigns = zeros(numel(leftBupsKeep) + numel(rightBupsKeep),1);
leftBupPosition = 1;
rightBupPosition = 1;
for numClick = 1:numel(clickTimes)
    if leftBupPosition <= numel(leftBupsKeep) && (rightBupPosition > numel(rightBupsKeep) || leftBupsKeep(leftBupPosition) <= rightBupsKeep(rightBupPosition))
        clickTimes(numClick) = leftBupsKeep(leftBupPosition);
        clickSigns(numClick) = -1;
        leftBupPosition = leftBupPosition + 1;
    else
        clickTimes(numClick) = rightBupsKeep(rightBupPosition);
        clickSigns(numClick) = +1;
        rightBupPosition = rightBupPosition + 1;
    end
end


%% Determine the spiking history (binary matrix) and the start frame number
%First, setsome important reference frame numbers
backwardBleedFrames = ceil(backwardBleedTime / frameLength);
%timeFromStartToEnd = obj.trialData{sessionIdx}(trialIdx).spoke_in - obj.trialData{sessionIdx}(trialIdx).stim_start;
timeFromStartToEnd = min(obj.trialData{sessionIdx}(trialIdx).cpoke_end + bleedTimeCenterFixation, obj.trialData{sessionIdx}(trialIdx).cpoke_out) - obj.trialData{sessionIdx}(trialIdx).stim_start;
numFramesWithStimStartUntilIncludeTrialEnd = max(ceil(timeFromStartToEnd / frameLength), ceil((clickTimes(end) + forwardBleedTime)/frameLength));  
%Above: (9/4/14) Note that if bleedTimeCenterFixation = forwardBleedTime and
%cpoke_end >= clickTimes(end) (which it should be), then timeFromStartToEnd
%>= clickTimes(end) + forwardBleedTime, so the second argument to the max
%is never larger.
startFrameNumber = backwardBleedFrames + 1;

%Construct spikingHistory: NOTE: THIS IS A BINARY VECTOR, NOT A COUNT
%VECTOR
%spikingHistoryBinary   = false(backwardBleedFrames + numFramesWithStimStartUntilIncludeTrialEnd,nNeurons);
spikingHistoryInt   = zeros(backwardBleedFrames + numFramesWithStimStartUntilIncludeTrialEnd,nNeurons);  %Perhaps this could be made a sparse matrix?
redundantSpikeCount = 0;
for cellIdx = 1:nNeurons
    for spikeIdx = 1:numel(obj.trialData{sessionIdx}(trialIdx).spikes{cellIdx})
        timePostStart = obj.trialData{sessionIdx}(trialIdx).spikes{cellIdx}(spikeIdx) - obj.trialData{sessionIdx}(trialIdx).stim_start;
        if timePostStart >= -backwardBleedTime  && timePostStart < timeFromStartToEnd
            %Above: Slight fuzziness here: OK if backwardBleedTime represents an integer number of frames
            %Determine which is the correct time index within the matrix
            appropriateTimeIndex = floor(timePostStart / frameLength + 1) + backwardBleedFrames;
            %23 June 2015: Altered this, such that it puts spikes on the
            %boundary between frames into the later frame, rather than the
            %earlier.  This is only really an issue with the synthetic
            %data.
            
            %Now, set spikingHistory(...) ++
            %if ~spikingHistoryBinary(appropriateTimeIndex,cellIdx)
            %    spikingHistoryBinary(appropriateTimeIndex,cellIdx) = true;
            %else
            %    redundantSpikeCount = redundantSpikeCount + 1;
            %end
            
            spikingHistoryInt(appropriateTimeIndex,cellIdx) = spikingHistoryInt(appropriateTimeIndex,cellIdx) + 1;

            %For sanity checking:
           % fprintf(['Spike ',num2str(spikeIdx),' occurred at ',num2str(timePostStart),' in frame number ',num2str(appropriateTimeIndex), 'cell ', num2str(cellIdx), '.\n\n']);
        end
    end
end

%For letting the user know if the binary spiking representation is
%problematic
if redundantSpikeCount > 0
    fprintf(['There were ',num2str(redundantSpikeCount),' spikes which fell in the same bin as a spike already\n',...
        'processed from the same cell! \n',...
        'However, because the spiking representation is now non-binary, this is not a problem.\n\n'])
end


