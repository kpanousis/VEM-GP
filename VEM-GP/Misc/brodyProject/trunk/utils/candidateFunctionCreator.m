function candidateFunction = candidateFunctionCreator(type,selfHistoryTime,outputSamplingRate)
%CANDIDATEFUNCTIONCREATOR Creates column vector candidate self weights.
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%14 July 2015
%
%Inputs: type: string, either 'positiveBump' or 'noBump'.
%   selfHistoryTime: the time for the self history weights to be in effect,
%       in seconds.
%   outputSamplingRate: the sampling or frame rate, in Hz, at which the 
%       GLM model will be operating.
%Outputs: candidateFunction: a column vector of weights from the past onto
%   the firing rate at this time.  These are very roughly normalized for
%   the new GLM firing rate, such that the GLM should generate
%   approximately the same number of spikes per real time.

%% Set some constants
referenceSamplingRate = 1e3;  %Reference sampling rate is 1 kHz.
nRefSelfHistFrames = selfHistoryTime * referenceSamplingRate;
assert(nRefSelfHistFrames == round(nRefSelfHistFrames),'The self history length must be an integer multiple of 1 milisecond.');  

nNewSelfHistFrames = selfHistoryTime * outputSamplingRate;
assert(nNewSelfHistFrames == round(nNewSelfHistFrames),'The self history length must also be an integer multiple of 1/outputSamplingRate.')

%% Check type and generate the reference candidate function 
switch type
    case 'positiveBump'
        endFirstNormal  = (nRefSelfHistFrames * 0.5 - 5);
        endSecondNormal = (nRefSelfHistFrames * 0.5 - 1);
        endThirdNormal  = (nRefSelfHistFrames * 0.5 - 3);
        
        candidateFunctionReference = 0.05 * (0.01 * normpdf(-5:0.5:endFirstNormal) + 0.01 * normpdf(-1:0.5:endSecondNormal) + 0.01 * normpdf(-3:0.5:endThirdNormal)  - 2.5 * poisspdf(0:nRefSelfHistFrames,0.5)) ;
        candidateFunctionReference = fliplr(candidateFunctionReference(1:end-1));
        candidateFunctionReference = candidateFunctionReference';
    case 'noBump'
        candidateFunctionReference = 0.05 * (                                                                               -poisspdf(0:nRefSelfHistFrames,0.5)) ; %NO BUMP UPWARD
        candidateFunctionReference = fliplr(candidateFunctionReference(1:end-1));
        candidateFunctionReference = candidateFunctionReference';
    otherwise
        error('Unrecognized candidate function type.')
end

%% Now, use the candidateFunctionReference to create the up or down-sampled version.

%Create a linear function for which the averages are the reference values:
waypoints = zeros(numel(candidateFunctionReference)+1,1);
waypoints(end) = 5/4 * candidateFunctionReference(end);
waypoints(end-1) = 3/4 * candidateFunctionReference(end);
for waypointIdx = (numel(waypoints) - 2):-1:1
    waypoints(waypointIdx) = 2 * candidateFunctionReference(waypointIdx) - waypoints(waypointIdx + 1);
end

assert(all( abs(waypoints(1:end-1) + waypoints(2:end) - 2 * candidateFunctionReference) < 1e-8),'The waypoints do not average to the reference values as desired!')

%Find the corresponding timestamps for these waypoints
waypointTimestamps = (-nRefSelfHistFrames / referenceSamplingRate : 1/referenceSamplingRate : 0 )';
assert(numel(waypointTimestamps) == numel(waypoints),'These two objects should be the same number of elements!')

newTimestamps = (-nNewSelfHistFrames / outputSamplingRate : 1/outputSamplingRate : 0)';  %This is the set of bin edges of the new bins

candidateFunction = zeros(nNewSelfHistFrames,1);

for candidateFunctionIdx = nNewSelfHistFrames:-1:1
    %Calculate the corresponding weight.
   
   %The weight is the mean of the corresponding section of the piecewise
   %linear function. There are two cases here: 1. one or more waypoints lie
   %inside of this open interval; 2. no waypoints lie inside this open
   %interval.
   
   endpoints = [newTimestamps(candidateFunctionIdx);newTimestamps(candidateFunctionIdx + 1)];  %The left entry is more negative than the right.
   waypointsInThisInterval = waypointTimestamps > endpoints(1) & waypointTimestamps < endpoints(2);  %Waypoints strictly in the open interval.
   if ~isempty(waypointsInThisInterval)
       %At least one waypoint falls in this range
       timestampsOfRelevance = [endpoints(1);waypointTimestamps(waypointsInThisInterval);endpoints(2)];
       
       [bracketingTimestampsMoreNegative, bracketingFValuesMoreNegative] = findBrackets(endpoints(1),waypointTimestamps,waypoints);
       [bracketingTimestampsMorePositive, bracketingFValuesMorePositive] = findBrackets(endpoints(2),waypointTimestamps,waypoints);
       
       fValuesAtThesePoints = [bracketingFValuesMoreNegative(1) + (bracketingFValuesMoreNegative(2) - bracketingFValuesMoreNegative(1)) * (endpoints(1) - bracketingTimestampsMoreNegative(1));...
           waypoints(waypointsInThisInterval);...
           bracketingFValuesMorePositive(1) + (bracketingFValuesMorePositive(2) - bracketingFValuesMorePositive(1)) * (endpoints(2) - bracketingTimestampsMorePositive(1))    ];
   else
       %No waypoints fall in the open interval
       timestampsOfRelevance = endpoints;
       
       [bracketingTimestampsMoreNegative, bracketingFValuesMoreNegative] = findBrackets(endpoints(1),waypointTimestamps,waypoints);
       [bracketingTimestampsMorePositive, bracketingFValuesMorePositive] = findBrackets(endpoints(2),waypointTimestamps,waypoints);
       
       
       fValuesAtThesePoints = [bracketingFValuesMoreNegative(1) + (bracketingFValuesMoreNegative(2) - bracketingFValuesMoreNegative(1)) * (endpoints(1) - bracketingTimestampsMoreNegative(1));...
           bracketingFValuesMorePositive(1) + (bracketingFValuesMorePositive(2) - bracketingFValuesMorePositive(1)) * (endpoints(2) - bracketingTimestampsMorePositive(1))    ];
   end
   
   candidateFunction(candidateFunctionIdx) = sum((diff(timestampsOfRelevance) .* (fValuesAtThesePoints(1:end-1) + fValuesAtThesePoints(2:end))/2 )) / (endpoints(2) - endpoints(1));
    
end


%% Apply an approximate rescaling to account for the new frame rate
%Note: The idea here is to have a roughly equal RATE of spike generation,
%in real time, after the rescaling; Unfortunately, since the GLM is
%non-linear, an exact rescaling is not possible.  However, we can use the
%following approximation. This is based on the idea that for lambda[t] =
%exp(C a[t] + d + Wh(y_old)), lambda approximately scales linearly with W
%when the argument is nearly zero.  When the firing rate is relatively
%small, this should be an OK approximation.
candidateFunction = (referenceSamplingRate / outputSamplingRate) * candidateFunction;  


end

function [bracketingTimestamps,bracketingFValues] = findBrackets(t,TRef,FVals)
%We want to find the index of the last TRef which is <= t, and the first
%TRef which is >= t.

assert(numel(TRef) == numel(FVals),'TRef and FVals must be the same number of elements.');
assert(numel(t) == 1,'t must be a scalar.')

% 
bracketingTimestamps = nan(2,1);
bracketingFValues = nan(2,1);
idxLeftBracket = find(TRef<= t,1,'last');
idxRightBracket = find(TRef >= t, 1, 'first');

%
bracketingFValues(1) = FVals(idxLeftBracket);
bracketingFValues(2) = FVals(idxRightBracket);
bracketingTimestamps(1) = TRef(idxLeftBracket);
bracketingTimestamps(2) = TRef(idxRightBracket);

end