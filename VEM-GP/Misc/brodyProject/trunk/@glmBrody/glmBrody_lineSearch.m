function [functionValues, lambdaValuesStore, varargout] = glmBrody_lineSearch(obj, modelNumber, proposedParameterStepStructure, fZero, gradientInformation, fittingParameters, gradNotNewtonLogical,varargin)
%GLMBRODY_LINESEARCH implements a line search in the model parameters
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%22 August, 2014
%
%Method which implements a line search in the parameters of a model.
%
%Inputs:
%   modelNumber: model index (in obj.modelTypes) of the model to be
%       searched.
%   proposedParameterStepStructure: A structure, possibly having
%       modelGlobalParameters, modelSessionParameters, and modelTrialParameters
%       fields, which describes the proposed step from some other method;
%       this should be either a gradient or Newton method.
%   fZero: The objective function value at the current parameters.
%   gradientInformation: 
%       EITHER: fPrimeZero: Value of the function gradient along the unit vector in the
%           proposed direction, starting at the original hyperparameters,
%           rescaled by the norm of the proposed step, s.t. this is the
%           derivative of f(lambda) as used here.
%       OR: a structure which can be used with proposedParameterStepStructure
%            to calculate fPrimeZero
%   
%   fittingParameters: the standard structure which controls fitting: has
%       fields alphaNewtonLS (Newton) and
%       lambdaGradLS_outwardRatio,lambdaGradLS_inwardRatio, 
%       lambdaStartGradLS (Gradient) which control elements of the line
%       search.
%       Also has fields:
%       maxNumEvalsLineSearch: Maximum number of evaluations of the model in a line
%       search; and
%       objectiveFunctionIdx: which determines which of the
%       returned function values the system uses to control the line search
%       process.
%   gradOrNewtonLogical: Logical which indicates whether the proposed step
%       is a Newton step or a gradient step: true if gradient step.
%Outputs:
%   functionValues: returns the logLikelihood and unnormalizedLogPosterior
%       values (assuming these are available) for the individual steps in
%       the optimization.
%   lambdaValuesStore: The corresponding values of Lambda.
%   varargout: if present, returns the successive structures of parameters
%       which correspond with the individual evaluations in the line
%       search, and if two are present, the above, along with the latents
%       found by the final proposed set of parameters.

%% Output handling
nExtraArgsOut = max(0,nargout - 2);

%% Input handling
nExtraArgsIn = max(0,nargin - 7);

if nExtraArgsIn == 3
    if islogical(varargin{1}) && numel(varargin{1}) == numel(obj.sessionIDCell)
        activeSessionsLogical = varargin{1};
    elseif ~islogical(varargin{1})
        activeSessionsLogical = any(repmat(1:numel(obj.sessionIDCell),numel(varargin{1}),1) == repmat(reshape(varargin{1},[],1),1,numel(obj.sessionIDCell)),1);
    end
    sessionIndices = find(activeSessionsLogical);
    inactiveSessionContributionsToLL = sum(varargin{2}(~activeSessionsLogical));
    inactiveSessionContributionsToAULP = sum(varargin{3}(~activeSessionsLogical));
elseif nExtraArgsIn ~= 0
    error('Improper number of inputs to glmBrodyLineSearch!')
end

%% Prepare for the major loop
maxNumEvals = fittingParameters.maxNumEvalsLineSearch;
parametersStructuresStore = cell(1,maxNumEvals);
lambdaValuesStore = nan(1,maxNumEvals);
%functionValues = nan(1,maxNumEvals); %Note; the number of rows used may increase 

%% Determine a value for fPrimeZero
if ~isstruct(gradientInformation)
    fPrimeZero = gradientInformation; %We've calculated the gradient wrt lambda above and are passing it down.
else
    %We need to calculate this locally.  
    
    %BELOW: REPLACE WITH
    fPrimeZeroCheck = paramsStructDotProduct(gradientInformation,proposedParameterStepStructure);
    
    %The idea is fPrimeZero = gradientInformation .*
    %proposedParamterStepStructure.  Unfortunately, since these are both
    %structures, this is a bit more complicated.
    fPrimeZero = 0;
    if isfield(proposedParameterStepStructure,'modelGlobalParameters')
        assert(numel(gradientInformation.modelGlobalParameters) == 1,'More than one model present in gradientInformation structure!')
        deltaFPrimeZero = fPrimeZeroSubStructParser(gradientInformation.modelGlobalParameters{1},proposedParameterStepStructure.modelGlobalParameters{1});
        fPrimeZero = fPrimeZero + deltaFPrimeZero;
    end
    
    clear deltaFPrimeZero;
    
    if isfield(proposedParameterStepStructure,'modelSessionParameters')
        for sessionIdx = 1:size(proposedParameterStepStructure.modelSessionParameters,1)
            deltaFPrimeZero = fPrimeZeroSubStructParser(gradientInformation.modelSessionParameters{sessionIdx,1},proposedParameterStepStructure.modelSessionParameters{sessionIdx,1});
            fPrimeZero = fPrimeZero + deltaFPrimeZero;
        end
    end
    
    clear deltaFPrimeZero;
    
    if isfield(proposedParameterStepStructure,'modelTrialParameters')
        for sessionIdx = 1:size(proposedParameterStepStructure.modelTrialParameters{sessionIdx,1},1)
            for trialIdx = 1:size(proposedParameterStepStructure.modelTrialParameters{sessionIdx,1},2)
                deltaFPrimeZero = fPrimeZeroSubStructParser(gradientInformation.modelTrialParameters{sessionIdx,1}{trialIdx},proposedParameterStepStructure.modelTrialParameters{sessionIdx,1}{trialIdx});
                fPrimeZero = fPrimeZero + deltaFPrimeZero;
            end
        end
        
    end
   
    assert(abs(fPrimeZero - fPrimeZeroCheck) < 1e-6,'These values should check equal!')
    
end

if fPrimeZero < 0
    fprintf('The derivative of the linear approximation to f is locally negative; this is a bad search direction!')
    fPrimeZero = 0;  
    %This seems like a terrible idea: it says, "search for ANY improvement in this direction, hoping that something else we've done is wrong or that we jump sufficiently far as to recover."
end
%% Set the originalParametersStruct for use below
originalParametersStruct.modelGlobalParameters  = obj.modelGlobalParameters(modelNumber);
originalParametersStruct.modelSessionParameters = obj.modelSessionParameters(:,modelNumber);
originalParametersStruct.modelTrialParameters   = obj.modelTrialParameters(:,modelNumber);
%% Begin the main loop
for evaluationNum = 1:maxNumEvals
    %Determine where we should go next
    if gradNotNewtonLogical
        error('Not yet implemented!')
        %if evaluationNum == 1
        %    outwardPhase = true;
        %    lambda = fittingParameters.lambdaStartGradLS;
        %elseif outwardPhase
        %    lambda = lambdaValuesStore(evaluationNum-1) * fittingParameters.lambdaGradLS_outwardRatio;
        %else  %We're in the inward phase: we need to know what was the last good evaluation and then what section we should be searching
        %    
        %end
    else  %We're doing Newton steps
        if evaluationNum == 1
            lambda = 1;  %Npte that this does NOT put norm limit on how far we're willing to go; if the Newton step is absurdly large, we're going to go REALLY far out before we come back.
        elseif evaluationNum == 2  %The Newton step was too far; we're going to have to try to find the maximum in the interval between lambda = 0 and 1.
            lambda = -fPrimeZero / (2 * (thisIterFxnValue - fZero - fPrimeZero));
            lambda = min(0.5,max(0.1,lambda));  %Enforce some constraints on the value of lambda, as recommended by Numerical Recipes (pp. 384-385)
        else  %We have at least two past evaluations, so use these in a quadratic fit to find a good estimate for lambda.
            fLambda1 = functionValues(fittingParameters.objectiveFunctionIdx,evaluationNum - 1);  %Most recent value for lambda
            fLambda2 = functionValues(fittingParameters.objectiveFunctionIdx,evaluationNum - 2);  %Second-most recent value for lambda
            lambda1 = lambdaValuesStore(evaluationNum-1);
            lambda2 = lambdaValuesStore(evaluationNum-2);
            invMatLambdas = 1/(lambda1-lambda2) * [lambda1^(-2), -lambda2^(-2); -lambda2 * lambda1^(-2), -lambda1 * lambda2^(-2)];
            coeffsVec = invMatLambdas * [fLambda1 - fZero - fPrimeZero * lambda1; fLambda2 - fZero - fPrimeZero * lambda2];
            lambda = (-coeffsVec(2) + sqrt(coeffsVec(2)^2 - 3 * coeffsVec(1) * fPrimeZero)) / (3 * coeffsVec(1));
            lambda = min(0.5 * lambda1,max(0.1 * lambda1, lambda));
        end
    end
    
    
    %Build the parameters structure to pass to the evaluation function:
    %This may have global, session, and trial parameters.
    assert(~exist('lambdaParametersStructure','var'), 'The lambdaParametersStructure still exists from the last iteration!')
    
    
    %New version:
    lambdaParametersStructure = paramsStructAdd(originalParametersStruct,paramsStructScalarMultiply(proposedParameterStepStructure,lambda));  %a + lambda * b
    
    %Old version:
%     %Initialize with the current values:
%     if isfield(proposedParameterStepStructure,'modelGlobalParameters')
%         assert(size(proposedParameterStepStructure.modelGlobalParameters,2) == 1,'Proposed step structure is wider than a single model (globals)!')
%         %Update using the proposed step and lambda
%         deltaThetaFieldsGlobal = fieldnames(proposedParameterStepStructure.modelGlobalParameters{1});
%         for thetaFieldIdx = 1:numel(deltaThetaFieldsGlobal)
%             eval(['lambdaParametersStructure.modelGlobalParameters{1}.',deltaThetaFieldsGlobal{thetaFieldIdx},' = obj.modelGlobalParameters{1,modelNumber}.',deltaThetaFieldsGlobal{thetaFieldIdx},' + lambda * proposedParameterStepStructure.modelGlobalParameters{1}.',deltaThetaFieldsGlobal{thetaFieldIdx},';']);
%         end
%         
%     end
%     if isfield(proposedParameterStepStructure,'modelSessionParameters')
%         assert(size(proposedParameterStepStructure.modelSessionParameters,2) == 1,'Proposed step structure is wider than a single model (sessions)!')
%         lambdaParametersStructure.modelSessionParameters = obj.modelSessionParameters(:,modelNumber);
%         for sessionIdx = 1:size(lambdaParametersStructure.modelSessionParameters,1)
%             %update using the proposed step and lambda
%             deltaThetaFieldsThisSession = fieldnames(proposedParameterStepStructure.modelSessionParameters{sessionIdx});  %Note that writing this in this fashion means that the session parameters which change may be different for different sessions.
%             for thetaFieldIdx = 1:numel(deltaThetaFieldsThisSession)
%                 %try
%                     eval(['lambdaParametersStructure.modelSessionParameters{sessionIdx}.',deltaThetaFieldsThisSession{thetaFieldIdx},' = obj.modelSessionParameters{sessionIdx,modelNumber}.',deltaThetaFieldsThisSession{thetaFieldIdx},' + lambda * proposedParameterStepStructure.modelSessionParameters{sessionIdx}.',deltaThetaFieldsThisSession{thetaFieldIdx},';']);
%                 %catch errorStruct
%                 %    keyboard
%                 %end
%             end
%         end
%     end
%     if isfield(proposedParameterStepStructure,'modelTrialParameters')
%         error('Have not debugged modelTrialParameters handling code!')
%         assert(size(proposedParameterStepStructure.modelTrialParameters,2) == 1,'Proposed step structure is wider than a single model (trials)!')
%         lambdaParametersStructure.modelTrialParameters = obj.modelTrialParameters(:,modelNumber);
%         for sessionIdx = 1:size(lambdaParametersStructure.modelTrialParameters,1)
%             assert(all(size(lambdaParametersStructure.modelTrialParameters{sessionIdx}) == size(obj.trialData{sessionIdx,1})),'The size of the modelTrialParameters{sessionIdx} cell array is not 1xnTrials!\n')
%             for trialIdx = 1:size(lambdaParametersStructure.modelTrialParameters{sessionIdx},2)
%                 %update using the proposed step and lambda
%                 deltaThetaFieldsThisTrial = fieldnames(proposedParameterStepStructure.modelTrailParameters{sessionIdx}{1,trialIdx});  %Note that writing this in this fashion means that the session parameters which change may be different for different sessions.
%                 for thetaFieldIdx = 1:numel(deltaThetaFieldsThisTrial)
%                     eval(['lambdaParametersStructure.modelTrialParameters{sessionIdx}{1,trialIdx}.',deltaThetaFieldsThisTrial{thetaFieldIdx},' = obj.modelSessionParameters{sessionIdx,modelNumber}{1,trialIdx}.',deltaThetaFieldsThisTrial{thetaFieldIdx},' + lambda * proposedParameterStepStructure.modelSessionParameters{sessionIdx}{1,trialIdx}.',deltaThetaFieldsThisTrial{thetaFieldIdx},';']);
%                 end
%             end
%         end
%     end
    
    
    
    
    %Pass this structure to a function which evaluates the objectives at the chosen values of the parameters
    
    %Noisy version:
    %functionValues(:,evaluationNum) = obj.evaluateObjectivesForProposedThetaStep(modelNumber,fittingParameters,lambdaParametersStructure);
    
    %Quiet version:
    if nExtraArgsOut <= 2
        [~, functionValues(:,evaluationNum),latentTrajectoriesFromEvalObjForProposedThetaStep] = evalc('obj.evaluateObjectivesForProposedThetaStep(modelNumber,fittingParameters,lambdaParametersStructure)');
    else
        [~, functionValues(:,evaluationNum),latentTrajectoriesFromEvalObjForProposedThetaStep,~,~,perSessionContributionsToLL,perSessionContributionstoAULP] = evalc('obj.evaluateObjectivesForProposedThetaStep(modelNumber,fittingParameters,lambdaParametersStructure,sessionIndices,inactiveSessionContributionsToLL,inactiveSessionContributionsToAULP )');
    end
    
%     if all(activeSessionsLogical)
%         thisIterFxnValue = functionValues(fittingParameters.objectiveFunctionIdx,evaluationNum);
%     else
%         functionValues(:,evaluationNum) = functionValues(:,evaluationNum) + [inactiveSessionContributionsToLL;inactiveSessionContributionsToAULP];
%         thisIterFxnValue = functionValues(fittingParameters.objectiveFunctionIdx,evaluationNum);
%     end
    thisIterFxnValue = functionValues(fittingParameters.objectiveFunctionIdx,evaluationNum);

        
    
    
    
    %Store the obtained values
    parametersStructuresStore{evaluationNum} = lambdaParametersStructure;
    lambdaValuesStore(evaluationNum) = lambda;
    
    
    %Check for shift from outward phase to inward phase
    %Check the chosen function value
    if gradNotNewtonLogical
        error('Not yet implemented!')
        %if outwardPhase && ((evaluationNum == 1 && thisIterFxnValue < fZero) || (evaluationNum > 1 && thisIterFxnValue < functionValues(fittingParameters.objectiveFunctionIdx,evaluationNum-1)))
        %    outwardPhase = false;  %We've found a downturn, so now we need to begin sectioning inward.
        %    %Anything else we need to do?
        %elseif  %
        %
        %end
    else  %We're doing Newton steps
        correctionVal = fittingParameters.alphaNewtonLS * lambda * fPrimeZero;
        if thisIterFxnValue > fZero + correctionVal
            terminationCondition = true; %The proposed value worked; accept.
        else
            terminationCondition = false;
        end
    end
    
    %Clear lambdaParametersStructure, such that it doesn't cause a weird collision next time:
    clear lambdaParametersStructure
    
    
    %Check termination condition
    if terminationCondition
        break
    end
    

    
end

%% Truncate lambdaValuesStore and parametersStructuresStore:
lambdaValuesStore = lambdaValuesStore(1,1:evaluationNum);
parametersStructuresStore = parametersStructuresStore(1,1:evaluationNum);

%% For diagnostic purposes
figure(1)
clf
hold all
scatter(lambdaValuesStore,functionValues(fittingParameters.objectiveFunctionIdx,:))
lambdasForPlotting = 0:0.01:1;
if gradNotNewtonLogical
    fLinForPlotting = fZero + fPrimeZero * lambdasForPlotting;
    plot(lambdasForPlotting,fLinForPlotting)
else
    fQuadForPlotting = fZero + fPrimeZero * (-1/2 * lambdasForPlotting.^2 + lambdasForPlotting);
    plot(lambdasForPlotting,fQuadForPlotting)
end
ylabel('f(x_0 + \lambda \Delta x)')
xlabel('\lambda')
title('f values for used \lambda')
pause(10)


%% Varargout handling
if nExtraArgsOut == 1
    varargout{1} = parametersStructuresStore;
elseif nExtraArgsOut == 2
    varargout{1} = parametersStructuresStore;
    varargout{2} = latentTrajectoriesFromEvalObjForProposedThetaStep;
elseif nExtraArgsOut == 4
    varargout{1} = parametersStructuresStore;
    varargout{2} = latentTrajectoriesFromEvalObjForProposedThetaStep;
    varargout{3} = perSessionContributionsToLL;
    varargout{4} = perSessionContributionstoAULP;
end


end

function deltaFPrimeZero = fPrimeZeroSubStructParser(substruct1,substruct2)

deltaFPrimeZero = 0;

fieldNamesSubstruct1 = fieldnames(substruct1);
fieldNamesSubstruct2 = fieldnames(substruct2);

assert(numel(fieldNamesSubstruct1) == numel(fieldNamesSubstruct2),'Unequal numbers of field names in the two structures!')
for fieldIdx = 1:numel(fieldNamesSubstruct2)
    %Verify that every field is present for both.
    assert(sum(strcmp(fieldNamesSubstruct2{fieldIdx},fieldNamesSubstruct1)) == 1,'There should be one and only one matching field!')
    
    %For each field, calculate the dot product and add this to
    %fPrimeZero.
    size1 = size(eval(['substruct1.',fieldNamesSubstruct2{fieldIdx}]));
    size2 = size(eval(['substruct2.',fieldNamesSubstruct2{fieldIdx}]));
    eval(['deltaFPrimeZero = deltaFPrimeZero + sum(sum(substruct1.',fieldNamesSubstruct2{fieldIdx},' .* substruct2.',fieldNamesSubstruct2{fieldIdx},'));'])
end


end

