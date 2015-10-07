function [magDeltaTheta, varargout] = updateModelParams(obj,modelNumber,proposedParamsStruct)
%UPDATEMODELPARAMS updates the parameter structure of a model
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%27 August, 2014
%
%Method of glmBrody class which updates the model parameters of a model
%contained within a glmBrody class object.
%Major code stolen from tightenLatentsForObjective.

%Output handling:
nExtraArgsOut = max(nargout-1,0);



%First, copy out all of the parameters.
backupParamsStruct.modelGlobalParameters = obj.modelGlobalParameters(modelNumber);
backupParamsStruct.modelSessionParameters = obj.modelSessionParameters(:,modelNumber);
if size(obj.modelTrialParameters,1) > 0
    backupParamsStruct.modelTrialParameters = obj.modelTrialParameters(:,modelNumber);
end
backupParamsStruct.modelAltDataRepresentations = obj.modelAltDataRepresentations(:,modelNumber);


try
    magSumSquaresDeltaTheta = 0;
    
    
    %First, set modelNumberInProposedParamsStruct
    possiblePropModelNums = nan(1,3);
    if isfield(proposedParamsStruct,'modelGlobalParameters')
        possiblePropModelNums(1) = size(proposedParamsStruct.modelGlobalParameters,2);
    end
    
    if isfield(proposedParamsStruct,'modelSessionParameters')
        possiblePropModelNums(2) = size(proposedParamsStruct.modelSessionParameters,2);
    end
    
    if isfield(proposedParamsStruct,'modelTrialParameters')
        possiblePropModelNums(3) = size(proposedParamsStruct.modelTrialParameters,2);
    end
    
    uniqueSizes = unique(possiblePropModelNums(~isnan(possiblePropModelNums)));
    assert(numel(uniqueSizes) == 1,'Non-matching sizes for proposed parameter structure fields!')
    
    if numel(obj.modelTypes) == uniqueSizes
        modelNumberInProposedParamsStruct = modelNumber;
    elseif uniqueSizes == 1
        modelNumberInProposedParamsStruct = 1;
    else
        error('Error: proposedParameters structure implies a number of models ~=1, and also ~= number of models present!')
    end
    
    
    %Now, set the values to the propsal: this requires walking over all of
    %the parameter structures which are present.
    if isfield(proposedParamsStruct,'modelGlobalParameters')
        existingGlobalParamNames = fieldnames(obj.modelGlobalParameters{modelNumber});
        proposedGlobalParamNames = fieldnames(proposedParamsStruct.modelGlobalParameters{modelNumberInProposedParamsStruct});
        for fieldIdx = 1:numel(proposedGlobalParamNames)
            assert(any(strcmp(proposedGlobalParamNames{fieldIdx},existingGlobalParamNames)),['Field modelGlobalParameters.',proposedGlobalParamNames{fieldIdx},'is not a valid field name\n for the proposed parameters structure!'])
        end
        %We now know that all of the proposed fields are in fact
        %appropriate fields of the model.
        for fieldIdx = 1:numel(proposedGlobalParamNames)
            eval(['magSumSquaresDeltaTheta = magSumSquaresDeltaTheta + sum(sum( (- obj.modelGlobalParameters{modelNumber}.',proposedGlobalParamNames{fieldIdx},' + proposedParamsStruct.modelGlobalParameters{modelNumberInProposedParamsStruct}.',proposedGlobalParamNames{fieldIdx},').^2));'])
            eval(['obj.modelGlobalParameters{modelNumber}.',proposedGlobalParamNames{fieldIdx},' = proposedParamsStruct.modelGlobalParameters{modelNumberInProposedParamsStruct}.',proposedGlobalParamNames{fieldIdx},';']);
        end
    end
    
    if isfield(proposedParamsStruct,'modelSessionParameters')
        for sessionIdx = 1:numel(obj.sessionIDCell)
            existingSessionParamNames = fieldnames(obj.modelSessionParameters{sessionIdx,modelNumber});
            proposedSessionParamNames = fieldnames(proposedParamsStruct.modelSessionParameters{sessionIdx,modelNumberInProposedParamsStruct});
            for fieldIdx = 1:numel(proposedSessionParamNames)
                assert(any(strcmp(proposedSessionParamNames{fieldIdx},existingSessionParamNames)),['Field modelSessionParameters.',proposedSessionParamNames{fieldIdx},'is not a valid field name\n for the proposed parameters structure!'])
            end
            %We now know that all of the proposed fields are in fact
            %appropriate fields of the model.
            for fieldIdx = 1:numel(proposedSessionParamNames)
                eval(['magSumSquaresDeltaTheta = magSumSquaresDeltaTheta + sum(sum( (- obj.modelSessionParameters{sessionIdx,modelNumber}.',proposedSessionParamNames{fieldIdx},' + proposedParamsStruct.modelSessionParameters{sessionIdx,modelNumberInProposedParamsStruct}.',proposedSessionParamNames{fieldIdx},').^2));'])
                eval(['obj.modelSessionParameters{sessionIdx,modelNumber}.',proposedSessionParamNames{fieldIdx},' = proposedParamsStruct.modelSessionParameters{sessionIdx,modelNumberInProposedParamsStruct}.',proposedSessionParamNames{fieldIdx},';']);
            end
        end
    end
    
    if isfield(proposedParamsStruct,'modelTrialParameters')
        for sessionIdx = 1:numel(obj.sessionIDCell)
            for trialIdx = 1:numel(obj.trialData{sessionIdx})
                existingTrialParamNames = fieldnames(obj.modelTrialParameters{sessionIdx,modelNumber}{trialIdx});
                proposedTrialParamNames = fieldnames(proposedParamsStruct.modelTrialParameters{sessionIdx,modelNumberInProposedParamsStruct}{trialIdx});
                for fieldIdx = 1:numel(proposedTrialParamNames)
                    assert(any(strcmp(proposedTrialParamNames{fieldIdx},existingTrialParamNames)),['Field modelTrialParameters.',proposedTrialParamNames{fieldIdx},'is not a valid field name\n for the proposed parameters structure!'])
                end
                %We now know that all of the proposed fields are in fact
                %appropriate fields of the model.
                for fieldIdx = 1:numel(proposedTrialParamNames)
                    eval(['magSumSquaresDeltaTheta = magSumSquaresDeltaTheta + sum(sum( (- obj.modelTrialParameters{sessionIdx,modelNumber}{trialIdx}.',proposedTrialParamNames{fieldIdx},' + proposedParamsStruct.modelTrialParameters{sessionIdx,modelNumberInProposedParamsStruct}{trialIdx}.',proposedTrialParamNames{fieldIdx},').^2));'])
                    eval(['obj.modelTrialParameters{sessionIdx,modelNumber}{trialIdx}.',proposedTrialParamNames{fieldIdx},' = proposedParamsStruct.modelTrialParameters{sessionIdx,modelNumberInProposedParamsStruct}{trialIdx}.',proposedTrialParamNames{fieldIdx},';']);
                end
            end
        end
    end

    magDeltaTheta = sqrt(magSumSquaresDeltaTheta);
    
    if nExtraArgsOut == 1
        varargout{1} = backupParamsStruct;
    end
     
   
    
catch errorStructOut
    %Restore the parameters we changed:
    obj.modelGlobalParameters(modelNumber)          = backupParamsStruct.modelGlobalParameters;
    obj.modelSessionParameters(:,modelNumber)       = backupParamsStruct.modelSessionParameters;
    if size(obj.modelTrialParameters,1) > 0
        obj.modelTrialParameters(:,modelNumber)     = backupParamsStruct.modelTrialParameters;
    end
    obj.modelAltDataRepresentations(:,modelNumber)  = backupParamsStruct.modelAltDataRepresentations;
    %Now, throw the error:
    throw(errorStructOut)
end


