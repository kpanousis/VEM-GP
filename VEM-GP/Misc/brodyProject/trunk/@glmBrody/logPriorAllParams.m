function [logPriorOverThetaGLMandGP, varargout] = logPriorAllParams(obj,sessionIndices,modelNumber)
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%11 May 2015
%Method which computes the prior, gradient (vector), and hessian (matrix) with respect to
%the hyperparameters.

nExtraArgsOut = max(0,nargout-1);
assert(nExtraArgsOut <=2, 'Too many outputs of method logPriorAllParams: should have at most three outputs (f, gradf,hessf).')
assert(isempty(setdiff(sessionIndices, unique(sessionIndices))), 'Non-unique sessionsIndices!')
assert(any(strcmp(obj.modelTypes{modelNumber},{'glmBrodyLaplace','glmBrodyLaplaceSelfOnly'})),...
    'The chosen model is not of type glmBrodyLaplace or glmBrodyLaplaceSelfOnly, and thus logPriorAllParams cannot be used.\n\n')


nSessionsIncludingUnused = numel(obj.sessionIDCell);

if nExtraArgsOut >= 1
    gradPriorThetas = cell(nSessionsIncludingUnused,1);
    if nExtraArgsOut == 2
        hessPriorThetas = cell(nSessionsIncludingUnused,1);
    end
end
splitsCell = cell(nSessionsIncludingUnused,1);


logPriorOverThetaGLMandGP = 0;

%Get the GLM contributions
for sessionIdx = sessionIndices
    
    %Put together the current vector of
    %hyperparameters in the correct indexing
    
    C = obj.modelSessionParameters{sessionIdx,modelNumber}.C;  %nNeurons x 1
    d = obj.modelSessionParameters{sessionIdx,modelNumber}.d;  %nNeurons x 1
    W = obj.modelSessionParameters{sessionIdx,modelNumber}.W;  %nNeurons x (nNeurons * nCompsPseudoBasis);
    nNeurons = numel(d);
    thetaVec = [reshape(C,[],1);reshape(d,[],1);reshape(W,[],1);obj.modelGlobalParameters{1,modelNumber}.gpHypers];
    splitsCell{sessionIdx} = [numel(C),nNeurons,numel(W),numel(obj.modelGlobalParameters{1,modelNumber}.gpHypers)]';

    
    
    meanC   = obj.modelGlobalParameters{modelNumber}.glmHypersPriorMean(1);  %Sets the mean for C_i, d_i, W_{i,j}: scalar valued for each
    meanD   = obj.modelGlobalParameters{modelNumber}.glmHypersPriorMean(2);
    meanW   = obj.modelGlobalParameters{modelNumber}.glmHypersPriorMean(3);
    varC    = obj.modelGlobalParameters{modelNumber}.glmHypersPriorVariance(1);
    varD    = obj.modelGlobalParameters{modelNumber}.glmHypersPriorVariance(2);
    varW    = obj.modelGlobalParameters{modelNumber}.glmHypersPriorVariance(3);
    meanThetaGLM = [ones(splitsCell{sessionIdx}(1),1) * meanC ; ones(splitsCell{sessionIdx}(2),1) * meanD ;ones(splitsCell{sessionIdx}(3),1) * meanW];
    varThetaGLM = diag([ones(splitsCell{sessionIdx}(1),1) * varC ; ones(splitsCell{sessionIdx}(2),1) * varD ;ones(splitsCell{sessionIdx}(3),1) * varW]);
    
    if sessionIdx == sessionIndices(end)
        
        meanThetaGP = obj.modelGlobalParameters{1,modelNumber}.gpHypersPriorMean;
        varThetaGP  = obj.modelGlobalParameters{1,modelNumber}.gpHypersPriorVariance;  %Stored as a matrix
        %Get the GP contributions
        if nExtraArgsOut == 0
            logPriorThetaGP                                         = logPThetaGPAndDerivs(thetaVec,splitsCell{sessionIdx},meanThetaGP,varThetaGP);
        elseif nExtraArgsOut == 1
            [logPriorThetaGP, gradPriorThetaGP]                     = logPThetaGPAndDerivs(thetaVec,splitsCell{sessionIdx},meanThetaGP,varThetaGP);
            %add to the appropriate entry in gradPriorThetas{sessionIndices(1)}
        elseif nExtraArgsOut == 2
            [logPriorThetaGP, gradPriorThetaGP, hessPriorThetaGP]   = logPThetaGPAndDerivs(thetaVec,splitsCell{sessionIdx},meanThetaGP,varThetaGP);
            %add to the appropriate entry in gradPriorThetas{sessionIndices(1)} and in
            %hessianPriorThetas{sessionIndices(1)}
        end
    end
    
    

    
    
    if nExtraArgsOut == 0
        [logPriorThetasTS]                                                              = logPThetaGLMAndDerivs(thetaVec,splitsCell{sessionIdx},meanThetaGLM,varThetaGLM);
    elseif nExtraArgsOut == 1
        [logPriorThetasTS, gradPriorThetas{sessionIdx}]                                 = logPThetaGLMAndDerivs(thetaVec,splitsCell{sessionIdx},meanThetaGLM,varThetaGLM);
    elseif nExtraArgsOut == 2
        [logPriorThetasTS, gradPriorThetas{sessionIdx}, hessPriorThetas{sessionIdx}]    = logPThetaGLMAndDerivs(thetaVec,splitsCell{sessionIdx},meanThetaGLM,varThetaGLM);
    end
    logPriorOverThetaGLMandGP = logPriorOverThetaGLMandGP + logPriorThetasTS;
    
    if sessionIdx == sessionIndices(end)
        %Add in the GP contributions
        logPriorOverThetaGLMandGP = logPriorOverThetaGLMandGP + logPriorThetaGP;
        clear logPriorThetaGP
        if nExtraArgsOut >=1
            gradPriorThetas{sessionIdx} = gradPriorThetas{sessionIdx} + gradPriorThetaGP;
            clear gradPriorThetaGP
            if nExtraArgsOut == 2
                hessPriorThetas{sessionIdx} = hessPriorThetas{sessionIdx} + hessPriorThetaGP;
                clear hessPriorThetaGP
            end
        end
    end
end



%% Assemble the gradient and Hessian
if nExtraArgsOut >= 1
    varargout{1}        = assembleGradOrHessian(gradPriorThetas, splitsCell);
    if nExtraArgsOut == 2
        varargout{2}    = assembleGradOrHessian(hessPriorThetas, splitsCell); 
    end
end