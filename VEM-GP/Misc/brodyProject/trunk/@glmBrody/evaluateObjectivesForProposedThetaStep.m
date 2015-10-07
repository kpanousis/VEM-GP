function [functionValues,varargout] = evaluateObjectivesForProposedThetaStep(obj,modelNumber,fittingParameters,lambdaParametersStructure,varargin)
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%26 August, 2014
%
%Evaluate the function values request

nExtraArgsOut = max(0,nargout - 1);
%assert(nExtraArgsOut <=1,'Too many outputs from evaluateObjectivesForProposedThetaStep!')

switch obj.modelTypes{modelNumber}
    case 'glmBrodyLaplaceSelfOnly'
        functionValues = zeros(2,1);
        [functionValues(1), functionValues(2),varargout{1:(nargout-1)}] = obj.tightenLatentsForObjective(modelNumber,fittingParameters,lambdaParametersStructure,varargin{:});
        %keyboard
        %function varargout = wrapper( varargin )
        %[varargout{1:nargout}] = someFunction( varargin{:} ); 
    otherwise
        error('Unsupported model class!')
end

%if nExtraArgsOut == 1
%    varargout{1} = latentTrajectoriesCellArray;
%end

end