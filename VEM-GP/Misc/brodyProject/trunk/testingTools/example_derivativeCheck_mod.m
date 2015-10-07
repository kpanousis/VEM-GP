function example_derivativeCheck_mod(functionHandle,wTest,varargin)
%Modified from the example code from Mark Schmidt's minFunc_2012 toolbox, available here:
%http://www.di.ens.fr/~mschmidt/Software/minFunc.html
%
%Modifications by:
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%3 March, 2014


%functionHandle = @cAndDerivsTester;
%nVars = 2;
%wTest = [exp(randn);rand];
%functionHandle = @kernelAndDerivsTester;
%nVars = 7;
%wTest = randn(nVars,1);

%% Input handling
if ischar(functionHandle)
    functionHandle = str2func(functionHandle);
end

% %% Return the value of the function at the specified parameters
 fOut = feval(functionHandle,wTest,varargin{:})
% 
%% Checks of the derivative
% fprintf('Testing gradient using forward-differencing...\n');
derivativeOrder = 1;
% derivativeCheck(functionHandle,wTest,derivativeOrder,1,varargin{:});
% 
fprintf('Testing gradient using central-differencing...\n');
derivativeCheck(functionHandle,wTest,derivativeOrder,2,varargin{:});

% % fprintf('Testing gradient using complex-step derivative...\n');
% % derivativeCheck(functionHandle,wTest,derivativeOrder,3);
% 
% fprintf('\n\n\n');
% pause
% 
% % fprintf('Testing Hessian using forward-differencing\n');
% % derivativeOrder = 2;
% % derivativeCheck(functionHandle,wTest,derivativeOrder,1);

% fprintf('Testing Hessian using central-differencing\n');
% derivativeOrder = 2;
% derivativeCheck(functionHandle,wTest,derivativeOrder,2,varargin{:});

% fprintf('Testing Hessian using complex-step derivative\n');
% derivativeOrder = 2;
% derivativeCheck(functionHandle,wTest,derivativeOrder,3);

% fprintf('\n\n\n');
% pause
% 
% fprintf('Testing gradient using fastDerivativeCheck...\n');
% derivativeOrder = 1;
% %fastDerivativeCheck(functionHandle,wTest,derivativeOrder,1);
% fastDerivativeCheck(functionHandle,wTest,derivativeOrder,2,varargin{:});
% %fastDerivativeCheck(functionHandle,wTest,derivativeOrder,3);
% 
% fprintf('\n\n\n');
% pause
% 
% fprintf('Testing Hessian using fastDerivativeCheck...\n');
% derivativeOrder = 2;
% %fastDerivativeCheck(functionHandle,wTest,derivativeOrder,1);
% fastDerivativeCheck(functionHandle,wTest,derivativeOrder,2,varargin{:});
% %fastDerivativeCheck(functionHandle,wTest,derivativeOrder,3);


end
