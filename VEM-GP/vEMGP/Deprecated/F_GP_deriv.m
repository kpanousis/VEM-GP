function [ deriv_f_gp ] = F_GP_deriv(hypers,expect,derivsExpectCell,condVar,derivsCond )
%F_GP_DERIV Calculate the derivative of the GP term of the free energy with
%respect to the parameter indicated by the param_handler function
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%16 June 2015


%% Iteratively calculate the derivative


% Get derivs from matrix
derivs=zeros(size(derivsCond));
for hyperIdx=1:numel(hypers)
    if ~isempty(derivsExpectCell{hyperIdx,1})
        derivs(:,hyperIdx)=derivsExpectCell{hyperIdx,1};
    end
end
% assignin('base','derivsCell',derivsCell);
% assignin('base','derivs',derivs);
% assignin('base','condDerivs',derivsCond);
if (any(abs(condVar-0)<eps))
    hypers
    error('Zero Conditional Variance');
end


condVar2=repmat(condVar,1,numel(hypers));
first_term=(1./condVar2).*derivsCond;


%second long term derivative of expectation*
%condVar-expectation*derivative_condVar over condVAr
expect=repmat(expect,1,numel(hypers));

second_term=(derivs.*condVar2-(expect.*derivsCond))./(condVar2.^2);

%% Adjust the sign
deriv_gp=sum(first_term+second_term,1);
deriv_f_gp=-0.5*deriv_gp';

end

% deriv_gp=zeros(1,numel(hypers));
% for hyperIdx=1:numel(hypers)
%     for t=1:numel(t)
%         first_term=(1/condVar(t))*derivsCond(t,hyperIdx);
%         second_term=(derivs(t,hyperIdx)*condVar(t)-expect(t)*derivsCond(t,hyperIdx))/(condVar(t)^2);
%         deriv_gp(hyperIdx)=deriv_gp(hyperIdx)+first_term+second_term;
%     end
% end
% deriv_gp