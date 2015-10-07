function [ f_gp,derivsOut ] = FgpAndDerivs( time,hypers,clickTimes,clickSigns,mu,Vsm,VVsm )
%FgpAndDerivatives Calculate the F_GP term and all the derivatives with respect to
%the parameters sigma_i, sigma_s, sigma_a and lambda
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%16 June 2015

%Inputs:  times: the times to consider
%         hypers: the hyperparameter vector
%         clickTimes: the times of the clicks
%         clickSigns: the signs of the clicks
%         mu: the posterior mean
%         Vsm: the posteriro variance
%         VVsm: the posterior covariance

%Outputs: f_gp: the calculated f_gp
%         derivsOut: the derivatives of the f_gp with respect to the 6
%         hyperparameters
%% Input checking
nExtraArgsOut = max(nargout-1,0);
if nExtraArgsOut >1
    error('Wrong number of outputs to F_GP_AndDerivs!')
end

%% Hyper Check
if (numel(hypers)~=6)
    error('Wrong Number of Hyperparameters');
end


%% Calculations
%Time convention
t=time(2:end);
t_prime=time(1:end-1);

%Get the parameters
phi=hypers(5);
log_tauPhi=hypers(6);

%Calculate c(t) and the derivatives
%This is a function of Dr. Thomas Desautels -not to be made public or
%redistribute-
[ciVals,ciDerivs]=cAndDerivs(clickTimes,0,exp(log_tauPhi),phi);

%% Calculate the E_{q(\alpha)}[(a_t-\mu_{t|t'})^2] and its derivatives
[expect,derivsExpectCell]=EqAndDerivatives(t,t_prime,ciVals,ciDerivs,clickTimes,clickSigns,mu',Vsm,VVsm,hypers);

%% Calculate the conditional Variance and its derivatives
[condVar,derivsCondVar]=condVarAndDerivs(hypers,t',t_prime',ciVals,ciDerivs,clickTimes);

%% Error Checking
 if (any(condVar==0))
        error(' Zero Conditional Covariance');
 end

%% Calculate the GP part of the free energy
F_GP=log(condVar)+bsxfun(@rdivide,expect,condVar);
f_gp=-0.5*sum(F_GP);

%Add the term for a_1
f_gp=f_gp-0.5*hypers(1)-0.5*bsxfun(@rdivide,(Vsm(1)+mu(1)^2),exp(hypers(1)));

%% Calculate the derivatives
if (nExtraArgsOut==1)
    
    derivs=zeros(size(derivsCondVar));
    
    % Get the derivs of Eq from the cell structure
    for hyperIdx=1:numel(hypers)
        if ~isempty(derivsExpectCell{hyperIdx,1})
            derivs(:,hyperIdx)=derivsExpectCell{hyperIdx,1};
        end
    end
        
    %calculate the (1/condVar)*(d\sigma_{t|t'})/d\Theta_GP
    first_term=bsxfun(@ldivide,condVar,derivsCondVar);
    
    %second long term derivative of expectation*
    %condVar-expectation*derivative_condVar over condVar
    derivsTimesCond=bsxfun(@times,condVar,derivs);
    expectTimesderivsCond=bsxfun(@times,expect,derivsCondVar);
    second_term=derivsTimesCond-expectTimesderivsCond;
    second_term=bsxfun(@rdivide,second_term,condVar.^2);
    
    %% Sum over values,normalize and Adjust the sign
    deriv_gp=sum(first_term+second_term,1);
    %add the terms for a_1
    deriv_gp(1)=deriv_gp(1)+1-bsxfun(@rdivide,(Vsm(1)+mu(1)^2),exp(hypers(1)));
    derivsOut=-0.5*deriv_gp';
    
end

end

