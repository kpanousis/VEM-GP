function[expect,varargout]=EqAndDerivatives(t,t_prime,ciVals,ciDerivs,clickTimes,clickSigns,mu,Vsm,VVsm,hypers)
%Calculate the derivatives of the expectation of (a_t-\mu_{t|t'})^2
%under the posterior distribution q
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%4 June 2015
%Using the equations derived, calculate the derivatives with respect to
%parameter lambda

%Inputs: t: the times considered
%        t_prime: the t_{i-1} vector
%        clickTimes: the times of the clicks
%        clickSigns: the signs of the clicks +/-1
%        mu: the posterior mean
%        Vsm: the Cov(x(t),x(t)|y)
%        VVsm: the Cov(x(t+1),x(t)|y)
%        hypers: the 6 hyperparameters of the GP

%Outputs: expect: the calculated expectation
%         varargout: the derivs (if needed)


%% Input/Output Checking
if (numel(hypers)~=6)
    error('Wrong Number of GP Hyperparameters!');
end

%% Calculate Derivs?
calculateDerivs=false;

nExtraArgsOut = max(nargout-1,0);

if nExtraArgsOut==1
    calculateDerivs=true;
elseif nExtraArgsOut>1
    error('Wrong Number Of Outputs to EqAndDerivatives');
end


%% Get Variables and the positions for the derivatives
lambda=hypers(4);
log_tauPhi=hypers(6);

loc_lambdaVal=4;
loc_phi=5;
loc_logTauPhi=6;


%% Resize-Transform
Vsm_prime=Vsm(1:end-1);
mu_prime=mu(1:end-1);
mu=mu(2:end);
Vsm=Vsm(2:end);

%% Revised Version
% Calculate the expectation with the new conditional mean definition and
% its derivatives

%Initialize expect vector
expect=zeros(numel(t),1);

%In case we need the derivatives, initialize the vectors used
if (calculateDerivs)
    derivs=cell(numel(hypers),1);
    for i=1:numel(hypers)
        derivs{i,1}=zeros(numel(t),1);
    end
    values_derivs=zeros(numel(t),3);
end

%% Main Loop

for i=1:numel(t)
    value=0;
    diff=t(i)-t_prime(i);
    c_t=exp(lambda*diff);
    %% Calculate the sum term S of the calculation
    for clickIter=1:numel(clickTimes)
        if (clickTimes(clickIter)>=t_prime(i) && (t(i)>clickTimes(clickIter)))
            value=value+clickSigns(clickIter)*ciVals(clickIter)*exp(lambda*(t(i)-clickTimes(clickIter)));
            
            %If we need the derivatives calculate dS/dtheta_i
            if (calculateDerivs)
                % Derivative with respect to lambda
                values_derivs(i,1)=values_derivs(i,1)+clickSigns(clickIter)*ciVals(clickIter)*exp(lambda*(t(i)-clickTimes(clickIter)))...
                    *(t(i)-clickTimes(clickIter));
                %S dS/d(\log tau Phi)
                values_derivs(i,2)=values_derivs(i,2)+clickSigns(clickIter)*exp(lambda*(t(i)-clickTimes(clickIter)))...
                    *ciDerivs(clickIter,1);
                %wrt to phi
                values_derivs(i,3)=values_derivs(i,3)+clickSigns(clickIter)*exp(lambda*(t(i)-clickTimes(clickIter)))...
                    *ciDerivs(clickIter,2);
            end
            
        end
    end
    %% Calculate the expectation at time t(i)
   expect(i)=Vsm(i)-2*c_t*VVsm(i)+c_t^2*Vsm_prime(i)+(mu(i)-c_t*mu_prime(i)-value)^2;

     
    %% Derivatives Final Calculation
    if (calculateDerivs)
        %%Derivative of lambda
        derivs{loc_lambdaVal,1}(i)=-2*((VVsm(i)+mu(i)*mu_prime(i))*c_t*(diff)...
            +mu(i)*values_derivs(i,1))...
            +2*(diff)*c_t^2*(Vsm_prime(i)+mu_prime(i)^2)...
            +2*(diff*mu_prime(i)*c_t*value+mu_prime(i)*c_t*values_derivs(i,1))...
            +2*value*values_derivs(i,1);
        
        %%Derivative of phi
        derivs{loc_phi,1}(i)=-2*mu(i)*values_derivs(i,3)...
            +2*(c_t*mu_prime(i))*(values_derivs(i,3))...
            +2*value*values_derivs(i,3);
        
        %%Derivative of log Tau phi
        derivs{loc_logTauPhi,1}(i)=exp(log_tauPhi)*(-2*mu(i)*values_derivs(i,2)...
            +2*(c_t*mu_prime(i))*(values_derivs(i,2))...
            +2*value*values_derivs(i,2));
    end
end

%% Assign varargout
if (calculateDerivs)
    varargout{1}=derivs;
end

end