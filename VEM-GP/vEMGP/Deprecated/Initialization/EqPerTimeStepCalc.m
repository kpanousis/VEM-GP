function [ expectation] = EqPerTimeStepCalc(t,t_prime,clickTimes,hypers,expect,clicksEnabled)
%EQPERTIMESTEPCALC Calculate Eq[(a_t-\mu_t)^2]/denominator per time step for \sigma_s
%or \sigma_a (depending on the value of clickEnable)
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%12 June 2015
%Inputs: times: the times considered
%        clickTimes: the times of the clicks
%        clickSigns: the signs of the clicks +/-1
%        mu: the posterior mean
%        Vsm: the Cov(x(t),x(t)|y)
%        VVsm: the Cov(x(t+1),x(t)|y)
%        hypers: the 6 hyperparameters of the GP
%        clicksEnabled: boolean to determine if we are looking for times
%        with clicks only or without clicks only
%Outputs: expectation: the calculated expectation normalized

%% Input checking

if (numel(hypers)~=6)
    error('Wrong Number of Hyperparameters');
end

% Get Hypers
lambda=hypers(4);
log_tauPhi=hypers(5);
phi=hypers(6);

%% Main loop calc
expectation=0;
for j=1:numel(t)
    % if clicksEnabled is false, skip times with clicks else consider only
    % time with clicks
    flag_continue=false;
    if (~clicksEnabled)
        
        %if a click is in [t',t), skip that time step
        for clickIter=1:numel(clickTimes)
            if (clickTimes(clickIter)>=t_prime(j) && clickTimes(clickIter)<t(j))
                flag_continue=true;
                break;
            end
        end
        
        if (flag_continue==true)
            flag_continue=false;
            continue;
        end
        
        if (lambda==0)
            denom=t(j)-t_prime(j);
        else
            denom=(expm1(2*lambda*(t(j)-t_prime(j))))/(2*lambda);
        end
         %add the normalized term to the current value of the expectation
         expectation=expectation+expect(j)/denom;
    else
        ciVals=cAndDerivs(clickTimes,0,exp(log_tauPhi),phi);
        for clickIter=1:numel(clickTimes)
            if (clickTimes(clickIter)>=t_prime(j) && clickTimes(clickIter)<t(j))
               %get c values
                %is it really a sum? here only one time step
                sum_term=ciVals(clickIter)^2*exp(2*lambda*(t(j)-clickTimes(clickIter)));
                denom=sum_term;
                expectation=expectation+expect(j)/denom;
            end
        end
    end
    
   
    
end

end

