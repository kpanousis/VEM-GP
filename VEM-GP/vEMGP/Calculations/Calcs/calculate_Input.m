function [ u ] = calculate_Input(seq_trial,hypers)
%CALCULATE_INPUT Calculate the input u 
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%8 August 2015
%Inputs: seq_trial: the seq struct for the specific trial
%        hypers: the current set of hyperparameters

%Outputs: u: the recalculated input

%Assert the correct number of hypers
assert(numel(hypers)==6 || numel(hypers(1)),'Wrong Number of Hyperparameters');


%get the click times, signs and times of interest for this trial
clickTimes=seq_trial.clickTimes;
clickSigns=seq_trial.clickSigns;

times=seq_trial.times;
t=times(2:end);
t_prime=times(1:end-1);

%assign hyperparameters
%the else statement is for when calculate input is called for recalculation
%of the input for the popSpike_dyn package
if (numel(hypers)==6)
    lambda=hypers(4);
    phi=hypers(5);
    log_tauPhi=hypers(6);
else
    lambda=hypers;
    phi=0.8;
    log_tauPhi=log(0.05);
end

%Initialize vector u and set the input at t=1 to zero
u=zeros(1,size(seq_trial.y,2));
u(:,1)=0;

%Calculate the ciVals
ciVals=cAndDerivs(clickTimes,0,exp(log_tauPhi),phi);

%Calculate the input
for ts=1:size(u,2)-1
    value=0;
    for clickIter=1:numel(clickTimes)
        if (clickTimes(clickIter)>=t_prime(ts) && t(ts)>clickTimes(clickIter))
            value=value+ciVals(clickIter)*clickSigns(clickIter)*exp(lambda*(t(ts)-clickTimes(clickIter)));
        end
    end
    u(1,ts+1)=value;
end


end

