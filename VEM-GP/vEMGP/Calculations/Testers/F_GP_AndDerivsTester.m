function [ f_gp,varargout] = F_GP_AndDerivsTester( chypers,sessionData,varargin)
%F_GP_ANDDERIVSTESTER Wrapper validating F_GP_AndDerivs
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%17 June 2015

%Inputs: chypers: the hypers to optimize; in case we optimize less than 6
%        parameters, it is necessary to provide the remaining through the varargin
%        argument
%        sessionData: the struct containing all the necessary info for each
%        session. This info includes C and d of the GLM and the seq struct
%        that contains times, clickTimes etcetera
if (length(varargin)>=1)
    hypers=[chypers; varargin{1}];
else
    hypers=chypers;
end

%Input checking
assert(numel(hypers)==6,'Something Wrong With the Number of Hyperparameters');

%Initialize the f_gp and the derivsto zero
f_gp=0;
derivsOut = zeros(size(hypers));

%For each session get the seq struct and for each trial in the seq struct
%calculate the f_gp and the derivatives and add them to the current value
%of the free energy and the derivatives
for i=1:numel(sessionData)
    
    seq=sessionData(i).seq;
    
    %for all the trials in seq
    for trialIdx=1:numel(seq)
        
        [f,derivs]=FgpAndDerivs(seq(trialIdx).times',hypers,seq(trialIdx).clickTimes,seq(trialIdx).clickSigns,...
            seq(trialIdx).posterior.xsm, seq(trialIdx).posterior.Vsm, seq(trialIdx).posterior.VVsm);
        
        f_gp=f_gp+f;
        derivsOut = derivsOut+derivs;
    end
end

%Assign the output/derivatives
varargout{1} = derivsOut(1:numel(chypers));

%Change sign for the minimization
f_gp=-f_gp;
varargout{1}=-varargout{1};

end

