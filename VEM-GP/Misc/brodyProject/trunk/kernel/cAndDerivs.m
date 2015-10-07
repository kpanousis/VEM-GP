function [ ciOut, varargout] = cAndDerivs( tiRel,cZero,tauPhi,phi)
%CANDDERIVS calculates the values and derivatives of the c latent process
%in the Brody evidence accumulator model.
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%25 February, 2014
%   Using the recursion relationship (derived 30 Aug, 2013, rederived in
%   somewhat different terms 25 Feb, 2014), calulates the ci values (and
%   derivatives) corresponding to the relative times given by tiRel.
%   Inputs:
%       tiRel: column vector of relative times at which we need ci and
%           its derivatives (the click times).
%       cZero: assumed value of c at time zero. Typically assumed to be
%           zero if t = 0 is the time of the double click arrival.
%       tauPhi: Time constant of the exponential decay of c(t) toward 1.
%       phi:  Scaling constant for c: each time a click arrives, c(tPlus) =
%       c(tMinus) * phi;
%   Outputs:
%       ciOut: Column vector of values of ci at the times tiRel.
%       varargout: If present, contains dcidvarsOut,d2cidvars2Out, which contain
%       respectively the first and second derivatives with respect to
%       tauPhi and phi, in |tiRel| x 2 and  |tiRel| x 3 arrays.  The column
%       order in these arrays is: [d_tauPhi, d_phi] and [d_tauPhi2,d_phi2,
%       d_tauPhi d_phi]

%% Input checking
nExtraArgsOut = max(nargout-1,0);
if nExtraArgsOut >2
    error('Wrong number of outputs to cAndDerivs!')
end

assert(size(tiRel,2) ==1,'tiRel must be a column vector!')


%% Initialize the recursive calculation
%Set initial values
tLast = 0;
cLast = cZero;
ciOut = zeros(size(tiRel));

%If we're calculating derivatives
if nExtraArgsOut >0
    dcdtauPhiLast = 0;
    dcdphiLast = 0;
    dcidvarsOut = zeros(size(tiRel,1),2);
end

%If we're calculating second derivatives
if nExtraArgsOut >1
    d2cdtauPhi2Last = 0;
    d2cdphi2Last = 0;
    d2cdphidtauPhiLast = 0;
    d2cidvars2Out = zeros(size(tiRel,1),3);
end
%% Do the main recursive calculation

%Main loop
for timeIdx = 1:size(tiRel,1)
    ciOut(timeIdx) = 1 + exp(-(tiRel(timeIdx) - tLast)/tauPhi)*(phi*cLast -1);
    if nExtraArgsOut > 0
        dcidvarsOut(timeIdx,1) = exp(-(tiRel(timeIdx) - tLast)/tauPhi) * (   ((tiRel(timeIdx) - tLast)/(tauPhi^2)) * (phi*cLast -1)   +   phi*dcdtauPhiLast); %dtauPhi
        dcidvarsOut(timeIdx,2) = exp(-(tiRel(timeIdx) - tLast)/tauPhi) * (cLast + phi * dcdphiLast); %dphi
    end
    
    if nExtraArgsOut > 1
        %2nd wrt tauPhi:
        nastyTerm = exp(-(tiRel(timeIdx) - tLast)/tauPhi) * ((-2*(tiRel(timeIdx) - tLast)/(tauPhi^3))* (phi*cLast -1) + ((tiRel(timeIdx) - tLast)/(tauPhi^2)) * (phi*dcdtauPhiLast) + phi * d2cdtauPhi2Last); %this term results fro the differentiation of the right term in dtauPhi with respect to tauPhi
        d2cidvars2Out(timeIdx,1) = ((tiRel(timeIdx) - tLast)/(tauPhi^2)) * dcidvarsOut(timeIdx,1)   + nastyTerm; %d_tauPhi2
        
        %2nd wrt phi
        d2cidvars2Out(timeIdx,2) = exp(-(tiRel(timeIdx) - tLast)/tauPhi) * (2*dcdphiLast + phi*d2cdphi2Last); %d_phi2
        
        %mixed 2nd
        d2cidvars2Out(timeIdx,3) = ((tiRel(timeIdx) - tLast)/(tauPhi^2)) * dcidvarsOut(timeIdx,2)   +     exp(-(tiRel(timeIdx) - tLast)/tauPhi) * (dcdtauPhiLast + phi * d2cdphidtauPhiLast); %d_tauPhi d_phi
    end
    
    %Cleanup and set values for next step
    tLast = tiRel(timeIdx);
    cLast = ciOut(timeIdx);
    if nExtraArgsOut > 0
        dcdtauPhiLast = dcidvarsOut(timeIdx,1); %dtauPhi
        dcdphiLast= dcidvarsOut(timeIdx,2); %dphi
    end
    
    if nExtraArgsOut > 1
        d2cdtauPhi2Last = d2cidvars2Out(timeIdx,1); %d_tauPhi2
        d2cdphi2Last = d2cidvars2Out(timeIdx,2); %d_phi2
        d2cdphidtauPhiLast = d2cidvars2Out(timeIdx,3); %d_tauPhi d_phi
    end
end

%% Output handling
if nExtraArgsOut >0
    varargout{1} = dcidvarsOut;
    if nExtraArgsOut > 1
        varargout{2} = d2cidvars2Out;
    end
end
end

