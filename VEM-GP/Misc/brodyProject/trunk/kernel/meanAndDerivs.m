function [mu,varargout] = meanAndDerivs(times,clickTimes,clickSigns,hypers)
%MEANANDDERIVS calculates the mean vector and its derivatives.
%meanAndDerivs.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%5 March, 2014
%Inputs:
%   times: The times of interest for which the kernel should be
%       calculated (ordered from smallest to largest) (note: this property
%       should be useful for efficiency purposes)
%   clickTimes: The times at which the clicks were delivered (ordered from smallest to largest)
%   clickSigns: The signs of the clicks (+/-1)
%   hypers: The column vector of kernel hyperparameters, which control the
%       model.  These are, in order:[logSigmaSquaredI; logSigmaSquaredA;
%       logSigmaSquaredS; logSigmaSquaredN; lambda; phi; logTauPhi]; Note
%       that the hyperparameters which have strictly non-negative values
%       are expressed as natural logarithms.  Note also that
%       hyperparameters 1-4 are not used by the mean function and will
%       consequently return 0 for the derivatives.
%Outputs:
%   mu: The mean vector evaluated at these times
%   varargout: 1 entry: 7x1 cell array of the vectors which comprise the
%       derivatives of the entries of mu with respect to the given
%       hyperparameters; 2 entries: the above, and then the 7x7 cell array
%       comprising the "Hessian" of the mu matrix, where each entry is
%       filled in appropriately; values which are identically zero are left
%       blank, as are values which are redundant above and below the
%       diagonal; thus, a blank cell below the diagonal is a zero matrix.

%% Input checking
nExtraArgsOut = max(nargout-1,0);
if nExtraArgsOut >2
    error('Wrong number of outputs to meanAndDerivs!')
end

assert(size(times,2) ==1,'times must be a column vector!')
assert(size(clickTimes,2) == 1,'clickTimes must be a column vector!')
assert(size(clickSigns,2) == 1,'clickSigns must be a column vector!')
assert(all(abs(clickSigns) == 1),'clickSigns must be composed of +/-1 entries!')
assert(all(size(hypers) == [7,1]) || all(size(hypers) == [6,1]),'Hyperparameter vector must be 7x1 or 6x1')


%% Check for a really dumb case:
if numel(clickTimes) == 0;
    fprintf('You''ve passed an example with no clicks. Check this!')
    mu = zeros(size(times));
    if nExtraArgsOut > 0
        muDerivs = cell(numel(hypers),1);
        for nHyper = 1:numel(hypers)
            muDerivs{nHyper} = zeros(size(times));
        end
        varargout{1} = muDerivs;
        if nExtraArgsOut > 1
            muSecondDerivs = cell(numel(hypers),numel(hypers));
            %The empty entries are interpreted as identically zero,
            %according to our convention.
            varargout{2} = muSecondDerivs;
        end
    end
    return %These values to the invoking function.
end
        


%% Set some key values:
cZero = 0;
%Note: the above assumes we have removed the simultaneous clicks at the
%start of each trial and set our times with respect to this instant.


%Determine M(t) (i.e., the number of clicks so far heard) for each element of times
m_t = mTGen(times,clickTimes);

% %For efficiency, we can replace these in-line, but for now we'll precompute
% %them.
% S = times(:,ones(1,size(times,1))) + times(:,ones(1,size(times,1)))';
% Min_Mt_ttau = min(m_t(:,ones(1,size(m_t,1))), m_t(:,ones(1,size(m_t,1)))');
% %Counts how many clicks have yet arrived: This term should ONLY appear in
% %the k_s terms.  Note that since we're going to use this in indexing, we'll
% %have to add a zero onto the beginning of the objects its indexing (to
% %account for no clicks having yet arrived) and then add one to each of
% %these to get the appropriate index into that vector.
% AbsDiff = abs(times(:,ones(1,size(times,1))) - times(:,ones(1,size(times,1)))');
% %Above: |t-tau|, so we can write directionless
% minTime = min(times(:,ones(1,size(times,1))), times(:,ones(1,size(times,1)))');

%% Unpack the hyperparameters
nHypers = numel(hypers);
if  nHypers == 7
    % logSigmaSquaredI = hypers(1);
    % logSigmaSquaredA = hypers(2);
    % logSigmaSquaredS = hypers(3);
    % logSigmaSquaredN = hypers(4);
    lambda = hypers(5);
    phi = hypers(6);
    logTauPhi = hypers(7);
    
    loc_lambdaVal = 5;
    loc_phi = 6;
    loc_logTauPhi = 7;
elseif nHypers == 6
    lambda = hypers(4);
    phi = hypers(5);
    logTauPhi = hypers(6);
    
    loc_lambdaVal = 4;
    loc_phi = 5;
    loc_logTauPhi = 6;
end

%% Get the values of c and its derivatives we need
% These correspond to the click times.
if nExtraArgsOut == 0
    ciVals = cAndDerivs(clickTimes,cZero,exp(logTauPhi),phi);
elseif nExtraArgsOut == 1
    %cDerivs columns are in order: [d_tauPhi, d_phi]
    [ ciVals, cDerivs] = cAndDerivs(clickTimes,cZero,exp(logTauPhi),phi);
elseif nExtraArgsOut == 2
    %cSecondDerivs column order: [d_tauPhi2,d_phi2,d_tauPhi_d_phi]
    [ ciVals, cDerivs, cSecondDerivs] = cAndDerivs(clickTimes,cZero,exp(logTauPhi),phi);
end

%% Initialize for calculating k and its derivatives
%Instantiate the key storage matrices
if nExtraArgsOut > 0
    muDerivs = cell(nHypers,1);
end
if nExtraArgsOut > 1
    muSecondDerivs = cell(nHypers,nHypers);
end

%Now, run through the procedure calculating the value of mu and the
%derivatives.  The only contribution to the mean function is from the
%incorporation of the individual clicks.

%Calculate the summations:
individualSumEntries_mu = clickSigns .* ciVals .* exp(-lambda * clickTimes);
cumSumProdTerms_muWithZero = [0;cumsum(individualSumEntries_mu)];
mu = exp(lambda * times) .* cumSumProdTerms_muWithZero(m_t + 1);

if nExtraArgsOut > 0
    
    %Calculate the first derivative terms
    individualSumEntries_2_dmudlambda = clickTimes .* individualSumEntries_mu;
    cumSumProdTerms_2_dmudlambdaWithZeros = [0;cumsum(individualSumEntries_2_dmudlambda)];
    
    individualSumEntries_dmudphi = clickSigns .* cDerivs(:,2) .* exp(-lambda * clickTimes);
    cumSumProdTerms_dmudphi = [0;cumsum(individualSumEntries_dmudphi)];
    
    individualSumEntries_dmudtauPhi = clickSigns .* cDerivs(:,1) .* exp(-lambda * clickTimes);
    cumSumProdTerms_dmudtauPhi = [0;cumsum(individualSumEntries_dmudtauPhi)];
    
    muDerivs{loc_lambdaVal,1} = times .* mu - exp(lambda*times) .* cumSumProdTerms_2_dmudlambdaWithZeros(m_t + 1);
    muDerivs{loc_phi,1} =                     exp(lambda * times) .* cumSumProdTerms_dmudphi(m_t + 1);
    muDerivs{loc_logTauPhi,1} = exp(logTauPhi) *    exp(lambda * times) .* cumSumProdTerms_dmudtauPhi(m_t + 1);
   if nExtraArgsOut > 1
      %Calculate the second derivative terms 
      individualSumEntries_d2mudlambda2             = clickTimes.^2 .* clickSigns .* ciVals .* exp(-lambda * clickTimes);
      individualSumEntries_d2mudphi2                = clickSigns .* cSecondDerivs(:,2) .* exp(-lambda * clickTimes);
      individualSumEntries_d2mudtauPhi2             = clickSigns .* cSecondDerivs(:,1) .* exp(-lambda * clickTimes);
      individualSumEntries_d2mudphidtauphi          = clickSigns .* cSecondDerivs(:,3) .* exp(-lambda * clickTimes);
      individualSumEntries_d2mudphidlambda          = clickSigns .* clickTimes .* cDerivs(:,2).* exp(-lambda * clickTimes);
      individualSumEntries_d2mudtauPhidlambda       = clickSigns .* clickTimes .* cDerivs(:,1).* exp(-lambda * clickTimes);
      
      cumSumProdTerms_d2mudlambda2                  = [0;cumsum(individualSumEntries_d2mudlambda2)];
      cumSumProdTerms_d2mudphi2                     = [0;cumsum(individualSumEntries_d2mudphi2)];
      cumSumProdTerms_d2mudtauPhi2                  = [0;cumsum(individualSumEntries_d2mudtauPhi2)];
      cumSumProdTerms_d2mudphidtauphi               = [0;cumsum(individualSumEntries_d2mudphidtauphi)];
      cumSumProdTerms_d2mudphidlambda               = [0;cumsum(individualSumEntries_d2mudphidlambda)];
      cumSumProdTerms_d2mudtauPhidlambda            = [0;cumsum(individualSumEntries_d2mudtauPhidlambda)];
      
      
      muSecondDerivs{loc_lambdaVal,loc_lambdaVal} = times.^2 .* mu  + exp(lambda * times) .* (-2 * times .* cumSumProdTerms_2_dmudlambdaWithZeros(m_t + 1)  + cumSumProdTerms_d2mudlambda2(m_t + 1))  ;
      muSecondDerivs{loc_phi,loc_phi} =                                    exp(lambda * times) .* cumSumProdTerms_d2mudphi2(m_t + 1);
      muSecondDerivs{loc_logTauPhi,loc_logTauPhi} = muDerivs{loc_logTauPhi,1} + exp(2*logTauPhi) * exp(lambda * times) .* cumSumProdTerms_d2mudtauPhi2(m_t + 1);  %Since we have the derivative wrt log tau_phi
      
      muSecondDerivs{loc_phi,loc_lambdaVal} = times .* muDerivs{loc_phi,1}                      - exp(lambda * times) .* cumSumProdTerms_d2mudphidlambda(m_t + 1);
      muSecondDerivs{loc_logTauPhi,loc_lambdaVal} = times .* muDerivs{loc_logTauPhi,1}     - exp(logTauPhi) * exp(lambda * times) .* cumSumProdTerms_d2mudtauPhidlambda(m_t + 1);
      muSecondDerivs{loc_logTauPhi,loc_phi} =                 exp(logTauPhi)   * exp(lambda * times) .* cumSumProdTerms_d2mudphidtauphi(m_t + 1); %Since we have the derivative wrt log tau_phi
      
      
      
   end
    
end

%% Output handling
if nExtraArgsOut > 0
    varargout{1} = muDerivs;
    if nExtraArgsOut > 1
        varargout{2} = muSecondDerivs;
    end
end

