function [k,varargout] = kernelAndDerivs(times,clickTimes,clickSigns,hypers)
%KERNELANDDERIVS calculates the kernel matrix and its derivatives.
%kernelAndDerivs.m
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%28 February, 2014
%Inputs:
%   times: The times of interest for which the kernel should be
%       calculated (ordered from smallest to largest) (note: this property
%       should be useful for efficiency purposes)
%   clickTimes: The times at which the clicks were delivered (ordered from smallest to largest)
%   clickSigns: The signs of the clicks (+/-1)
%   hypers: The column vector of kernel hyperparameters, which control the
%       model.  These are, in order:[logSigmaSquaredI; logSigmaSquaredA;
%       logSigmaSquaredS; logSigmaSquaredN; lambdaVal; phi; logTauPhi]; Note
%       that the hyperparameters which have strictly non-negative values
%       are expressed as natural logarithms.
%Outputs:
%   k: The kernel matrix evaluated at these times
%   varargout: 1 entry: 7x1 cell array of the matrices which comprise the
%       derivatives of the entries of k with respect to the given
%       hyperparameters; 2 entries: the above, and then the 7x7 cell array
%       comprising the "Hessian" of the k matrix, where each entry is
%       filled in appropriately; values which are identically zero are left
%       blank, as are values which are redundant above and below the
%       diagonal; thus, a blank cell below the diagonal is a zero matrix.


%% Input checking
nExtraArgsOut = max(nargout-1,0);
if nExtraArgsOut >2
    error('Wrong number of outputs to kernelAndDerivs!')
end

assert(size(times,2) ==1,'times must be a column vector!')
assert(size(clickTimes,2) == 1,'clickTimes must be a column vector!')
if ~isempty(clickSigns)
    assert(size(clickSigns,2) == 1,'clickSigns must be a column vector!')
    assert(all(abs(clickSigns) == 1),'clickSigns must be composed of +/-1 entries!')
end
assert(all(size(hypers) == [7,1]) || all(size(hypers) == [6,1]),'Hyperparameter vector must be either 7x1 (including a white noise term) or 6x1')

%% Set some key values:
cZero = 0;
%Note: the above assumes we have removed the simultaneous clicks at the
%start of each trial and set our times with respect to this instant.


%Determine M(t) for each element of times
m_t = mTGen(times,clickTimes);

%For efficiency, we can replace these in-line, but for now we'll precompute
%them.
S = times(:,ones(1,size(times,1))) + times(:,ones(1,size(times,1)))';
Min_Mt_ttau = min(m_t(:,ones(1,size(m_t,1))), m_t(:,ones(1,size(m_t,1)))');  
%Counts how many clicks have yet arrived: This term should ONLY appear in
%the k_s terms.  Note that since we're going to use this in indexing, we'll
%have to add a zero onto the beginning of the objects its indexing (to
%account for no clicks having yet arrived) and then add one to each of
%these to get the appropriate index into that vector.
AbsDiff = abs(times(:,ones(1,size(times,1))) - times(:,ones(1,size(times,1)))');
%Above: |t-tau|, so we can write directionless
minTime = min(times(:,ones(1,size(times,1))), times(:,ones(1,size(times,1)))');


%% Unpack the hyperparameters
nHypers = numel(hypers);
if  nHypers == 7
    %Determine which locations the appropriate values are within the
    %hyperparamters vector
    loc_logSSI = 1;
    loc_logSSA = 2;
    loc_logSSS = 3;
    loc_logSSN = 4;
    loc_lambdaVal = 5;
    loc_phi = 6;
    loc_logTauPhi = 7;
    
    % Unpack
    logSigmaSquaredI = hypers(loc_logSSI);
    logSigmaSquaredA = hypers(loc_logSSA);
    logSigmaSquaredS = hypers(loc_logSSS);
    logSigmaSquaredN = hypers(loc_logSSN);
    lambdaVal = hypers(loc_lambdaVal);
    phi = hypers(loc_phi);
    logTauPhi = hypers(loc_logTauPhi);
elseif nHypers == 6
    %Determine which locations the appropriate values are within the
    %hyperparamters vector
    loc_logSSI = 1;
    loc_logSSA = 2;
    loc_logSSS = 3;
    %loc_logSSN = 4;
    loc_lambdaVal = 4;
    loc_phi = 5;
    loc_logTauPhi = 6;
    
    % Unpack
    logSigmaSquaredI = hypers(loc_logSSI);
    logSigmaSquaredA = hypers(loc_logSSA);
    logSigmaSquaredS = hypers(loc_logSSS);
    %logSigmaSquaredN = ;
    lambdaVal = hypers(loc_lambdaVal);
    phi = hypers(loc_phi);
    logTauPhi = hypers(loc_logTauPhi);
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
k = zeros(size(times,1));
if nExtraArgsOut > 0
    kDerivs = cell(nHypers,1);
end
if nExtraArgsOut > 1
    kSecondDerivs = cell(nHypers,nHypers);
end

%Now, run through the procedure calculating the value of k and the
%derivatives.  Clear what can be cleared after each section.

%% Components with respect to sigma^2 I
if nHypers == 7
    k = k + exp(logSigmaSquaredN) * eye(size(k));
    if nExtraArgsOut > 0
        kDerivs{loc_logSSN,1} = k;
    end
    if nExtraArgsOut > 1
        kSecondDerivs{loc_logSSN,loc_logSSN} = k;
    end
end

%% Components with respect to K_i
if nExtraArgsOut == 0
    k = k + exp(logSigmaSquaredI) * exp(lambdaVal * times) * exp(lambdaVal * times');
elseif nExtraArgsOut > 0
    kDerivs{loc_logSSI,1} = exp(logSigmaSquaredI) * exp(lambdaVal * times) * exp(lambdaVal * times');
    kDerivs{loc_lambdaVal,1} = kDerivs{loc_logSSI,1} .* (S); %k_i * (t + tau)  %This should be possible to do more efficiently
    if nExtraArgsOut > 1
        kSecondDerivs{loc_logSSI,loc_logSSI} = kDerivs{loc_logSSI,1};
        kSecondDerivs{loc_lambdaVal,loc_logSSI} = kDerivs{loc_lambdaVal,1}; % k_i * (t + tau)^2
        kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kDerivs{loc_lambdaVal,1} .* (S); % k_i * (t + tau)^2
    end
    k = k + kDerivs{loc_logSSI,1};
end

%% Components with respect to K_a

%Below: in calculating k_a or dkadlambda for lambdaVal very small but not
%identically zero, we have a potential problem with numerical precision, I
%think.

if nExtraArgsOut == 0
    if lambdaVal == 0
        k = k + exp(logSigmaSquaredA) * minTime;
    else
        k = k +         exp(logSigmaSquaredA) * exp(lambdaVal *(AbsDiff)) .* (  expm1(2 * lambdaVal * minTime) ) / ( 2 * lambdaVal);
    end
elseif nExtraArgsOut > 0
    if lambdaVal == 0
        kDerivs{loc_logSSA,1} = exp(logSigmaSquaredA) * minTime;
        if nExtraArgsOut > 1
            kSecondDerivs{loc_lambdaVal,loc_logSSA} = exp(logSigmaSquaredA) * (times * times');
            kDerivs{loc_lambdaVal,1} = kDerivs{loc_lambdaVal,1} + kSecondDerivs{loc_lambdaVal,loc_logSSA};
        else
            kDerivs{loc_lambdaVal,1} = kDerivs{loc_lambdaVal,1} + exp(logSigmaSquaredA) * (times * times');  % derivative wrt other terms so far + dKa/dlambda
        end
    else
        kDerivs{loc_logSSA,1} =  exp(logSigmaSquaredA) * exp(lambdaVal *(AbsDiff)) .* (  expm1(2 * lambdaVal * minTime) ) / ( 2 * lambdaVal);
        
        %Line below is a problem: we have a potential divide by zero.
        %Somehow, we need to rearrange such that for t = 0, we get rid of
        %this problem, while still maintaining the quick formulation given
        %here.
        minTimeEqZero = minTime == 0;
        intermediateMatrix = zeros(size(AbsDiff));
        intermediateMatrix( ~minTimeEqZero) = AbsDiff(~minTimeEqZero) - 1/lambdaVal + 2 * minTime(~minTimeEqZero) .* (1 + 1 ./ (expm1(2 * lambdaVal * minTime(~minTimeEqZero))));
        intermediateMatrix(  minTimeEqZero) = AbsDiff( minTimeEqZero);  %This is because the limit as t --> 0 of 2t(exp(2 lambda t))/(exp(2 lambda t) - 1) = 1/lambda, which zeroes out.
        if nExtraArgsOut > 1
            kSecondDerivs{loc_lambdaVal,loc_logSSA} = kDerivs{loc_logSSA,1} .* intermediateMatrix;
            kDerivs{loc_lambdaVal,1} = kDerivs{loc_lambdaVal,1} + kSecondDerivs{loc_lambdaVal,loc_logSSA};  %dKa/dlambda
        else
            kDerivs{loc_lambdaVal,1} = kDerivs{loc_lambdaVal,1} + kDerivs{loc_logSSA,1} .* intermediateMatrix;  % derivative wrt other terms so far + dKa/dlambda
        end
        clear intermediateMatrix %will get rid of that matrix in the future.
        clear minTimeEqZero
    end
    
    
    if nExtraArgsOut > 1
        kSecondDerivs{loc_logSSA,loc_logSSA} = kDerivs{loc_logSSA,1};  %This is k_a
        if lambdaVal == 0
            % timesSquared = times.^2;
            % timesCubed = times .* timesSquared;
            % kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} + exp(logSigmaSquaredA) / 3 * ( min(timesCubed(:,ones(1,size(times,1))), timesCubed(:,ones(1,size(times,1)))') + 3 * minTime .* max(timesSquared(:,ones(1,size(times,1))),timesSquared(:,ones(1,size(times,1)))'));
              
            %timesSquared = times.^2;
            %timesCubed = times .* timesSquared;
            kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} + (exp(logSigmaSquaredA) / 3) * ( minTime.^3 + 3 * (max(times(:,ones(1,size(times,1))), times(:,ones(1,size(times,1)))').^2).* minTime );
            
            
            %Above: derivative wrt other terms so far + d2Ka/dlambda2 = 1/3
            %sigma_a^2 * (t^3 + t * tau^2)
        else
            % intermediateMatrixA = (times(:,ones(1,size(times,1))) - times(:,ones(1,size(times,1)))').^2;
            % intermediateMatrixB = - 2/lambdaVal * AbsDiff;
            % intermediateMatrixC = (1 - exp(-2 * lambdaVal * minTime)).^(-1) .* (4 * (times * times') - 4/ lambdaVal * minTime);
            % intermediateMatrix = intermediateMatrixA + intermediateMatrixB   + intermediateMatrixC + 2 * lambdaVal^(-2);
            % kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} + kDerivs{loc_logSSA,1}.* intermediateMatrix;
            
            % % intermediateMatrixA = ;
            % % intermediateMatrixB = ;
            % % intermediateMatrixC = ;
            % % intermediateMatrix = ;
            %kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} + kDerivs{loc_logSSA,1}.* ((times(:,ones(1,size(times,1))) - times(:,ones(1,size(times,1)))').^2 + - 2/lambdaVal * AbsDiff   + (1 - exp(-2 * lambdaVal * minTime)).^(-1) .* (4 * (times * times') - 4/ lambdaVal * minTime) + 2 * lambdaVal^(-2));
            
            kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} +  exp(lambdaVal * AbsDiff) .* ((-AbsDiff.^2 ./lambdaVal) + exp(2 * lambdaVal * minTime) .* (lambdaVal^(-1) * S.^2  - 2 *lambdaVal^(-2) * S) + (2 * lambdaVal^(-2) * AbsDiff) + 2 * lambdaVal^(-3) * expm1(2 * lambdaVal * minTime)) * exp(logSigmaSquaredA)/2;
            
            
            %kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} + exp(logSigmaSquaredA)/2 * exp(lambdaVal * AbsDiff) .* (  (expm1(2*lambdaVal*minTime) .* ((lambdaVal^(-1) * AbsDiff.^2) - (2 * lambdaVal^(-2) * AbsDiff) + 2 *lambdaVal^(-3))) + exp(2*lambdaVal * minTime) .* (lambdaVal^(-1) * 4 * (times * times') - 4 *lambdaVal^(-2) *minTime));
            
            
            
            %Above: derivative wrt other terms so far + d2Ka/dlambda2
            clear('intermediateMatrix','intermediateMatrixA', 'intermediateMatrixB', 'intermediateMatrixC')
        end
    end
    
    %Add k_a onto k
    k = k + kDerivs{loc_logSSA,1};
end

%So far: have computed k_n , k_i, k_a, and added them into k_soFar.  Have
%also computed first and second derivatives of k with respect to
%log(sigma_n^2), log(sigma_i^2),log(sigma_a^2); these are kDerivs{loc_logSSN,1},
%kDerivs{loc_logSSI,1}, kDerivs{loc_logSSA,1} and kSecondDerivs{loc_logSSN,loc_logSSN}, kSecondDerivs{loc_logSSI,loc_logSSI},
%kSecondDerivs{loc_logSSA,loc_logSSA}.  Have also computed the mixed derivatives of these and
%lambdaVal, which are kSecondDerivs{loc_lambdaVal,loc_logSSN} (empty, since zero),
%kSecondDerivs{loc_lambdaVal,loc_logSSI}, and kSecondDerivs{loc_lambdaVal,loc_logSSA}.  Finally, have computed part
%of kDerivs{loc_lambdaVal,1} and kSecondDerivs{loc_lambdaVal,loc_lambdaVal}, the first and second derivatives
%with respect to lambdaVal.


%% Components with respect to K_s
%Note: this is the only place the terms with respect to phi, log(tauPhi)
%come in.

%Precompute some key values:
expNegTwoLambdaTi = exp(-2 * lambdaVal * clickTimes);
individualProdTerms_Ks = ciVals.^2 .* expNegTwoLambdaTi;
%Create a version with zero at the beginning so we have a "no effects so
%far since no clicks have arrived" object to index to.
cumSumProdTerms_Ks_WithZero = [0; cumsum(individualProdTerms_Ks)];




if nExtraArgsOut == 0
    k = k + exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* cumSumProdTerms_Ks_WithZero(Min_Mt_ttau + 1); %add on k_s
elseif nExtraArgsOut >0
    %Generate the sequences of products we'll need to sum across
    individualProdTerms_dksdlac = clickTimes .* individualProdTerms_Ks;
    individualProdTerms_dksdphi = ciVals .* expNegTwoLambdaTi .*cDerivs(:,2);
    individualProdTerms_dksdtauphi = ciVals .* expNegTwoLambdaTi .*cDerivs(:,1);
    
    %Again, generate cumulative sums with a zero entry at the beginning
    cumSumProdTerms_dksdlac_WithZero = [0;cumsum(individualProdTerms_dksdlac)];
    cumSumProdTerms_dksdphi_WithZero = [0;cumsum(individualProdTerms_dksdphi)];
    cumSumProdTerms_dksdtauphi_WithZero = [0;cumsum(individualProdTerms_dksdtauphi)];
    
    
    kDerivs{loc_logSSS,1} = exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* cumSumProdTerms_Ks_WithZero(Min_Mt_ttau + 1); % k_s
    kDerivs{loc_phi,1} = 2 * exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* cumSumProdTerms_dksdphi_WithZero(Min_Mt_ttau + 1);
    kDerivs{loc_logTauPhi,1} = 2 * exp(logTauPhi) * exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* cumSumProdTerms_dksdtauphi_WithZero(Min_Mt_ttau + 1);

    %Add onto dkdlac
    if nExtraArgsOut == 2
        kSecondDerivs{loc_lambdaVal,loc_logSSS} =            exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* (  (S) .* cumSumProdTerms_Ks_WithZero(Min_Mt_ttau + 1) - 2 *cumSumProdTerms_dksdlac_WithZero(Min_Mt_ttau + 1)) ;
        kDerivs{loc_lambdaVal,1} = kDerivs{loc_lambdaVal,1} + kSecondDerivs{loc_lambdaVal,loc_logSSS};
    else
        kDerivs{loc_lambdaVal,1} = kDerivs{loc_lambdaVal,1} +   exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* (  (S) .* cumSumProdTerms_Ks_WithZero(Min_Mt_ttau + 1) - 2 *cumSumProdTerms_dksdlac_WithZero(Min_Mt_ttau + 1)) ; %dksdlambda
    end
    
    if nExtraArgsOut == 2
        %Generate the sequences of products we'll need to sum across
        individualProdTerms_d2ksdlac2 = clickTimes.^2 .* ciVals.^2 .* expNegTwoLambdaTi;
        individualProdTerms_d2ksdphi2 = expNegTwoLambdaTi.* (cDerivs(:,2).^2 + ciVals .* cSecondDerivs(:,2));
        individualProdTerms_d2ksdtauphi2 = expNegTwoLambdaTi.* (cDerivs(:,1).^2 + ciVals .* cSecondDerivs(:,1));
        individualProdTerms_d2ksdtauphidphi = expNegTwoLambdaTi.* (cDerivs(:,1).*cDerivs(:,2) + ciVals .* cSecondDerivs(:,3));
        individualProdTerms_d2ksdlacdphi = clickTimes.*individualProdTerms_dksdphi;  %we'll need another term from dks_dphi
        individualProdTerms_d2ksdlacdtauphi = clickTimes.*individualProdTerms_dksdtauphi; %and here another term from dks_dtauphi
        
        %Once more, generate versions with a zero entry at the beginning
        cumSumProdTerms_d2ksdlac2_WithZero         = [0;cumsum(individualProdTerms_d2ksdlac2)];
        cumSumProdTerms_d2ksdphi2_WithZero         = [0;cumsum(individualProdTerms_d2ksdphi2)];
        cumSumProdTerms_d2ksdtauphi2_WithZero      = [0;cumsum(individualProdTerms_d2ksdtauphi2)];
        cumSumProdTerms_d2ksdtauphidphi_WithZero   = [0;cumsum(individualProdTerms_d2ksdtauphidphi)];
        cumSumProdTerms_d2ksdlacdphi_WithZero      = [0;cumsum(individualProdTerms_d2ksdlacdphi)];
        cumSumProdTerms_d2ksdlacdtauphi_WithZero   = [0;cumsum(individualProdTerms_d2ksdlacdtauphi)];
        
        %Do the calculations and store the results.
        kSecondDerivs{loc_logSSS,loc_logSSS} = kDerivs{loc_logSSS,1};
        kSecondDerivs{loc_phi,loc_lambdaVal} =                    exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .*    (2 * S .*cumSumProdTerms_dksdphi_WithZero(Min_Mt_ttau + 1)         -4 * cumSumProdTerms_d2ksdlacdphi_WithZero(Min_Mt_ttau + 1));
        kSecondDerivs{loc_phi,loc_phi} = 2 * exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* cumSumProdTerms_d2ksdphi2_WithZero(Min_Mt_ttau + 1);
        kSecondDerivs{loc_logTauPhi,loc_lambdaVal} = exp(logTauPhi) *   exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .*    (2 * S .*cumSumProdTerms_dksdtauphi_WithZero(Min_Mt_ttau + 1)      -4 * cumSumProdTerms_d2ksdlacdtauphi_WithZero(Min_Mt_ttau + 1));
        kSecondDerivs{loc_logTauPhi,loc_phi} = exp(logTauPhi) * 2 * exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .* cumSumProdTerms_d2ksdtauphidphi_WithZero(Min_Mt_ttau + 1);
        kSecondDerivs{loc_logTauPhi,loc_logTauPhi} = kDerivs{loc_logTauPhi,1} + exp(2* logTauPhi) * (2 * exp(logSigmaSquaredS) * exp(lambdaVal * (S)) .*cumSumProdTerms_d2ksdtauphi2_WithZero(Min_Mt_ttau + 1));  %Recall that this will be for LOG(TAUPHI^2) and so has a more complex form
        
        kSecondDerivs{loc_phi,loc_logSSS} = kDerivs{loc_phi,1};
        kSecondDerivs{loc_logTauPhi,loc_logSSS} = kDerivs{loc_logTauPhi,1};
        
        %Add on to d2ksdlac2:
        summand1 = (S).^2 .* cumSumProdTerms_Ks_WithZero(Min_Mt_ttau + 1);
        summand2 =  -4 * (S) .*   cumSumProdTerms_dksdlac_WithZero(Min_Mt_ttau + 1);
        summand3 = 4 *  cumSumProdTerms_d2ksdlac2_WithZero(Min_Mt_ttau + 1);
        kSecondDerivs{loc_lambdaVal,loc_lambdaVal} = kSecondDerivs{loc_lambdaVal,loc_lambdaVal} +  exp(logSigmaSquaredS) * exp(lambdaVal * (S)).* ( summand1 + summand2 + summand3);
    end
    
    %Update k:
    k = k + kDerivs{loc_logSSS,1}; %add on k_s
end 

%% Output handling
if nExtraArgsOut > 0
    varargout{1} = kDerivs;
    if nExtraArgsOut > 1
        varargout{2} = kSecondDerivs;
    end
end

end




