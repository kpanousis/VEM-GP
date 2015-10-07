function testerForThetaUpdateTerms(obj, functionHandles, modelNumber, sessionIdx, times, clickTimes, clickSigns, aMAP, spikeNumbers)
%testerForThetaUpdateTerms
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%University College London
%18 March, 2014
%
%Method which implements some derivative and Hessian consistency checks on
%the terms in the expression maximized to find thetaHat.
%Used to be a script


%Each of the functions below calculates one of the components of that
%expression (per trial).  We'll need to establish some values first.
%Currently, we won't be changing theta_GP, nor will we be changing a[t];
%we'll thus get theta_GP --> K, dK, d2K, K^-1, mu[t] and these will be
%fixed.  Additionally, we'll get a[t], y[t] fixed.  The changing inputs to
%the functions will be theta_GLM (within thetaVec, the vertcat of theta_GLM
%and theta_GP) and lambda (which we'll unfortunately need to calculate
%somehow in the wrapper, or pass down).

%functionHandles = {'@logPThetaGLMAndDerivs_Tester',...
%'@logPThetaGPAndDerivs_Tester',...
%'@logPYCondAThetaAndDerivs_noAdep_Tester',...
%'@logPACondThetaAndDerivs_noAdep_Tester',...
%'@logDetSigmaInvLambdaAndDerivs_noAdep_Tester'};


%Each one of the above functions should pass up f, gradf_theta_GLM,
%Hessianf_theta_GLM
%% Input check
%if ~islogical(spikeNumbers)
%    fprintf('In testerForThetatUpdateTerms.m, spikeNumbers isn''t a logical!  This code shouldn''t have aproblem with this directly, but hasn''t been debugged.\n')
%    fprintf('Ensure that all of logPYCondAThetaAndDerivs_noAdep_Tester, logDetSigmaInvLambdaAndDerivs_noAdep_Tester, \n and logDetTerms_noAdep_Tester have been debugged for this case and disable this message.\n\n')
%end
    

%% Shall we fit theta_GP?
fitThetaGP = true;

%% Shall we calculate the determinant terms in log(p(A | theta))?
doDetTermCalculations = false;


%% Pull down key values from the object:
C = obj.modelSessionParameters{sessionIdx,modelNumber}.C;  %nNeurons x 1
d = obj.modelSessionParameters{sessionIdx,modelNumber}.d;  %nNeurons x 1
W = obj.modelSessionParameters{sessionIdx,modelNumber}.W;  %nNeurons x (nNeurons * nCompsPseudoBasis);


%% Establsh Key Constants
gpHypers = obj.modelGlobalParameters{modelNumber}.gpHypers;

% nNeurons = 1;
% nHistoryFeatures = ;
% nCs = nNeurons;
% nDs = nNeurons;
% nWs = nNeurons * nHistoryFeatures; %could also be thought of as nNeurons,nNeurons * nCompsPseudoBasis
% splits = [nCs,nDs,nWs,nThetaGP];  %Tells us how to split up the hyperparameter vector into the different components.

%Key question: how are we going to order the W elements into the
%hyperparameter vector?  Should be consistent with however we're unfolding
%them other places.
thetaGLM = [reshape(C,[],1);reshape(d,[],1);reshape(W,[],1)];  %This is how we implemented this in calculateFiringRates_BrodyGLM.
splits = [numel(C),numel(d),numel(W),numel(gpHypers)]';  %Column vector

%assumedBasalFiringRate = 10; %Hz
%deltaT = 0.001; %Seconds

%% Establish Priors over theta_GLM, theta_GP
%meanThetaGLM = [zeros(numel(C),1),log(assumedBasalFiringRate * deltaT) *ones(numel(d),1),zeros(numel(W),1)];  %If all of our weights are zero, the cells should fire at (more or less) the basal firing rate with this value chosen.
%covThetaGLM = diag(ones(numel(splits),1));  %This is a bit silly, but what the hell, we're just testing.  This SHOULD scale with deltaT somehow.

%%%Won't need the below for now, since we're just doing the theta_GP terms.
%%meanThetaGP = hypersGP;
%%covThetaGP = diag(ones(size(hypersGP),1));  %Again, a bit silly, but OK for testing.  We expect this to be diagonal, anyway.  %This SHOULD NOT scale at all with deltaT.

%% Do key preparatory calculations:
[K,dK,d2K]      = kernelAndDerivs(times,clickTimes,clickSigns,gpHypers);
%[Kinv]          = kAISInverse(times,clickTimes,gpHypers);
[Kinv]          = kAISInverse2(times,clickTimes,gpHypers);
[mu,dMu,d2Mu]   = meanAndDerivs(times,clickTimes,clickSigns,gpHypers);
%aMAP            = ;
%firingRates     = ;

%% Run the tests


%Arguments needed:
%For logPThetaGLMAndDerivs_Wrapper: (thetaVec,splits,meanThetaGLM,covThetaGLM)
%   logPYCondAThetaAndDerivs_noAdep_Wrapper: (thetaVec,splits,spikeLogicals,firingRates,aMAP,featuresHistoryRepresentation)
%   logPACondThetaAndDerivs_noAdep_Wrapper:  (thetaVec,splits,aMAP,K,dK,d2K,mu,dMu,d2Mu)
%   logDetSigmaInvLambdaAndDerivs_noAdep_Wrapper:


%NOTE: PASSING THE WHOLE OBJECT AS AN ARGUMENT TO THESE FUNCTIONS.  THIS IS
%HIDEOUSLY INEFFICIENT, BUT MAY BE NECESSARY TO GET THINGS TO WORK IN THE
%DERIVATIVE TESTER.

for functionIdx = 1:numel(functionHandles)
    switch functionHandles{functionIdx}
        case '@logPThetaGPAndDerivs_Tester'
            extraArgsCell = {obj, modelNumber, splits};
        
        case '@logPThetaGLMAndDerivs_Tester'
            %Arguments needed:
            %For logPThetaGLMAndDerivs_Wrapper: (thetaVec,splits,meanThetaGLM,covThetaGLM)
            %Extra inputs are not different if we want to fit thetaGP or not.
            extraArgsCell = {obj, modelNumber, splits};
        case '@logPYCondAThetaAndDerivs_noAdep_Tester'
            %   logPYCondAThetaAndDerivs_noAdep_Wrapper: (thetaVec,splits,spikeLogicals,firingRates,aMAP,featuresHistoryRepresentation)
            if fitThetaGP
                extraArgsCell = {obj,modelNumber, sessionIdx,           splits,spikeNumbers,aMAP};
            else
                extraArgsCell = {obj,modelNumber, sessionIdx, gpHypers, splits,spikeNumbers,aMAP};
            end
        case '@logPACondThetaAndDerivs_noAdep_Tester'
            %   logPACondThetaAndDerivs_noAdep_Wrapper:  (thetaVec,splits,aMAP,K,dK,d2K,mu,dMu,d2Mu)
            if fitThetaGP
                extraArgsCell = {         splits,aMAP,                          times, clickTimes,clickSigns, doDetTermCalculations};
            else
                extraArgsCell = {gpHypers,splits,aMAP,K,dK,d2K,Kinv,mu,dMu,d2Mu, doDetTermCalculations};
            end
        case '@logDetSigmaInvLambdaAndDerivs_noAdep_Tester'
            %   logDetSigmaInvLambdaAndDerivs_noAdep_Tester:
            if fitThetaGP
                extraArgsCell = {obj,modelNumber, sessionIdx,           splits, spikeNumbers, aMAP, times, clickTimes, clickSigns};
            else
                extraArgsCell = {obj,modelNumber, sessionIdx, gpHypers, splits, spikeNumbers, aMAP, K, dK, d2K, Kinv};
            end
        case '@logDetTerms_noAdep_Tester'
            if fitThetaGP
                extraArgsCell = {obj,modelNumber, sessionIdx,           splits, spikeNumbers, aMAP, times, clickTimes, clickSigns};
            else
                extraArgsCell = {obj,modelNumber, sessionIdx, gpHypers, splits, spikeNumbers, aMAP, K, dK, d2K, Kinv};
            end
            
    end
    
    fprintf(['Beginning to test function ',functionHandles{functionIdx},' via the derivative and Hessian checker.\n\n'])
    %     switch functionHandles{functionIdx}
    %         case '@logDetSigmaInvLambdaAndDerivs_noAdep_Tester'
    %             fprintf('Testing the hessian of the logDet term');
    %             derivativeCheck(functionHandles{functionIdx},thetaGLM,2,2,extraArgsCell{:});
    %         otherwise
    if fitThetaGP
        example_derivativeCheck_mod(functionHandles{functionIdx},   [thetaGLM;gpHypers],     extraArgsCell);
    else
        example_derivativeCheck_mod(functionHandles{functionIdx},   thetaGLM,               extraArgsCell);
    end
    %     end
end
