function [curr_params, sessionDataOut, varBound, EStepTimes, MStepTimes ] = vEM_GP(params,sessionDataIn,nParams)
%VEM_GP Main File for the implementation of the modified variational EM
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%8 June 2015

%Inputs: SessionDataIn: structure that contains the parameters C and d,
%        for each session along with the seq struct
%              seq: Containing fields
%              seq.y: NxT for a spike raster from N neurons
%              seq.T: the 2nd dimension of seq.Y
%        params:  initialized by pop_spike_dyn and added fields for gp
%        Parameters
%        nParams: number of params to optimize

%Outputs: curr_params: the updated struct
%         sessionDataOut: the updated struct
%         curr_params: the updated parameters
%         varBound: the variational Bound
%         EStepTimes: the times of the e step
%         MStepTimes: the times of the m step

%% Initialize Parameters
maxIter=params.opts.algorithmic.EMIterations.maxIter;
progTolvarBound = params.opts.algorithmic.EMIterations.progTolvarBound;
maxCPUTime      = params.opts.algorithmic.EMIterations.maxCPUTime;

ParamPenalizerHandle = params.model.ParamPenalizerHandle;
InferenceMethod      = params.model.inferenceHandle;

EStepTimes      = nan(maxIter,1);
MStepTimes      = nan(maxIter+1,1);
varBoundSessions        = nan(maxIter,numel(sessionDataIn));
varBound = nan(maxIter,1);

sessionDataOut=sessionDataIn;
previous_params=params;
curr_params=params;

varBoundMax=-inf;

begin_time=cputime;

% get the total number of trials through all the sessions
T_sum=0;
for j=1:numel(sessionDataIn)
    T_sum=T_sum+sum([sessionDataIn(j).seq.T]);
end

%% EM LOOP
for iter=1:maxIter
    
    curr_params.state.EMiter=iter;
    fprintf('\nEM Current Iteration: %d\n',iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%E-STEP%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %record time
    eTimeBegin=cputime;
    curr_params.opt.EMiter=iter;
    
    varBound(iter)=0;
    
    %Perform inference for each session
    %The inference method is in VariationalInferenceDualLDS function, that
    %is originally of popSpike_dyn Package, but some modifications were
    %performed to recalculate and transform some stuff
    for i=1:numel(sessionDataOut)
        curr_params.model.C=sessionDataOut(i).C;
        curr_params.model.d=sessionDataOut(i).d;
        [sessionDataOut(i).seq , varBoundSessions(iter,i)]=InferenceMethod(curr_params,sessionDataOut(i).seq);
        varBound(iter)=varBound(iter)+varBoundSessions(iter,i);
    end
    
    %end of E step, record time again
    eTimeEnd=cputime;
    EStepTimes(iter)=eTimeEnd-eTimeBegin;
    
    %add regularizer cost to varBound (see how that works because we
    %calculate different things, unless it depends on sigma and mu
    varBound(iter)=varBound(iter)-ParamPenalizerHandle(curr_params);
    varBound(iter)
    
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%TERMINATION CRITERIA%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if params.opts.algorithmic.EMIterations.abortDecresingVarBound && (varBound(iter)<varBoundMax)    % check if varBound is increasing!
        curr_params = previous_params;	   % parameter backtracking
        fprintf('\n ');
        warning('Variational lower bound in vEM is decreasing, aborting EM & backtracking');
        break;
    end
    
    if params.opts.algorithmic.EMIterations.abortDecresingVarBound && ((abs(varBound(iter)-varBoundMax)/T_sum)<progTolvarBound)
        fprintf('\nReached progTolvarBound for EM, aborting');
        break;
    end
    
    if (cputime-begin_time)>maxCPUTime
        fprintf('\nReached maxCPUTime for EM, aborting');
        break;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%Assign Params%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    previous_params=curr_params;
    varBoundMax=varBound(iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%M-STEP%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    mStepTimeBegin=cputime;
    [curr_params,sessionDataOut]=variationalMStep(curr_params,sessionDataOut,nParams);
    mStepTimeEnd=cputime;
    MStepTimes(iter)=mStepTimeEnd-mStepTimeBegin;
    
    
%%%End Iteration%%%    
end

%% End/ Finish Message
fprintf('\n----------------------------------------------------------------------------\n')
disp('Variational EM Fnished!');

end