classdef glmBrody < handle
    properties
        sessionIDCell  %Cell array of session number strings (nSessions x 1).
        
        cellsIDCell
        %nSessions x 1 cell array: each cell contains a cell array of
        %strings, which are the designators of the active cells in that
        %session.
        
        log
        %Cell array of increasing size which records the actions taken by
        %the object; automatic logging at the conclusion of each method.
        
        trialData
        %Cell array of structures, nSessions x 1.  Each structure (one for
        %each trial) has the following fields:
        %bup_diff: difference between number of right_bups and left_bups (
        %   negative values mean more left bups, positive values mean more
        %   right bups).
        %pokedR: Logical value saying if the animal poked right.
        %left_bups: times of left bups.
        %right_bups: Times of right bups.
        %spikes: Times of spikes in the cells indicated.
        %state_0_exits: Time of leaving state 0 (definition?) in Brody
        %   group's finite state model of the trial. IN GLOBAL TIME
        %spoke_in: time of entering the side poke chosen (?)
        %stim_start: time of stimulus start
        
        prunedTrialData
        %Cell array of structures, nSessions x 1.  Holds data we decide
        %(during acquireData) should be pruned for some reason.
        %Has the same fields as trialData, with the additional field:
        %originalTrialIdx: Index into the original order of the trials.
        %This enables re-interleaving if desired.
        %NOTE THAT THE SEGREGATION OF DATA MEANS THAT TRIALDATA IS NO
        %LONGER THE COMPLETE RECORD OF THE SESSION.
        
        modelTypes
        %Cell array of strings: 1 x nModels. Each string is the name of a
        %model type we'll fit.
        
        modelGlobalParameters
        %Cell array of structures, 1 x nModels.  Each structure includes
        %globally important information regarding the model, its
        %parameters, priors, etc.
        
        modelSessionParameters
        %Cell array of structures, nSessions x nModels.  Each structure
        %contains the parameters pertaining to this model for this session
        
        modelTrialParameters
        %Cell array of cell arrays, nSessions x nModels.  Each cell in the
        %cell array contains a structure which gives the parameters for
        %this model for this trial.
        
        modelAltDataRepresentations
        %Cell array of (possibly) cell arrays, nSessions x nModels.  Each
        %cell in the cell array contains a structure giving local versions
        %of the important information, e.g., a feature representation of
        %the spiking history.
        
        
    end
    properties (SetAccess = immutable)
        animalDesignator  %String which designates the animal's unique identifier.
    end
    
    methods
        %Class constructor method
        function obj = glmBrody(animalDesignator,varargin)
            nExtraArgsIn = max(0,nargin -1);
            if nExtraArgsIn == 2
                soonToBeSessionIDCell = varargin{1};
                soonToBeCellsIDCell = varargin{2};
                filePrefixStr = 'packdata_cell';
            elseif nExtraArgsIn == 3
                soonToBeSessionIDCell = varargin{1};
                soonToBeCellsIDCell = varargin{2};
                filePrefixStr = varargin{3};
            else
                assert(nExtraArgsIn == 0, 'Wrong number of inputs to glmBrody class creator.')
            end
            %Create the log cell array:
            obj.log = cell(0,3);
            
            %Create the key data storage cell arrays.
            obj.sessionIDCell = cell(0,1);
            obj.cellsIDCell = cell(0,1);
            obj.trialData = cell(0,1);
            obj.prunedTrialData = cell(0,1);
            obj.modelTypes = cell(1,0);
            obj.modelGlobalParameters = cell(1,0);
            obj.modelSessionParameters = cell(0,0);
            obj.modelTrialParameters = cell(0,0);
            obj.modelAltDataRepresentations = cell(0,0);
            
            
            %Initialize the animal designator.
            obj.animalDesignator = animalDesignator;
            obj.logEntry(['Created a new glmBrody object for animal ',animalDesignator,'.'],animalDesignator)
            
            if any(nExtraArgsIn == [2,3])
                %Acquire the first data to be used with this animal.
                obj.acquireData(soonToBeSessionIDCell, soonToBeCellsIDCell,filePrefixStr);
            end
        end
        
        
        % Acquire data
        function acquireData(obj,sessionIDCell, cellsIDCell,filePrefixStr)
            assert(numel(sessionIDCell) == numel(unique(sessionIDCell)),'Non-unique session ID numbers are included!')
            assert(numel(sessionIDCell) == numel(cellsIDCell),'The selected session ID cell array and cell ID cell array do not have the same number of sessions!\n')
            
            sessionErrFlag = false(size(sessionIDCell));  %only those sessions which pass will be included in the data loaded.
            allSessionsStr = cell(numel(sessionIDCell),1);
            prunedTrialDataFromSessionsStr = cell(numel(sessionIDCell),1);
            
            %For each session:
            for sessionIdx = 1:numel(sessionIDCell)
                trialCountMismatchFlag = false;
                sessionErrFlag(sessionIdx) = false;
                if numel(unique(cellsIDCell{sessionIdx})) ~= numel(cellsIDCell{sessionIdx}) %i.e., there are redundant cell listings
                    sessionErrFlag(sessionIdx) = true;
                    fprintf(['Error occurred in session ',sessionIDCell{sessionIdx},'. Check this!'])
                end
                
                
                
                %Try to load each designated cell: require all to work to
                %keep data at all.
                for cellIdx = 1:numel(cellsIDCell{sessionIdx})
                    if ~sessionErrFlag(sessionIdx)
                        %errFlag = false;
                        try
                            allData = load([filePrefixStr,cellsIDCell{sessionIdx}{cellIdx},'.mat']);
                            assert(strcmp(allData.ratname,obj.animalDesignator),'The animal this data object represents and the one in the file are not the same!')
                            assert(strcmp(sessionIDCell{sessionIdx},num2str(round(allData.sessid))), 'The session ID given and the session ID in this data file do not match!')
                            
                            %For each trial, pull out the relevant data;
                            %this should be the vec_data.bup_diff value,
                            %array_data.left_bups, array_data.spikes,
                            %array_data.right_bups.
                            if cellIdx == 1
                                
                                %Can we initialize this structure somehow?
                                
                                for trialIdx = 1:numel(allData.vec_data.bup_diff);
                                    %Initialize the common, per-trial, cross-cell values.
                                    thisSessionDataStr(trialIdx).bup_diff = allData.vec_data.bup_diff(trialIdx);
                                    thisSessionDataStr(trialIdx).pokedR = allData.vec_data.pokedR(trialIdx);
                                    thisSessionDataStr(trialIdx).state_0_exits = allData.vec_data.state_0_exits(trialIdx);
                                    thisSessionDataStr(trialIdx).cpoke_end = allData.vec_data.cpoke_end(trialIdx);  %When the center fixation cue LED turns off.
                                    thisSessionDataStr(trialIdx).cpoke_out = allData.vec_data.cpoke_out(trialIdx);  %When the animal withdraws from the center port.
                                    %thisSessionDataStr(trialIdx).spoke_in = allData.vec_data.spoke_in(trialIdx);  %When the animal pokes into the side port registered on the trial (?).
                                    thisSessionDataStr(trialIdx).stim_start = allData.vec_data.stim_start(trialIdx);  %Start of stimulus (marked by double-click)
                                    thisSessionDataStr(trialIdx).gamma = allData.vec_data.gamma(trialIdx);  %variable describing rate discrepancy
                                    
                                    thisSessionDataStr(trialIdx).left_bups = allData.array_data(trialIdx).left_bups;
                                    thisSessionDataStr(trialIdx).right_bups = allData.array_data(trialIdx).right_bups;
                                    %Save the spiking data.
                                    thisSessionDataStr(trialIdx).spikes = cell(1,numel(cellsIDCell{sessionIdx}));
                                    thisSessionDataStr(trialIdx).spikes(cellIdx) = {allData.array_data(trialIdx).spikes};  %Check how this is done in debugging.
                                end
                            else
                                for trialIdx = 1:max(numel(allData.vec_data.bup_diff),numel(thisSessionDataStr));
                                    if trialIdx <= numel(allData.vec_data.bup_diff) && trialIdx <= numel(thisSessionDataStr)
                                        %Check against the common, per-trial,
                                        %cross-cell values; Throw an error if
                                        %any of these checks fail.
                                        assert(all(thisSessionDataStr(trialIdx).bup_diff == allData.vec_data.bup_diff(trialIdx)),'Bups differences between cells do not match!')
                                        assert(all(thisSessionDataStr(trialIdx).pokedR == allData.vec_data.pokedR(trialIdx)),'Behavioral action differs between cells!')
                                        assert(all(thisSessionDataStr(trialIdx).state_0_exits == allData.vec_data.state_0_exits(trialIdx)),'State 0 exit time differs between cells!')
                                        assert(all(thisSessionDataStr(trialIdx).cpoke_end == allData.vec_data.cpoke_end(trialIdx)),'Center poke LED end (cpoke_end) time differs between cells!')
                                        assert(all(thisSessionDataStr(trialIdx).cpoke_out == allData.vec_data.cpoke_out(trialIdx)) || (isnan(thisSessionDataStr(trialIdx).cpoke_out) && isnan(allData.vec_data.cpoke_out(trialIdx))),'Center poke withdraw (cpoke_out) time differs between cells!')
                                        %assert(all(thisSessionDataStr(trialIdx).spoke_in == allData.vec_data.spoke_in(trialIdx)),'Side poke in time differs between cells!')
                                        assert(all(thisSessionDataStr(trialIdx).stim_start == allData.vec_data.stim_start(trialIdx)),'Stimulus start time differs between cells!')
                                        assert(all(thisSessionDataStr(trialIdx).gamma == allData.vec_data.gamma(trialIdx)),'Value of gamma differs between cells!')
                                        
                                        assert(all(thisSessionDataStr(trialIdx).left_bups == allData.array_data(trialIdx).left_bups),'Left bups do not match between cells!')
                                        assert(all(thisSessionDataStr(trialIdx).right_bups == allData.array_data(trialIdx).right_bups),'Eight bups do not match between cells!')
                                        
                                        %Save the spiking data from this cell.
                                        thisSessionDataStr(trialIdx).spikes(cellIdx) = {allData.array_data(trialIdx).spikes};
                                    elseif trialIdx > numel(allData.vec_data.bup_diff)
                                        thisSessionDataStr(trialIdx).spikes(cellIdx) = {};
                                        if ~trialCountMismatchFlag
                                            fprintf('SHOULD NOT HAPPEN!!!!! CHECK !!!! In session %s of animal %s, cell %s has fewer trials recorded than the others from the same session!\n',sessionIDCell{sessionIdx},obj.animalDesignator,cellsIDCell{sessionIdx}{cellIdx})
                                            trialCountMismatchFlag = true;
                                        end
                                    else  %This cell has more trials than others in the same session!
                                        %Expand the fields of the trial
                                        %data structure for the preceeding
                                        %cells.
                                        thisSessionDataStr(trialIdx).bup_diff = allData.vec_data.bup_diff(trialIdx);
                                        thisSessionDataStr(trialIdx).pokedR = allData.vec_data.pokedR(trialIdx);
                                        thisSessionDataStr(trialIdx).state_0_exits = allData.vec_data.state_0_exits(trialIdx);
                                        thisSessionDataStr(trialIdx).cpoke_end = allData.vec_data.cpoke_end(trialIdx);  %When the center fixation cue LED turns off.
                                        thisSessionDataStr(trialIdx).cpoke_out = allData.vec_data.cpoke_out(trialIdx);  %When the animal withdraws from the center port.
                                        %thisSessionDataStr(trialIdx).spoke_in = allData.vec_data.spoke_in(trialIdx);  %When the animal pokes into the side port registered on the trial (?).
                                        thisSessionDataStr(trialIdx).stim_start = allData.vec_data.stim_start(trialIdx);  %Start of stimulus (marked by double-click)
                                        thisSessionDataStr(trialIdx).gamma = allData.vec_data.gamma(trialIdx);  %variable describing rate discrepancy
                                        
                                        thisSessionDataStr(trialIdx).left_bups = allData.array_data(trialIdx).left_bups;
                                        thisSessionDataStr(trialIdx).right_bups = allData.array_data(trialIdx).right_bups;
                                        %Save the spiking data.
                                        thisSessionDataStr(trialIdx).spikes = cell(1,numel(cellsIDCell{sessionIdx}));
                                        
                                        for cellIdxFillIn = 1:(cellIdx-1)
                                            thisSessionDataStr(trialIdx).spikes(cellIdxFillIn) = {};
                                        end
                                        
                                        thisSessionDataStr(trialIdx).spikes(cellIdx) = {allData.array_data(trialIdx).spikes};  %Check how this is done in debugging.
                                        
                                        
                                        
                                        if ~trialCountMismatchFlag
                                            fprintf('SHOULD NOT HAPPEN!!!!! CHECK !!!! In session %s of animal %s, cell %s has MORE trials recorded than the others from the same session!\n',sessionIDCell{sessionIdx},obj.animalDesignator,cellsIDCell{sessionIdx}{cellIdx})
                                            trialCountMismatchFlag = true;
                                        end
                                    end
                                    
                                    
                                end
                            end
                            
                        catch errStruct
                            %errFlag = true;
                            sessionErrFlag(sessionIdx) = true;
                            fprintf(['Error occurred in session ',sessionIDCell{sessionIdx},' file packdata_cell',cellsIDCell{sessionIdx}{cellIdx},'.mat . Check this!\n'])
                            obj.logEntry(['Error occurred in acquiring session ',sessionIDCell{sessionIdx},' file packdata_cell',cellsIDCell{sessionIdx}{cellIdx},'.mat .'],errStruct)
                            if exist('thisSessionDataStr','var')
                                clear thisSessionDataStr
                            end
                        end
                    end
                    
                end
                if exist('thisSessionDataStr','var')
                    %Run a check to prune the acquired trials of any which
                    %would be troublesome
                    keepInThisSessionDataStr = true(size(thisSessionDataStr));
                    pruningIndices = [];
                    for trialIdx = 1:numel(thisSessionDataStr)
                        %Run check
                        if isnan(thisSessionDataStr(trialIdx).cpoke_out) || (~any(thisSessionDataStr(trialIdx).left_bups ~= thisSessionDataStr(trialIdx).stim_start) && ~any(thisSessionDataStr(trialIdx).right_bups ~= thisSessionDataStr(trialIdx).stim_start) )
                            %the number of left bups and right bups post-start are zero
                            keepInThisSessionDataStr(trialIdx) = false;
                            pruningIndices = horzcat(pruningIndices,trialIdx);
                            if isnan(thisSessionDataStr(trialIdx).cpoke_out)
                                fprintf('Pruning due to NaN timestamp in cpoke_out of trial %d of session %s!\n\n',trialIdx,sessionIDCell{sessionIdx})
                            end
                        end
                    end
                    
                    if any(~keepInThisSessionDataStr)
                        %Create the appropriate pruned data structure
                        prunedSessionDataStr = thisSessionDataStr(~keepInThisSessionDataStr);
                        for prunedIdx = 1:numel(pruningIndices)
                            prunedSessionDataStr(prunedIdx).originalTrialIndex = pruningIndices(prunedIdx);
                        end
                        
                        %Excise the pruned data and keep the rest.
                        thisSessionDataStr = thisSessionDataStr(keepInThisSessionDataStr);
                        fprintf(['Pruned some data from session ',sessionIDCell{sessionIdx},', check that this is satisfactory.\n\n'])
                    else
                        prunedSessionDataStr = [];
                    end
                    
                    allSessionsStr{sessionIdx} = thisSessionDataStr;
                    prunedTrialDataFromSessionsStr{sessionIdx} = prunedSessionDataStr;
                    
                end
                clear thisSessionDataStr
                %Above: CRITICAL, but hadn't done before.
                
            end
            
            %Add data to structure and log its (attempted and/or successful) acquisition.
            obj.logEntry('Attempted to acquire data from sessions and cells.',horzcat(sessionIDCell,cellsIDCell))
            if any(~sessionErrFlag)
                %Concatenate these onto the data arrays:
                obj.sessionIDCell = vertcat(obj.sessionIDCell, sessionIDCell(~sessionErrFlag) );
                obj.cellsIDCell = vertcat(obj.cellsIDCell, cellsIDCell(~sessionErrFlag));
                obj.trialData = vertcat(obj.trialData,allSessionsStr(~sessionErrFlag));  %Figure this bit out.
                obj.prunedTrialData = vertcat(obj.prunedTrialData,prunedTrialDataFromSessionsStr(~sessionErrFlag));  %Figure this bit out.
                
                obj.modelSessionParameters = vertcat(obj.modelSessionParameters,cell(sum(~sessionErrFlag),size(obj.modelSessionParameters,2)));
                obj.modelTrialParameters = vertcat(obj.modelTrialParameters,cell(sum(~sessionErrFlag),size(obj.modelTrialParameters,2)));
                obj.modelAltDataRepresentations = vertcat(obj.modelAltDataRepresentations,cell(sum(~sessionErrFlag),size(obj.modelAltDataRepresentations,2)));
                %Should we populate these somehow?
                
                obj.logEntry('Successfully acquired data from all cells in sessions:',sessionIDCell(~sessionErrFlag))
                if any(sessionErrFlag)
                    obj.logEntry('Unsuccessful loading from sessions:',sessionIDCell(sessionErrFlag))
                end
            else
                obj.logEntry('No sessions successfully loaded.')
            end
        end
        
        
        %Make a log entry:
        function logEntry(obj,message,varargin)
            nExtraArgs = max(0,nargin - 2); %Check this is correct:
            if nExtraArgs > 1
                error('Too many arguments to logEntry!')
            elseif nExtraArgs == 1
                objectAppended = varargin{1};
                %if ~iscell(cellToAdd)
                %    cellToAdd = {cellToAdd};
                %end
            else
                objectAppended = [];
            end
            
            %Concatenate this onto the log.
            obj.log = vertcat(obj.log,{datestr(now,'ddmmyyyy_HHMM'),message,objectAppended});
            
            
        end
        
        %Add a model type
        function addModel(obj,modelTypeStr)
            switch modelTypeStr
                case 'glmBrodyLaplace'
                    %Create a glmBrody type model which uses Laplace
                    %approximations.
                    
                case 'glmBrodyLaplaceSelfOnly'
                    %Create a glmBrody type model which uses Laplace
                    %approximations, with the restriction that the glm
                    %history weights be zero from the other cells' history.
                    
                case 'glmBrodyFixedSelfOnly'
                    %Create a glmBrody type model which uses a fixed latent
                    %equal to the GP mean, where glm cross-weights are
                    %fixed to be zero.
                    
                case 'glmBrodyVariational'
                    %Create a glmBrody type model which uses a global
                    %variational approximation to the latent.
                    
                case 'glmRecurrenceOnly'
                    %Create a model of spiking behavior on a
                    %session-by-session basis, explained only by the mutual
                    %firing activity of the visible units.
                    
                case 'glmWithBupHistory'
                    %Create a GLM model of spiking behavior which uses both
                    %the firing of the visible units and the stimulus
                    %history to explain the firing activity.
                    
                otherwise
                    error('Model type string not recognized; none added.')
            end
            
            %Expand and initialize the parameters of the model
            obj.modelTypes = horzcat(obj.modelTypes, {modelTypeStr});
            obj.modelGlobalParameters = horzcat(obj.modelGlobalParameters,{});
            obj.modelSessionParameters = horzcat(obj.modelSessionParameters,cell(size(obj.modelSessionParameters,1),1));
            obj.modelTrialParameters = horzcat(obj.modelTrialParameters,cell(size(obj.modelSessionParameters,1),1));
            obj.modelAltDataRepresentations = horzcat(obj.modelAltDataRepresentations,cell(size(obj.modelSessionParameters,1),1));
            obj.logEntry(['Added model %d, of type %s.',numel(obj.modelTypes),modelTypeStr])
        end
        
        function initializeParameters(obj,modelNumber,globalParametersSpecification,varargin)
            %initialize any parameters used by the model
            switch obj.modelTypes{modelNumber}
                case 'glmBrodyLaplace'
                    
                case 'glmBrodyVariational'
                    
                case 'glmRecurrenceOnly'
                    %Initialize parameters for a recurrent-only GLM model
                    
                case 'glmWithBupHistory'
                    %Initialize parameters for a GLM model which uses
                    %recurrent connections between the visible units and
                    %dependence on the history.
                    
                case 'glmBrodyFixedSelfOnly'
                    %Input checking
                    assert(any(nargin-3 == [0,1]),'Incorrect number of input arguments!')
                    
                    %Initialize the global parameters
                    obj.modelGlobalParameters{1,modelNumber} = globalParametersSpecification;
                    if nargin-3 == 1
                        sessionParametersSpecification = varargin{1};
                        error('Initialization of session parameters not yet implemented!')
                        %for sessionIdx = 1:size(obj.trialData,1)
                        %    obj.modelSessionParameters{sessionIdx,modelNumber} = ;
                        %end
                    end
                    
                    
                    
                otherwise
                    error('Model type string not recognized; no model parameters initialized.')
            end
            obj.logEntry(['Initialized parameters of model number ',num2str(modelNumber),' which is of type ',obj.modelTypes{modelNumber},'.'])
        end
        
        function createBasis(obj, modelNumber, basisTypeString, selfOrStim, basisParametersStructure)
            assert(numel(obj.modelTypes)>= modelNumber, 'You have not yet created this model number!')
            
            %Create an appropriate basis
            switch basisTypeString
                case 'hybridPoissonGaussian'
                    %This basis divides up the region in the recent past
                    %into poisson-shaped components (the most-recent) and
                    %Gaussian-bump shaped components (after the Poisson
                    %bumps begin to become Gaussian in shape, we hope).
                    
                    %sampleRate = obj.modelGlobalParameters(modelNumber).sampleRate; %Set the sample rate to be used in computing the basis
                    switch selfOrStim
                        case 'self'
                            relevantSamples = obj.modelGlobalParameters{modelNumber}.relevantSelfHistory * obj.modelGlobalParameters{modelNumber}.samplingRate;
                        case 'stim'
                            relevantSamples = obj.modelGlobalParameters{modelNumber}.relevantStimHistory * obj.modelGlobalParameters{modelNumber}.samplingRate;
                    end
                    nPoissonCompsBasis = basisParametersStructure.nPoissonCompsBasis;
                    assert(nPoissonCompsBasis >=1, 'At least one Poisson bump must be present for the hybridPoissonGaussian basis choice.')
                    
                    %Place the centers evenly
                    if nPoissonCompsBasis > 1
                        poissonCenters = [0,1,2:2:2*(nPoissonCompsBasis-2)];
                    elseif nPoissonCompsBasis == 1
                        poissonCenters = 0;
                    end
                    
                    %Match the width of the Gaussians to the last poisson
                    sigmaGaussians = sqrt(poissonCenters(end));
                    
                    if exist('basisParametersStructure','var') && isfield(basisParametersStructure,'nGaussianCompsBasis')
                        nGaussianCompsBasis = basisParametersStructure.nGaussianCompsBasis;
                        fprintf('Prespecifying the number of Gaussian pseudo-basis components may lead to poor spacing/coverage!\n')
                        gaussianCenters = poissonCenters(end) + (1:nGaussianCompsBasis) * ((relevantSamples-1) - poissonCenters(end))/nGaussianCompsBasis;
                    else
                        nGaussianCompsBasis = floor((relevantSamples - poissonCenters(end))/(2*sigmaGaussians));
                        gaussianCenters = poissonCenters(end) + (1:nGaussianCompsBasis) * (2*sigmaGaussians);
                    end
                    
                    nCompsPseudoBasis = nPoissonCompsBasis + nGaussianCompsBasis;
                    centers = [poissonCenters,gaussianCenters];
                    
                    
                    pseudoBasis = zeros(relevantSamples,nCompsPseudoBasis);
                    for idx = 1:nCompsPseudoBasis
                        if idx <= nPoissonCompsBasis
                            pseudoBasis(:,idx) = poisspdf(0:relevantSamples-1,centers(idx));
                        else
                            pseudoBasis(:,idx) = normpdf(0:relevantSamples-1,centers(idx),sigmaGaussians);    %Tiles remainder of interval with normal PDFs
                            renormalizationConst = 1/sum(normpdf(centers(idx) - 100 : centers(idx) + 100, centers(idx), sigmaGaussians));
                            pseudoBasis(:,idx) = renormalizationConst * pseudoBasis(:,idx);
                        end
                    end
                    
                    
                case 'gaussian'
                    switch selfOrStim
                        case 'self'
                            relevantSamples = obj.modelGlobalParameters{modelNumber}.relevantSelfHistory * obj.modelGlobalParameters{modelNumber}.samplingRate;
                        case 'stim'
                            relevantSamples = obj.modelGlobalParameters{modelNumber}.relevantStimHistory * obj.modelGlobalParameters{modelNumber}.samplingRate;
                    end
                    nGaussianCompsBasis = basisParametersStructure.nGaussianCompsBasis;
                    nCompsPseudoBasis = nGaussianCompsBasis;
                    sigmaGaussians = 1/2 * (relevantSamples-1)/(nGaussianCompsBasis -1);
                    centers = 1:(2*sigmaGaussians):relevantSamples;
                    
                    pseudoBasis = zeros(relevantSamples,nCompsPseudoBasis);
                    for idx = 1:nGaussianCompsBasis;
                        pseudoBasis(:,idx) = normpdf(1:relevantSamples,centers(idx),sigmaGaussians);    %Tiles remainder of interval with normal PDFs
                        renormalizationConst = 1/sum(normpdf(centers(idx) - 100 : centers(idx) + 100, centers(idx), sigmaGaussians));
                        pseudoBasis(:,idx) = renormalizationConst * pseudoBasis(:,idx);
                    end
                otherwise
                    error('Basis type string not recognized; none created.')
            end
            
            pseudoBasis = flipud(pseudoBasis);
            
            %Note: The current basis is parameterized in time such that the
            %most recent components are at the bottom rows, and oldest
            %timestamps are at the top.  The basis functions are also
            %ordered such that the most recent are in the left columns and
            %the oldest are in the right columns.
            
            
            %Now, store the created basis in the appropriate entry
            switch selfOrStim
                case 'self'
                    fprintf(['For model ', num2str(modelNumber),', setting the newly-created pseudo-basis to parameterize recurrence among visible units.\n'])
                    obj.modelGlobalParameters{modelNumber}.pseudoBasisSelf = pseudoBasis;
                    obj.logEntry('Self pseudo-basis set to:',pseudoBasis)
                case 'stim'
                    fprintf(['For model ', num2str(modelNumber),', setting the newly-created pseudo-basis to parameterize the stimulus influence.\n'])
                    obj.modelGlobalParameters{modelNumber}.pseudoBasisStim = pseudoBasis;
                    obj.logEntry('Stimulus pseudo-basis set to:',pseudoBasis)
            end
            
            
        end
        
        function varargout = fitModel(obj,modelNumber,fittingParameters)
            nExtraArgsOut = nargout;
            assert(modelNumber <= numel(obj.modelTypes))
            switch obj.modelTypes{modelNumber}
                case 'glmBrodyLaplace'
                    %Perform the fitting procedure
                    if isfield(fittingParameters,'fittingMethod')
                        if strcmp(fittingParameters.fittingMethod,'fit_glmBrodyLaplace_GLMandGP')
                            fittingStr = 'y';
                        elseif strcmp(fittingParameters.fittingMethod,'fit_glmBrodyLaplace_GLMonly')
                            fittingStr = 'n';
                        else
                            error('Unrecognized fittingMethod field of fittingParameters input!')
                        end
                        fittingParameters = rmfield(fittingParameters,'fittingMethod');
                    else
                        fittingStr = input('Do you want to fit all parameters, or just the GLM weights?\n\n','s');
                    end
                    
                    switch fittingStr
                        case 'y'
                            %error('Model not yet created!')
                            [deltaThetas] = obj.fit_glmBrodyLaplace_GLMandGP(modelNumber,fittingParameters);
                            disp(deltaThetas)
                            if nExtraArgsOut == 1
                                varargout{1} = deltaThetas;
                            end
                        case 'n'
                            [deltaThetas] = obj.fit_glmBrodyLaplace_GLMonly(modelNumber,fittingParameters);
                            disp(deltaThetas)
                            if nExtraArgsOut == 1
                                varargout{1} = deltaThetas;
                            end
                    end
                    
                case 'glmBrodyLaplaceSelfOnly'
                    %Perform the fitting procedure
                    %%Old Method:
                    %switch input('Do you want to fit all parameters, or just the GLM weights?\n\n','s')
                    
                    if isfield(fittingParameters,'fittingMethod')
                        if strcmp(fittingParameters.fittingMethod,'fit_glmBrodyLaplace_GLMandGP')
                            fittingStr = 'y';
                        elseif strcmp(fittingParameters.fittingMethod,'fit_glmBrodyLaplace_GLMonly')
                            fittingStr = 'n';
                        else
                            error('Unrecognized fittingMethod field of fittingParameters input!')
                        end
                        fittingParameters = rmfield(fittingParameters,'fittingMethod');
                    else
                        %New method, using waiting and timing.
                        %Core code follows
                        %http://wiki.stdout.org/matlabcookbook/Collecting%20data/Getting%20time-sensitive%20keyboard%20input/,
                        %example #2, accessed 23 June, 2014.
                        if ~strcmp(which('KbCheck'),'')  %The psychtoolbox is installed and on the path
                            fprintf('Do you want to fit all parameters (y) or just the GLM weights (n)? Waiting 10 seconds will also fit all weights.\n\n')
                            
                            yCode = KbName('y');
                            nCode = KbName('n');
                            
                            startTimestamp = GetSecs();
                            while true
                                %Listen at the keyboard
                                [~, ~, keycodes] = KbCheck();
                                
                                %Check if y or n has been pushed, or if we've
                                %exceeded 10 seconds since the start of the
                                %process.
                                if keycodes(yCode) || GetSecs() - startTimestamp > 10
                                    fprintf('Fitting using all hyperparameters!\n')
                                    fittingStr = 'y';
                                    break
                                elseif keycodes(nCode)
                                    fprintf('Fitting using only the GLM hyperparameters!\n')
                                    fittingStr = 'n';
                                    break
                                end
                            end
                        else
                            fprintf('Unable to query user for parameters to fit; defaulting to fitting GP and GLM parameters.\n\n')
                            fittingStr = 'y';
                        end
                    end
                    
                    
                    %Now use the switching structure, as before.
                    switch fittingStr
                        case 'y'
                            [deltaThetas,logLikelihood,unnormalizedLogPosteriorOverTheta] = obj.fit_glmBrodyLaplace_GLMandGP(modelNumber,fittingParameters);
                            fprintf('The deltaThetas, logLikelihood, and unnormalized logPosteriorOverTheta (each is a row) were:\n\n')
                            disp([deltaThetas;logLikelihood;unnormalizedLogPosteriorOverTheta])
                            %disp(logLikelihood)
                            %disp(unnormalizedLogPosteriorOverTheta)
                            obj.logEntry('Obtained the following delta theta and objective function values for the optimization (deltaTheta; logLikelihood; unNLogPosterior):',[deltaThetas;logLikelihood;unnormalizedLogPosteriorOverTheta])
                            if nExtraArgsOut == 3
                                varargout{1} = deltaThetas;
                                varargout{2} = logLikelihood;
                                varargout{3} = unnormalizedLogPosteriorOverTheta;
                            end
                        case 'n'
                            [deltaThetas,logLikelihood,unnormalizedLogPosteriorOverTheta] = obj.fit_glmBrodyLaplace_GLMonly(modelNumber,fittingParameters);
                            fprintf('The deltaThetas were:\n\n')
                            disp(deltaThetas)
                            obj.logEntry('Obtained the following delta theta values for the optimization:',deltaThetas)
                            if nExtraArgsOut == 1
                                varargout{1} = deltaThetas;
                            elseif nExtraArgsOut == 3
                                varargout{1} = deltaThetas;
                                varargout{2} = logLikelihood;
                                varargout{3} = unnormalizedLogPosteriorOverTheta;
                            end
                    end
                case 'glmBrodyFixedSelfOnly' 
                    %Perform the fitting procedure
                    [logLikelihood, unnormalizedLogPosteriorOverTheta] = obj.fit_glmBrodyFixed(modelNumber,fittingParameters);
                    %The resulting fitting procedure will recover some
                    %optimized values for the GLM parameters.
                    fprintf('The logLikelihood (log(p(y|a,thetaGP,thetaGLM,B)) )  and unnormalized log posterior ( log(p(y|a,thetaGP,thetaGLM,B)) + log(p(thetaGLM)) ) were:\n\n')
                    disp([logLikelihood;unnormalizedLogPosteriorOverTheta])
                    obj.logEntry('Obtained the following objective function values for the optimization (logLikelihood, unNLogPosterior):',[logLikelihood;unnormalizedLogPosteriorOverTheta])
                    if nExtraArgsOut == 2
                        varargout{1} = logLikelihood;
                        varargout{2} = unnormalizedLogPosteriorOverTheta;
                    end
                    
                    
                case 'glmBrodyVariational'
                    %Perform the fitting procedure
                    %obj.fit_glmBrodyVariational(modelNumber,fittingParameters)
                    error('Model not yet created!')
                    
                case 'glmRecurrenceOnly'
                    %Perform the fitting procedure
                    obj.fit_glmRecurrenceOnly(modelNumber,fittingParameters)
                    
                case 'glmWithBupHistory'
                    %Perform the fitting procedure
                    obj.fit_glmWithBupHistory(modelNumber,fittingParameters)
                    
                otherwise
                    fprintf('Model type string not recognized; none fit.')
            end
            obj.logEntry(['Fit parameters of model number ',num2str(modelNumber),' which is of type ',obj.modelTypes{modelNumber},'.'])
        end
    end
end