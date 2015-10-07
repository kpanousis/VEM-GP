function makeTopLevelSlurmScript(highLevelScriptName, highLevelScriptPath, midLevelSlurmScriptPathsCell,varargin)
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%12 August, 2014
%
%Function which generates a simple, top-level script to execute a set of
%mid-level slurm scripts.

%% Input handling
if ischar(midLevelSlurmScriptPathsCell)
    midLevelSlurmScriptPathsCell = {midLevelSlurmScriptPathsCell};
else
    midLevelSlurmScriptPathsCell = reshape(midLevelSlurmScriptPathsCell,[],1);
end

nExtraArgsIn = max(0,nargin-3);
assert(nExtraArgsIn <= 2,'Too many inputs to makeTopLevelSlurmScript!')
if nExtraArgsIn == 1
    sleepLength = varargin{1};
    if ~ischar(sleepLength)
        sleepLength = sprintf('%0.4f',sleepLength);
    end
elseif nExtraArgsIn == 2
    if ~isempty(varargin{1})
        assert(iscell(varargin{1}),'Input intended to be prependCell is not a cell array!')
        prependCell = varargin{1};
    end
    if ~isempty(varargin{2})
        assert(iscell(varargin{2}),'Input intended to be appendCell is not a cell array!')
        appendCell = varargin{2};
    end
end
        

%% Generate the text of the script

if exist('sleepLength','var')
    sizeScriptCell = 1 + max(0,(2 * numel(midLevelSlurmScriptPathsCell) - 1));
else
    sizeScriptCell = 1 + numel(midLevelSlurmScriptPathsCell);
end

scriptCell = cell(sizeScriptCell,1);
scriptCell{1} = '#!/bin/bash\n';

rowIdxInMidLevelPaths = 1;
for rowIdx = 2:sizeScriptCell
    if ~exist('sleepLength','var') || mod(rowIdx,2) == 0  %so if we're either not sleeping between submissions or this is an even-numbered row
        %Print the next mid-level slurm script to run.
        scriptCell{rowIdx} = ['sbatch ',midLevelSlurmScriptPathsCell{rowIdxInMidLevelPaths},' &\n'];
        rowIdxInMidLevelPaths = rowIdxInMidLevelPaths + 1;
    else %sleepLength exists and this is an odd numbered row
        scriptCell{rowIdx} = ['sleep ',sleepLength,'\n'];  %Note: think it's really important to have no & here, such that this DOESN'T run in the background.
    end
end

%Append and or prepend
if exist('prependCell', 'var')
    scriptCell = [scriptCell(1);prependCell;scriptCell(2:end)];  %Since we need the '#!/bin/bash\n' line to be the first.
end

if exist('appendCell', 'var')
    scriptCell = [scriptCell;appendCell];
end



% Create the file:
fid = fopen(fullfile(highLevelScriptPath,highLevelScriptName),'w');

%Write the script to the file
for rowIdx = 1:size(scriptCell,1)
    fprintf(fid,scriptCell{rowIdx});
end

%Close the file
fclose(fid);

%Set the permissions to allow execution
system(['chmod 744 ',fullfile(highLevelScriptPath,highLevelScriptName)]);