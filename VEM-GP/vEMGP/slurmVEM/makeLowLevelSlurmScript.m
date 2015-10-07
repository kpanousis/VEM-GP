function  makeLowLevelSlurmScript(evalString, outputScriptPath, outputScriptName, varargin)
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%12 August, 2014
%
%Function which generates a SLURM script to execute a MATLAB command.

%% Template
%A typical Script should read as follows (minus the left-margin %'s).
%
%#!/bin/bash
%
%cd /nfs/nhome/live/tdesautels/slurmTesting/
%/opt/matlab-R2012b/bin/matlab -nodisplay -nosplash -singleCompThread -r "fprintf('Hello World!'); exit;"

%% Input and output handling
nExtraArgsIn = max(nargin-3,0);
if nExtraArgsIn > 0
    matlabFlags = varargin{1};
    if nExtraArgsIn > 1
        matlabLocation = varargin{2};
    else
        matlabLocation = '/opt/matlab-R2012b/bin/matlab';
    end
else
    matlabLocation = '/opt/matlab-R2012b/bin/matlab';
    matlabFlags = '-nodisplay -nosplash -singleCompThread';
end

%% Create the file
fid = fopen(fullfile(outputScriptPath,outputScriptName),'w');

%Create the script
if ischar(evalString)
    scriptCell = {'#!/bin/bash\n';'\n';['cd ',outputScriptPath,'\n'];[matlabLocation,' ',matlabFlags,' -r "',evalString,'"']};
elseif iscellstr(evalString)
    scriptCell = cell(3+numel(evalString),1);
    scriptCell{1,1} = '#!/bin/bash\n';
    scriptCell{2,1} = '\n';
    scriptCell{3,1} = ['cd ',outputScriptPath,'\n'];
    for evalStringLineNum = 1:numel(evalString)
        scriptCell{3+evalStringLineNum} = [matlabLocation,' ',matlabFlags,' -r "',evalString{evalStringLineNum},'"\n'];
    end
end
%Write the script to the file
for rowIdx = 1:size(scriptCell,1)
    fprintf(fid,scriptCell{rowIdx});
end

%Close the file
fclose(fid);

%Set the permissions to allow execution
system(['chmod 744 ',fullfile(outputScriptPath,outputScriptName)]);