function makeMidLevelSlurmScript(midLevelScriptName, midLevelScriptPath, targetScriptFullPath, timeEstimateString,qosString, varargin)
%Thomas Desautels
%Gatsby Computational Neuroscience Unit
%12 August, 2014
%
%Function which generates a SLURM script to execute a lower-level slurm script.

%Typical script is of this form (minus leading %'s):
%
%#!/bin/bash
%#SBATCH --time=0-0:03:00   %%%%Note that this is in days-hours:minutes:seconds
%#SBATCH --output job.%j.out
%#SBATCH --error job.%j.err
%
%
%srun -n 1 /nfs/nhome/live/tdesautels/slurmTesting/script1bTom


%% Input and output handling
nExtraArgsIn = max(nargin-5,0);
if nExtraArgsIn > 0
    nProcessors = varargin{1};
else
    nProcessors = 1;
end

%% Create the file
[fid]= fopen(fullfile(midLevelScriptPath,midLevelScriptName),'w');

%assert(~strcmp('~',midLevelScriptPath(1)),'The midLevelScriptPath specified is a relative path!  It must be specified as an absolute path!')
if strcmp('~',midLevelScriptPath(1))
    currDir = pwd;
    cd(midLevelScriptPath)
    midLevelScriptPath = pwd;
    cd(currDir)
end
%Create the script
scriptCell = {'#!/bin/bash\n';['#SBATCH --time=',timeEstimateString,'\n'];...
    ['#SBATCH --qos=',qosString,'\n'];...
    ['#SBATCH --output ', fullfile(midLevelScriptPath, 'job.%%j.out'),'\n'];...
    ['#SBATCH --error ',  fullfile(midLevelScriptPath, 'job.%%j.err'),'\n'];...
    '\n';...
    '\n';...
    ['srun -n ',num2str(nProcessors),' ',targetScriptFullPath]};
%Above: In a change from the previous version, the .out and .err files are
%now to be created in the folder where the midLevelScript itself is being
%created.  Note that for this to work (or for the SLURM script to run at
%all, actually) the midLevelScriptPath must be an ABSOLUTE path.



%Write the script to the file
for rowIdx = 1:size(scriptCell,1)
    fprintf(fid,scriptCell{rowIdx});
end

%Close the file
fclose(fid);

%Set the permissions to allow execution
system(['chmod 744 ',fullfile(midLevelScriptPath,midLevelScriptName)]);