%varBound_ExractionFromRuns.m
%Simple script that loads the varBounds for each slurm run in a folder and
%then plots the varBound for each run (in one plot) and the varBound of
%VEM-GP against the VEM and lEM of Macke et al.
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%27 August 2015
clc;
clear all;
close all;


run('set_path.m');
addpath 'matlab2tikz/';
runFolder='vEMGP/slurmVEM/Runs/';
current_slurmRun='Run_s5_2000tr_10ms_Neurons3';

%get the number of subfolders
d=dir(fullfile(runFolder,current_slurmRun));
numRuns=numel(d)-2;

%Initialize struct for the bounds
varBound(numRuns).vEM=[];
varBound(numRuns).PopSpikeEM=[];
varBound(numRuns).PopSpikeEM_LaPlace=[];

%Runs iterator
max_Bound_vEMGP=-inf;
max_Bound_PopSpikevEM=-inf;
max_Bound_PopSpikelEM=-inf;
index_maxBound_vEMGP=-1;
index_maxBound_PopSpikevEM=-1;
index_maxBound_PopSpikelEM=-1;

for run=1:numRuns
   currentRun=['vEMGP_Run_',num2str(run)];
   structDir=what(fullfile(runFolder,current_slurmRun,currentRun));
   fileToLoad=fullfile(runFolder,current_slurmRun,currentRun,char(structDir.mat));
   varBound(run).vEM=load(fileToLoad,'varBound_vEM'); 
   varBound(run).PopSpikeEM=load(fileToLoad,'varBound_PopSpikeEM'); 
   varBound(run).PopSpikeEM_LaPlace=load(fileToLoad,'varBound_PopSpikeEM_LaPlace'); 
  
   %Find the max and the run that it occured
   if (max(varBound(run).vEM.varBound_vEM)>max_Bound_vEMGP)
       max_Bound_vEMGP=max(varBound(run).vEM.varBound_vEM);
       index_maxBound_vEMGP=run;
   end
   if (max(varBound(run).PopSpikeEM.varBound_PopSpikeEM)>max_Bound_PopSpikevEM)
       max_Bound_PopSpikevEM=max(varBound(run).PopSpikeEM.varBound_PopSpikeEM);
       index_maxBound_PopSpikevEM=run;
   end
   if (max(varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace)>max_Bound_PopSpikelEM)
       max_Bound_PopSpikelEM=max(varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace);
       index_maxBound_PopSpikelEM=run;
   end
   

end

%Get min and max for each iteration for error bars
minVarBound=nan(numel(varBound(1).vEM.varBound_vEM),1);
maxVarBound=nan(numel(varBound(1).vEM.varBound_vEM),1);

minVarBoundvEM=nan(numel(varBound(1).vEM.varBound_vEM),1);
maxVarBoundvEM=nan(numel(varBound(1).vEM.varBound_vEM),1);

minVarBoundlEM=nan(numel(varBound(1).vEM.varBound_vEM),1);
maxVarBoundlEM=nan(numel(varBound(1).vEM.varBound_vEM),1);


for index=1:numel(varBound(1).vEM.varBound_vEM)
    
    for run=1:numRuns
        
        %%Find min and max for our Bound
        if isnan(varBound(run).vEM.varBound_vEM(index))
            continue;
        end
        
        %find min bound amongst runs
        if isnan(minVarBound(index,1))
            minVarBound(index)=varBound(run).vEM.varBound_vEM(index);
        elseif (minVarBound(index)>varBound(run).vEM.varBound_vEM(index))
                minVarBound(index)=varBound(run).vEM.varBound_vEM(index);
        end
        
        %find max bound
        if isnan(maxVarBound(index,1))
            maxVarBound(index)=varBound(run).vEM.varBound_vEM(index);
        elseif (maxVarBound(index)<varBound(run).vEM.varBound_vEM(index))
                maxVarBound(index)=varBound(run).vEM.varBound_vEM(index);
        end
        
        
        %%Find min and max for VEM POP SPIKE
        if isnan(varBound(run).PopSpikeEM.varBound_PopSpikeEM(index))
            continue;
        end
        
        %find min bound amongst runs
        if isnan(minVarBoundvEM(index,1))
            minVarBoundvEM(index)=varBound(run).PopSpikeEM.varBound_PopSpikeEM(index);
        elseif (minVarBoundvEM(index)>varBound(run).PopSpikeEM.varBound_PopSpikeEM(index))
                minVarBoundvEM(index)=varBound(run).PopSpikeEM.varBound_PopSpikeEM(index);
        end
        
        %find max bound
        if isnan(maxVarBoundvEM(index,1))
            maxVarBoundvEM(index)=varBound(run).PopSpikeEM.varBound_PopSpikeEM(index);
        elseif (maxVarBoundvEM(index)<varBound(run).PopSpikeEM.varBound_PopSpikeEM(index))
                maxVarBoundvEM(index)=varBound(run).PopSpikeEM.varBound_PopSpikeEM(index);
        end
        
        %%Find min and max LAPLACE 
        if isnan(varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index))
            continue;
        end
        
        %find min bound amongst runs
        if isnan(minVarBoundlEM(index,1))
            minVarBoundlEM(index)=varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index);
        elseif (minVarBoundlEM(index)>varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index))
                minVarBoundlEM(index)=varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index);
        end
        
        %find max bound
        if isnan(maxVarBoundlEM(index,1))
            maxVarBoundlEM(index)=varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index);
        elseif (maxVarBoundlEM(index)<varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index))
                maxVarBoundlEM(index)=varBound(run).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(index);
        end
        
        
    end
    
    
end

%find the upper and lower bounds for the bounds
Upper=nan(numel(maxVarBound),1);
Lower=nan(numel(minVarBound),1);

UppervEM=nan(numel(maxVarBound),1);
LowervEM=nan(numel(minVarBound),1);

UpperlEM=nan(numel(maxVarBound),1);
LowerlEM=nan(numel(minVarBound),1);

numRunToPresent=2;
numRunPopSpike=5;
for i=1:numel(maxVarBound)
    if ~isnan(varBound(numRunToPresent).vEM.varBound_vEM(i))
        if (mod(i,5)==0)
            continue;
        else
        Upper(i,1)=maxVarBound(i)-varBound(numRunToPresent).vEM.varBound_vEM(i);
        Lower(i,1)=varBound(numRunToPresent).vEM.varBound_vEM(i)-minVarBound(i);
        end
    end
       
    %POP SPIKE vEM UPPER LOWER
     if ~isnan(varBound(numRunPopSpike).PopSpikeEM.varBound_PopSpikeEM(i))
        if (mod(i,5)==0)
            continue;
        else
        UppervEM(i,1)=maxVarBoundvEM(i)-varBound(numRunPopSpike).PopSpikeEM.varBound_PopSpikeEM(i);
        LowervEM(i,1)=varBound(numRunPopSpike).PopSpikeEM.varBound_PopSpikeEM(i)-minVarBoundvEM(i);
        end
     end
     
     %POP SPIKE lEM UPPER LOWER
     if ~isnan(varBound(numRunPopSpike).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(i))
        if (mod(i,5)==0)
            continue;
        else
        UpperlEM(i,1)=maxVarBoundlEM(i)-varBound(numRunPopSpike).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(i);
        LowerlEM(i,1)=varBound(numRunPopSpike).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(i)-minVarBoundlEM(i);
        end
    end
    
end

x=1:100;
h=figure('Color',[0.8 0.8 0.8]);
errorbar(x,varBound(numRunToPresent).vEM.varBound_vEM(1:100),Lower(1:100),Upper(1:100),'--', 'color', [0.3 0.3 0.3],'LineWidth',0.8);
hold on;
errorbar(x,varBound(numRunToPresent).PopSpikeEM.varBound_PopSpikeEM(1:100),LowervEM(1:100),UppervEM(1:100),'--', 'color', [0 34/255 139/255],'LineWidth',0.8);
errorbar(x,varBound(numRunToPresent).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace(1:100),LowerlEM(1:100),UpperlEM(1:100),'--', 'color', [0.5 0 0.5],'LineWidth',0.8);
hold off;
xlabel('Iterations','Interpreter','latex');
ylabel('Variational Bound','Interpreter','latex');
ylim([min(minVarBoundlEM)-max(LowerlEM), max(varBound(index_maxBound_vEMGP).vEM.varBound_vEM)+1000]);
xlim([0 30])
set(gca,'xscale','log');
set(gca,'yscale','log');
 box off
 grid off
color = get(h,'Color');
set(gca,'TickDir','out')
leg=legend('VEM-GP','vEM','lEM','location','SouthEast');
legend('boxoff');
title(['\indent Methods Comparison',char(10),'Sessions: 1, Trials: 2000, Neurons: 1',char(10),'6 Hyperparameters Optimization'],'Interpreter','latex')
%print -depsc graphs/s1_tr2000_n1_6_hypers
% matlab2tikz();

%plot(x,varBound(index_maxBound_vEMGP).vEM.varBound_vEM(1:200),'--k','LineWidth',1);

%Plot our method against the other in the best runs for each
% figure();
% semilogy(varBound(index_maxBound_vEMGP).vEM.varBound_vEM,'--b','LineWidth',1);
% xlabel('Number of Iterations');
% ylabel('Variational Bound');
% title(sprintf('VEM-GP/vEM/lEM Comparison\nSessions: 1 Trials: 2000 Neuron: 1'));
% hold on;
% semilogy(varBound(index_maxBound_PopSpikevEM).PopSpikeEM.varBound_PopSpikeEM,'--k','LineWidth',1);
% semilogy(varBound(index_maxBound_PopSpikelEM).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace,'--g','LineWidth',1);
% legend('VEM-GP','vEM','lEM','location','East');
% hold off;
% 
%  figure();
% %  subplot(3,1,1)
%  loglog(varBound(1).vEM.varBound_vEM(1:100));
%  xlabel('Number of vEM Iterations');
%  ylabel('Variational Bound');
%  title(sprintf('Variational Bound: VEM-GP\nSessions: 1 Trials: 2000 Neuron: 1'));
%  hold on;
%  for i=2:numRuns
%      loglog(varBound(i).vEM.varBound_vEM(1:100));
%  end
%  hold off;
% %  matlab2tikz
%  
%  figure();
% %  subplot(3,1,2)
%  semilogy(varBound(1).PopSpikeEM.varBound_PopSpikeEM,':');
%  xlabel('Number of PopSpike\_vEM Iterations');
%  ylabel('Variational Bound');
%  title(sprintf('Variational Bound: vEM\nSessions: 1 Trials: 2000 Neuron: 1'));
%  hold on;
%  for i=2:numRuns
%      semilogy(varBound(i).PopSpikeEM.varBound_PopSpikeEM,':');
%  end
%  hold off;
%  matlab2tikz
 
%  figure()
% %  subplot(3,1,3)
%  plot(varBound(1).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace,':');
%  xlabel('Number of PopSpike_lEM Iterations');
%  title('Variational Bound/Number of Iterations');
%  hold on;
%  for i=2:numRuns
%      plot(varBound(i).PopSpikeEM_LaPlace.varBound_PopSpikeEM_LaPlace);
%  end
%  hold off;
 
 