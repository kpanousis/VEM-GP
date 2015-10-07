function [ F_GP ] = F_GP_calc(hypers,expect,condVar )
%F_GP_CALC Calculate the free energy function for the Gaussian Process term
%for the PLDS Model
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%15 June 2015

%% Input Checking
if (numel(hypers)~=6)
    error('Something wrong with the number of hyperparameters');
end


% F_GP=0;
% for t=1:numel(times)
%     F_GP=F_GP+log(condVar(t))+expect(t)/condVar(t);
% end
% F_GP=-(0.5)*F_GP;


end

