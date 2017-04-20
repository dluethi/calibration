function errpred=ctrlpred(metamodel,parameters,datamatrix)

% Evaluate if the model data on which the metamodel is estimated is
% able to reproduce the data with zero error
% NAME 
%   ctrlpred
% PURPOSE 
%   Predict modeldata on which the metamodel has been estimated to.
% INPUTS 
%   The structure metamodel, parameters and datamatrix 
%   are used for the inpute of neelin_p. Addionally the parameter
%   matrix is read 
%
%   parameters.experiments:
%
%           Parameter matrix of experiments on which metamodel is
%           estimated
%
% OUTUTS 
%   Plot: Histogram of the difference between predicted and actual
%   model data 
% HISTORY 
% First version: 11.10.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)



%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

N=length(parameters)
ds=2*N+N*(N-1)/2; % Number of fitted simulations in determined design
pmatrix=parameters(1).experiments(1:ds,:)

%--------------------------------------------------------------------
% COMPUTE Data for each parameter experiment
%--------------------------------------------------------------------

% Help variable to select all matrix dimensions if no score used
dd=ndims(datamatrix.refdata);
if dd>2
  for i=1:dd
    indd{i}=':';
  end
else
  indd=':';
end

ctrl=NaN([size(datamatrix.refdata),length(pmatrix)]); % Allocate data

for i=1:length(pmatrix)
  ctrl(indd{:},i)=neelin_p(metamodel,parameters,datamatrix,pmatrix(i,:));
end
    	
errpred=ctrl-datamatrix.moddata(indd{:},1:length(pmatrix));

figure;
hist(errpred(:))
set(gca,'Fontsize',12)
ylabel('Number of counts','Fontsize',14)
xlabel('Error of prediction','Fontsize',14)
title('Verification on design data','Fontsize',14)




