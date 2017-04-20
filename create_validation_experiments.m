function [xexp]=create_validation_experiments(parameters,lhacc)

% Optimise model parameters using a latin hypercube sampling 
% NAME 
%   create_validation_experiments
% PURPOSE 
%   Create a sample of parameters using a latin hypercube design
%   and predict the model performance of the sample using the
%   metamodel.
% INPUTS 
%   From the structure parameters the following fields are
%   processed (mind the same naming in the input)
%   
%   parameters.range:
%
%            Range of values for each paramter to normalize the
%            scale.
%
%   datamatrix.reffdata:
%
%            Modeldata of when using default parameter values to 
%            to center the datamatrix
%
%   lhacc:   number of experiments for which parameter sets are generated
%            
% OUTUTS 
%   dmatrix: Predicted data for parameter experiment
% HISTORY 
% First version: 11.10.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)
% NOTE 
% Currently routine does only allow to compute one experiment at
% the time, could possibly changed by adapting matrix operations.



%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

N=length(parameters); % Number of model parameters
range={parameters.range}; % Parameter ranges

%--------------------------------------------------------------------
% CREATE Latin Hypercube design
%--------------------------------------------------------------------

for i=1:N
  UB(i)=range{i}(2);
  LB(i)=range{i}(1);
end

lh=lhsdesign(lhacc,N,'criterion','maximin','iterations',20);
xexp=repmat(LB,[lhacc,1])+lh.*repmat((UB-LB),[lhacc,1]);


