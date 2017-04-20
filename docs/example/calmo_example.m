% Master file for CALMO calibration suite
% NAME 
%   calmo.m
% PURPOSE 
%   Definitions for calibration, I/O of data, calling routines
% HISTORY 
%   First version: 11.10.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)

% Three data structures are defined:
%
% - Parameters: Contains all default parameter information 
%               and settings of experiments
%
% - Datamatrix: Contains all data of experiments and validation
%               runs. Defines if a score is computed out of the
%               data or if directly a score is read
%
% - Metamodel:  Contains the metamodel parameters which are fitted

%--------------------------------------------------------------------
% DEFINE Calibration structures (Parameters,Datamatrix,Metamodel)
%--------------------------------------------------------------------

clear all; close all;

% Define Parameters

paramn={'rlam_heat','entr_sc','qi0','uc1','rootdp','tkhmin','radfac','soilhyd'}; % No TEX interpreter
paramnt={'rlam\_heat','entr\_sc','qi0','uc1','rootdp','tkhmin','radfac','soilhyd'}; % Parameter names

range={[0.1 10]; % Parameter ranges (min,default,max)
	[log(3E-5) log(3E-3)];
	[0 0.0001];
	[0 1.6];
	[0.5 1.5];
        [0.01 2];
        [0.3 0.9];
        [1 6]}

default={[1 roundn(log(3E-4),-4) 0.00000 0.8 1 1 0.5 1]};


expval=NaN % Read-in: Experiment values to fit metamodel
valcon=NaN % Read-in: Experiment values to additionally constrain the metamodel
valval=NaN % Read-in: Experiment values to validate the metamodel

parameters=struct('name',paramn,'range',range','default',default,'experiments', ...
		  expval,'constrain',valcon,'validation',valval,'name_tex',paramnt)

% For the example a already computed parameter structure is loaded 

load('./data/parameters_cal_ea')

% Define Datamatrix

moddata=[]; 
valdata=[];
obsdata=[];	
refdata=[];
variables={4,'T2M','PR','CLCT'} % Index of variables in moddata,
                                  % and name of each model variable
				  % if no score data used

scoren='ps'; % If scoren [], no computation of score assumed
                  % and data values corresponding to score values

datamatrix=struct('moddata',moddata,'refdata',refdata,'valdata',valdata,'obsdata',...
obsdata,'score',scoren)

% For the example a already computed datamatrix structure is loaded 

load('./data/datamatrix_cal_ea')


%-----------------------------------------------------------------
% FIT METAMODEL
%-----------------------------------------------------------------

% For each data vector an independent metamodel has to be fit.
% If a model score (one value per simulation) is estimated only one
% metamodel needs to be fitted. Otherwise you need to loop over the 
% whole data matrix.

metamodel=neelin_e(parameters,datamatrix);
metamodel=neelin_c(parameters,datamatrix,metamodel);
save('./data/metamodel_cal_ea','metamodel')

%-----------------------------------------------------------------
% VALIDATION METRICS
%-----------------------------------------------------------------

% The fitted metamodel is analysed in terms of accuracy and non-linearity
% Different functions can be used for this purpose as described in 
% Bellprat et al. (2012) JGR.

% (1) Check if fitted data is reproduced with 0 errors

ctrlpred(metamodel,parameters,datamatrix)

% (2) Estimate the error of independent simulations and plot scattorplots

[errstd]=errmeta(metamodel,parameters,datamatrix)

% (3) Visualize mean metamodel parameters for linear, quadratic and
%     interaction terms

metaparam(metamodel,parameters,datamatrix);

% (4) Test quadratic approximation along a speicifc parameter dimension

% squarmetamodel=neelin_c(parameters,datamatrix,metamodel);
%squreefit(metamodel,parameters,datamatrix) % Not yet finished

% (5) Visualize performance landscape for each parameter pair
%     between all experiments

planes(metamodel,parameters,datamatrix)


%-----------------------------------------------------------------
% OPTIMIZATION OF PARAMETERS
%-----------------------------------------------------------------

% The validated metamodel can now be used to find optimal parameter
% values and to perform a perfect model experiment.

% (1) Find optimal model parameters using a latin hypercube
% optimisation

lhacc=1000000; % Number of experiments to sample parameter space

% lhscore: Modelscore for all experiments; 
% lhexp: Latin hypercube parameter experiments
% popt: Parameter setting with highest score

[lhscore lhexp popt]=lhopt(metamodel,parameters,datamatrix,lhacc);

% (2) Plot performance range covered  by the metamodel and compare to
% reference simulation

histplot(lhscore,datamatrix)

% (3) Plot optimised parameter distributions
% errm: error of metamodel

optparam(parameters,lhscore,lhexp,popt,0.025)



