function metamodel=neelin_e_lsq(parameters, datamatrix)

% Quadratic regression metamodel as described in Neelin et al. (2010) PNAS
% NAME 
%   neelin_f
% PURPOSE 
%   Estimate a mutlivariate quadratic metamodel which estimates quadratic
%   regressions in each parameter dimensions and computes interaction
%   terms for all pair of parameter experiments
% INPUTS 
%   From the structure parameters and datamatrix the following fields are
%   processed (mind the same naming in the input)
%  
%   parameters.experiments:
%
%            Parameter values for each experiment with the
%            dimension of [N, 2*N+N*(N-1)/2]
%            The structure NEEDS to be as follows
%            Example for 2 parameters (p1,p2)
%           
%            [p1_l dp2 ] ! Low parameter value for p1 default dp2
%            [p1_h dp2 ] ! High parameter value for p1 default dp2
%            [dp1  p2_l] ! Low parameter value for p2 default dp1
%            [dp1  p2_h] ! Hihg parameter value for p2 default dp1
%            [p1_l p2_h] ! Experiments with interaction (no default)
%                        ! Additional experiments used to
%                        constrain interaction terms
%   parameters.range:
%
%            Range of values for each paramter to normalize the
%            scale.
%
%   parameters.default:
%
%            Default values of parameters to center the scale
% 
%   datamatrix.moddata:
%           
%            Modeldata corresponding to the dimenoins of
%            parameter.experiments
%    
%   datamatrix.reffdata:
%
%            Modeldata of when using default parameter values to 
%            to center the datamatrix
%            fitted. Not needed if metamodel fits score data directly
% OUTUTS 
%   structure metamodel.
%   a: Metamodel parameter for linear terms [N,1]
%   B: Metamodel parameter for quadratic and interaction terms
%      [N,N]. Quadratic terms in the diagonal, interaction terms
%      in the off-diagonal. Matrix symetric, B(i,j)=B(j,i).
% HISTORY 
% First version: 11.10.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)

%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

N=length(parameters); % Number of model parameters
refp=parameters(1).default; % Default modelparameters
range={parameters.range}; % Parameter ranges
pmatrix=parameters(1).experiments;  % Parameter values 
sd=size(datamatrix.moddata);
nd=prod(sd(1:end-1)); % Number of datapoints                                  
				    

% Liearize all dimensions 
refd=datamatrix.refdata(:);
dvector=reshape(datamatrix.moddata,[nd,sd(end)]);

dm=2*N;
ds=2*N+N*(N-1)/2; %Number of experiments required to estimate the metamodel
di=N*(N-1)/2; %Number of all possible pairs
dp=size(pmatrix);

% Normalize parameter values by the total range and center around
% default value

for i=1:N
  varp(i)=abs(diff(range{i}));
end

pmatrix=roundn((pmatrix-repmat(refp,[dp(1),1]))./repmat(varp,[dp(1),1]),-3);
dvector=dvector-repmat(refd,[1,sd(end)]);


%--------------------------------------------------------------------
% ALLOCATE Output variables
%--------------------------------------------------------------------

a=zeros(dp(1),ds);
coeff=zeros([sd(1:end-1),ds]);

%--------------------------------------------------------------------
% DEFINE Additional needed vectors
%--------------------------------------------------------------------

% Compute index vector for all possible pairs
for i=1:N
	a(:,i) = pmatrix(:,i);
end
pqn=allcomb(1:N,1:N);
cnt=1;
for i=1:length(pqn)
  if pqn(i,1)>=pqn(i,2)
	a(:,N+cnt) = pmatrix(:,pqn(i,1)).*pmatrix(:,pqn(i,2));
    cnt=cnt+1;
  end
end

dvector = permute(dvector,[2,1]);

mm_coef = a \ dvector;

mm_coef = permute(mm_coef,[2,1]);

metamodel.coeff = reshape(mm_coef,[sd(1:end-1),ds]);



