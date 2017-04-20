function metamodel=neelin_c(parameters, datamatrix, metamodel)

% Constrain linear and quadratic metamodel terms with additinal simulations
% NAME 
%   neelin_c
% PURPOSE 
%   Use additional simulations to narrow uncertainty of the
%   metamodel paramters, particulary if strong unequal distancies
%   between default and min/max values.
% INPUTS 
%   From the structure parameters and datamatrix the following fields are
%   processed (mind the same naming in the input)
%  
%   parameters.constrain:
%
%            Parameter values for each experiment for
%            additional parameter samplilng
%
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
%            fitted. Not needed if metamodel fits score data
%            directly
%
%   datamatrix.constrain
%
%            Simulation data of experiments used to further sample
%            parameter ranges
%
% OUTUTS 
%   Updatated metamodel structure
% HISTORY 
% First version: 4.11.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)

%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

N=length(parameters); % Number of model parameters
refp=parameters(1).default; % Default modelparameters
range={parameters.range}; % Parameter ranges
pmatrix=parameters(1).experiments;  % Parameter values of neelin experiments
pmatrixc=parameters(1).constrain;

sd=size(datamatrix.moddata);
sc=size(datamatrix.constrain);
nd=prod(sd(1:end-1)); % Number of datapoints                                     
a=reshape(metamodel.a,[nd,N]);
c=zeros(size(a));
B=reshape(metamodel.B,[nd,N,N]);


% Liearize all dimensions 
refd=datamatrix.refdata(:);
dvector=reshape(datamatrix.moddata,[prod(sd(1:end-1)),sd(end)]);
if min(size(pmatrixc))>1
  dvectorc=reshape(datamatrix.constrain,[prod(sc(1:end-1)),sc(end)]);
else
  dvectorc=reshape(datamatrix.constrain,[prod(sc) 1]);
end
dm=2*N;
ds=2*N+N*(N-1)/2; %Number of experiments required to estimate the metamodel
di=N*(N-1)/2; %Number of all possible pairs
dp=size(pmatrix);
dc=size(pmatrixc);
rmsest=false; % Least-square estimation of inter-action terms
intest=false; % Determination of inter-action terms 

% Find parameters that are further constrained

ptmp=pmatrixc-repmat(refp,[dc(1),1]);

if min(size(pmatrixc))>1
  pind=find(round(sum(ptmp))>0);
else
  pind=find(round(ptmp)>0);
end

% Normalize parameter values by the total range and center around
% default value

for i=1:N
  varp(i)=abs(diff(range{i}));
end

pmatrix=(pmatrix-repmat(refp,[dp(1),1]))./repmat(varp,[dp(1),1]);
pmatrixc=(pmatrixc-repmat(refp,[dc(1),1]))./repmat(varp,[dc(1),1]);

dvector=dvector-repmat(refd,[1,sd(end)]);

if min(size(pmatrixc))>1
  dvectorc=dvectorc-repmat(refd,[1,sc(end)]);
else
  dvectorc=dvectorc-refd;
end

%--------------------------------------------------------------------
% CONSTRAIN Neelin parameters
%--------------------------------------------------------------------

for j=1:length(pind)
  pall=[pmatrixc(ptmp(:,pind(j))~=0,pind(j))' pmatrix(1+(pind(j)-1)*2,pind(j)) ...
	pmatrix(pind(j)*2,pind(j)) 0];
  dall=[dvectorc(:,ptmp(:,pind(j))~=0) dvector(:,1+(pind(j)-1)) ...
	dvector(:,pind(j)*2) zeros(size(dvector(:,pind(j)*2)))];
  for i=1:nd
    xv=pall;
    yv=dall(i,:);
    abtemp=polyfit(xv,yv,2); % Second order polynomial regression
    a(i,pind(j))=abtemp(2); %Write into output variables
    B(i,pind(j),pind(j))=abtemp(1);   
    c(i,pind(j))=abtemp(3);
  end
end

% reshape a and B to original data structure
a=reshape(a,[sd(1:end-1),N]);c=reshape(c,[sd(1:end-1),N]);
B=reshape(B,[sd(1:end-1),N,N]);

metamodel.a=a; metamodel.B=B; metamodel.c=c;



