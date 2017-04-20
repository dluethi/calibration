function [dmatrix]=neelin_p_lsq(metamodel,parameters,datamatrix,pvector)

% Forecast using regression metamodel as described in Neelin et al. (2010) PNAS
% NAME 
%   neelin_p
% PURPOSE 
%   Predict data using the metamodel for a parameter matrix
% INPUTS 
%   From the structure metamodel, parameters and datamatrix the following fields are
%   processed (mind the same naming in the input)
%   
%   metamodel.coeff:
%
%          Vector of coefficients for the metamodel [...,N(2+(N-1)/2)] 
%          additional data dimensions possible [Regions,Variables,Time]
%          data dimensions possible (ex:a~[Regions,Variables,Time,N,1])
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
%   datamatrix.refdata:
%
%            Modeldata of when using default parameter values to 
%            to center the datamatrix
%
%   pvector: Parameter values for one experiment with the
%            dimension of [N,1] N=Number parameters
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

coeff=metamodel.coeff;
N=length(parameters); % Number of model parameters
refp=parameters(1).default; % Default modelparameters
range={parameters.range}; % Parameter ranges
refd=datamatrix.refdata; % Reference data
ds=2*N+N*(N-1)/2; %Number of coefficients in the metamodel


%--------------------------------------------------------------------
% CHECK Input consistency
%--------------------------------------------------------------------
if ~exist('coeff','var')
  error('Please give metamodel parameter "coeff" as input')
end

if ~exist('pvector','var')
  error('Please give parameter matrix as input')
end

da=size(coeff);
nd=prod(da(1:end-1));
nda=length(da); %Number of dimensiones in a
dp=size(pvector); 
dd=da(1:end-1); %Dimension of the data
cvector = reshape(coeff,[nd,da(end)]);

if length(pvector) > N
   error(['Parameter vector wrongly specified, only one set of' ...
	  ' parameters allowed'])
elseif da(end)~=ds
  error('number of coefficients in metamodel disagrees with parameter count')
elseif length(pvector)<N
  error(['Number of parameters in parameter matrix does not agree' ...
	 ' with coeff'])
end

%Normalize parameter values by the total range

for i=1:N
  varp(i)=abs(diff(range{i}));
end

indent=1.5;
pvector=roundn((pvector-refp)./varp,-3);

%--------------------------------------------------------------------
% ALLOCATE Output variables
%--------------------------------------------------------------------

dmatrix=NaN(dd);

%--------------------------------------------------------------------
% COMPUTE Data for each parameter experiment
%--------------------------------------------------------------------

% Compute index vector for all possible pairs
for i=1:N
	a(i) = pvector(i);
end
pqn=allcomb(1:N,1:N);
cnt=1;
for i=1:length(pqn)
  if pqn(i,1)>=pqn(i,2)
	a(N+cnt) = pvector(pqn(i,1)).*pvector(pqn(i,2));
    cnt=cnt+1;
  end
end

cmatrix = cvector * a';
dmatrix = reshape(cmatrix,[dd]);

dmatrix = dmatrix + refd;
  
if isfield(datamatrix,'limits')
  for i=1:length(datamatrix.variables)-1
    indd{datamatrix.variables{1}}=i;
    % Find values in datamatrix which are small than lower limit
    % and set to lower limit
    tmpdata=dmatrix(indd{:});
    llim=datamatrix.limits{i}(1);
    tmpdata(tmpdata < llim)= llim; 
    dmatrix(indd{:})=tmpdata;

    % Find values in datamatrix which are larger than upper limit
    % and set to upper limit
    tmpdata=dmatrix(indd{:});
    hlim=datamatrix.limits{i}(2);
    tmpdata(tmpdata > hlim)= hlim;
    dmatrix(indd{:})=tmpdata;
  end
end                                            
