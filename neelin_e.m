function metamodel=neelin_e(parameters, datamatrix)

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
dvector=reshape(datamatrix.moddata,[prod(sd(1:end-1)),sd(end)]);


%--------------------------------------------------------------------
% CHECK Input consistency
%--------------------------------------------------------------------

dm=2*N;
ds=2*N+N*(N-1)/2; %Number of experiments required to estimate the metamodel
di=N*(N-1)/2; %Number of all possible pairs
dp=size(pmatrix);
rmsest=false; % Least-square estimation of inter-action terms
intest=true; % Determination of inter-action terms 


if dp(2)~=N
  error('Dimension of pmatrix does not correspond to number of parameters')
elseif dp(1)<dm
  error('Linear set of equations is underdetermined, parameter matrix too short')
elseif dp(1)>=dm && dp(1)<ds
  display(['Not enough experiments specified, interaction parameters' ...
	   'not determined'])
  intest=false;
elseif dp(1)>ds
  display(['Linear system of equations overdetermined, addtional experiments'...
           'used to constrain interaction terms using least-squares'])
 rmsest=true;
end

% Check if default value is in the center of the parameter ranges
lwb=false(1,N);upb=false(1,N);

for i=1:N
  if sum(find(pmatrix(:,i)<refp(i)))==0
    lwb(i)=true;
    display(['Default of parameter ' parameters(i).name ...
	     ' is taken at the lower bound'])
  end
  if sum(find(pmatrix(:,i)>refp(i)))==0
    upb(i)=true;
    display(['Default of parameter ' parameters(i).name ...
             ' is taken at the upper bound'])
  end
end

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

a=zeros(nd,N,1); B=zeros(nd,N,N);

%--------------------------------------------------------------------
% DEFINE Additional needed vectors
%--------------------------------------------------------------------

% Compute index vector for all possible pairs
pqn=allcomb(1:N,1:N);
cnt=1;
for i=1:length(pqn)
  if pqn(i,1)>=pqn(i,2)
   cind(cnt)=i;
   cnt=cnt+1;
  end
end
pqn(cind,:)=[];

%--------------------------------------------------------------------
% DETERMINE PARAMETERS FOR MULTIVARIATE QUADRATIC MODEL
%--------------------------------------------------------------------

for i=1:nd % Estimate metamodel for each datapoint
  for n=1:N % Loop over number of parameters
    ne=1+(n-1)*2; %Index of single parameter experiments 
    xv=[pmatrix(ne,n),0,pmatrix(ne+1,n)];
    yv=[dvector(i,ne),0,dvector(i,ne+1)];
    if lwb(n)
      xv=[0,pmatrix(ne,n),pmatrix(ne+1,n)];
      yv=[0,dvector(i,ne),dvector(i,ne+1)];
    end
    if upb(n)
      xv=[pmatrix(ne,n),pmatrix(ne+1,n),0];
      yv=[dvector(i,ne),dvector(i,ne+1),0];
    end
    abtemp=polyfit(xv,yv,2); % Second order polynomial regression
    a(i,n)=abtemp(2); %Write into output variables
    B(i,n,n)=abtemp(1);
  end
  
  % Estimate interaction terms
  if intest
    for n=1:di % Loop over all possible combinations of pairs
      i1=pqn(n,1); i2=pqn(n,2); ne=n+2*N; % Indices of parameters for interactions
	intcont=(dvector(i,ne) ...
            ... % Substract linear contribution ...
        -a(i,i1)*pmatrix(ne,i1)-a(i,i2)*pmatrix(ne,i2) ...
            ... % Substract quadratic contribution
        -B(i,i1,i1)*pmatrix(ne,i1).^2-B(i,i2,i2)*pmatrix(ne,i2).^2);
	intcont=dvector(i,ne)-dvector(i,(i1-1)*2+1)-dvector(i,(i2-1)*2+1);
	B(i,i1,i2)=intcont/(2*pmatrix(ne,i1)*pmatrix(ne,i2));
    end
  end
  
  % Use additional simulations to constrain interaction parameters
  % using a least square minimization
end

if rmsest
  intind=pmatrix./pmatrix;
  rgv=[100,2];
  for k=1:length(rgv) % Loops of error reduction
    for n=1:di % Number of interactions
      i1=pqn(n,1); i2=pqn(n,2); ne=n+2*N; rg=rgv(k);acc=500; 
      % Indices of parameters for interactions
      % Search for interaction experiments
      expi=find(sum(intind(:,[pqn(n,1),pqn(n,2)]),2)==2);  
      % Create a vector of interaction parameters with space rg and
      % acuracy acc centered around the original estimated
      % interaction term
      
      if nd>1
	for p=1:nd
	  bint(p,:)=linspace(B(p,i1,i2)-rg*B(p,i1,i2),B(p,i1,i2)+rg*B(p,i1,i2),acc);
	end
      else
	bint=linspace(B(i1,i2)-rg*B(i1,i2),B(i1,i2)+rg*B(i1,i2),acc);
      end
      for j=1:length(expi)
	for i=1:acc
	  Btmp=B;
	  x=reshape(pmatrix(expi(j),:),[1 N]);
	  xa=squeeze(repmat(x,[nd 1]));
	  xh1=reshape(pmatrix(expi(j),:),[1 1 N]);
	  xh2=reshape(pmatrix(expi(j),:),[1 N 1]);
	  xb1=squeeze(repmat(xh1,[nd,N,1]));
	  xb2=squeeze(repmat(xh2,[nd,1,N]));
	  xb1(:,i1,i1)=pmatrix(expi(j),i1);xb1(:,i2,i2)=pmatrix(expi(j),i2);
          xb2(:,i1,i1)=pmatrix(expi(j),i1);xb2(:,i2,i2)=pmatrix(expi(j),i2);
	  if nd>1
	    Btmp(:,i1,i2)=bint(:,i); Btmp(:,i2,i1)=bint(:,i);            
	    dtmp(i,j,:)=sum(xa'.*a')+sum(sum(xb2.*xb1.*Btmp,3),2)';
	  end
	end
	 if nd>1
	   dcmp(j,:)=dvector(:,expi(j));
	 else 
	   dcmp(j)=dvector(expi(j));
	 end
      end 
      display(['Use least-square esitmation for interaction parameter' ...
	       ' B' num2str(pqn(n,1)) num2str(pqn(n,2))]);
      if nd>1
	[m im]=min(squeeze(sum((dtmp-repmat(reshape(dcmp,[1 ...
		    size(dcmp)]),[acc,1,1])).^2,2)));
	[mo imo]=min(squeeze((dtmp-repmat(reshape(dcmp,[1 ...
                    size(dcmp)]),[acc,1,1])).^2));
	imo=squeeze(imo);
	size(imo);
	vsum(k,n)=mean(var(imo));
	for p=1:nd
	  B(p,i1,i2)=bint(p,im(p)); % B(p,i2,i1)=bint(p,im(p));
	end
      else
      end
    end % for n
  end % for k
end % if rmest

% reshape a and B to original data structure
a=reshape(a,[sd(1:end-1),N]);B=reshape(B,[sd(1:end-1),N,N]);

metamodel.a=a; metamodel.B=B;



