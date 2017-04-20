function [errstd errps]=errmeta(metamodel,parameters,datamatrix)

% Estimate the error of the metamodel to predict indpendent model information
% NAME 
%   errmeta
% PURPOSE 
%   Predict modeldata on of independt simulations and estimate the
%   standard error of the metamodel
%
% INPUTS 
%   The structure metamodel, parameters and datamatrix 
%   are used for the inpute of neelin_p. Addionally the parameter
%   matrix is read 
%
%   parameters.experiments:
%
%           Parameter matrix of experiments on which metamodel is
%           estimated
%   phyd:
%
%           Index in the model data along which the prediction
%           error is physically seperated (ex: ind of different
%           model variabiles (temperature,precipitation,clouds))
% OUTUTS 
%   Plot: Scatter plot for each defined variable of simulated and
%   predicted points
%   errstd: Standard error to predict the model data
%   errps: Standard error to predict the model score
% HISTORY 
% First version: 11.10.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)



%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

N=length(parameters);
pmatrix=parameters(1).validation % Parameter values of validation experiments
sd=size(datamatrix.moddata);
nd=prod(sd(1:end-1)); % Number of datapoints                                  
numvar=length(datamatrix.variables)-1;
indvar=datamatrix.variables{1};
refd=datamatrix.refdata;
sr=size(refd);

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

predval=NaN([size(datamatrix.refdata),length(pmatrix)]); % Allocate data

for i=1:length(pmatrix)
  predval(indd{:},i)=neelin_p(metamodel,parameters,datamatrix,pmatrix(i,:));
end
    	
errpred=predval-datamatrix.valdata;

figure;
for i=1:numvar
  mysubplot(1,numvar,i)
  indd{indvar}=i;
  tmpr=datamatrix.obsdata(indd{:});
  tmps=squeeze(datamatrix.valdata(indd{:},:))-repmat(tmpr,[size(tmpr)./size(tmpr) 10]);
  tmpp=squeeze(predval(indd{:},:))-repmat(tmpr,[size(tmpr)./size(tmpr) 10]);
  lwb=(min(tmpp(:))-abs(min(tmpp(:)))*.15);
  upb=(max(tmpp(:))+abs(min(tmpp(:)))*.15);
  [hs xs]=hist(tmps(:),linspace(lwb,upb,50));
  [hp xp]=hist(tmpp(:),linspace(lwb,upb,50));
  hold on
  errs(i,:)=prctile(tmps(:)-tmpp(:),[5 95]);
  refln=linspace(lwb,upb,1000);
  f = [refln-errs(i,1)/2;flipdim(refln-errs(i,2)/2,1)];
  plot(tmps(:),tmpp(:),'.','MarkerSize',.1,'color','k');
  fiv=patch([refln; flipdim(refln,1)], f, [5 5 5]/8, 'EdgeColor',[5 5 5]/8);
  hr=plot(refln,refln,'Linewidth',1.5,'color','k');
  ylabel(['Predicted ' datamatrix.variables{i+1}],'Fontsize',14)
  xlabel(['Simulated ' datamatrix.variables{i+1}],'Fontsize',14)
  set(gca,'Fontsize',14,'Ylim',[lwb upb],'Xlim',[lwb upb],'Box','on')
  hpdf=(upb-lwb)/5; scalef=hpdf/max(hs);
  plot(smooth(hs,10)*scalef+lwb,xs,'k','Linewidth',1.5)
  plot(xp,smooth(hp,10)*scalef+lwb,'k','Linewidth',1.5)
  rstat=regstats(tmps(:),tmpp(:));
  if i==1
    hl=legend([fiv(1),hr],['95% range'],['R^2=' num2str(roundn(rstat.rsquare,-2))],2);
    set(hl,'Box','off')
  else
    hl=legend([hr],['R^2=' num2str(roundn(rstat.rsquare,-2))],2);
    set(hl,'Box','off')
  end
end


errstd=nanstd(errpred,0,length(sd)); %Standard error

% Compute resulting error in model performance by adding
% "nsam" times Gaussian whithe noise to the uncertainty to the reference  
% simulation and computing multiple times the effect on the
% model score

nsam=1000; % Number of samples drawn

if isfield(datamatrix,'score')
   sampref=repmat(refd,[sr./sr nsam]);
   linsamp=reshape(sampref,[nd nsam]);
   linerr=reshape(errstd,[1 nd]);
   for n=1:nd
       linsamp(n,:)=linsamp(n,:)+normrnd(0,linerr(n),1,nsam);
   end
   sampref=reshape(linsamp,[sr nsam]);
   if strcmp(datamatrix.score,'ps')
     load('/home/omarb/CLMEVAL/emulate/NEELIN/calmo/data/stddata')
     [pi ps]=pscalc(sampref,datamatrix.obsdata,stddata);
   elseif strcmp(datamatrix.score,'psnam')
     load('/home/omarb/CLMEVAL/emulate/NEELIN/calmo/data/stddata_nam')
     [pi ps]=pscalc(sampref,datamatrix.obsdata,stddata);
   end
   errps=std(ps);
else
  errps=NaN;
end






