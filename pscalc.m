function [pi ps]=pscalc(moddata,obsdata,stddata)

% Function that computes the Performance Score (PS) introduced in 
% Bellprat et al. (2012) JClimate for a number of simulations
% NAME 
%   pscalc
% PURPOSE 
%
%   Compute a PS, an approximation of the Gaussian Likelyhood
%   without the consideration of covariances, scaled by the
%   interannual variability, observational uncertainty and internal
%   variability. 
%
% INPUTS 
%
%   moddata: Data of the model output simulations, last dimension
%            needs to correspond to different simulations if any
%            existent
%
%   obsdata: Data of the observations. Dimensions need to agree
%            with moddata, except of one possible additional
%            dimension in moddata, which is assumed to correspond
%            to multple simulations
%
%   stddata: Standard devation of uncertainties to scale the least
%            square errors.
%
% OUTUTS 
%  pi:      Performance index for all data points 
%  ps:      PS value for each simulation
% HISTORY 
% First version: 11.10.2013
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)


%--------------------------------------------------------------------
% CHECK Input consistency
%--------------------------------------------------------------------

sd=size(moddata);
so=size(obsdata);
ss=size(stddata);
numsim=1;

if (ndims(moddata)-ndims(obsdata))==1
 numsim=sd(end);

  if so~=sd(1:end-1)
    error('Dimensions of model and observation matrix does not agree')
  end

elseif sd~=so
  error('Dimensions of model and observation matrix does not agree')
end

if ss~=so
    error('Dimensions of observation and error matrix does not agree')
end

% Expand observation and standard error data to dimenions of the
% data matrix

if numsim>1
  obsdata=repmat(reshape(obsdata,[1, so]),[numsim ones(1,length(so))]);
  obsdata=permute(obsdata,[2:length(so)+1 1]);
  stddata=repmat(reshape(stddata,[1, so]),[numsim ones(1,length(so))]);
  stddata=permute(stddata,[2:length(so)+1 1]);
end

% Selection vector for variable multidimensional data

for i=1:ndims(moddata)-1
    indd{i}=':';
end

%--------------------------------------------------------------------
% COMPUTE Data for each parameter experiment
%--------------------------------------------------------------------
% Compute erformance index
err=moddata-obsdata;
pi=((moddata-obsdata)./stddata).^2;

%varweight=squeeze(mean(mean(mean(pi,1),2),3)./sum(mean(mean(mean(pi,1),2),3)))*100;
%egweight=squeeze(mean(mean(mean(pi,1),2),4)./sum(mean(mean(mean(pi,1),2),4)))*100;

% Compute performance score
if numsim>1
  for i=1:numsim
    pitemp=pi(indd{:},i);
    ps(i)=exp(-.5*nanmean(pitemp(:)));
  end
else
  ps=exp(-.5*nanmean(pi(:)));
end




