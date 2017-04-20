function metaparam(metamodel,parameters,datamatrix)

% Visualize fitted metamodel parameters
% NAME 
%   metaparam
% PURPOSE 
%   Show normalized linear,quadratic and interaction terms for each 
%   parameter and interaction
% INPUTS 
%   The structure metamodel, parameters,datamatrix are used
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
refd=parameters(1).default;
pmatrix=parameters(1).experiments(1:ds,:);
range={parameters.range};
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

% Select all matrix dimensions if no score used
dd=ndims(datamatrix.refdata);
if dd>2
  for i=1:dd
    indd{i}=':';
  end
else
  indd=':';
end

%Normalize parameter values by the total range

for i=1:N
  varp(i)=abs(diff(range{i}));
end

pmatrix=(pmatrix-repmat(refd,[ds,1]))./repmat(varp,[ds,1]);

% Standard deviation of the data along all experiments

stddata=std(datamatrix.moddata,0,dd+1);

indent=1.5;
%--------------------------------------------------------------------
% PLOT Metamodel parameters
%--------------------------------------------------------------------

figure;

ax1=axes('Position',[.35 .2 0.5 .65])
t=eye(N+1,N+1);
t2=zeros(1,N+1);

for i=1:length(pqn)
   t(pqn(i,1),pqn(i,2))=0.33;
end

for i=1:N
  tmp=abs(squeeze(metamodel.a(indd{:},i))*((parameters(i).range(2)-parameters(1).default(i))./varp(i)))./stddata;
  lc(i)=mean(tmp(:));
  tmp=abs(squeeze(metamodel.B(indd{:},i,i))*((parameters(i).range(2)-parameters(1).default(i))./varp(i))^2)./stddata;
  lq(i,i)=mean(tmp(:));
end

for p=1:length(pqn)
  ne=2*N+p;
  tmp=abs((pmatrix(ne,pqn(p,1))+indent)*(pmatrix(ne,pqn(p,2))+indent)...
	  *squeeze(metamodel.B(indd{:},pqn(p,1),pqn(p,2))))./stddata;
  lq(pqn(p,1),pqn(p,2))=mean(tmp(:));
end

pcolor(t');
colormap([1 1 1 ; 0.9 0.9 0.9;  0.7 0.7 0.7; 0.5 0.5 0.5])
caxis([0 1]);

for p=1:length(pqn)
  text(pqn(p,1)+0.15,pqn(p,2)+0.3,num2str(roundn(lq(pqn(p,1),pqn(p,2)),-2)),'Fontsize',12,'Color','k','Linewidth',1.5)
end

for i=1:N
  text(i+0.15,i+0.3,num2str(roundn(lq(i,i),-2)),'Fontsize',12,'Color','k','Linewidth',1.5)
end

pnames={parameters.name};
pnames_tex={parameters.name_tex};
set(gca,'XTick',[1.5:N+.5],'XTickLabel',pnames)
ticks=get(gca,'xtick');
set(gca,'xticklabel',[],'yticklabel',[]);
text(ticks,ones(N,1)*0.8,pnames_tex,'rotation',45,'Fontsize',14,'HorizontalAlignment', 'Right')
textLabels = findall(gca,'Tag','XTickLabel');
set(textLabels, 'HorizontalAlignment', 'Right')
set(gca,'Fontsize',14)
ax2=axes('Position',[.2 .2 0.10 .65])
te=ones(N+1,2)*0.66;
pcolor(te)
caxis([0 1]);
set(gca,'Xtick',[],'YTick',[1.5:N+1.5],'YTickLabel',pnames)
set(gca,'Fontsize',14)

for i=1:N
  text(1.1,i+0.5,num2str(roundn(lc(i),-2)),'Fontsize',12,'Color','k','Linewidth',1.5)
end

axes('Visible','off')
text(0.59,.95,'B','Fontsize',16)
text(0.13,.95,'a','Fontsize',16)
set(gcf,'Paperposition',[0 0 8 6])




