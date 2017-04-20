function mdata=read_cal_f(startyr,endyr,vars,sims,indir)

% Input routine which reads CCLM simulations based on monthly means
% postprocessed by post-processing script of E.Zubler.
%
% NAME 
%   read_cal_f
%
% PURPOSE 
%   Read data of a number of experiments for specified period and
%   variables and computes spatial averages for PRUDENCE regions
%
% INPUTS 
%   startyr: Start year of data window [integer]
%   endyr:   End year of data window [integer]
%   vars:    Vector of variables that are read [cell array!!]
%   sims:    Cell array of simulation names [cell array]
%   indir:   Directory where simulations are stored [string]
%   prud:    Switch if averaged to PRUDENCE regions [logic]
%
% OUTUTS 
%   mdata:   Data matrix with dimensions
%   [Year,Month,Variable,Simulation] 
% 
% HISTORY 
% First version: 11.10.2013
%
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)


%--------------------------------------------------------------------
% DEFINITIONS Of variable metadata
%--------------------------------------------------------------------

[varname,ofact1,ofact2,mfact1,mfact2,unit,clat,clon,gridd,obsmask]=var_meta(vars);

%--------------------------------------------------------------------
% DEFINE Indices
%--------------------------------------------------------------------

% Month indices	 
months={'01','02','03','04','05','06','07','08','09','10','11', '12'};

% Prudence regions	 
prudence(1,:)=[50 59 -10  2];    % BI: British Isles
prudence(2,:)=[36 44 -10  3];    % IP: Iberian Peninsula
prudence(3,:)=[44 50  -5  5];    % FR: France
prudence(4,:)=[48 55  -2 16];    % ME: Mid-Europe
prudence(5,:)=[55 70   5 30];    % SC: Scandinavia
prudence(6,:)=[44 48   5 15];    % AL: Alps
prudence(7,:)=[36 44   3 25];    % MD: Mediterranean
prudence(8,:)=[44 55  16 30];    % EA: Eastern Europe

nyears=endyr-startyr+1;
nvar=length(vars);
nsim=length(sims);

% Allocate output variable
mdata=NaN(nyears,12,8,nvar,nsim); % Years, Months, Regions, Variables, Simulations

%--------------------------------------------------------------------
% READ simulation data
%--------------------------------------------------------------------


for n=1:nsim
  display(['Reading data of simulation ' sims{n} ' for period ' ...
	  num2str(startyr) '-' num2str(endyr)])
  for i=1:nvar

    % Check observations mask
    if obsmask{i}
      ncid = netcdf.open([obsmask{i}],'NC_NOWRITE'); 
      varid = netcdf.inqVarID(ncid,'mask');
      mask =  netcdf.getVar(ncid,varid);
      mask(mask==0)=NaN;
      netcdf.close(ncid);
    else
      mask=false;
    end

    % Read grid data for each variable
    ncid = netcdf.open([gridd{i}],'NC_NOWRITE'); 
    varid = netcdf.inqVarID(ncid,char(clat{i}));
    vlat =  netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,char(clon{i}));
    vlon =  netcdf.getVar(ncid,varid);
    netcdf.close(ncid);

    for y=1:nyears
      % Directory of file
      
      simdir=[indir '/' char(sims{n}) '/clm_out/monmean/' ...
	      num2str(startyr+y-1)];
      for k=1:12
	% File name 
	fname=[char(sims{n}) '_' vars{i} '_rm_monmean_' ...
	       num2str(startyr+y-1) char(months{k}) '.nc'];
	ncid = netcdf.open([simdir '/' fname],'NC_NOWRITE'); 
	varid = netcdf.inqVarID(ncid,varname{i});
	data  = netcdf.getVar(ncid,varid);
	m_fv=netcdf.getAtt(ncid,varid,'_FillValue');
	data(data==m_fv)=NaN;
	data(isnan(mask))=NaN; % Mask points where no observations available
	data=data*mfact1(i) + mfact2(i);
	netcdf.close(ncid);
	for r=1:8
	  [dx,dy]=find(vlat>=prudence(r,1) & vlat<=prudence(r,2) & ...
		       vlon>=prudence(r,3) & vlon<=prudence(r,4));
	  datatmp=data(min(dx):max(dx),min(dy):max(dy));
	  mdata(y,k,r,i,n)=nanmean(nanmean(datatmp(:)));
	  clear datatmp
	end
      end
    end
  end
end




