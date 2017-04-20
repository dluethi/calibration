function [varname,ofact1,ofact2,mfact1,mfact2,unit,clat,clon,gridd,obsmask]=var_meta(v)

% Retrieve meta information of model variables that a read in
% read_cal_f and read_obs
% NAME 
%   var_meta
% PURPOSE 
%   Allocate meta information for netcdf I/O and 
%
% INPUTS 
%   v: Cell array of variables names that are used
%
% OUTUTS 
%   varname:  Variable name in netcdf files
%   ofact1:   Factor for observation input
%   ofcat2:   Addition term for observation input
%   ofact1:   Factor for model input
%   ofcat2:   Addition term for model input
%   unit  :   Unit of variable
%   clat/lon: Latitude/Longitude name of variables
%
% HISTORY 
% First version: 11.10.2013
%
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)


nvar=length(v);

for l=1:nvar
  if strcmp(v(l),'t2m')
    varname{l}='T_2M';
    ofact1(l) = 0.01;
    ofact2(l) = 0;
    mfact1(l) = 1;
    mfact2(l) = -273.15;
    unit{l} = 'degC';
    clat{l} = 'Actual_latitude';
    clon{l} = 'Actual_longitude';
    gridd{l}   = '/lhome/omarb/observations/grid_eobs.nc';
    obsmask{l} = '/lhome/omarb/observations/eobs_mask.nc';
    bn = [-30:0.5:30] ;
  elseif strcmp(v(l),'pr')
    varname{l}='pr';
    ofact1(l) = 1;
    ofact2(l) = 0;
    mfact1(l) = 8;
    mfact2(l) = 0;
    unit{l} = 'mm/d';
    clat{l} = 'Actual_latitude';
    clon{l} = 'Actual_longitude';
    gridd{l}   = '/lhome/omarb/observations/grid_eobs.nc';
    obsmask{l} =  '/lhome/omarb/observations/eobs_mask.nc';
    bn = [1:10] ;
  elseif strcmp(v(l),'clct') 
    varname{l}='CLCT';
    ofact1(l) = 1;
    ofact2(l) = 0.0;
    mfact1(l) = 100.0;
    mfact2(l) = 0.0;
    unit{l} = '%';
    obsset{l} = 'CRU_TS_2.1';
    clat{l} = 'latitude';
    clon{l} = 'longitude';
    gridd{l}   = '/lhome/omarb/observations/grid_cru.nc';
    obsmask{l} =   '/lhome/omarb/observations/cru_mask.nc';
    bn = [0:5:100];
  elseif strcmp(v(l),'clcl') 
    varname{l}='CLCL';
    ofact1(l) = 1;
    ofact2(l) = 0.0;
    mfact1(l) = 100.0;
    mfact2(l) = 0.0;
    unit{l} = '%';
    obsset{l} = 'CRU_TS_2.1';
    clat{l} = 'latitude';
    clon{l} = 'longitude';
    gridd{l}   = '/lhome/omarb/observations/grid_cru.nc';
    obsmask{l} =  '/lhome/omarb/observations/cru_mask.nc';
    bn = [0:5:100];  
  elseif strcmp(v(l),'srad')
    varname{l}='ASOB_S';
    ofact1(l) =4.62963e-05;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0.0;
    unit{l} = 'W/m^2';
    obsset{l} = 'EMCWF';
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'trad')
    varname{l}='ATHB_S';
    ofact1(l) = 4.62963e-05;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0.0;
    unit{l} = 'W/m^2';
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'aod')
    varname{l}='AOD550';
    ofact1(l) = 0.01;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = -273.15;
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'ang')
    varname{l}='ANGSTR';
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'abs')
    varidobs(l)=1;
    varname{l}='AODABS';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
  elseif strcmp(v(l),'geop500')
    varname{l}='FI';
    ofact1(l) = 0.010197162;
    ofact2(l) = 0.0;
    mfact1(l) = 0.010197162;
    mfact2(l) = 0;
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'geop700')
    varname{l}='FI';
    ofact1(l) =0.010197162;
    ofact2(l) =0.0;
    mfact1(l) =0.010197162;
    mfact2(l) =0.0;
    unit{l} ='dm';
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'geop850')
    varname{l}='FI';
    ofact1(l) = 0.010197162;
    ofact2(l) = 0.0;
    mfact1(l) = 0.010197162;
    mfact2(l) = 0;
    unit{l} = 'degC';
    clat{l} = 'olat';
    clon{l} = 'olon';
  elseif strcmp(v(l),'toaswd')
    varname{l}='ASOB_T';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    unit{l} = 'w/m^2';
    clat{l} = 'latitude';
    clon{l} = 'longitude';
  elseif strcmp(v(l),'toalwu')
    varname{l}='ATHB_T';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    unit{l} = 'w/m^2';
    clat{l} = 'latitude';
    clon{l} = 'longitude';
  elseif strcmp(v(l),'dtr')
    varname{l}='DTR';
    ofact1(l) = 0.01;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    bn = [-5:0.5:30] ;
    unit{l} = 'K';
    clat{l} = 'Actual_latitude';
    clon{l} = 'Actual_longitude';
  elseif strcmp(v(l),'tminmax')
    varname{l}='TMAX_2M';
    ofact1(l) = 0.01;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    bn = [0.5:0.5:25] ;
    unit{l} = 'K';
    clat{l} = 'lat';
    clon{l} = 'lon';
  elseif strcmp(v(l),'soil')
    varname{l}='W_SO';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    unit{l} = 'm';
    clat{l} = 'lat';
    clon{l} = 'lon';
  elseif strcmp(v(l),'ws')
    varname{l}='W_SO';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    unit{l} = 'm';
    clat{l} = 'Actual_latitude';
    clon{l} = 'Actual_longitude';
  elseif strcmp(v(l),'soil')
    varname{l}='W_SO';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    unit{l} = 'm';
    clat{l} = 'Actual_latitude';
    clon{l} = 'Actual_longitude';
  elseif strcmp(v(l),'hlat')
    varname{l}='ALHFL_S';
    ofact1(l) = 1.0;
    ofact2(l) = 0.0;
    mfact1(l) = 1.0;
    mfact2(l) = 0;
    unit{l} = 'm';
    clat{l} = 'lat';
    clon{l} = 'lon';
  end;
end;