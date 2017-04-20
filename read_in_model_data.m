function mod_matrix = read_in_model_data(filepath)

var={'T_2M','TOT_PREC','CLCT','DTR'};

mod_matrix = NaN(5,12,8,4);

for m=1:size(var,2)
  aux0=ncread(filepath,var{m});
  aux1=permute(reshape(aux0,[12,5,8]),[2 1 3]);
  mod_matrix(:,:,:,m)=aux1;
end
mod_matrix(:,:,:,1)=mod_matrix(:,:,:,1) - 273.15;
mod_matrix(:,:,:,3)=mod_matrix(:,:,:,3) * 100.0;

clear var aux0 aux1 m;

    

    
