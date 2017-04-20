
var={'tg','rr','cld','tx'};

obsmatrix = NaN(5,12,8,4);

    filepath='obs_dat.nc';
    for m=1:size(var,2)
      aux0=ncread(filepath,var{m});
      aux1=permute(reshape(aux0,[12,5,8]),[2 1 3]);
      obsmatrix(:,:,:,m)=aux1;
    end


    

    
