function data_on_grid = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_to_regular_grid(...
    ecco_xc,ecco_yc,ecco_data)

% ecco_xc,ecco_yc,ecco_data should all have the same size and only depend on space, e.g.
% ecco_data(lon,lat,tile)

% i = 1;
% bfr_d = reshape(data_in.data,data_in.size_source);
% bfr_x = data_in.xc;
% bfr_y = data_in.yc;
% % select
% bfr_d = bfr_y(:,:,:,i);
% bfr_x = bfr_x(:,:,:,i);
% bfr_y = bfr_y(:,:,:,i);

if min(ecco_xc(:))<0 || max(ecco_xc(:))>=380
    stop2check
else
    [Xnew,Ynew] = ndgrid(20.5:379.5,-89.5:89.5);
end

bfr_F = scatteredInterpolant(ecco_xc(:),ecco_yc(:),ecco_data(:),'nearest');

data_on_grid.data = reshape(bfr_F(Xnew(:),Ynew(:)),size(Xnew));
data_on_grid.X    = Xnew;
data_on_grid.Y    = Ynew;
return