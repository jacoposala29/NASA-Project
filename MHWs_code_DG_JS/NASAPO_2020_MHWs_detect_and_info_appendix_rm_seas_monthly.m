function data3d_noseas = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_monthly(data3d)

% data3d is monthly data and a function of (dim1,dim2,dim3): dim3 is time and
% the timestep is monthly

data3d_noseas = data3d;
for mm=1:12
    data3d_noseas(:,:,mm:12:end) = data3d_noseas(:,:,mm:12:end) - ...
        repmat(mean(data3d(:,:,mm:12:end),3),[1 1 size(data3d_noseas(:,:,mm:12:end),3)]);
end

end