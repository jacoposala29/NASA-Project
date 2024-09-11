function data3d_noseas = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_daily(data3d,data3d_datenum,delta_days)
    
% data3d is daily data and a function of (dim1,dim2,dim3): dim3 is time and
% the timestep is daily

data3d_noseas = data3d;

mmdd = unique([month(data3d_datenum) day(data3d_datenum)],'row');
for i=1:size(mmdd,1)
    clear msk msk2
    
    msk = month(data3d_datenum)==mmdd(i,1) & ...
        day(data3d_datenum)==mmdd(i,2);
    msk2 = find_days(i, data3d_datenum, delta_days);

%     
%     if mmdd(i,1)==2 && mmdd(i,2)==29 
        % for "Feb 29", let's use the mean over all the "Feb 28"
        %         msk2 = month(data3d_datenum)==2 & ...
        %         day(data3d_datenum)==28;
        % for "Feb 29", let's average between "Feb 28" and "Mar 1" - updated
        % June 4
%         msk2 = (month(data3d_datenum)==2 & ...
%             day(data3d_datenum)==28) | ...
%             (month(data3d_datenum)==3 & ...
%             day(data3d_datenum)==1);
%     else
%         msk2 = msk;
    
%     end
    data3d_noseas(:,:,msk) = data3d_noseas(:,:,msk) - ...
        repmat(mean(data3d(:,:,msk2),3),[1 1 size(data3d(:,:,msk),3)]);
end

end