function data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
    start_tstep_msk,end_tstep)
% INPUT:
% data2ave is data2ave(dim1,dim2,time) is the data to average
%
% start_tstep_msk is start_tstep_msk(dim1,dim2,time) is a mask indicating
% the initial timesteps for the average: the average is stored at timesteps
% where start_tstep_msk is true
%
% end_tstep is end_tstep(dim1,dim2,time) and stores the ending timestep for
% the everage
%
% PLEASE NOTE that if start_tstep_msk is equal to find_MHWs_info.start_tstep_msk
% in main and  end_tstep is equal to find_MHWs_info.end_tstep in main THEN
% THE AVERAGE IS DONE OVER EACH EVENT and stored at the timestep that
% corresponds to the start of the event
%


% now show how to do an average of a variable (e.g. input data_in.data) over each
% event "relevant" timesteps
% 
data_select_event_ave = nan(size(data2ave));
tstep3d =  permute(repmat([1:size(data2ave,3)]',[1 size(data2ave,1) size(data2ave,2)]),[2 3 1]);
% disp('Running NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs.m')
for i = 1:size(data2ave,3)
    clear bfr*
    bfr_end =  end_tstep(:,:,i);
    bfr_start_msk = start_tstep_msk(:,:,i);
    
%     %%%% tests and edge cases for when we tried onset without peak (no used
%     %%%% in the end, hence the if statement should always be false)
%     if min(bfr_end(:))<i
%        disp(['Timestep ' num2str(i) ': at ' num2str(sum(bfr_end(:)<i)) ...
%            ' locations, the peak is at 1st MHW timestep: onset will be the 1st MHW timestep as a fix'])
%        
%        if min(bfr_end(:))-i<-1
%            thisshouldnothappen
%        else
%            bfr_end(bfr_end<i) = bfr_end(bfr_end<i)+1;
%        end
%     end
    %%%% tests and edge cases for the decline phase to not include the peak
    %%%% time step (this could have been triggered also for the test to not
    %%%% include peak in onset phase, which we did not use)
    if min(bfr_end(:))<i
        disp('This statement should only appear for the decline phase or (maybe?) for monthly data')
        disp(['Timestep ' num2str(i) ': at ' num2str(sum(bfr_end(:)<i)) ...
            ' locations, the peak is at the LAST MHW timestep: ' ...
            'there is no decline phase for these events!!!'])
        
        if min(bfr_end(:))-i<-1
            thisshouldnothappen
        else
            bfr_start_msk(bfr_end<i) = 0;
            bfr_end(bfr_end<i) = nan;
        end
    end
    %
    bfr_max_end_tstep = max(bfr_end(:));
    %
    if ~isnan(bfr_max_end_tstep)
        bfr_sel = data2ave(:,:,i: bfr_max_end_tstep);
        bfr_tstep_sel=tstep3d (:,:,i: bfr_max_end_tstep);
        bfr_sel(bfr_tstep_sel>repmat(bfr_end,[1,1,size(bfr_tstep_sel,3)]) | ~(repmat(bfr_start_msk,[1,1,size(bfr_tstep_sel,3)]))) = nan;
        data_select_event_ave(:,:,i) = nanmean(bfr_sel,3);
    end
end