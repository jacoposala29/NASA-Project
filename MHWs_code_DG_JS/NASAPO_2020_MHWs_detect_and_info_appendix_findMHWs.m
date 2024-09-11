function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_findMHWs(data_in,...
    delta_tstep,percentile,flag_remove_trend,fig_path)
%
% data_in.data
% delta_tstep is the minimum duration of MHW events of interest
% percentile is the percentile to use to define MHW events

% data_in.data is data_in.data(dim1,dim2,time) hence:
dim=3; % (this param should not change)

data_out.delta_tstep = delta_tstep;
data_out.percentile  = percentile;
data_out.flag_remove_trend = flag_remove_trend;

% create timestep
data_in.tstep        = [1:length(data_in.data_datenum)]';

%>>> case for daily data (include check in case Feb has always 28 days)

%%%%%%%%%%

%%%%%%%%%% remove seasonal cycle
data_in.data_backup = data_in.data;

disp('---> Remove the seasonal cycle')
tic;
% data_with_seas = data_in.data;
if data_in.flag_monthly_timeseries %isfield('flag_monthly_timeseries', data_in) && 
    % remove seasonal cycle for a monthly timeseries
    data_in.data = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_monthly(data_in.data);
    
elseif data_in.flag_daily_timeseries %isfield('flag_daily_timeseries', data_in) && 
    % remove seasonal cycle for a daily timeseries
    data_in.data = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_daily(data_in.data,data_in.data_datenum,...
        data_in.delta_days_4_percentile);
    
else
    fields(data_in)
    writemorehere
end
toc;

% set to nan x,y locations with constant timeseries
msk_constant_tseries = min(data_in.data,[],3)==max(data_in.data,[],3);
data_in.data(repmat(msk_constant_tseries,[1 1 size(data_in.data,3)])) = nan;

%%%%%%%%%%

data_in.data_seasonal = data_in.data_backup - data_in.data;

% this line should be eventually commented to reduce file size
data_out.data_backup = data_in.data_backup;

%%%%%%%%%% detrend
% data_with_trend = data_in.data;

disp('---> Detrend')
tic;
if flag_remove_trend
    data_in.data = reshape(detrend(reshape(data_in.data,size(data_in.data,1)*size(data_in.data,2),[])')',...
        size(data_in.data,1),size(data_in.data,2),[]);
end
toc;
%%%%%%%%%%
% calculate percentile (in the following, the percentile is calculated
% differently than https://psl.noaa.gov/marine-heatwaves/; the calculation
% here relies on the fact that the seasonal cycle has already been removed)

% data_percentile = prctile(data_in.data,percentile,dim);
% data_percentile3d = repmat(data_percentile,[1,1,length(data_in.tstep)]);
disp('---> Calculate percentile')
tic;
if isfield(data_in, 'flag_monthly_timeseries') && data_in.flag_monthly_timeseries
    % trying to replicate the method by noaa
    
    % Adjusted below so that it works also for daily data
    disp('monthly.....................')
    data_percentile3d = nan(size(data_in.data));
    mmm = {[12 1 2] [1 2 3] [2 3 4] [3 4 5] [4 5 6] [5 6 7] [6 7 8] [7 8 9] ...
        [8 9 10] [9 10 11] [10 11 12] [11 12 1]};
    
    % for mm=1:12
    %     msk_mm = false(size([1:size(data_in.data,3)]'));
    %     msk_mm(mmm{mm}(1):12:end) = true;
    %     msk_mm(mmm{mm}(2):12:end) = true;
    %     msk_mm(mmm{mm}(3):12:end) = true;
    %
    %     data_percentile3d(:,:,mm:12:end) = repmat(...
    %         prctile(data_in.data(:,:,msk_mm),percentile,dim),...
    %         [1 1 size(data_percentile3d(:,:,mm:12:end),3)]);
    % end
    for mm=1:12
        msk_mm = false(size([1:size(data_in.data,3)]'));
        msk_mm(month(data_in.data_datenum)==mmm{mm}(1)) = true;
        msk_mm(month(data_in.data_datenum)==mmm{mm}(2)) = true;
        msk_mm(month(data_in.data_datenum)==mmm{mm}(3)) = true;
        
        data_percentile3d(:,:,month(data_in.data_datenum)==mm) = repmat(...
            prctile(data_in.data(:,:,msk_mm),percentile,dim),...
            [1 1 size(data_percentile3d(:,:,month(data_in.data_datenum)==mm),3)]);
    end
elseif isfield(data_in, 'flag_daily_timeseries') && data_in.flag_daily_timeseries
    
    %%% replicating Hobday with the option to use a different +/- delta
    %%% time, i.e. data_in.delta_days_4_percentile
    %     data_in.delta_days_4_percentile = 5
    data_percentile3d = nan(size(data_in.data));
    
    %     doy_month = month(data_in.data_datenum);
    %     doy_day = day(data_in.data_datenum);
    %     doy_year = year(data_in.data_datenum);
    %     doy_year_unique = unique(doy_year);
    %     data_in.datetime = datetime(data_in.data_datenum, 'ConvertFrom', 'datenum');
    %     doy = day(data_in.data_datenum,'dayofyear');
    for iday = 1:length(data_in.data_datenum)
        bfr_data_percentile3d_test = data_percentile3d(:,:,iday);
        if sum(~isnan(bfr_data_percentile3d_test(:))) == 0
            %         if ~(doy_day == 29 && doy_month == 2)
            clear bfr_day* bfr_month* bfr_check_leap* bfr_day_msk*
            bfr_day_msk = find_days(iday, data_in.data_datenum, data_in.delta_days_4_percentile);
            
            data_percentile3d(:,:,iday) = repmat(...
                prctile(data_in.data(:,:,bfr_day_msk),percentile,dim),...
                [1 1 size(data_percentile3d(:,:,iday),3)]);
            bfr_msk_same_day = day(data_in.data_datenum) == day(data_in.data_datenum(iday)) & ...
                month(data_in.data_datenum) == month(data_in.data_datenum(iday));
            data_percentile3d(:,:,bfr_msk_same_day) = repmat(...
                data_percentile3d(:,:,iday), ...
                [1 1 sum(bfr_msk_same_day)]);
            
        end
    end

    % the calculation of data_percentile3d_smooth is not efficient at this
    % time and needs revisting
    
    if isfield(data_in, 'flag_daily_timeseries_smoothed_percentile') && data_in.flag_daily_timeseries_smoothed_percentile
        
        sz_win_left_side = 45;%30;%15; % the window size is sz_win_left_side*2+1
        sz_win = sz_win_left_side*2+1;
        disp(['For daily data smooth the percentile using a ' num2str(sz_win) ' day window'])
        
        
        data_percentile3d_smooth = nan(size(data_percentile3d));
        for ix=1:size(data_percentile3d,1)
            for iy=1:size(data_percentile3d,2)
                
                data_percentile3d_smooth(ix,iy,sz_win_left_side+1:end-sz_win_left_side) = ...
                    NASAPO_smoothField_TheaS(squeeze(data_percentile3d(ix,iy,:)),'movmean',sz_win);
                
                % fill in the nans at the start and end
                inds_start = find(day(data_in.data_datenum)==1 & month(data_in.data_datenum)==1,2);
                data_percentile3d_smooth(ix,iy,1:sz_win_left_side) = ...
                    data_percentile3d_smooth(ix,iy,inds_start(2):inds_start(2)+sz_win_left_side-1);
                
%                 inds_end   = find(day(data_in.data_datenum)==31-sz_win_left_side & month(data_in.data_datenum)==12,1);
                inds_end   = find(day(data_in.data_datenum)==31 & month(data_in.data_datenum)==12);
                data_percentile3d_smooth(ix,iy,end-sz_win_left_side+1:end) = ...
                    data_percentile3d_smooth(ix,iy,inds_end(end-1):inds_end(end-1)+sz_win_left_side-1);
                
            end
        end
    end
end
toc;
close all

if isfield(data_in, 'flag_daily_timeseries_smoothed_percentile') && data_in.flag_daily_timeseries_smoothed_percentile
    data_percentile3d_not_smooth = data_percentile3d;
    data_percentile3d            = data_percentile3d_smooth;
end

% test1 = data_percentile3d(:,:,day(data_in.data_datenum) == day(data_in.data_datenum(iday)) & ...
%                 month(data_in.data_datenum) == month(data_in.data_datenum(iday)) & ...
%                 year(data_in.data_datenum) == 1999);
% 
% test2 = data_percentile3d(:,:,day(data_in.data_datenum) == day(data_in.data_datenum(iday)) & ...
%                 month(data_in.data_datenum) == month(data_in.data_datenum(iday)) & ...
%                 year(data_in.data_datenum) == 2003);
% 
% test_diff = test1 - test2;
% 
% unique(test_diff(:))
% pcolor(test1 - test2); shading flat; colorbar;

% ix_test = 10; %800
% iy_test = 5; %500
% % to visualize data backup (original)
% nexttile();plot(data_in.data_datenum, squeeze(data_in.data_backup(ix_test,iy_test,:)))
% 
% % to visualize data without seasonal cycle
% nexttile();plot(data_in.data_datenum, squeeze(data_in.data(ix_test,iy_test,:)))
% 
% % to visualize just the seasonal cycle
% nexttile();plot(data_in.data_datenum, squeeze(data_in.data_backup(ix_test,iy_test,:)) - squeeze(data_in.data(ix_test,iy_test,:)))
% 
% nexttile();plot(squeeze(data_in.data_backup(ix_test,iy_test,end-366:end)) - ...
%     squeeze(data_in.data(ix_test,iy_test,end-366:end)))

% % to visualize just the mask
% nexttile();plot(data_in.data_datenum, squeeze(data_in.data(ix_test,iy_test,:)))
% hold on;plot(data_in.data_datenum(bfr_day_msk), squeeze(data_in.data(ix_test,iy_test,bfr_day_msk)), '.', 'linestyle', 'none')
%
% % to visualize the mask in a different way
% doy_day = day(data_in.data_datenum);
% ind_test = 1:length(data_in.data_datenum);
% nexttile();plot(ind_test, doy_day)
% hold on;plot(ind_test(bfr_day_msk), doy_day(bfr_day_msk), '.', 'linestyle', 'none')

% day_test = 28;
% month_test = 2;
% unique(data_in.data(ix_test,iy_test,day(data_in.data_datenum)==day_test & ...
%     month(data_in.data_datenum)==month_test))
% save('./test_30may_data_percentile3d.mat', 'data_percentile3d', '-v7.3')

% ciao

% >>> remove for daily data
if data_in.flag_dim_lon_lat_time && isfield(data_in, 'flag_daily_timeseries') && ~data_in.flag_daily_timeseries
    figure;set(gcf,'color','w','position',[0       0      1800        1200]);
    for mm=1:12
        nexttile()
        pcolor(data_in.x,data_in.y,data_percentile3d(:,:,mm)');shading flat;colorbar
        title([num2str(percentile) 'th percentile for month ' num2str(mm)])
        if ~isempty(data_in.cax_percentile)
            caxis(data_in.cax_percentile)
        end
        set(gca,'fontsize',18)
        
    end
    if ~isempty(fig_path)
        print('-dpng',[fig_path data_in.case_tag '.png'])
    end
end
%%%%%%%%%%
% find events, i.e.
% 1. create a 3d variable (start_tstep_msk) that indicates which timesteps
%    are the start of a MHW event at any location
% 2. create a 3d variable with the last tstep for the MHW at each location
%
% this should work for both delta_tstep>1 and delta_tstep == 1 (algorithm
% not yet checked/tested)
disp('>>> Finding MHW events and event peaks takes these number of seconds:')
tic;
%%%%%%%%%%%%%%%%%%%%%%
% mask for the start of the event
start_tstep_msk = false(size(data_in.data));
% mask for the event (true if a timestep is part of an event, regardless of
% the anomaly (at that timestep) being larger than the selected percentile)
tstep_part_of_event_msk = false(size(data_in.data));
% these below are values stored at the start tstep for MHW
end_tstep       = nan(size(data_in.data));
peak_tstep      = nan(size(data_in.data));
peak_value      = nan(size(data_in.data));
%%%%%%%%%%%%%%%%%%%%%%
% mask for the peak of the event
peak_tstep_msk  = false(size(data_in.data));
% these below are values stored at the peak of the event
end_tstep_stored_at_peak = reshape(nan(size(data_in.data)),[],1);
%end_tstep_stored_at_peak = nan(size(data_in.data));


ind3d           = permute(repmat(1:size(start_tstep_msk,3),...
    [size(start_tstep_msk,1) 1 size(start_tstep_msk,2)]),...
    [1 3 2]);

for i = 1:length(data_in.tstep)-delta_tstep+1
    % for tstep i we look for MHWs that start at tstep i
    for j = i+delta_tstep-1:length(data_in.tstep)
        clear bfr*
        % let's first check if there is a continous event
        bfr_msk= ...
            (sum(data_in.data(:,:,i:j)>data_percentile3d(:,:,i:j),3)==length(i:j)); %

        if j == i+delta_tstep-1
            if sum(bfr_msk(:)) == 0
                break
            end
        end
        % allow for gaps of data_in.num_gap_tsteps_in_event tsteps
        if j > i+delta_tstep-1 && data_in.num_gap_tsteps_in_event>0

            if j+data_in.num_gap_tsteps_in_event<=length(data_in.tstep)
                jjjend = j+data_in.num_gap_tsteps_in_event;
            else
                jjjend = length(data_in.tstep);
            end

            % bfr_ind3d = ind3d(:,:,i+delta_tstep:jjjend);
            % bfr_ind3d(~(bfr_ind3d-repmat(end_tstep(:,:,i+delta_tstep:jjjend),...
            %     [1 1 size(bfr_ind3d,3)])>=0)) = nan;
            % bfr_ind3d(bfr_ind3d>0)= 1;

            % bfr_ind3d = ind3d(:,:,i+delta_tstep:jjjend);
            % bfr_ind3d(~(bfr_ind3d-repmat(end_tstep(:,:,i),...
            %     [1 1 size(bfr_ind3d,3)])>=0)) = nan;

            % check if there is an event at all at that i AND
            % check that the time step of interest in the iteration (j) is
            % not too far from the last time step of the event (what "too
            % far" means depends on whether time step j is greater or not
            % than the percentile)
            bfr_msk_ievent = start_tstep_msk(:,:,i)==1;

            % bfr_j_gr_prc_j_end_tstep = data_in.data(:,:,j)>data_percentile3d(:,:,j) & ...
            %     j - end_tstep(:,:,i) <= data_in.num_gap_tsteps_in_event+1;
            % 
            % bfr_j_nogr_prc_j_end_tstep = ~(data_in.data(:,:,j)>data_percentile3d(:,:,j)) & ...
            %     j - end_tstep(:,:,i) <= data_in.num_gap_tsteps_in_event;


            clear bfr_ind3d
            bfr_ind3d = ind3d(:,:,j:jjjend);
            bfr_ind3d(~(data_in.data(:,:,j:jjjend)>data_percentile3d(:,:,j:jjjend))) = nan;
            
            bfr_msk = bfr_msk ...
                | ...
                ( ...
                bfr_msk_ievent & ...
                min(bfr_ind3d,[],3)-end_tstep(:,:,i) <= data_in.num_gap_tsteps_in_event+1 ...
                );


            % bfr_ind3d(~(bfr_ind3d-repmat(end_tstep(:,:,i+delta_tstep:jjjend),...
            %     [1 1 size(bfr_ind3d,3)])>=0)) = nan;
            % bfr_ind3d(bfr_ind3d>0)= 1;
            clear bfr_ind3d


            % bfr_ind3d(bfr_ind3d>0)= 1;

            % bfr_msk = bfr_msk ...
            %     | ...
            %     ( ...
            %     bfr_msk_ievent & ...
            %     bfr_j_gr_prc_j_end_tstep ...
            %     ) ...
            %     | ...
            %     ( ...
            %     bfr_msk_ievent & ...
            %     bfr_j_nogr_prc_j_end_tstep & ...
            %     sum(data_in.data(:,:,j:jjjend)>...
            %     data_percentile3d(:,:,j:jjjend),3)>=1 ...
            %     );

            % sum(data_in.data(:,:,i:j)>data_percentile3d(:,:,i:j),3)>= ...
            %     length(i:j) - data_in.num_gap_tsteps_in_event & ...

            % if j+data_in.num_gap_tsteps_in_event+1<=length(data_in.tstep)
            % 
            %     bfr_msk = bfr_msk ...
            %         | ...
            %         ( ...
            %         start_tstep_msk(:,:,i)==1 & ...
            %         ~(data_in.data(:,:,j)>data_percentile3d(:,:,j)) & ...
            %         sum(data_in.data(:,:,j:j+data_in.num_gap_tsteps_in_event)>...
            %         data_percentile3d(:,:,j:j+data_in.num_gap_tsteps_in_event),3)>=1 ...
            %         );
            % else
            %     bfr_msk = bfr_msk ...
            %         | ...
            %         ( ...
            %         start_tstep_msk(:,:,i)==1 & ...
            %         ~(data_in.data(:,:,j)>data_percentile3d(:,:,j)) & ...
            %         sum(data_in.data(:,:,j:end)>...
            %         data_percentile3d(:,:,j:end),3)>=1 ...
            %         );
            % end
        end
            
        % this should exclude tsteps from previous events (to not count the
        % same event more than once)
        if i>1 
            bfr_msk_before = ~tstep_part_of_event_msk(:,:,i-1);
            
            bfr_msk        = bfr_msk + bfr_msk_before == 2;
        end
        
        if sum(bfr_msk(:)) == 0
            break
        end
        
        if j == i+delta_tstep-1
            start_tstep_msk(:,:,i)= bfr_msk;
            % initialize end tstep
            bfr_end = nan(size(data_in.data(:,:,1)));
            
        else
            bfr_end = end_tstep(:,:,i);
        end

        bfr_end(bfr_msk) = j;
        
        end_tstep(:,:,i) = bfr_end;
        
    end
    clear bfr_ind3d
    bfr_ind3d = ind3d;
    bfr_ind3d(bfr_ind3d<i) = nan;
    bfr_ind3d(~(bfr_ind3d-repmat(end_tstep(:,:,i),[1 1 size(ind3d,3)])<=0)) = nan;
    bfr_ind3d(bfr_ind3d>0)= 1;
    
    [bfrMax,bfrMaxI] = max(data_in.data.*bfr_ind3d,[],3);
    bfrMaxI(isnan(bfrMax)) = nan;

    % fill in the mask where we have an event
    tstep_part_of_event_msk(bfr_ind3d==1) = true;
    
    peak_tstep(:,:,i) = bfrMaxI;
    peak_value(:,:,i) = bfrMax;
    
    % let's now create peak_tstep_msk
    peak_tstep_msk(ind3d==repmat(bfrMaxI,[1 1 size(ind3d,3)])) = true;
    
    %     %%%% attempt #1 to create end_tstep_stored_at_peak
    %     bfrbfr = end_tstep(:,:,i);
    %     end_tstep_stored_at_peak(ind3d==repmat(bfrMaxI,[1 1 size(ind3d,3)])) = ...
    %         bfrbfr(~isnan(bfrbfr));
    %%%% attempt #2 to create end_tstep_stored_at_peak
    bfrbfr1 = reshape(ind3d==repmat(bfrMaxI,[1 1 size(ind3d,3)]),[],1);
    bfrbfr2 = reshape(repmat(end_tstep(:,:,i),[1 1 size(ind3d,3)]),[],1);
    end_tstep_stored_at_peak(bfrbfr1) = bfrbfr2(bfrbfr1);
    
    
    if mod(i,round(length(data_in.tstep)-delta_tstep+1.*.2))==0
        disp([num2str(i/length(data_in.tstep)-delta_tstep+1.*100) '% completed'])
    end
end
toc;

data_out.peak_tstep_msk = peak_tstep_msk;
data_out.end_tstep_stored_at_peak = reshape(end_tstep_stored_at_peak,size(ind3d));
data_out.peak_tstep = peak_tstep;
data_out.peak_value = peak_value;

data_out.end_tstep = end_tstep;
data_out.start_tstep_msk = start_tstep_msk;
data_out.tstep_part_of_event_msk = tstep_part_of_event_msk;
data_out.data_percentile3d = data_percentile3d;

data_out.data_datenum = data_in.data_datenum;

if isfield(data_in, 'flag_daily_timeseries_smoothed_percentile') && data_in.flag_daily_timeseries_smoothed_percentile
    data_out.data_percentile3d_not_smooth = data_percentile3d_not_smooth;
end

data_out.data_used4MHWs = data_in.data;
data_out.data_seasonal = data_in.data_seasonal;

%%%%%%%%%% find more info about MHW events

% let's also save a mask that will allow a comparison with maps at https://psl.noaa.gov/marine-heatwaves/
% This mask is true if data are above the percentile selected at a certain
% (x,y) point and tstep
data_out.data_mhw_tstep_msk = data_in.data>=data_percentile3d;

% timesteps when a MHW event starts (created directly from tstep and
% masking out tsteps when a MHW event does not start)
data_out.start_tstep        = data_in.tstep;
data_out.start_tstep        = permute(repmat(data_out.start_tstep,[1, size(data_in.data,1), size(data_in.data,2)]),[2 3 1]);

data_out.start_tstep(isnan(end_tstep)) = nan;

% duration of each event (in tsteps): for each (x,y) point, we store the
% duration at each timestep when a MHW starts
data_out.events_duration_in_tsteps = end_tstep - data_out.start_tstep + 1;

% onset duration
data_out.onset_duration_in_tsteps = peak_tstep - data_out.start_tstep + 1;

% number of events over the all time period (at each point in space)
data_out.events_number             = nansum(start_tstep_msk,3);
return