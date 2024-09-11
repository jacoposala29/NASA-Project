function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily_box(input_params)%input_params years,tag_file)

% load params (this needs to be included in all the scripts that use input
% params
% load('./NASAPO_2020_MHWs_input_params.mat','input_params');

% fpath = input_params.fpath_oisst_daily_box; % /scratch/alpine/jasa1084/inputs_MHW_project_blanca/oisst/

data_out.case_tag = input_params.case_tag;%'oisst_v2';
data_out.flag_dim_lon_lat_time = input_params.flag_dim_lon_lat_time;%false;

% fname = 'OISST_1993_2017';

data_out.flag_monthly_timeseries = input_params.flag_monthly_timeseries;%false;
data_out.flag_daily_timeseries   = input_params.flag_daily_timeseries;%true;
data_out.flag_daily_timeseries_smoothed_percentile = input_params.flag_daily_timeseries_smoothed_percentile;
data_out.num_gap_tsteps_in_event = input_params.num_gap_tsteps_in_event;


years_files = min(input_params.years):max(input_params.years);

bfr = read_vars_in_ncfile([input_params.fpath_oisst_daily_box input_params.fname input_params.tag_file '.nc']);
bfr.SST(bfr.SST==bfr.SST__FillValue) = nan;
    
data_in.sst(1,1,:) = bfr.SST;
data_in.time = bfr.time;

if strcmp(bfr.time_units,'days since 1992-01-01 00:00:00')
    time_datenum = data_in.time+datenum(1992,1,1);
else
    check_time_units
end

msk_years    = year(time_datenum)>=min(input_params.years) & year(time_datenum)<=max(input_params.years);

% select data in the years of interest
data_out.data         = data_in.sst(:,:,msk_years);
data_out.data_datenum = time_datenum(msk_years);

% data_out.x = squeeze(data_in.lon);
% data_out.y = squeeze(data_in.lat);
% 
% % convert longitudes 20-380
% bfr_x = data_out.x;
% bfr_x(bfr_x<20) = bfr_x(bfr_x<20) + 360;
% [data_out.x,I] = sort(bfr_x);
% 
% % arrange based on the new longitudes
% data_out.data = data(I,:,:);
% 
% data_out.maplonlimit = [min(data_out.x) max(data_out.x)];
% 
% [data_out.X,data_out.Y] = ndgrid(data_out.x,data_out.y);

% define indices to look at in the timeseries plot at the end
data_out.i_all= input_params.i_all;%[1]; %[180 170];
data_out.j_all=input_params.j_all;%[1]; %[110 120];

% % caxis for percentile plot
% data_out.cax_percentile = [0 6];
% delta days before/after to calculate percentile - should be the same for
% ECCO daily
data_out.delta_days_4_percentile = input_params.delta_days_4_percentile;%5;

data_out.data2ave_ALL_from_data_in = {};

% variables to average over MHW event-relevant timesteps (e.g. over all the
% timesteps for an event, over the final ones only...)
data_out.data2ave_ALL_from_find_MHWs_info = {'data_used4MHWs'};


return
