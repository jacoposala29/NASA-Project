function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily(input_params)
%July29
% each time fpath is used, it will be replaced with input_params.fpath
% fpath = '/scratch/alpine/jasa1084/inputs_MHW_project_blanca/oisst/'; % /scratch/alpine/jasa1084/inputs_MHW_project_blanca/oisst/

%July29
% each time you store a var in data_out, take it directly from
% input_params

data_out.case_tag = input_params.case_tag;%'oisst_v2'; 
data_out.flag_dim_lon_lat_time = input_params.flag_dim_lon_lat_time;%true;

% fname = 'sst.day.mean.';

data_out.flag_monthly_timeseries = input_params.flag_monthly_timeseries;%false;
data_out.flag_daily_timeseries   = input_params.flag_daily_timeseries;%true;
data_out.num_gap_tsteps_in_event = input_params.num_gap_tsteps_in_event;
years_files = min(input_params.years):max(input_params.years);

for iyr=1:length(years_files)
    clear bfr
%     bfr = read_vars_in_ncfile([fpath fname num2str(years_files(iyr)) '.nc']);
%     bfr = read_vars_in_ncfile([fpath fname num2str(years_files(iyr)) '_subset_lonmin209.5_lonmax225.5_latmin_39.5_latmax50.5.nc']);
%     bfr = read_vars_in_ncfile([fpath fname num2str(years_files(iyr)) '_subset_lonmin209.5_lonmax215.5_latmin_39.5_latmax43.5.nc']);
    bfr = read_vars_in_ncfile([input_params.fpath_oisst_daily input_params.fname num2str(years_files(iyr)) '_rcon.nc']);
    bfr.sst(bfr.sst==bfr.sst_missing_value) = nan;
    if iyr ==1
        data_in = bfr;
    else
        data_in.sst  = cat(3,data_in.sst,bfr.sst);
        data_in.time = cat(1,data_in.time,bfr.time);
    end
end
    
time_datenum = data_in.time+datenum(1800,1,1);

% msk_years    = year(time_datenum)>=min(years) & year(time_datenum)<=max(years);
msk_years    = year(time_datenum)>=min(input_params.years) & ...
    year(time_datenum)<=max(input_params.years);

% select data in the years of interest
data         = data_in.sst(:,:,msk_years);
data(data<=ceil(data_in.sst_missing_value)) = nan;
data_out.data_datenum = time_datenum(msk_years);

data_out.x = squeeze(data_in.lon);
data_out.y = squeeze(data_in.lat);

% convert longitudes 20-380
bfr_x = data_out.x;
bfr_x(bfr_x<20) = bfr_x(bfr_x<20) + 360;
[data_out.x,I] = sort(bfr_x);

% arrange based on the new longitudes
data_out.data = data(I,:,:);

data_out.maplonlimit = [min(data_out.x) max(data_out.x)];

[data_out.X,data_out.Y] = ndgrid(data_out.x,data_out.y);

% define indices to look at in the timeseries plot at the end
data_out.i_all=[818]; %[180 170];
data_out.j_all=[218]; %[110 120];

% caxis for percentile plot
data_out.cax_percentile = [0 6];
% delta days before/after to calculate percentile - should be the same for
% ECCO daily
data_out.delta_days_4_percentile = 5;

% variables to average over MHW event-relevant timesteps (e.g. over all the
% timesteps for an event, over the final ones only...)
data_out.data2ave_ALL_from_find_MHWs_info = {'data_used4MHWs'}; % keep it here


return
