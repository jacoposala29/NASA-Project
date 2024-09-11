function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_argo(input_params)

switch input_params.fname_prefix
    case 'intTempFullFieldSpaceTimeTrend'
        data_out.case_tag = ['argo_ohc' input_params.plev_tag];
        cp0                      = input_params.cp0;%3989.244; % as in McDougall 2003
        rho0                     = input_params.rho0;%1030;%1025;%

        factor                   = cp0*rho0;
    otherwise
        checkandwrite
end
        
data_out.flag_dim_lon_lat_time = input_params.flag_dim_lon_lat_time;%true;
data_out.num_gap_tsteps_in_event = input_params.num_gap_tsteps_in_event;

%fname = 'sst.mon.mean.nc';

data_out.flag_monthly_timeseries = input_params.flag_monthly_timeseries;%true;
data_out.flag_daily_timeseries   = input_params.flag_daily_timeseries;%false;

% data_in      = read_vars_in_ncfile([fpath fname]);
% time_datenum = data_in.time+datenum(1800,1,1);
% 
% msk_years    = year(time_datenum)>=min(years) & year(time_datenum)<=max(years);

data_out.x = input_params.xgrid_input;%[20.5:379.5]';
data_out.y = input_params.ygrid_input;%[-89.5:89.5]';

[data_out.X,data_out.Y] = ndgrid(data_out.x,data_out.y);

data_out.maplonlimit = [min(data_out.x) max(data_out.x)];

data_out.data = nan(size(data_out.X,1),size(data_out.X,2),length(input_params.years)*12);

n = 0;
for iyr=1:length(input_params.years)
    for mm=1:12
        
        clear bfr
        bfr = load([input_params.fpath input_params.fname_prefix '_' input_params.plev_tag '_' num2str(mm,'%02d') ...
            '_' num2str(input_params.years(iyr))]);
        n = n + 1;
        data_out.data(:,:,n) = bfr.fullFieldGrid.*factor;

        data_out.data_datenum(n,1) = datenum(input_params.years(iyr),mm,15);
        
    end
end

% define indices to look at in the timeseries plot at the end
data_out.i_all=input_params.i_all;%[184]; %[180 170];
data_out.j_all=input_params.j_all;%[55]; %[110 120];

% caxis for percentile plot
data_out.cax_percentile = input_params.cax_percentile;%[0 6].*1e8;%[0 6];

% variables to average over MHW event-relevant timesteps (e.g. over all the
% timesteps for an event, over the final ones only...)
data_out.data2ave_ALL_from_find_MHWs_info = {'data_used4MHWs'};


return
