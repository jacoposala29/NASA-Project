function params = set_params_NASAPO_2020_MHWs_OISST_MONTHLY(addpath_matlab_code,...
    num_gap_tsteps_in_event, years, user_path_OISST, output_dir)
% Initialize structure
params = struct();
% Set path where the code is
params.addpath_matlab_code = addpath_matlab_code;%'/Users/jacoposala/Downloads/';
% Set the range of years to analyze
params.years = years;%1992:2017;
% Num of tsteps allowed for a gap within an event
params.num_gap_tsteps_in_event = num_gap_tsteps_in_event;
% Create the fig_path dynamically using the extracted years
params.fig_path = sprintf([output_dir 'OISST_monthly_%d_%d/'], years(1), years(end));
% Flags for different data sources and options
params.flag_noaa_oisst = true;
% params.flag_argo = false;
% params.flag_ecco_zint = false;
% params.flag_ecco = false;
% params.flag_ecco_daily = false;
% params.flag_noaa_oisst_daily = false;
% params.flag_ecco_daily_box = false;
% params.flag_noaa_oisst_daily_box = false;
% params.flag_fake_data = false;
% Time step delta
params.delta_tstep = 1;
% Percentile for detection
params.percentile = 90;
% Flags for data processing options
params.flag_remove_trend = true;
params.flag_save_MHW_info = true;
params.flag_save_MHW_ave_from_info = true;
% File name
params.fname = 'sst.mon.mean.nc';
% Case tag
params.case_tag = 'oisst_v2';
% User path for OISST data
params.user_path_OISST = user_path_OISST;%'/Users/jacoposala/Desktop/CU/3.RESEARCH/NASA_project/OISSTv2/DATA/';
% Flag for monthly time series
params.flag_monthly_timeseries = true;
% Flag for daily time series
params.flag_daily_timeseries = false;
% Delta days for percentile
params.delta_days_4_percentile = 5;
% Color axis limits for percentile
params.cax_percentile = [0 6];
% Color axis limits for number of events
params.cax_num_events = [0 30];
% Color axis limits for average duration
params.cax_aveDuration = [0 18];
% Color axis limits for data
params.cax_data = [-6 6];
% Flag for dimension (longitude, latitude, time)
params.flag_dim_lon_lat_time = true;
% Indices for example cases
params.i_all = [818];
params.j_all = [218];
% Flags for map summary and maps at each time step
params.flag_maps_summary = true;
params.flag_maps_at_each_tstep = false;
% Call the main function with the parameters
NASAPO_2020_MHWs_detect_and_info_main(params);

