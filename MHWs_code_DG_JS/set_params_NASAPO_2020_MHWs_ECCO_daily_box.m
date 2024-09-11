function params = set_params_NASAPO_2020_MHWs_ECCO_daily_box(addpath_matlab_code,...
    num_gap_tsteps_in_event,years, tag_file,tag_years_filename,user_path_ecco,output_dir)
% addpath_matlab_code     = '/Users/dogi7244/Downloads/MHWs_code_DG/';
% num_gap_tsteps_in_event = 2;
% years                   = 1992:2017;
% tag_file                = '_avg_box_NEP_OHC';
% tag_years_filename      = '1992_2017';
% user_path_ecco          = '/Users/dogi7244/Downloads/MHWs_data/ECCO_daily_NEP_avg_box/';
% output_dir              = '/Users/dogi7244/Downloads/MHWs_data_ouput/';

% Initialize structure
params = struct();
% Set path where the code is
params.addpath_matlab_code = addpath_matlab_code; % DG9Aug
% Num of tsteps allowed for a gap within an event
params.num_gap_tsteps_in_event = num_gap_tsteps_in_event; % DG9Aug  % 0, 2
% Set the range of years to analyze
params.years = years; % DG9Aug 1992:2017;
% Tag file 
params.tag_file = tag_file; % DG9Aug  '_avg_box_NEP_OHC';
% DG9Aug (see next line)
params.fig_path = sprintf([output_dir 'ECCO_daily%s_%d_%d/'], ...
    params.tag_file, params.years(1), params.years(end));
% Set the z-level
params.zlev = [];
% Flags for different data sources and options
% params.flag_noaa_oisst = false;
% params.flag_argo = false;
% params.flag_ecco_zint = false;
% params.flag_ecco = false;
% params.flag_ecco_daily = false;
% params.flag_noaa_oisst_daily = false;
params.flag_ecco_daily_box = true;
% params.flag_noaa_oisst_daily_box = false;
% params.flag_fake_data = false;
% Time step delta
params.delta_tstep = 5;
% Percentile for detection
params.percentile = 90;
% Flags for data processing options
params.flag_remove_trend = true;
params.flag_save_MHW_info = true;
params.flag_save_MHW_ave_from_info = true;
% User path for ECCO data
params.user_path_ecco = user_path_ecco; % DG9Aug
% params.user_path_ecco = '/Users/jacoposala/Downloads/ECCO_daily_region_test_DG/';
% Filename years
params.tag_years_filename = tag_years_filename; % DG9Aug '1992_2017';
% ECCO data set name
params.ecco_name = 'ECCOv4r4';
% Flag for daily time series
params.flag_daily_timeseries = true;
% Flag for monthly time series
params.flag_monthly_timeseries = false;
% Flag for percentile
params.flag_daily_timeseries_smoothed_percentile = true;
% % Number of tiles
% params.num_of_tiles = 13;
% % Number of points in tiles
% params.ecco_num_of_points_in_tiles = 90;
% % Number of levels in monthly ECCO data
% params.ecco_monthly_num_of_levs = 1;
% Field tag
params.tag_field = 'heat';
% Delta days for percentile
params.delta_days_4_percentile = 5;
% Color axis limits for percentile
params.cax_percentile = [0 6];
% Color axis limits for number of events
params.cax_num_events = [0 30];
% Color axis limits for average duration
params.cax_aveDuration = [0 18];
% Color axis limits for data
params.cax_data = [-3 3];
% Flag for dimension (longitude, latitude, time)
params.flag_dim_lon_lat_time = false;
% Indices for example cases % 156.5W (around 204 east), 35.1757S
params.i_all = [1 ];
params.j_all = [1 ];
% Flags for map summary and maps at each time step
params.flag_maps_summary = true;
params.flag_maps_at_each_tstep = false;
% Call the main function with the parameters
% NASAPO_2020_MHWs_detect_and_info_main(params); % DG9Aug

