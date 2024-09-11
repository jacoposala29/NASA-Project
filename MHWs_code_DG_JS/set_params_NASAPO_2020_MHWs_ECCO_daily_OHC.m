function params = set_params_NASAPO_2020_MHWs_ECCO_daily(addpath_matlab_code,...
    num_gap_tsteps_in_event, years, tag_file, year_start_fname, year_end_fname, ...
    num_days, user_path_ecco, output_dir)

% Initialize structure
params = struct();
% Set path where the code is
params.addpath_matlab_code = addpath_matlab_code;%'/Users/jacoposala/Downloads/';
% Tag for days between each event
params.num_gap_tsteps_in_event = num_gap_tsteps_in_event;%0; % 0, 2
% Set the range of years to analyze
params.year_start_fname = year_start_fname;% 1992;
params.year_end_fname = year_end_fname; %1994;
params.years = years;%1992:1994;
% Set the number of days
params.num_days = num_days; %730;
% Create the fig_path dynamically using the extracted years
params.fig_path = sprintf([output_dir 'ECCO_daily_ohc_k0_k5_%d_%d/'], ...
    params.year_start_fname, params.year_end_fname);
% Set the z-level
params.zlev = [];
% Flags for different data sources and options
% params.flag_noaa_oisst = false;
% params.flag_argo = false;
% params.flag_ecco_zint = false;
% params.flag_ecco = false;
params.flag_ecco_daily = true;
% params.flag_noaa_oisst_daily = false;
% params.flag_ecco_daily_box = false;
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
params.user_path_ecco = user_path_ecco; %'/Users/jacoposala/Downloads/ecco_daily_input_files_1992_1994/';
% Filename years
% params.tag_years_filename = '1992_1994';
% Metadata years
params.tag_years_metadata = '1993_2017';
% tag_file
params.tag_file = tag_file;%'_cut';
% ECCO data set name
params.ecco_name = 'ECCOv4r4';
% Flag for daily time series
params.flag_daily_timeseries = true;
% Flag for monthly time series
params.flag_monthly_timeseries = false;
% Number of tiles
params.num_of_tiles = 13;
% Number of points in tiles
params.ecco_num_of_points_in_tiles = 90;
% Number of levels in monthly ECCO data
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
params.i_all = [6586 5548];
params.j_all = [9 9];
% Flags for map summary and maps at each time step
params.flag_maps_summary = true; %true
params.flag_maps_at_each_tstep = false;
% Call the main function with the parameters
NASAPO_2020_MHWs_detect_and_info_main(params);


