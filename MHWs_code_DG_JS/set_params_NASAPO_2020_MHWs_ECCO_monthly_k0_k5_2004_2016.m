function params = set_params_NASAPO_2020_MHWs_ECCO_monthly_k0_k5_2004_2016(addpath_matlab_code,...
    num_gap_tsteps_in_event,output_dir)

% _ohc_to50m_k0_k5_1992_2017

% addpath_matlab_code     = '/Users/dogi7244/Downloads/MHWs_code_DG/';
% num_gap_tsteps_in_event = 0;
% output_dir              = '/Users/dogi7244/Downloads/MHWs_data_ouput/';



% Initialize structure
params = struct();
% Flag use only DATA
params.flag_use_only_DATA = true;
% Set path where the code is
params.addpath_matlab_code = addpath_matlab_code; % DG9Aug '/Users/jacoposala/Downloads/';
% Num of tsteps allowed for a gap within an event
params.num_gap_tsteps_in_event = num_gap_tsteps_in_event; % DG9Aug 2; % 0, 2

% Set the range of years to analyze
params.years = 2004:2016;
% Set the path for the figure folder
% params.fig_path = '/Users/jacoposala/Downloads/test_ECCO_monthly_Aug1/';
% Create the fig_path dynamically using the extracted years
% DG9Aug (see next line)
params.fig_path = sprintf([output_dir 'ECCO_monthly_heat_ohc_to50m_zint_k0_k5_%d_%d/'],...
    params.years(1), params.years(end));
% Set the z-level
params.zlev = 1;%[];
% Flags for different data sources and options
% params.flag_noaa_oisst = false;
% params.flag_argo = false;
params.flag_ecco_zint = true;
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
% Tag file 
params.tag_file = '_zint_k0_k5';
% User path for ECCO data
params.user_path_ecco = '/Users/jacoposala/Desktop/CU/3.RESEARCH/NASA_project/heatBudgetECCO_jupyter/nc_files_zlev_or_zint_monthly/ohc_k0_k5/';
% Filename years
params.year_start_filename = 1993;
params.year_end_filename = 2016;
% Years ACTUALLY in the file
params.year_start_in_file = 1993;
params.year_end_in_file = 2016;
% ECCO data set name
params.ecco_name = 'ECCOv4r4';
% Metadata years
params.tag_years_metadata = '1993_2017';
% Flag for daily time series
params.flag_daily_timeseries = false;
% Flag for monthly time series
params.flag_monthly_timeseries = true;
% Number of tiles
params.num_of_tiles = 13;
% Number of points in tiles
params.ecco_num_of_points_in_tiles = 90;
% Number of levels in monthly ECCO data
params.ecco_monthly_num_of_levs = 1;
% ECCO monthly tag field
params.ECCO_monthly_tag_field = '';
% Field tag
params.tag_field = 'ohc_to50m';
% Color axis limits for percentile
params.cax_percentile = [0 6];
% Color axis limits for number of events
params.cax_num_events = [0 30];
% Color axis limits for average duration
params.cax_aveDuration = [0 18];
% Color axis limits for data
params.cax_data = [-6 6];
% Flag for dimension (longitude, latitude, time)
params.flag_dim_lon_lat_time = false;
% Indices for example cases
params.i_all = [5548];
params.j_all = [9];
% Flags for map summary and maps at each time step
params.flag_maps_summary = true;
params.flag_maps_at_each_tstep = false;
% Call the main function with the parameters
% NASAPO_2020_MHWs_detect_and_info_main(params); % DG9Aug
return

