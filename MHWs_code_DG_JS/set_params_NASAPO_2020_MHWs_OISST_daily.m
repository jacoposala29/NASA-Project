function params = set_params_NASAPO_2020_MHWs_OISST_daily(addpath_matlab_code,...
    num_gap_tsteps_in_event, years, tag_file,input_fname,fpath_oisst_daily,output_dir)
% Initialize structure
params = struct();
% Set path where the code is
params.addpath_matlab_code = addpath_matlab_code; % DG9Aug '/Users/dogi7244/Downloads/MHWs_code_DG/';
% Num of tsteps allowed for a gap within an event
params.num_gap_tsteps_in_event = num_gap_tsteps_in_event; % DG9Aug 2; % 0, 2
% Set the range of years to analyze
params.years = years; % DG9Aug 1993:2016;
% Tag file 
params.tag_file = tag_file; % DG9Aug ''
% cleaned_tag_file = strrep(params.tag_file, '_', '');
% Create the fig_path dynamically using the extracted years
% DG9Aug (see next line)
params.fig_path = sprintf([output_dir 'OISST_daily%s_%d_%d/'], ...
    params.tag_file, params.years(1), params.years(end));
disp(params.fig_path)
% Flags for different data sources and options
% params.flag_noaa_oisst_daily_box = false;
% params.flag_noaa_oisst = false;
% params.flag_argo = false;
% params.flag_ecco_zint = false;
% params.flag_ecco = false;
% params.flag_ecco_daily = false;
params.flag_noaa_oisst_daily = true;
% params.flag_ecco_daily_box = false;
% params.flag_fake_data = false;
% Time step delta (minimum lenght of MHWs in tsteps; to not be confused with 
% delta_days_4_percentile)
params.delta_tstep = 5;
% Percentile for detection
params.percentile = 90;
% Flags for data processing options
params.flag_remove_trend = true;
params.flag_save_MHW_info = true;
% flag_save_MHW_ave_from_info is False only to read a file previously saved
% instead of finding the MHWs using the code
params.flag_save_MHW_ave_from_info = true;
% zlev
params.zlev = [];
% File name
params.fname = input_fname; % DG9Aug 'sst.day.mean.';
% Case tag
params.case_tag = 'oisst_v2';
% User path for OISST data
params.fpath_oisst_daily = fpath_oisst_daily; % DG9Aug '/Users/dogi7244/Downloads/MHWs_data/OISST/';
% Flag for monthly time series
params.flag_monthly_timeseries = false;
% Flag for daily time series
params.flag_daily_timeseries = true;
% Flag for smoothed percentile (only for daily)
params.flag_daily_timeseries_smoothed_percentile = true;
% Delta days for percentile (used for calculation of seas cycle and
% percentile (to not be confused with delta_tstep above))
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
params.i_all = [1];
params.j_all = [1];
% Flags for map summary and maps at each time step
params.flag_maps_summary = true;
params.flag_maps_at_each_tstep = false;
% Call the main function with the parameters
% NASAPO_2020_MHWs_detect_and_info_main(params); % DG9Aug
return

