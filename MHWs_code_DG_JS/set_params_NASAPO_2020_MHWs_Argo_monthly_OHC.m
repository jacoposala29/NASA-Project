function params = set_params_NASAPO_2020_MHWs_Argo_monthly_OHC(addpath_matlab_code,...
    num_gap_tsteps_in_event, years, output_dir)
% Initialize structure
params = struct();
% Set path where the code is
params.addpath_matlab_code = addpath_matlab_code;%'/Users/jacoposala/Downloads/';
% Set the range of years to analyze
params.years = years;%2004:2016;
% Create the fig_path dynamically using the extracted years
params.fig_path = sprintf([output_dir 'Argo_monthly_ohc_%d_%d/'], ...
    params.years(1), params.years(end));
% Flags for different data sources and options
% params.flag_noaa_oisst = false;
params.flag_argo = true;
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
% Num of tsteps allowed for a gap within an event
params.num_gap_tsteps_in_event = num_gap_tsteps_in_event;
% Flags for data processing options
params.flag_remove_trend = true;
params.flag_save_MHW_info = true;
params.flag_save_MHW_ave_from_info = true;
% lev tag
params.plev_tag = '15_50';
% Case tag
params.case_tag = ['argo_ohc' params.plev_tag];
% User path for Argo data
params.fpath = '/Users/jacoposala/Desktop/CU/3.RESEARCH/NASA_project/NEW_heatBudgetECCO/data/GCOS_Global_2004_2022_15_50/Results/FullField/';
% Argo filename prefix
params.fname_prefix = 'intTempFullFieldSpaceTimeTrend';
% Grid in input - only needed in Argo as lon/lat are not saved in file
params.xgrid_input = [20.5:379.5]'; 
params.ygrid_input = [-89.5:89.5]';
% Params for OHC conversion - do not change, unless really needed
params.cp0 = 3989.244; % as in McDougall 2003
params.rho0 = 1030; %1025;%
% Flag for monthly time series
params.flag_monthly_timeseries = true;
% Flag for daily time series
params.flag_daily_timeseries = false;
% Color axis limits for percentile
params.cax_percentile = [0 6].*1e8;
% Color axis limits for number of events
params.cax_num_events = [0 20];
% Color axis limits for average duration
params.cax_aveDuration = [0 18];
% Color axis limits for data
params.cax_data = [-6 6].*1e8;
% Flag for dimension (longitude, latitude, time)
params.flag_dim_lon_lat_time = true;
% Indices for example cases
params.i_all = [184];
params.j_all = [55];
% Flags for map summary and maps at each time step
params.flag_maps_summary = true;
params.flag_maps_at_each_tstep = false;
% Call the main function with the parameters
% NASAPO_2020_MHWs_detect_and_info_main(params);

