function NASAPO_2020_MHWs_detect_and_info_main(input_params_at_start)
%
% README on how to run the code:
% if you are OK with the default parameters in
% NASAPO_2020_MHWs_detect_and_info_main_input.m (these were set for the
% ECCO daily box case), you could techinically run
% the code with just NASAPO_2020_MHWs_detect_and_info_main(struct([]))
%
% if you want to run the code with params of interest, you need to specify
% all the params needed and then call this function using input_params as input:
% NASAPO_2020_MHWs_detect_and_info_main(input_params);

% for an example on how to creat input_params, please see
% NASAPO_2020_MHWs_detect_and_info_main_input.m




% CHECK case_tag line 313



% clear all
close all
set(0, 'DefaultFigureVisible', 'off')

%%%
input_params = NASAPO_2020_MHWs_detect_and_info_main_input(input_params_at_start);

% Create automatically the folder - if not already created
if ~exist(input_params.fig_path, 'dir')
    mkdir(input_params.fig_path)
end
disp(['output will be in ' input_params.fig_path])

% load params (this needs to be included in all the scripts that use input
% params
% load('./NASAPO_2020_MHWs_input_params.mat','input_params');
%%%

addpath(input_params.addpath_matlab_code)
% fig_path = '/Users/dgiglio/Work/projects/NASAPO2020_OHC/Figures/MHWs_events/';
%%% July29
%%% we do not define any intermediate variable in this script; we use
%%% directly what we have in input_params, i.e. we remove the definition of
%%% fig_path and use input_params.fig_path
% fig_path = input_params.fig_path;


% select years for which all the months are available
% years = 2004:2016; % all (ecco daily goes to 2017)
%%% July29 1993:2016 will be input_params.years and it will be given in
%%% input by the user
% years = 1993:2016; % all (ecco daily goes to 2017)

%years = 2004:2016; % all
%years = 1993:2016; % sst, ecco
%years = 2004:2022; % sst, argo
%years = 1993:2022; % sst only


% flag_noaa_oisst= false; % true; %
% flag_argo      = false;
% flag_ecco_zint = false;
% flag_ecco      = false; % false; %
% flag_ecco_daily= true;
% flag_noaa_oisst= false;
% flag_argo      = false;
% flag_ecco_zint = false;
% flag_ecco      = false;
% flag_ecco_daily= false;
% flag_noaa_oisst_daily= false;
%
% flag_ecco_daily_box= false; % created bash 17 Jun
% flag_noaa_oisst_daily_box = true; % created bash 19 Jun
%
% % flag_fake_data = false;
% flag_fake_data = false;

% set minimum duration of MHW event and selected percentile
% delta_tstep = 1;
% delta_tstep = 1;
% if flag_ecco_daily
%    delta_tstep = 5;
% end
% if flag_ecco_daily || flag_noaa_oisst_daily || flag_noaa_oisst_daily
%     delta_tstep = 5;
% end
% percentile  = 90;
% percentile = 90;
%%%%% flags to run some parts of the code or not
% % flag to save .mat files (now there are two of them, just to isolate the
% % first step while testing... this could be changed back in the future)
% flag_save_MHW_info = true; % false;
% flag_save_MHW_ave_from_info = true; % false;
%
% % flag for summary maps
% flag_maps_summary  = true;

% flag_save_MHW_info = true;
% flag_save_MHW_info = true;
% flag_save_MHW_ave_from_info = true;
% flag_save_MHW_ave_from_info = true;
% flag_maps_summary  = true;
% flag_maps_summary  = true;

% flag for monthly maps
% flag_maps_at_each_tstep = false; %true;%
% flag_maps_at_each_tstep = false; %true;%

% ECCO_monthly_tag_field = 'bash_ECCO_monthly_tag_field'; %%%%%%% CHECK

%if flag_ecco_daily
%   fig_path = [fig_path  'daily/'];
%end


% set param for plots
% cax_aveDuration = [0 18];
% if flag_ecco_daily || flag_noaa_oisst_daily || flag_ecco_daily_box || flag_noaa_oisst_daily_box
%     cax_aveDuration = [0 60];
% end
%
% switch num2str([years(1) years(end)])
%     case {num2str([2004 2016]),num2str([2004 2022])}
%         cax_num_events  = [0 20];
%         if flag_ecco_daily || flag_noaa_oisst_daily || flag_ecco_daily_box || flag_noaa_oisst_daily_box
%             cax_num_events = [0 60];
%         end
%     case {num2str([1993 2016]),num2str([1993 2022])}
%         cax_num_events  = [0 30];
%     otherwise
%         cax_num_events  = [0 20];
% end

%%%%%

% flag to remove trend
% flag_remove_trend = true;

%%%% set other param in the following if statement
% July29 change the if statement for each product
if isfield(input_params, 'flag_noaa_oisst') && input_params.flag_noaa_oisst
    %%% July29: to each function that load the input data, we give only
    %%% input_params
    data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst(input_params);
    %     data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst(years);
    %     data_in.cax_data        = [-6 6];%[0 6];
elseif isfield(input_params, 'flag_noaa_oisst_daily') && input_params.flag_noaa_oisst_daily
    data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily(input_params);
    %     data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily(years);
    %     data_in.cax_data        = [-6 6];%[0 6];
    
elseif isfield(input_params, 'flag_argo') && input_params.flag_argo
    % index to select along vertical dim
    %     plev_tag = '15_50';% layers available are '15_50' '50_150' '150_300'
    % if plev_tag changes, data_in.cax_data needs changing
    
    data_in  = NASAPO_2020_MHWs_detect_and_info_appendix_argo(input_params);
    %     data_in  = NASAPO_2020_MHWs_detect_and_info_appendix_argo(years,plev_tag);
    %     data_in.cax_data        = [-6 6].*1e8;%[0 6].*1e8;
    
    %data_in.case_tag = [data_in.case_tag '_plev' plev_tag];
    
    % elseif flag_ecco_daily
    %     zlev    = 0;
    %     tag_file= '_cut';
    %     data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily(years,zlev,'heat',tag_file);
    %
elseif isfield(input_params, 'flag_ecco') && input_params.flag_ecco || ...
        isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily || ...
        isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box || ...
        isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box || ...
        isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint
    %     zlev = input_params.zlev; % should be [] for the files with vertical integral
    %     tag_file = input_params.tag_file; % '_zint_to200m' for version with vertical integral _zint_k1_k5
    
    if isfield(input_params, 'flag_ecco') && input_params.flag_ecco
        data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco(input_params); %years,zlev,ECCO_monthly_tag_field,tag_file);
    end
    
    if isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily
        data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily(input_params);
    end
    
    if isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box
        data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily_box(input_params);
    end
    
    if isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box
        data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily_box(input_params);
    end
    
    if isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint
        data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco(input_params);
    end
    
    
    
    
    % elseif flag_ecco || flag_ecco_daily || flag_ecco_daily_box || flag_noaa_oisst_daily_box
    %     if flag_ecco
    %         % index to select along vertical dim
    %         %         zlev = 1;% should be [] for the files with vertical integral
    %         zlev = 1;% should be [] for the files with vertical integral
    %         tag_file = '';% '_zint_to200m' for version with vertical intergral
    %         data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco(years,zlev,ECCO_monthly_tag_field,tag_file);
    %     elseif flag_ecco_daily
    %         %         %%%%%%%% for ONE level
    %         %         zlev    = 0;
    %         %         tag_file= '_cut';
    %         %%%%%%%% for OHC
    %         %         zlev = [];
    %         zlev = 0;
    %         %         tag_file= '_ohc_k0_k5';
    %         tag_file= '_cut';
    %         data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily(years,zlev,'heat',tag_file);
    %     elseif flag_ecco_daily_box
    %         zlev = [];
    %         tag_file = '_avg_box_NEP'; % create bash 17 Jun
    %         data_in = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily_box(years,'heat',tag_file);
    %     elseif flag_noaa_oisst_daily_box
    %         zlev = [];
    %         tag_file = '_avg_box_TASMAN'; % create bash 19 Jun
    %         %%% July29
    %         data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily_box(years,tag_file);
    % %         data_in = NASAPO_2020_MHWs_detect_and_info_appendix_oisst_daily_box(years,tag_file);
    %     end
    
    
    %     data_in.cax_data        = [-6 6]; %[0 6];
    
    if isfield(input_params, 'zlev')
        zlev_cbar_buffer = input_params.zlev;
    end
   
    if ~isfield(input_params, 'zlev') || isempty(input_params.zlev)
        % caxis to use for ohc with daily data are for now the same as for lev 1 for daily
        % data
        zlev_cbar_buffer = 1;
    end
    % July29 ***
    switch zlev_cbar_buffer
        case {0,1}
            data_in.cax_data_heat_budget_terms = {...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-5 ...
                [-4.5 4.5].*1e-5 [-4.5 4.5].*1e-5 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-7 };
            
            if (isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily) || ...
                    (isfield(input_params, 'flag_noaa_oisst_daily') && input_params.flag_noaa_oisst_daily)
                data_in.cax_data_heat_budget_terms = {...
                    [-4.5 4.5].*1e-7 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-7 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6...
                    };%[-4.5 4.5].*1e-6 [-4.5 4.5].*1e-7
            else
                data_in.cax_data_heat_budget_terms_monthly = {[-6 6] ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                    };
            end
        case {5,9} % this should be the zlev corresponding to 95 m
            %
            % edit values here for the summary plots:
            data_in.cax_data_heat_budget_terms = {...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-5 ...
                [-4.5 4.5].*1e-5 [-4.5 4.5].*1e-5 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-7 };
            %
            % edit values here for the monthly plots:
            data_in.cax_data_heat_budget_terms_monthly = {[-6 6] ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                [-4.5 4.5].*1e-6 [-4.5 4.5].*1e-6 ...
                };
    end
    %%% *******************************************************************
    %data_in.case_tag = [data_in.case_tag '_zlev' num2str(zlev,'%02d')];
    
    
    %     data_in.cax_data        = [-180 180];
    %%% *******************************************************************
    %%%% COMMENTED ON July 31 - TO CHECK
    %     data_in.cax_data_heat_budget_terms = { ...
    %         [-6 6].*1e-5 [-6 6].*1e-4 ...
    %         [-6 6].*1e-4 [-6 6].*1e-4 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5...
    %         };
    %
    %     data_in.cax_data_heat_budget_terms_monthly = {[-180 180] ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         [-6 6].*1e-5 [-6 6].*1e-5 ...
    %         };
    %%% *******************************************************************
    %data_in.case_tag = [data_in.case_tag tag_file];
    % elseif flag_fake_data
    %     data_in.data_datenum = [1:12*20]';
    %     data_in.data  = randn(4,3,length(data_in.data_datenum));
    %     data_in.flag_monthly_timeseries = true;
    %     data_in.cax_percentile = [0 6];
end

% set some param to plot
data_in.cax_aveDuration = input_params.cax_aveDuration; %[0 18];
data_in.cax_num_events  = input_params.cax_num_events; %[0 30];
data_in.cax_data        = input_params.cax_data; %[-6 6];

% case tag additions
data_in.case_tag = [data_in.case_tag '_' ...
    num2str(input_params.years(1)) '_' num2str(input_params.years(end)) '_' ...
    'prcnt' num2str(input_params.percentile)];


if isfield(data_in, 'flag_daily_timeseries_smoothed_percentile') && data_in.flag_daily_timeseries_smoothed_percentile
    data_in.case_tag = [data_in.case_tag '_smooth'];
end

if input_params.flag_remove_trend
    data_in.case_tag = [data_in.case_tag '_noTrend'];
end

data_in.case_tag = [data_in.case_tag '_minLen_' num2str(input_params.delta_tstep) 'tsteps'];


if isfield(data_in, 'num_gap_tsteps_in_event')
    data_in.case_tag = [data_in.case_tag '_maxGap_' num2str(data_in.num_gap_tsteps_in_event) 'tsteps'];
end

input_params.fname_mat = data_in.case_tag;


% check that the time period selected does not include parial years
if isfield(data_in, 'flag_monthly_timeseries') && data_in.flag_monthly_timeseries && ...
        mod(size(data_in.data,3),12)~=0
    disp('One of the selected years does not have all the months!')
    ciao
elseif isfield(data_in, 'flag_monthly_timeseries') && ~data_in.flag_monthly_timeseries && ...
        isfield(data_in, 'flag_daily_timeseries') && ~data_in.flag_daily_timeseries
    writewhatneedschecking
end

%%
%%%%%%%%%% find info about MHW events (including when they start and at
%%%%%%%%%% which tstep they end)

if input_params.flag_save_MHW_info
    find_MHWs_info = NASAPO_2020_MHWs_detect_and_info_appendix_findMHWs(data_in,...
        input_params.delta_tstep,input_params.percentile,input_params.flag_remove_trend,input_params.fig_path);
    
    find_MHWs_info.years = input_params.years;
    
    disp('>>> Saving find_MHWs_info takes this number of seconds:')
    tic;
    save('-v7.3',[input_params.fig_path input_params.fname_mat '.mat'],'find_MHWs_info');
    toc;
    
elseif input_params.flag_save_MHW_ave_from_info
    load([input_params.fig_path input_params.fname_mat '.mat'])
end

if input_params.flag_save_MHW_ave_from_info
    
    disp('>>> Computing MHW event/onset averages takes these number of seconds:')
    tic;
    %%%% average selected variable in find_MHWs_info over events
    for i=1:length(data_in.data2ave_ALL_from_find_MHWs_info)
        clear data2ave data_select_event_ave
        data2ave = eval(['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{i}]);
        data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
            find_MHWs_info.start_tstep_msk,find_MHWs_info.end_tstep);
        eval(['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{i} '_eventAve = data_select_event_ave;'])
    end
    if input_params.flag_daily_timeseries && ~input_params.flag_monthly_timeseries
        %%%% average selected variable in find_MHWs_info over onset phase
        for i=1:length(data_in.data2ave_ALL_from_find_MHWs_info)
            clear data2ave data_select_event_ave
            data2ave = eval(['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{i}]);
            data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
                find_MHWs_info.start_tstep_msk,find_MHWs_info.peak_tstep);
            eval(['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{i} '_onsetAve = data_select_event_ave;'])
        end
        %%%% average selected variable in find_MHWs_info over decline phase
        for i=1:length(data_in.data2ave_ALL_from_find_MHWs_info)
            clear data2ave data_select_event_ave msk4ave endtstep4ave
            
            msk4ave      = find_MHWs_info.peak_tstep_msk;
            endtstep4ave = find_MHWs_info.end_tstep_stored_at_peak;
            
            % we want the decline phase to start at the timestep just after the peak timestep
            % if the peak timestep is the last timestep of the event, then
            % there is no decline phase in
            % NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs.m
            % Please note that the avg during the decline phase is saved at the
            % first timestep after the peak as the decline phase starts at
            % that timestep
            msk4ave(:,:,2:end) = msk4ave(:,:,1:end-1);
            msk4ave(:,:,1)     = 0;
            
            endtstep4ave(:,:,2:end) = endtstep4ave(:,:,1:end-1);
            endtstep4ave(:,:,1)     = nan;
            
            data2ave = eval(['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{i}]);
            data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
                msk4ave,endtstep4ave);
            eval(['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{i} '_declineAve = data_select_event_ave;'])
        end
    end
    
    %%%% average selected variable in data_in over events, before event, at
    %%%% the end of event
    if isfield(input_params, 'flag_ecco') && input_params.flag_ecco || ...
            isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint || ...
            isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily || ...
            isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box || ...
            isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box
        if isfield(data_in,'data2ave_ALL_from_data_in')
            for i=1:length(data_in.data2ave_ALL_from_data_in)
                clear data2ave data_select_event_ave* bfrbfr

                data2ave = eval(['data_in.' data_in.data2ave_ALL_from_data_in{i}]);

                %%%%%%%%%%%%%
                % average selected variable in data_in OVER events
                data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
                    find_MHWs_info.start_tstep_msk,find_MHWs_info.end_tstep);
                eval(['find_MHWs_info.' data_in.data2ave_ALL_from_data_in{i} '_eventAve = data_select_event_ave;'])

                if input_params.flag_daily_timeseries && ~input_params.flag_monthly_timeseries
                    %%%%%%%%%%%%%
                    % average selected variable in data_in during event onset
                    data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
                        find_MHWs_info.start_tstep_msk,find_MHWs_info.peak_tstep);
                    eval(['find_MHWs_info.' data_in.data2ave_ALL_from_data_in{i} '_onsetAve = data_select_event_ave;'])

                    %%%%%%%%%%%%%
                    % average selected variable in data_in during event decline
                    msk4ave      = find_MHWs_info.peak_tstep_msk;
                    endtstep4ave = find_MHWs_info.end_tstep_stored_at_peak;

                    % we want the decline phase to start at the timestep just after the peak timestep
                    % if the peak timestep is the last timestep of the event, then
                    % there is no decline phase in
                    % NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs.m
                    msk4ave(:,:,2:end) = msk4ave(:,:,1:end-1);
                    msk4ave(:,:,1)     = 0;

                    endtstep4ave(:,:,2:end) = endtstep4ave(:,:,1:end-1);
                    endtstep4ave(:,:,1)     = nan;

                    data_select_event_ave = NASAPO_2020_MHWs_detect_and_info_appendix_ave_overMHWs(data2ave,...
                        msk4ave,endtstep4ave);
                    eval(['find_MHWs_info.' data_in.data2ave_ALL_from_data_in{i} '_declineAve = data_select_event_ave;'])
                end
            end
        end
    end
    toc;
    
    
    disp('>>> Saving find_MHWs_info_withAVE takes this number of seconds:')
    tic;
    save('-v7.3',[input_params.fig_path input_params.fname_mat '_withAVE.mat'],'find_MHWs_info','input_params')
    toc;
    
    eval(['!rm ' input_params.fig_path input_params.fname_mat '.mat'])
else
    disp('>>> Upload the file with all the MHW info computed previously')
    tic;
    load([input_params.fig_path input_params.fname_mat '_withAVE.mat'])
    toc;
end

%%
%%%%%%%%%% let's make plots
close all

warning('off','all')

%%%%% maps
if input_params.flag_maps_summary && ...
        ~(isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) && ...
        ~(isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
    clear maps_to_plot
    %%%%%%%%% prepare input for the plotting functions
    % NOTE: we will use data_in.X, data_in.Y only for non-ECCO data
    maps_to_plot{1} = {'data_in.X','data_in.Y',...
        'nanmean(find_MHWs_info.events_duration_in_tsteps,3)',data_in.maplonlimit,data_in.cax_aveDuration,...
        ['Average duration of MHW events (' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
        ' tsteps) during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
        [input_params.fname_mat '_' ...
        'aveDuration'],input_params.fig_path};
    
    maps_to_plot{2} = {'data_in.X','data_in.Y',...
        'find_MHWs_info.events_number',data_in.maplonlimit,data_in.cax_num_events,...
        ['Number of MHW events (' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
        ' tsteps) during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
        [input_params.fname_mat '_' ...
        'numOFevents'],input_params.fig_path};
    
    jcasetags = {''};
    %     if ~flag_ecco_daily
    %         jcasetags = {'' '_at_tstep_before'};
    %     end
    
    % ave over event
    for jcase=1:1
        
        n = 0;
        
        for i=length(maps_to_plot)+1:length(maps_to_plot)+1+length(data_in.data2ave_ALL_from_find_MHWs_info)-1
            
            n = n + 1;
            maps_to_plot{i} = {'data_in.X','data_in.Y',...
                ['nanmean(find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{n} ...
                '_eventAve'  jcasetags{jcase} ',3)'],...
                data_in.maplonlimit,data_in.cax_data,...
                ['Mean ' data_in.data2ave_ALL_from_find_MHWs_info{n} ' over MHW events' ...
                strrep(jcasetags{jcase},'_',' ') ...
                ' (' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
                ' tsteps): average during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
                [input_params.fname_mat '_' ...
                'MEAN' data_in.data2ave_ALL_from_find_MHWs_info{n} '_ave' jcasetags{jcase}],input_params.fig_path};
        end
    end
    
    if input_params.flag_daily_timeseries && ~input_params.flag_monthly_timeseries
        % onset ave
        for jcase=1:1
            n = 0;
            
            for i=length(maps_to_plot)+1:length(maps_to_plot)+1+length(data_in.data2ave_ALL_from_find_MHWs_info)-1
                
                n = n + 1;
                maps_to_plot{i} = {'data_in.X','data_in.Y',...
                    ['nanmean(find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{n} ...
                    '_onsetAve'  jcasetags{jcase} ',3)'],...
                    data_in.maplonlimit,data_in.cax_data,...
                    ['Mean ' data_in.data2ave_ALL_from_find_MHWs_info{n} ' over MHW onset' ...
                    strrep(jcasetags{jcase},'_',' ') ...
                    ' (' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
                    ' tsteps): average during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
                    [input_params.fname_mat '_' ...
                    'MEAN' data_in.data2ave_ALL_from_find_MHWs_info{n} '_onset' jcasetags{jcase}],input_params.fig_path};
            end
        end
        
        % decline ave
        for jcase=1:1
            n = 0;
            
            for i=length(maps_to_plot)+1:length(maps_to_plot)+1+length(data_in.data2ave_ALL_from_find_MHWs_info)-1
                
                n = n + 1;
                maps_to_plot{i} = {'data_in.X','data_in.Y',...
                    ['nanmean(find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{n} ...
                    '_declineAve'  jcasetags{jcase} ',3)'],...
                    data_in.maplonlimit,data_in.cax_data,...
                    ['Mean ' data_in.data2ave_ALL_from_find_MHWs_info{n} ' over MHW decline' ...
                    strrep(jcasetags{jcase},'_',' ') ...
                    ' (' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
                    ' tsteps): average during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
                    [input_params.fname_mat '_' ...
                    'MEAN' data_in.data2ave_ALL_from_find_MHWs_info{n} '_decline' jcasetags{jcase}],input_params.fig_path};
            end
        end
    end
    
    % ave over event
    if isfield(input_params, 'flag_ecco') && input_params.flag_ecco || ...
            isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint || ...
            isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily || ...
            isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box || ...
            isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box
        for jcase=1:length(jcasetags)
            n = 0;
            if isfield(data_in,'data2ave_ALL_from_data_in')
                for i=length(maps_to_plot)+1:length(maps_to_plot)+1+length(data_in.data2ave_ALL_from_data_in)-1
                    n = n + 1;
                    maps_to_plot{i} = {'data_in.X','data_in.Y',...
                        ['nanmean(find_MHWs_info.' data_in.data2ave_ALL_from_data_in{n} ...
                        '_eventAve'  jcasetags{jcase} ',3)'],...
                        data_in.maplonlimit,data_in.cax_data_heat_budget_terms{n},...
                        ['Mean ' data_in.data2ave_ALL_from_data_in{n} ' over MHW events' ...
                        strrep(jcasetags{jcase},'_',' ') ...
                        '(' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
                        ' tsteps): average during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
                        [input_params.fname_mat '_' ...
                        'MEAN' data_in.data2ave_ALL_from_data_in{n} '_ave' jcasetags{jcase}],input_params.fig_path};
                end
            end
        end
    end
    
    if input_params.flag_daily_timeseries && ~input_params.flag_monthly_timeseries
        % ave over onset
        if isfield(input_params, 'flag_ecco') && input_params.flag_ecco || ...
                isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint || ...
                isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily || ...
                isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box || ...
                isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box
            for jcase=1:length(jcasetags)
                n = 0;
                if isfield(data_in,'data2ave_ALL_from_data_in')
                    for i=length(maps_to_plot)+1:length(maps_to_plot)+1+length(data_in.data2ave_ALL_from_data_in)-1
                        n = n + 1;
                        maps_to_plot{i} = {'data_in.X','data_in.Y',...
                            ['nanmean(find_MHWs_info.' data_in.data2ave_ALL_from_data_in{n} ...
                            '_onsetAve'  jcasetags{jcase} ',3)'],...
                            data_in.maplonlimit,data_in.cax_data_heat_budget_terms{n},...
                            ['Mean ' data_in.data2ave_ALL_from_data_in{n} ' over MHW onset' ...
                            strrep(jcasetags{jcase},'_',' ') ...
                            '(' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
                            ' tsteps): average during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
                            [input_params.fname_mat '_' ...
                            'MEAN' data_in.data2ave_ALL_from_data_in{n} '_onset' jcasetags{jcase}],input_params.fig_path};
                    end
                end
            end
        end
        
        % ave over decline
        if isfield(input_params, 'flag_ecco') && input_params.flag_ecco || ...
                isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint || ...
                isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily || ...
                isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box || ...
                isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box
            for jcase=1:length(jcasetags)
                n = 0;
                if isfield(data_in,'data2ave_ALL_from_data_in')
                    for i=length(maps_to_plot)+1:length(maps_to_plot)+1+length(data_in.data2ave_ALL_from_data_in)-1
                        n = n + 1;
                        maps_to_plot{i} = {'data_in.X','data_in.Y',...
                            ['nanmean(find_MHWs_info.' data_in.data2ave_ALL_from_data_in{n} ...
                            '_declineAve'  jcasetags{jcase} ',3)'],...
                            data_in.maplonlimit,data_in.cax_data_heat_budget_terms{n},...
                            ['Mean ' data_in.data2ave_ALL_from_data_in{n} ' over MHW decline' ...
                            strrep(jcasetags{jcase},'_',' ') ...
                            '(' num2str(input_params.percentile) '%, \Deltat=' num2str(input_params.delta_tstep)...
                            ' tsteps): average during ' num2str(input_params.years(1)) '-' num2str(input_params.years(end))],...
                            [input_params.fname_mat '_' ...
                            'MEAN' data_in.data2ave_ALL_from_data_in{n} '_decline' jcasetags{jcase}],input_params.fig_path};
                    end
                end
            end
        end
    end
    %%%%%%%%%
    
    %%%%%%%%% loop over all the plots to make
    for i=1:length(maps_to_plot)
        close all
        figure;set(gcf,'color','w','position',[0       0      1800        1200]);
        
        if ~(isfield(data_in, 'flag_dim_lon_lat_time') && data_in.flag_dim_lon_lat_time)
            clear bfr*
            
            % ecco_xc,ecco_yc,ecco_data should all have the same
            bfr_ecco_data   = reshape(eval(maps_to_plot{i}{3}),data_in.size_source(1:3));
            bfr_ecco_xc     = data_in.xc(:,:,:,1);
            bfr_ecco_yc     = data_in.yc(:,:,:,1);
            
            % put on a regular grid
            data_on_grid = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_to_regular_grid(...
                bfr_ecco_xc,bfr_ecco_yc,bfr_ecco_data);
            
            % plot maps
            NASAPO_2020_MHWs_detect_and_info_appendix_plot_map(...
                data_on_grid.X,...
                data_on_grid.Y,data_on_grid.data,...
                maps_to_plot{i}{4},...
                maps_to_plot{i}{5},maps_to_plot{i}{6},...
                maps_to_plot{i}{7},maps_to_plot{i}{8});
        else
            NASAPO_2020_MHWs_detect_and_info_appendix_plot_map(...
                eval(maps_to_plot{i}{1}),eval(maps_to_plot{i}{2}),...
                eval(maps_to_plot{i}{3}),maps_to_plot{i}{4},...
                maps_to_plot{i}{5},maps_to_plot{i}{6},...
                maps_to_plot{i}{7},maps_to_plot{i}{8});
        end
    end
end

%% lets plot maps at each timestep
%%%%%% let's plot info for heat budget terms
if input_params.flag_maps_at_each_tstep && ...
        ((isfield(input_params, 'flag_ecco') && input_params.flag_ecco) || ...
        (isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint))
    clear bfr*
    bfr_msk=find_MHWs_info.data_used4MHWs>=find_MHWs_info.data_percentile3d;
    
    data_in.bfr1           = data_in.G_total;
    data_in.bfr1(~bfr_msk) = nan;
    data_in.bfr2           = data_in.G_forcing;
    data_in.bfr2(~bfr_msk) = nan;
    
    data_in.bfr3           = data_in.G_diffusion_conv;
    data_in.bfr3(~bfr_msk) = nan;
    data_in.bfr3v          = data_in.dif_vConv;
    data_in.bfr3v(~bfr_msk)= nan;
    data_in.bfr3h          = data_in.G_diffusion_conv-data_in.dif_vConv;
    %data_in.bfr3h(~bfr_msk)= nan;
    
    data_in.bfr4           = data_in.G_advection_conv;
    data_in.bfr4(~bfr_msk) = nan;
    data_in.bfr4v          = data_in.adv_vConv;
    data_in.bfr4v(~bfr_msk)= nan;
    data_in.bfr4h          = data_in.G_advection_conv-data_in.adv_vConv;
    %data_in.bfr4h(~bfr_msk)= nan;
    
    vars2plot_cax_ecco_heat = data_in.cax_data_heat_budget_terms_monthly;
    
    vars2plot = {'find_MHWs_info.data_used4MHWs' ...
        'data_in.G_total' 'data_in.bfr1' ...
        'data_in.G_forcing' 'data_in.bfr2' ...
        'data_in.G_diffusion_conv' 'data_in.bfr3' ...
        'data_in.dif_vConv' 'data_in.bfr3v' ...
        'data_in.bfr3h' ...
        'data_in.G_advection_conv' 'data_in.bfr4' ...
        'data_in.adv_vConv' 'data_in.bfr4v' ...
        'data_in.bfr4h' ...
        };
    
    vars2plot_title = {'data_used_for_mhw_def' ...
        'G_total' 'G_total only at MHW' ...
        'G_forcing' 'G_forcing only at MHW' ...
        'G_diffusion_conv' 'G_diffusion_conv only at MHW' ...
        'dif_vConv' 'dif_vH only at MHW' ...
        'dif_hConv'  ...
        'G_advection_conv' 'G_advection_conv only at MHW' ...
        'adv_vConv' 'adv_vH only at MHW' ...
        'adv_hConv'  ...
        };
    
    vars2plot_cax = {data_in.cax_data data_in.cax_data data_in.cax_data []};
    
    NASAPO_2020_MHWs_detect_and_info_appendix_maps_at_tsteps(data_in,find_MHWs_info,vars2plot,...
        vars2plot_cax_ecco_heat,vars2plot_title,'_heatBudget',input_params.fig_path)
    
    fields = {'bfr1','bfr2','bfr3','bfr4'};
    data_in = rmfield(data_in,fields);
end

%%%% let's plot info that is available for all the products
if input_params.flag_maps_at_each_tstep
    data_in.bfr1=find_MHWs_info.data_used4MHWs;
    data_in.bfr1(~find_MHWs_info.data_mhw_tstep_msk) = nan;
    
    vars2plot = {'find_MHWs_info.data_used4MHWs' 'data_in.bfr1' ...
        'find_MHWs_info.data_used4MHWs_eventAve' ...
        'find_MHWs_info.events_duration_in_tsteps'};
    
    vars2plot_title = {'data_used_for_mhw_def' ['data_used_for_mhw_def (only ge ' ...
        num2str(input_params.percentile) 'th prcnt)'] ...
        'data_used_for_mhw_def: ave over event' ...
        'event duration (plot at event start)'};
    
    vars2plot_cax = {data_in.cax_data data_in.cax_data data_in.cax_data []};
    
    NASAPO_2020_MHWs_detect_and_info_appendix_maps_at_tsteps(data_in,find_MHWs_info,vars2plot,...
        vars2plot_cax,vars2plot_title,'',input_params.fig_path)
    fields = {'bfr1'};
    data_in = rmfield(data_in,fields);
end
%% LAST SECTION OF THE CODE: bar plots
%tiledlayout(length(data_in.i_all),1)

%%%%%%%%%% CHANGE here to change i_all, j_all without having to run the
%%%%%%%%%% script from the beginning
%
% uncomment here to set other i_all, j_all of interest (for now only the
% usual two points are included, you can include more)
%data_in.i_all = [6586 5548]; % 156.5W (around 204 east), 35.1757S
%data_in.j_all = [9 9];
%%%%%%%%%%

% %++++++++% CHANGE here to zoom in on a specific time period
% %++++++++%
% % in x_lim_select_ALL, you can include [] to plot the full time series
% % or you can include e.g. [datenum(2015,3,1) datenum(2016,6,30)] to plot only a specific
% % time period
%

x_lim_select_ALL = {[] ...
    [datenum(2015,3,1) datenum(2016,6,30)]...
    [datenum(2008,3,1) datenum(2011,3,31)]...
    [datenum(2012,3,1) datenum(2016,8,31)]...
    };

if (isfield(input_params, 'flag_ecco') && input_params.flag_ecco) || ...
        (isfield(input_params, 'flag_ecco_zint') && input_params.flag_ecco_zint) || ...
        (isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily) || ...
        (isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) || ...
        (isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
    
    for i_x_lim_select_ALL=1:length(x_lim_select_ALL)
        clear x_lim_select
        x_lim_select = x_lim_select_ALL{i_x_lim_select_ALL};
        
        for nij=1:length(data_in.i_all)
            i=data_in.i_all(nij);
            j=data_in.j_all(nij);
            close all
            figure;set(gcf,'color','w','position',[0       0      3000        1500]);
            clear d2pl_all
            
            n = 0;
            
            %     n = 1;
            %     %d2pl_all{n} = 'data_in.data-repmat(mean(data_in.data,3),[1 1 size(data_in.data,3)])';
            %     %d2pl_all_title{n} = 'Source DATA after removing time mean';
            %     d2pl_all{n} = 'data_in.data-273.15';
            %     d2pl_all_title{n} = 'Source DATA in degC';
            
            for ivar=1:length(data_in.data2ave_ALL_from_find_MHWs_info)
                n = n+1;
                d2pl_all{n} = ['find_MHWs_info.' data_in.data2ave_ALL_from_find_MHWs_info{ivar}];
                d2pl_all_title{n} = data_in.data2ave_ALL_from_find_MHWs_info{ivar};
            end
            
            if isfield(data_in,'data2ave_ALL_from_data_in')
                for ivar=1:length(data_in.data2ave_ALL_from_data_in)
                    if ~contains(data_in.data2ave_ALL_from_data_in{ivar},'vConv')
                        n = n+1;
                        d2pl_all{n} = ['data_in.' data_in.data2ave_ALL_from_data_in{ivar}];
                        d2pl_all_title{n} = data_in.data2ave_ALL_from_data_in{ivar};
                    end
                end
            end
            
            tiledlayout(length(d2pl_all),1)
            
            fcol = {[.75 .75 .75],[255 178 102]./255,'none'};% [0 .9 0]
            ecol = {[.5 .5 .5] 'k' 'r'};
            
            if (isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily) || ...
                    (isfield(input_params, 'flag_noaa_oisst_daily') && input_params.flag_noaa_oisst_daily) || ...
                    (isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) || ...
                    (isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
                
                fcol = {[.75 .75 .75],[255 178 102]./255,'none', 'none', 'none'};% [204 255 204]./255
                ecol = {'none' 'none' 'r' 'b' 'm'};
            end
            
            for ivar=1:length(d2pl_all)
                nexttile()
                clear d2pl0
                d2pl0 = eval(d2pl_all{ivar});
                bar(data_in.data_datenum,squeeze(d2pl0(i,j,:)),'facecolor',fcol{1},'edgecolor',ecol{1})
                hold on
                
                % let's indicate which tsteps are above the selected percentile
                d2pl = d2pl0;
                d2pl(~find_MHWs_info.data_mhw_tstep_msk) = nan;
                bar(data_in.data_datenum,squeeze(d2pl(i,j,:)),'facecolor',fcol{2},'edgecolor',ecol{2})
                
                d2pl = d2pl0;
                d2pl(~find_MHWs_info.start_tstep_msk) = nan;
                bar(data_in.data_datenum,squeeze(d2pl(i,j,:)),'facecolor',fcol{3},'edgecolor',ecol{3},...
                    'linewidth',1.5)
                
                % disp('Code to complete here and include percentile, peak, ...')
                
                if (isfield(input_params, 'flag_ecco_daily') && input_params.flag_ecco_daily) || ...
                        (isfield(input_params, 'flag_noaa_oisst_daily') && input_params.flag_noaa_oisst_daily) || ...
                        (isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) || ...
                        (isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
                    
                    clear d2pl ind2pl
                    d2pl = d2pl0;
                    d2pl(~find_MHWs_info.peak_tstep_msk) = nan;
                    bar(data_in.data_datenum,squeeze(d2pl(i,j,:)),'facecolor',fcol{5},'edgecolor',ecol{5},...
                        'linewidth',1.5,'linestyle',':')
                    
                    clear d2pl ind2pl
                    d2pl = d2pl0;
                    ind2pl = reshape(find_MHWs_info.peak_tstep(i,j,~isnan(find_MHWs_info.peak_tstep(i,j,:))),[],1);
                    %             bar(data_in.data_datenum(ind2pl),squeeze(d2pl(i,j,ind2pl)),'facecolor',fcol{4},'edgecolor',ecol{4},...
                    %                 'linewidth',1.5,'linestyle',':')
                    plot(data_in.data_datenum(ind2pl),squeeze(d2pl(i,j,ind2pl)),'linestyle','none','marker','+','color',ecol{4})
                    
                    %                 data_out.peak_tstep_msk = peak_tstep_msk;
                    %                 data_out.end_tstep_stored_at_peak = end_tstep_stored_at_peak;
                    
                    clear d2pl ind2pl
                    d2pl = d2pl0;
                    ind2pl = reshape(find_MHWs_info.end_tstep_stored_at_peak(i,j,~isnan(find_MHWs_info.end_tstep_stored_at_peak(i,j,:))),[],1);
                    %             bar(data_in.data_datenum(ind2pl),squeeze(d2pl(i,j,ind2pl)),'facecolor',fcol{4},'edgecolor',ecol{4},...
                    %                 'linewidth',1.5,'linestyle',':')
                    plot(data_in.data_datenum(ind2pl),squeeze(d2pl(i,j,ind2pl)),'linestyle','none','marker','o','color',ecol{4})
                    
                    clear d2pl ind2pl
                    d2pl = d2pl0;
                    ind2pl = reshape(find_MHWs_info.end_tstep(i,j,~isnan(find_MHWs_info.end_tstep(i,j,:))),[],1);
                    %             bar(data_in.data_datenum(ind2pl),squeeze(d2pl(i,j,ind2pl)),'facecolor',fcol{4},'edgecolor',ecol{4},...
                    %                 'linewidth',1.5,'linestyle',':')
                    plot(data_in.data_datenum(ind2pl),squeeze(d2pl(i,j,ind2pl)),'linestyle','none','marker','^','color',ecol{4})
                    
                end
                
                %hline data_percentile3d
                if ivar==1
                    plot(data_in.data_datenum,squeeze(find_MHWs_info.data_percentile3d(i,j,:)),'--k')
                end
                
                if ~data_in.flag_dim_lon_lat_time && ...
                    ~(isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) && ...
                        ~(isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
                    bfrx = reshape(data_in.xc,size(d2pl0));
                    bfry = reshape(data_in.yc,size(d2pl0));
                    bfrx = bfrx(i,j,1);
                    bfry = bfry(i,j,1);
                    
                    %             % to find i,j of interest (see next figure as it
                    %             is more helpful)
                    %             % m = 9;plot(bfrx(:,m,1),bfry(:,m,1),'*')
                    %             % find(bfrx(:,m,1)==215.5)
                    %             % plot(find(bfrx(:,m,1)==215.5),bfry(find(bfrx(:,m,1)==215.5),m,1),'*')
                    %
                    %             %%%% to find i_all, j_all of interest: maps
                    %             figure;set(gcf,'color','w','position',[0       0      1800        1200]);
                    %             for iiii=1:size(bfrx,2)
                    %                 nexttile()
                    %                 scatter(bfrx(:,iiii,1),bfry(:,iiii,1),10,1:length(bfrx(:,iiii,1)),'filled')
                    %                 colorbar
                    %                 title(['j_all = ' num2str(iiii) ', i_all in color'],...
                    %                     'interpreter','none')
                    %                 xlabel('Longitudes')
                    %                 ylabel('Latitudes')
                    %                 set(gca,'fontsize',14)
                    %                 nexttile()
                    %                 scatter(bfrx(:,iiii,1),bfry(:,iiii,1),10,1:length(bfrx(:,iiii,1)),'filled')
                    %                 colorbar
                    %                 title(['j_all = ' num2str(iiii) ', i_all in color'],...
                    %                     'interpreter','none')
                    %                 xlabel('Longitudes')
                    %                 ylabel('Latitudes')
                    %                 set(gca,'fontsize',14)
                    %             end
                    
                elseif ~(isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) && ...
                        ~(isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
                    
                    bfrx = data_in.x(i);
                    bfry = data_in.y(j);
                else
                    bfrx = i;
                    bfry = j;
                end
                
                if ~(isfield(input_params, 'flag_ecco_daily_box') && input_params.flag_ecco_daily_box) && ...
                        ~(isfield(input_params, 'flag_noaa_oisst_daily_box') && input_params.flag_noaa_oisst_daily_box)
                    
                    title([d2pl_all_title{ivar} ', Lon : ' num2str(bfrx) ', Lat : ' num2str(bfry) ...
                        ', percentile : ' num2str(input_params.percentile)],'interpreter','none')
                else
                    title([d2pl_all_title{ivar} ', index dim1: ' num2str(bfrx) ', index dim2: ' num2str(bfry) ...
                        ', percentile : ' num2str(input_params.percentile)],'interpreter','none')
                end
                set(gca,'fontsize',24)
                
                x_lim_select_tag = '';
                if ~isempty(x_lim_select)
                    % uncomment here to set different limits for the x-axis
                    set(gca,'xlim',x_lim_select)
                    x_lim_select_tag = ['_' num2str(year(min(x_lim_select))) ...
                        num2str(month(min(x_lim_select)),'%02d') ...
                        num2str(day(min(x_lim_select)),'%02d') '_' ...
                        num2str(year(max(x_lim_select))) ...
                        num2str(month(max(x_lim_select)),'%02d') ...
                        num2str(day(max(x_lim_select)),'%02d') ...
                        ];
                end
                
                
                datetick('x','dd-mm-yy','keepticks')
                
                
            end
            
            %             fname_mat = [data_in.case_tag bfr_smoothed_percentile '_minLen_' num2str(input_params.delta_tstep) 'tsteps'];
            
            print('-dpng',[input_params.fig_path input_params.fname_mat '_i' num2str(i) '_j' num2str(j) x_lim_select_tag '.png'])
        end
    end
end

% !rm ./NASAPO_2020_MHWs_input_params.mat
return



