function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily_box(input_params)

% input_vars.ecco_name= 'ECCOv4r4';
input_vars.ecco_name = input_params.ecco_name; %'ECCOv4r4';

% % input_vars.num_of_tiles = 13;
% input_vars.num_of_tiles = 13;

% delta days before/after to calculate percentile - should be the same for
% OISST daily
data_out.delta_days_4_percentile = input_params.delta_days_4_percentile; %5;
data_out.num_gap_tsteps_in_event = input_params.num_gap_tsteps_in_event;

% % expected size of the variables (for some reason, G_forcing has different
% % dimensions from other variables)
% size1 = [90    90    input_vars.num_of_tiles input_vars.num_days];
%
% size2 = [];

%%% adjust for OHC
% if ~isempty(zlev)
%     data_out.case_tag = [input_vars.ecco_name '_' tag_field '_zlev' num2str(zlev,'%02d') '_daily'];
% else
%     data_out.case_tag = [input_vars.ecco_name '_' tag_field '_daily'];
% end

data_out.case_tag = [input_vars.ecco_name '_' input_params.tag_field '_daily_box'];

input_vars.fpath     = input_params.user_path_ecco;

switch input_params.tag_field
    case 'heat'
        %         input_vars.var_names = {'DATA' ...
        %             'G_total' 'G_forcing' ...
        %             'G_advection_conv' 'adv_vConv' ...
        %             'G_advection_conv_mer' 'G_advection_conv_zon' ...
        %             'G_diffusion_conv' 'dif_vConv' ...
        %              }; % 'adv_meridH'
        
        input_vars.var_names = {'DATA' ...
            'G_total' 'G_forcing' ...
            'G_advection' 'adv_vConv' ...
            'G_diffusion' 'dif_vConv' ...
            }; % 'adv_meridH'
        % input_vars.var_names = {'DATA' ...
        %     'G_total' 'G_forcing' ...
        %     'G_advection' ...
        %     'G_diffusion' ...
        %      }; % 'adv_meridH'
        data_out.cax_percentile = input_params.cax_percentile; %[0 6];
        
    otherwise
        to complete
end
% %input_vars.ecco_other_vars = {'area' 'XC_lon' 'YC_lat' 'Z_depth' 'vol'};
% input_vars.ecco_other_vars = {'area' 'XC_lon' 'YC_lat' 'Z_depth'};

data_out.var_names = input_vars.var_names;
%
data_out.flag_dim_lon_lat_time = input_params.flag_dim_lon_lat_time;
data_out.flag_monthly_timeseries = input_params.flag_monthly_timeseries;
data_out.flag_daily_timeseries   = input_params.flag_daily_timeseries;
data_out.flag_daily_timeseries_smoothed_percentile = input_params.flag_daily_timeseries_smoothed_percentile;

%

% load other vars
% vars_all = input_vars.ecco_other_vars;
% vars_all{end+1} = input_vars.var_names{1};

% vars_all = horzcat(input_vars.ecco_other_vars,input_vars.var_names);

vars_all = input_vars.var_names;

time_flag = true;
for i=1:length(vars_all)
    clear bfr
    bfr = read_vars_in_ncfile([input_vars.fpath input_vars.ecco_name ...
        '_' vars_all{i} '_' ...
        input_params.tag_years_filename input_params.tag_file '.nc']);
    
    % 17 June: for the box data, I suggest to create files with the main
    % variable named as the variable in the file name... when that happens,
    % the command at line 108 can be used instead of that on line 107
    clear data_bfr
    %     data_bfr = bfr.xarray_dataarray_variable;
    eval(['data_bfr = bfr.' vars_all{i} ';'])
    %
    % NOTE: data needs to be 3d... since the input file includes a 2d
    % variable, we make data 3d here:
    if length(size(data_bfr))==2
        if size(data_bfr, 2)==1
            data_bfr = data_bfr';
        end
        disp('Let''s move the 2d input data to 3d')
        data(:,1,:) = data_bfr;
    else
        checkwhatwrong
    end
    
    % save time variable
    if time_flag
        % time
        switch bfr.time_units
            case 'hours since 2004-01-01 12:00:00'
                time_datenum = datenum(2004,1,1,12,0,0) + bfr.time./24;
            case 'hours since 1992-01-02 12:00:00'
                time_datenum = datenum(1992,1,2,12,0,0) + bfr.time./24;
            case 'hours since 1992-01-02T12:00:00'  % Handle the format with 'T'
                time_datenum = datenum(1992,1,2,12,0,0) + bfr.time./24;
            otherwise
                case2code
        end
        msk_years    = year(time_datenum)>=min(input_params.years) & ...
            year(time_datenum)<=max(input_params.years);
        
        data_out.data_datenum = time_datenum(msk_years);
        
        time_flag = false;
    end
    
    % select data in the years of interest
    data                    = data(:,:,msk_years);
    
    
    % except for DATA, remove seasonal cycle and trend
    if ~strcmp(vars_all{i},'DATA')
        % remove seasonal cycle for a daily timeseries
        data = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_daily(data, data_out.data_datenum, ...
            data_out.delta_days_4_percentile);
        
        % detrend
        data = reshape(detrend(reshape(data,size(data,1)*size(data,2),[])')',...
            size(data,1),size(data,2),[]);
    end
    
    %
    eval([vars_all{i} ' = data;'])
    clear data
end

data_out.size_source(4) = sum(msk_years);

% % data_out.z_source    = bfr.depth;


% data_out.xc          = repmat(XC_lon,[1 1 1 data_out.size_source(4)]);
% % convert longitudes to be from 20 to less than 380
% data_out.xc(data_out.xc<20) = data_out.xc(data_out.xc<20) + 360;
%
% data_out.yc          = repmat(YC_lat,[1 1 1 data_out.size_source(4)]);
% data_out.area4d      = repmat(area_bfr,[1 1 1 data_out.size_source(4)]);
% if ~contains(tag_file,'zint') %isempty(tag_file)
%     data_out.Z_depth     = Z_depth(zlev);
%     data_out.zlev        = zlev;
% end

% data_out.Z_depth     = Z_depth(zlev);
% data_out.zlev        = zlev;

data_out.data = DATA;

n = 0;
for i=1:length(input_vars.var_names)
    if ~strcmp(input_vars.var_names{i},'DATA')
        %         eval(['data_out.' input_vars.var_names{i} ' = ' ...
        %             input_vars.var_names{i}  '(:,:,msk_years);'])
        eval(['data_out.' input_vars.var_names{i} ' = ' ...
            input_vars.var_names{i}  ';'])
        n = n + 1;
        data_out.data2ave_ALL_from_data_in{n} = input_vars.var_names{i};
    end
end



% % data_out.x = squeeze(data_in.lon);
% % data_out.y = squeeze(data_in.lat);
% data_out.maplonlimit = [min(data_out.xc(:)) max(data_out.xc(:))];

% [data_out.X,data_out.Y] = ndgrid(data_out.x,data_out.y);

% define indices to look at in the timeseries plot at the end
data_out.i_all = input_params.i_all; %[1 2]; % 156.5W (around 204 east), 35.1757S
data_out.j_all = input_params.j_all; %[1 1];

% caxis for percentile plot


% variables to average over MHW event-relevant timesteps (e.g. over all the
% timesteps for an event, over the final ones only...)
% NOTE: see also data_out.data2ave_ALL_from_data_in above
data_out.data2ave_ALL_from_find_MHWs_info = {'data_used4MHWs'};

return
